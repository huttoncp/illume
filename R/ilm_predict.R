## illume: prediction.
##
## Smooths: smooth2random() reparameterises the basis, so a new-data design must
## be pushed through the SAME transform.  From
##   b_original = trans.U %*% (trans.D * c(b_random, b_fixed))
## (verified numerically to 1e-12; trans.D is not optional) it follows that
##   M = PredictMat(sm, newdata) %*% trans.U %*% diag(trans.D)
## splits column-wise into the random blocks then the fixed null-space block.
##
## WHICH TERMS ARE HELD vs INTEGRATED.  After smooth2random a spline looks like
## a random effect but is not a population to average over -- it IS the mean
## structure, so basis terms are always evaluated.  Grouping factors and AR terms
## are the populations: zeroed (conditional) or integrated out (marginal).

#' Evaluate a smooth's design at new covariate values
#'
#' A smooth is stored in a reparameterised form (see [ilm_smooth()]), so
#' predicting at new values means pushing the new design through the *same*
#' transformation the fit used. `mgcv::PredictMat()` rebuilds the original basis;
#' multiplying by the stored transformation recovers the fixed and random blocks
#' in the order the model expects.
#'
#' @param sob A stored smooth, an element of `fit$smooths`.
#' @param newdata A data frame containing the smooth's variables.
#' @return A list with `Xr` (a list of random blocks) and `Xf` (fixed
#'   null-space columns).
#' @keywords internal
#' @noRd
ilm_smooth_design <- function(sob, newdata) {
  ## a smooth of a column that needs backticks was built on a stand-in name,
  ## and has to be evaluated on one too
  newdata <- ilm_add_standins(newdata, sob$name_map)
  X0 <- mgcv::PredictMat(sob$sm, newdata)
  D  <- sob$re$trans.D
  M  <- X0 %*% sob$re$trans.U %*% diag(D, length(D))
  nrk <- vapply(sob$rand, ncol, 1L); nr <- sum(nrk); off <- c(0L, cumsum(nrk))
  list(Xr = lapply(seq_along(nrk), function(k) M[, (off[k] + 1L):off[k + 1L], drop = FALSE]),
       Xf = if (ncol(sob$Xf)) M[, nr + seq_len(ncol(sob$Xf)), drop = FALSE] else NULL)
}

#' Fitted random coefficients for one term
#'
#' @param object A fitted `"ilm_model"` object.
#' @param k Integer index of the random term.
#' @param bvec Optional replacement for the estimated random effects, used when
#'   simulating to propagate their uncertainty.
#' @return A matrix of random coefficients.
#' @keywords internal
#' @noRd
ilm_Bhat_term <- function(object, k, bvec = NULL) {
  v <- if (is.null(bvec)) object$sdr$par.random else bvec
  matrix(v[object$b_idx[[k]]], object$nlk[k] * object$dk[k], object$wk[k])
}

#' Build the fixed design for new data
#'
#' Reconstructs the design matrix in exactly the column order the fit used,
#' including any smooth null-space columns, using the stored `terms`, `xlev` and
#' `contrasts` so that factor levels and contrasts match the original fit.
#'
#' @param object A fitted `"ilm_model"` object.
#' @param newdata A data frame.
#' @return A list with `X` and the per-smooth designs.
#' @keywords internal
#' @noRd
ilm_newX <- function(object, newdata) {
  mt <- stats::delete.response(object$terms)
  mf <- stats::model.frame(mt, newdata, xlev = object$xlev)
  X  <- ilm_drop_intercept(
    stats::model.matrix(mt, mf, contrasts.arg = object$contrasts), object)
  sd_list <- list()
  for (lab in names(object$smooths)) {
    sd_list[[lab]] <- ilm_smooth_design(object$smooths[[lab]], newdata)
    if (!is.null(sd_list[[lab]]$Xf)) X <- cbind(X, sd_list[[lab]]$Xf)
  }
  ## A flexible parametric baseline puts spline columns at the front of the
  ## fitted design. They are a function of the response, not of the covariates,
  ## so new data cannot supply them and the comparison is against what is left.
  nrp <- if (is.null(object$rp)) 0L else length(object$rp$cols)
  if (ncol(X) != ncol(object$X) - nrp)
    stop("new-data design has ", ncol(X), " columns but the fit has ",
         ncol(object$X) - nrp,
         if (nrp) paste0(" besides its ", nrp, " baseline spline columns") else "",
         call. = FALSE)
  list(X = X, smooths = sd_list)
}

#' Linear predictor on the sum-to-zero scale
#'
#' Smooth terms are always evaluated; grouping and AR terms are left at zero.
#' That asymmetry is intentional. After reparameterisation a smooth *looks* like
#' a random effect, but it is not a population to average over -- it is part of
#' the mean structure. Grouping factors are the populations.
#'
#' @param object A fitted `"ilm_model"` object.
#' @param nd Output of [ilm_newX()].
#' @param beta Optional replacement fixed effects.
#' @param bvec Optional replacement random effects.
#' @return A matrix with `N` rows and `C` columns.
#' @keywords internal
#' @noRd
ilm_eta <- function(object, nd, beta = NULL, bvec = NULL) {
  if (is.null(beta)) beta <- object$beta
  eta <- nd$X %*% beta
  for (k in seq_along(object$re)) {
    e <- object$re[[k]]
    if (e$kind != "basis") next                        # populations: left at zero
    lab <- sub("\\.[0-9]+$", "", names(object$re)[k])
    sdl <- nd$smooths[[lab]]
    if (is.null(sdl)) sdl <- nd$smooths[[names(object$re)[k]]]
    if (is.null(sdl)) stop("no stored smooth for term '", names(object$re)[k], "'")
    blk <- if (length(sdl$Xr) == 1L) sdl$Xr[[1]] else {
      i <- suppressWarnings(as.integer(sub(".*\\.", "", names(object$re)[k])))
      sdl$Xr[[if (is.na(i)) 1L else i]]
    }
    ctb <- blk %*% ilm_Bhat_term(object, k, bvec)
    if (identical(object$re_struct[[k]]$type, "rr")) ctb <- ctb %*% t(object$Lambda[[k]])
    eta <- eta + ctb
  }
  eta
}

## ---- drawing a random effect the way the model actually stores one ---------
##
## WHY THIS EXISTS: a random effect is not one number per group. For a term
## with `dk` within-group dimensions -- `(1 + time | id)` has two -- and `C`
## linear-predictor dimensions, one group's effect is a `dk x C` MATRIX, and
## the objective gives it row covariance `Sigma_d` and column covariance
## `Sigma`. That is a matrix normal, so
##
##     U = A %*% Z %*% B      Z iid N(0, 1), A A' = Sigma_d, B' B = Sigma
##
## and the contribution at a row is `z_row %*% U`, where `z_row` is that row's
## value of the bar's left-hand side.
##
## Treating the whole thing as a per-group SHIFT -- drawing `Sigma` alone and
## adding it to every row -- is right if and only if `dk == 1`. With a random
## slope it drops the slope variance entirely, and it also makes the shift
## constant across rows when the true spread GROWS with distance from where
## the slope is centred. Both matter: the marginal mean of a nonlinear link
## depends on how much is being integrated over.
##
## When `dk == 1` these functions reproduce the old draw exactly, including
## the random number stream, so nothing that was already right moves.

## The two factors, for term k.
#' @keywords internal
#' @noRd
ilm_re_factors <- function(object, k) {
  d  <- if (is.null(object$dk)) 1L else object$dk[k]
  Sd <- if (d > 1L) object$Sigma_d[[names(object$re)[k]]] else NULL
  ## ilm_msqrt() rather than a Cholesky: an intercept-slope correlation can
  ## reach +/-1, where the matrix is singular and chol() stops
  A <- if (is.null(Sd)) matrix(1, 1L, 1L) else ilm_msqrt(Sd)
  B <- if (identical(object$re_struct[[k]]$type, "rr"))
    t(object$Lambda[[k]]) else ilm_msqrt(object$Sigma[[k]])
  list(A = A, B = B, d = if (is.null(Sd)) 1L else nrow(Sd))
}

## The bar's left-hand side, evaluated on whatever rows are being predicted.
## `object$re[[k]]$Z` is the design for the rows the model was FITTED to, so
## new data needs the bar evaluated again. Smooth bases are prepended to the
## random-effect list while `bars` holds only the bars, so the index into
## `bars` is the position among the NON-basis terms.
#' @keywords internal
#' @noRd
ilm_re_design <- function(object, k, data, d) {
  if (d == 1L) return(matrix(1, nrow(data), 1L))
  gk <- which(vapply(object$re, function(e) !identical(e$kind, "basis"), TRUE))
  j <- match(k, gk)
  b <- if (is.na(j)) NULL else object$bars[[j]]
  if (is.null(b)) return(NULL)
  ## from the expression, not pasted text: a lone backticked slope variable
  ## deparses bare, fails to parse, and the tryCatch below used to turn that
  ## into a quiet NULL -- a random slope dropped from the draws without a word
  env <- environment(object$formula)
  if (is.null(env)) env <- parent.frame()
  Z <- tryCatch(
    stats::model.matrix(ilm_one_sided(b[[2L]], env), data),
    error = function(e) NULL)
  if (is.null(Z) || ncol(Z) != d) NULL else Z
}

## `ndraw` draws of one group's effect, as a list of dk x C matrices.
##
## The normals come from one column-major block of the same shape the old code
## used, so when dk == 1 the stream is identical and every number that was
## already correct stays exactly where it was.
#' @keywords internal
#' @noRd
ilm_re_draws <- function(A, B, ndraw, C) {
  d <- nrow(A)
  Z <- matrix(stats::rnorm(ndraw * d * C), ndraw, d * C)
  lapply(seq_len(ndraw), function(m) A %*% matrix(Z[m, ], d, C) %*% B)
}

## Gauss-Hermite nodes and weights for E[f(Z)], Z ~ N(0, 1), by the
## Golub-Welsch eigenvalue method: E[f(Z)] ~ sum(w * f(z)). With one linear
## predictor, a row's whole latent contribution is a single normal, so its
## average over the random effects is a one-dimensional integral, and 40 points
## do it to the digits a probability is ever quoted to.
#' @keywords internal
#' @noRd
ilm_gh <- function(n = 40L) {
  i <- seq_len(n - 1L)
  J <- matrix(0, n, n)
  J[cbind(i, i + 1L)] <- J[cbind(i + 1L, i)] <- sqrt(i / 2)
  e <- eigen(J, symmetric = TRUE)
  list(z = sqrt(2) * e$values, w = e$vectors[1L, ]^2)
}

#' Convert linear predictors to category probabilities
#'
#' Applies the softmax (multinomial logistic) transform so each row gives
#' probabilities summing to one.
#'
#' @param eta Matrix of linear predictors on the sum-to-zero scale.
#' @param Tc Sum-to-zero contrast matrix.
#' @return A matrix of probabilities with `J` columns.
#' @keywords internal
#' @noRd
ilm_softmax_J <- function(eta, Tc) { P <- exp(eta %*% t(Tc)); P / rowSums(P) }

#' Draw fixed and random parameters together
#'
#' Samples from the joint distribution of all parameters using the joint
#' precision matrix, which requires a fit made with `joint = TRUE`.
#'
#' Drawing jointly matters for smooths. A smooth's unpenalised (null-space) term
#' and its penalised coefficients are strongly correlated -- they trade off
#' against each other -- so varying one while holding the other fixed produces
#' far more variability in the fitted curve than the model actually implies.
#' Joint draws respect that correlation.
#'
#' @param object A fitted `"ilm_model"` object with `jointPrecision`.
#' @param nsim Integer. Number of draws.
#' @param seed Integer. Random seed, for reproducibility.
#' @return A list with `draws` (parameters by draws) and `which` (which block
#'   each row belongs to), or `NULL` if unavailable.
#' @keywords internal
#' @noRd
ilm_joint_draws <- function(object, nsim, seed) {
  Q <- object$jointPrecision
  if (is.null(Q)) return(NULL)
  R <- tryCatch(Matrix::chol(Q), error = function(e) NULL)
  if (is.null(R)) return(NULL)
  ## The joint precision is in the template's own parameter order, random
  ## blocks where they are declared -- not the fixed parameters and then the
  ## random ones. Centred as c(opt$par, par.random), every fixed parameter
  ## declared after the random block (a dispersion, cut points, a zero part,
  ## an AR term) moved the draws off their estimates: for a gaussian y ~ s(x),
  ## a fit of 0.245 had a "95% interval" of -3.30 to -2.98. Each block is
  ## centred by name instead; under REML beta is one of the random blocks.
  rn <- rownames(Q)
  pf <- object$opt$par; pr <- object$sdr$par.random
  mu <- rep(NA_real_, length(rn))
  for (nm in unique(rn)) {
    src <- if (nm %in% names(pr)) pr[names(pr) == nm] else pf[names(pf) == nm]
    at <- which(rn == nm)
    if (length(src) != length(at)) return(NULL)
    mu[at] <- src
  }
  set.seed(seed)
  Z <- matrix(rnorm(nrow(Q) * nsim), nrow(Q), nsim)
  D <- as.matrix(Matrix::solve(R, Z))            # columns ~ N(0, Q^-1)
  list(draws = sweep(D, 1L, mu, "+"), which = rn)
}

#' Predictions from a fitted model
#'
#' Returns the fitted mean on the response scale -- one probability per
#' category for a multinomial or ordinal outcome, as
#' `nnet::multinom(type = "probs")` does -- or the linear predictor, optionally
#' with standard errors and intervals.
#'
#' @section Conditional versus population-averaged:
#' This is the choice that matters most, and there is no safe default that suits
#' everyone.
#'
#' With `marginal = FALSE` the random effects are set to zero, giving the
#' probabilities for a **typical** group -- one exactly at the population
#' average. With `marginal = TRUE` the prediction is averaged over the
#' distribution of random effects, giving the probabilities for the
#' **population as a whole**.
#'
#' These differ, sometimes substantially, because averaging and the softmax
#' transform do not commute: the average of the transformed values is not the
#' transform of the average. The population-averaged probabilities are pulled
#' toward being more even across categories. Which you want depends on the
#' question -- "what do I expect for an average subject?" or "what proportion of
#' the population falls in each category?"
#'
#' Under an **identity link** the two coincide exactly, because the random
#' effects have mean zero and nothing nonlinear stands between. `marginal` is
#' then answered in closed form rather than by simulation, so the result does
#' not depend on `ndraw` and carries no Monte Carlo noise.
#'
#' A **random slope** is averaged over as a slope. The amount being integrated
#' over then depends on the row -- it grows with distance from wherever the
#' slope is centred -- so the marginal and conditional curves separate by more
#' at the ends of the range than in the middle. Averaging such a term as if it
#' were an intercept understates that, and the error grows with the slope
#' variance and with distance from centre.
#'
#' An **AR(1) or CAR(1) term** is averaged over too: at any one row its latent
#' value has the stationary distribution, whatever the time.
#'
#' With **one linear predictor** -- every family but the multinomial -- a row's
#' whole latent contribution is a single normal, and the average is taken by
#' Gauss-Hermite quadrature: exact for any purpose, the same on every call, and
#' free of `ndraw`. A multinomial outcome has one dimension per category, and
#' is averaged over `ndraw` draws with common random numbers.
#'
#' @section Uncertainty:
#' Standard errors and intervals come from simulation rather than a formula,
#' because the softmax makes the quantity nonlinear in the parameters. Intervals
#' are **percentile** intervals from the simulated draws, so they always lie
#' within 0 and 1; a symmetric interval on the probability scale would not.
#'
#' If the fit was made with `joint = TRUE` the draws include the penalised smooth
#' coefficients. Without it only the fixed effects vary, which breaks the
#' correlation described in `ilm_joint_draws()` and distorts intervals around
#' smooths; a warning says so. [ilm_model()] enables it automatically when the model
#' contains smooths.
#'
#' Random-effect draws are held fixed across rows and across parameter draws
#' ("common random numbers"). Without that, Monte Carlo noise would swamp
#' comparisons between grid points, and `marginaleffects` would be unable to
#' compute stable numerical derivatives.
#'
#' @param object A fitted `"ilm_model"` object.
#' @param newdata Optional data frame. If omitted, predictions are for the data
#'   used to fit the model.
#' @param type `"response"` for probabilities (the default), `"link"` for linear
#'   predictors, or `"class"` for the most likely category.
#' @param marginal Logical. Average over the random-effect distribution
#'   (population-averaged) rather than setting it to zero (conditional).
#' @param se.fit Logical. Return standard errors.
#' @param interval `"none"` or `"confidence"`.
#' @param level Numeric. Interval coverage, default 0.95.
#' @param nsim Integer. Parameter draws used for uncertainty.
#' @param ndraw Integer. Random-effect draws used when `marginal = TRUE` for a
#'   multinomial outcome; a single linear predictor is averaged by quadrature
#'   and does not use it.
#' @param seed Integer. Random seed, so results are reproducible.
#' @param ... Unused.
#'
#' @return A matrix of probabilities (or linear predictors), or a factor for
#'   `type = "class"`. When standard errors or intervals are requested, a list
#'   with `fit`, `se.fit`, `lower`, `upper`, `level` and `joint`.
#'
#' @references
#' Skrondal, A., & Rabe-Hesketh, S. (2009). Prediction in multilevel generalized
#' linear models. *Journal of the Royal Statistical Society, Series A*, 172(3),
#' 659--687. (On the distinction between conditional and marginal prediction.)
#'
#' @seealso [ilm_model()], [ilm_fitted()].
#' @export
predict.ilm_model <- function(object, newdata = NULL,
                         type = c("response", "link", "class"),
                         marginal = FALSE, se.fit = FALSE,
                         interval = c("none", "confidence"), level = 0.95,
                         nsim = 200L, ndraw = 200L, seed = 1L, ...) {
  type <- match.arg(type); interval <- match.arg(interval)
  ## A flexible parametric model's linear predictor depends on TIME through the
  ## spline, so there is no fitted value for a covariate pattern alone. Asking
  ## for one is a question about the survival curve.
  if (!is.null(object$rp) && !is.null(newdata))
    stop("a flexible parametric model has no fitted value for a covariate ",
         "pattern on its own: its linear predictor depends on time through ",
         "the baseline spline. Use ilm_survival(object, newdata, times) for ",
         "the survival curve, or ilm_plot_survival() to see it.",
         call. = FALSE)
  want_unc <- isTRUE(se.fit) || interval != "none"
  if (marginal && type == "link")
    stop("marginal = TRUE applies on the response scale; use type = \"response\"")
  ## needed before point() closes over it
  multinom0 <- object$C > 1L
  Tc <- contr.sum(object$J)
  nd <- if (is.null(newdata)) list(X = object$X, smooths = NULL) else ilm_newX(object, newdata)
  if (is.null(newdata) && length(object$smooths))
    nd$smooths <- lapply(object$smooths, ilm_smooth_design, newdata = object$model)

  gk <- which(vapply(object$re, function(e) e$kind != "basis", TRUE))
  ## An AR or CAR term is a population to average over too. Its latent value at
  ## any one row is N(0, Sigma$ar) -- the stationary covariance, whatever the
  ## time -- and it used to be left out: a Poisson model with a stationary AR
  ## variance of 0.59 averaged 1.40 where its own simulations averaged 1.86.
  has_ar <- !is.null(object$ar) && !is.null(object$Sigma[["ar"]])
  ## Under an identity link the average over the random effects IS the
  ## conditional value: E[eta + z'u] = eta, because the random effects have
  ## mean zero and nothing nonlinear stands between. Simulating it instead
  ## returns a noisy estimate of a number already known exactly -- with 200
  ## draws and a random slope that noise reached 0.14 on the response scale.
  ## A zero part scales the mean by a constant, so it stays exact too.
  if (marginal && !isTRUE(object$ordinal) && object$C == 1L &&
      !is.null(object$family) && identical(object$family$link, "identity"))
    marginal <- FALSE
  integ <- marginal && (length(gk) > 0L || has_ar)
  draws <- NULL; vrow <- NULL
  if (integ) {
    pdat <- if (is.null(newdata)) object$model else newdata
    tms <- lapply(gk, function(k) {
      f  <- ilm_re_factors(object, k)
      Zb <- ilm_re_design(object, k, pdat, f$d)
      if (is.null(Zb)) {
        ## the bar varies over something the prediction rows do not carry, so
        ## the only honest option left is the intercept part -- said out loud,
        ## because a marginal average over less than the whole term is a
        ## different quantity from the one that was asked for
        warning("the random-effect term '", names(object$re)[k], "' varies ",
                "within a group over a column the prediction data does not ",
                "have, so the marginal average integrates its intercept ",
                "only. Supply that column to average over the whole term.",
                call. = FALSE)
        f$A <- matrix(1, 1L, 1L)
        Zb  <- matrix(1, nrow(pdat), 1L)
      }
      list(Zb = Zb, A = f$A, B = f$B)
    })
    if (object$C > 1L) {
      ## C linear predictors: draws, with common random numbers, in the same
      ## stream as before for the grouping terms, and the AR term after them
      set.seed(seed)
      draws <- lapply(tms, function(tk)
        list(Zb = tk$Zb, U = ilm_re_draws(tk$A, tk$B, ndraw, object$C)))
      if (has_ar)
        draws[[length(draws) + 1L]] <- list(
          Zb = matrix(1, nrow(pdat), 1L),
          U = ilm_re_draws(matrix(1, 1L, 1L), ilm_msqrt(object$Sigma[["ar"]]),
                           ndraw, object$C))
    } else {
      ## One linear predictor: a row's latent contribution is a single normal,
      ## with variance z' Sigma_d z times the term's scale for each term, plus
      ## the AR variance. Averaged by quadrature rather than by draws, whose
      ## noise grows with that variance and moves finite-difference effects
      ## with the seed.
      vrow <- numeric(nrow(pdat))
      for (tk in tms) vrow <- vrow + rowSums((tk$Zb %*% tk$A)^2) * sum(tk$B^2)
      if (has_ar) vrow <- vrow + object$Sigma[["ar"]][1L, 1L]
    }
  }
  ## A univariate family has one linear predictor and its own inverse link; the
  ## multinomial has C dimensions that the softmax maps onto J probabilities.
  multinom <- object$C > 1L
  ## An ordered response has one linear predictor and J categories, so it is
  ## neither of the two shapes the rest of this function knows about.
  ordinal <- isTRUE(object$ordinal)
  linkinv <- if (!is.null(object$family)) object$family$linkinv else identity
  ## A zero part makes the response scale something other than the inverse
  ## link of the linear predictor: that is the mean of the COUNT process, and
  ## the mean of the RESPONSE has to carry the zeros too. The design is built
  ## for the prediction rows, since the probability depends on their own
  ## covariates rather than on the ones the model was fitted to.
  Zp <- if (!is.null(object$Zzi))
    ilm_zi_design(object$zi_formula,
                  if (is.null(newdata)) object$model else newdata,
                  colnames(object$Zzi)) else NULL
  zi_adj <- function(mu) if (is.null(Zp)) mu else ilm_zi_mean(object, mu, Zp)
  ## The shift for draw m, as an n x C matrix rather than one C-vector added to
  ## every row. With a random slope the amount being integrated over depends on
  ## the row -- it grows with distance from where the slope is centred -- and a
  ## constant shift cannot represent that. With a random intercept alone `Zb`
  ## is a column of ones and this is the constant shift it always was.
  shift <- function(m, n, C) {
    S <- matrix(0, n, C)
    for (q in seq_along(draws)) S <- S + draws[[q]]$Zb %*% draws[[q]]$U[[m]]
    S
  }
  point <- function(beta, bvec = NULL) {
    eta <- ilm_eta(object, nd, beta, bvec)
    if (type == "link") return(if (multinom) eta %*% t(Tc) else eta[, 1, drop = FALSE])
    if (ordinal) {
      if (!integ)
        return(ilm_ord_probs(eta[, 1], object$zeta, object$family$pfun))
      gh <- ilm_gh(); P <- 0
      for (q in seq_along(gh$z))
        P <- P + gh$w[q] * ilm_ord_probs(eta[, 1] + sqrt(vrow) * gh$z[q],
                                         object$zeta, object$family$pfun)
      return(P)
    }
    if (!integ)
      return(if (multinom) ilm_softmax_J(eta, Tc) else
               matrix(zi_adj(linkinv(eta[, 1])), ncol = 1L))
    if (!multinom) {
      ## the zero part inside the integral: a hurdle's mean is not linear in
      ## the count mean, so it cannot be applied to the average afterwards
      gh <- ilm_gh(); P <- 0
      for (q in seq_along(gh$z))
        P <- P + gh$w[q] * zi_adj(linkinv(eta[, 1] + sqrt(vrow) * gh$z[q]))
      return(matrix(P, ncol = 1L))
    }
    P <- matrix(0, nrow(eta), object$J)
    for (m in seq_len(ndraw)) {
      sh <- shift(m, nrow(eta), ncol(eta))
      P <- P + ilm_softmax_J(eta + sh, Tc)
    }
    P / ndraw
  }
  est <- point(object$beta)
  colnames(est) <- if (multinom || (ordinal && type != "link")) object$ylevels
    else if (type == "link") "link" else "response"
  if (type == "class") {
    if (!multinom && !ordinal)
      stop("type = \"class\" applies only to a categorical response ",
           "(multinomial or ordinal)", call. = FALSE)
    cls <- factor(object$ylevels[max.col(est, ties.method = "first")],
                  levels = object$ylevels, ordered = ordinal)
    if (!want_unc) return(cls)
    warning("se.fit / interval are not defined for type = \"class\"", call. = FALSE)
    return(cls)
  }
  if (!want_unc) return(est)

  ## ---- uncertainty by simulation -----------------------------------------
  p <- ncol(object$X); C <- object$C
  jd <- if (!is.null(object$jointPrecision)) ilm_joint_draws(object, nsim, seed + 1L) else NULL
  if (is.null(jd)) {
    ## NOT "too narrow" -- measured, it is the opposite.  Drawing beta from its
    ## marginal covariance while HOLDING the penalised coefficients fixed breaks
    ## the strong negative correlation between a smooth's null-space term and its
    ## penalised part (they trade off against each other).  Ignoring that
    ## compensation inflates the fitted smooth's variability: measured mean
    ## interval width 0.358 without joint draws versus 0.123 with them, and the
    ## non-joint intervals ran outside [0, 1].
    if (length(object$smooths))
      warning("fit was made without joint = TRUE: the penalised smooth ",
              "coefficients are held fixed, which breaks their correlation with ",
              "the null-space term and makes smooth intervals far TOO WIDE. ",
              "Refit with joint = TRUE.", call. = FALSE)
    V <- suppressWarnings(vcov(object))
    R <- tryCatch(chol(V + diag(1e-12, ncol(V))), error = function(e) NULL)
    if (is.null(R)) { warning("covariance unusable; no intervals"); return(est) }
    set.seed(seed + 1L); bb <- as.vector(object$beta)
    acc <- array(0, c(nrow(est), ncol(est), nsim))
    for (s in seq_len(nsim))
      acc[, , s] <- point(matrix(bb + as.vector(rnorm(length(bb)) %*% R), p, C))
  } else {
    isb <- jd$which == "beta"; isr <- jd$which == "bvec"
    acc <- array(0, c(nrow(est), ncol(est), nsim))
    for (s in seq_len(nsim))
      acc[, , s] <- point(matrix(jd$draws[isb, s], p, C), jd$draws[isr, s])
  }
  se <- apply(acc, 1:2, sd)
  a <- (1 - level) / 2
  lo <- apply(acc, 1:2, quantile, probs = a,     names = FALSE)
  hi <- apply(acc, 1:2, quantile, probs = 1 - a, names = FALSE)
  dimnames(se) <- dimnames(lo) <- dimnames(hi) <- dimnames(est)
  list(fit = est, se.fit = se, lower = lo, upper = hi, level = level,
       joint = !is.null(jd))
}
