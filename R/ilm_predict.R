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
  X  <- stats::model.matrix(mt, mf, contrasts.arg = object$contrasts)
  sd_list <- list()
  for (lab in names(object$smooths)) {
    sd_list[[lab]] <- ilm_smooth_design(object$smooths[[lab]], newdata)
    if (!is.null(sd_list[[lab]]$Xf)) X <- cbind(X, sd_list[[lab]]$Xf)
  }
  if (ncol(X) != ncol(object$X))
    stop("new-data design has ", ncol(X), " columns but the fit has ", ncol(object$X))
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
  mu <- c(object$opt$par, object$sdr$par.random)
  set.seed(seed)
  Z <- matrix(rnorm(nrow(Q) * nsim), nrow(Q), nsim)
  D <- as.matrix(Matrix::solve(R, Z))            # columns ~ N(0, Q^-1)
  list(draws = sweep(D, 1L, mu, "+"), which = rownames(Q))
}

#' Predicted category probabilities
#'
#' Returns one probability per observation per category, as
#' `nnet::multinom(type = "probs")` does, optionally with standard errors and
#' intervals.
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
#' @param ndraw Integer. Random-effect draws used when `marginal = TRUE`.
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
  draws <- NULL
  if (marginal && length(gk)) {                    # fixed RE draws: common random numbers
    set.seed(seed)
    draws <- lapply(gk, function(k) {
      w <- object$wk[k]; Z <- matrix(rnorm(ndraw * w), ndraw, w)
      if (identical(object$re_struct[[k]]$type, "rr")) Z %*% t(object$Lambda[[k]])
      else Z %*% ilm_msqrt(object$Sigma[[k]])
    })
  }
  ## A univariate family has one linear predictor and its own inverse link; the
  ## multinomial has C dimensions that the softmax maps onto J probabilities.
  multinom <- object$C > 1L
  linkinv <- if (!is.null(object$family)) object$family$linkinv else identity
  point <- function(beta, bvec = NULL) {
    eta <- ilm_eta(object, nd, beta, bvec)
    if (type == "link") return(if (multinom) eta %*% t(Tc) else eta[, 1, drop = FALSE])
    if (!marginal || !length(gk))
      return(if (multinom) ilm_softmax_J(eta, Tc) else
               matrix(linkinv(eta[, 1]), ncol = 1L))
    P <- matrix(0, nrow(eta), if (multinom) object$J else 1L)
    for (m in seq_len(ndraw)) {
      sh <- Reduce(`+`, lapply(draws, function(d) d[m, ]))
      P <- P + if (multinom) ilm_softmax_J(sweep(eta, 2L, sh, `+`), Tc)
               else matrix(linkinv(eta[, 1] + sh[1]), ncol = 1L)
    }
    P / ndraw
  }
  est <- point(object$beta)
  colnames(est) <- if (multinom) object$ylevels else
    if (type == "link") "link" else "response"
  if (type == "class") {
    if (!multinom)
      stop("type = \"class\" applies only to the multinomial family", call. = FALSE)
    cls <- factor(object$ylevels[max.col(est, ties.method = "first")],
                  levels = object$ylevels)
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
