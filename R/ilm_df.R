## ---------------------------------------------------------------------------
## Denominator degrees of freedom for a mixed model.
##
## A Wald statistic treats the variance components as known. They are not, and
## with few clusters the resulting chi-square reference is anti-conservative --
## intervals too narrow, p-values too small. Satterthwaite and Kenward-Roger
## both answer the same question: how much is a contrast's standard error
## itself worth trusting?
##
## SATTERTHWAITE matches a scaled chi-square to the estimated variance of the
## contrast. For a contrast l,
##
##     nu  =  2 (l' Vb l)^2  /  ( g' A g ),    g = d(l' Vb l)/d theta
##
## where Vb is the fixed-effect covariance with theta HELD FIXED, and A is the
## sampling covariance of theta_hat. Both pieces are already available here:
## A is the theta block of sdreport's cov.fixed, and Vb(theta) is the inverse
## of the beta block of the Hessian of the Laplace objective, which can be
## recovered by differencing obj$gr(). That route is worth spelling out,
## because it is why this works for EVERY structure the package fits rather
## than only the ones someone wrote a formula for: it asks the objective, and
## the objective already knows about AR(1), smooths, dispersion models and the
## rest. Checked against vcov() it reproduces it to 8e-11.
##
## KENWARD-ROGER goes further and inflates Vb itself for the uncertainty in
## theta, then computes a df for the inflated version. It needs
##
##     Q_ij = X' V^-1 V_i V^-1 V_j V^-1 X     and     R_ij = X' V^-1 V_ij V^-1 X
##
## SEPARATELY, where V is the MARGINAL covariance and V_i its derivative.
## Differentiating the beta block of the Hessian twice gives only the
## combination Q_ij + Q_ji - R_ij, so the objective route does not reach it --
## except when V is LINEAR in the variance parameters, which makes V_ij zero,
## R_ij vanish, and the second derivative exactly 2*sym(Q_ij). That is the case
## for random intercepts and slopes plus a residual, which is where Kenward and
## Roger derived it and where repeated-measures designs live. It is NOT the
## case once a correlation parameter enters -- AR(1), CAR(1), or an
## unstructured term with off-diagonal covariances -- and there ilm_denom_df()
## refuses rather than returning a number it cannot stand behind.
##
## Both are defined for LINEAR mixed models. For a non-gaussian family the
## small-sample reference is a different problem and neither applies; the
## remedy there is ilm_pb_lrt(), which simulates the null instead of
## approximating it.
## ---------------------------------------------------------------------------

#' Stop a likelihood comparison that REML makes meaningless
#'
#' A restricted likelihood is the likelihood of contrasts orthogonal to `X`.
#' Change `X` and it is a likelihood of different data, so two REML fits with
#' different fixed effects have nothing comparable about them -- their
#' difference is not a likelihood ratio and its reference is not chi-square.
#' The same applies to AIC and to any bootstrap built on the ratio.
#'
#' This is not a numerical nicety. The comparison runs perfectly happily and
#' returns a number, which is why it has to be refused rather than warned
#' about.
#'
#' @param object A fitted `"ilm_model"` object.
#' @param what What the caller was trying to do, named in the message.
#' @return `TRUE`, invisibly, or an error.
#' @keywords internal
#' @noRd
ilm_stop_reml_lrt <- function(object, what = "a likelihood-ratio test") {
  if (!isTRUE(object$reml)) return(invisible(TRUE))
  stop(what, " compares likelihoods across different fixed-effect structures, ",
       "and this model was fitted by REML. A restricted likelihood belongs to ",
       "contrasts orthogonal to the design matrix, so changing the fixed ",
       "effects changes which data it is the likelihood of and the two are ",
       "not comparable. Refit with `reml = FALSE` to compare fixed effects, ",
       "then switch back once the structure is settled -- or use a Wald test, ",
       "which is valid under REML.", call. = FALSE)
}

#' Which parameters are fixed effects, and which are variance components
#'
#' The optimiser works on one vector holding both. Everything here needs to
#' tell them apart, and `names(opt$par)` is how the fit records it.
#'
#' @param object A fitted `"ilm_model"` object.
#' @return A list with integer vectors `beta` and `theta`.
#' @keywords internal
#' @noRd
ilm_par_blocks <- function(object) {
  p <- object$opt$par
  nm <- names(p)
  if (is.null(nm))
    stop("the fit does not record parameter names, so the fixed effects ",
         "cannot be separated from the variance components", call. = FALSE)
  list(beta = which(nm == "beta"), theta = which(nm != "beta"))
}

#' The fixed-effect covariance with the variance components held fixed
#'
#' `vcov()` reports the block of the inverse joint Hessian, which carries the
#' uncertainty in `theta` with it. Satterthwaite needs the other thing: the
#' covariance you would have if `theta` were known, as a FUNCTION of `theta`,
#' so that it can be differentiated. That is the inverse of the beta block of
#' the Hessian, and the Hessian is recovered by central differences on the
#' gradient -- `obj$he()` is not available once anything is integrated out.
#'
#' @param object A fitted `"ilm_model"` object.
#' @param par Parameter vector to evaluate at; defaults to the fit's own.
#' @param h Step size for the central difference.
#' @return The beta-by-beta covariance matrix.
#' @keywords internal
#' @noRd
ilm_vbeta_at <- function(object, par = NULL, h = 1e-5) {
  ix <- ilm_par_blocks(object)
  p <- if (is.null(par)) object$opt$par else par
  ## Under REML `beta` was integrated out, so object$obj no longer takes it as
  ## an argument and cannot be differentiated over it. The fit keeps an
  ## ML-SHAPED twin for exactly this: V_beta(theta) is a function of theta and
  ## the data, not of how theta was estimated, so evaluating the twin at the
  ## REML estimates gives the right derivatives. Verified against
  ## vcov(lmer(REML = TRUE)) to 2.4e-08.
  obj <- if (!is.null(object$obj_ml)) object$obj_ml else object$obj
  nb <- length(ix$beta)
  H <- matrix(0, nb, nb)
  for (k in seq_len(nb)) {
    up <- p; up[ix$beta[k]] <- up[ix$beta[k]] + h
    dn <- p; dn[ix$beta[k]] <- dn[ix$beta[k]] - h
    H[, k] <- (obj$gr(up)[ix$beta] - obj$gr(dn)[ix$beta]) / (2 * h)
  }
  H <- (H + t(H)) / 2                       # symmetrise the differencing noise
  V <- try(solve(H), silent = TRUE)
  if (inherits(V, "try-error")) return(NULL)
  V
}

#' Is the marginal covariance linear in the variance parameters?
#'
#' Kenward-Roger needs this and Satterthwaite does not. A random intercept or
#' slope contributes a variance that enters `V` linearly once it is put on the
#' variance scale; a correlation does not. `V` also stops being linear when the
#' residual variance is itself modelled, because then it varies by row in a way
#' the variance components do not control.
#'
#' @param object A fitted `"ilm_model"` object.
#' @return `TRUE` when Kenward-Roger is applicable, otherwise a character
#'   string saying what makes it inapplicable.
#' @keywords internal
#' @noRd
ilm_kr_applicable <- function(object) {
  if (!is.null(object$ar))
    return("the model has an AR(1) or CAR(1) term, whose correlation parameter makes the marginal covariance non-linear in the variance components")
  if (!is.null(object$disp_formula))
    return("the model has a dispersion formula, so the residual variance is not a single variance component")
  st <- object$re_struct
  if (length(st)) {
    for (k in seq_along(st)) {
      ty <- st[[k]]$type
      if (identical(ty, "us") && isTRUE(object$C > 1L))
        return("an unstructured category covariance carries correlations, which enter the marginal covariance non-linearly")
      if (identical(ty, "rr"))
        return("a reduced-rank term is a product of loadings, not a set of variances")
    }
  }
  if (length(object$re) && any(object$wk > 1L))
    return("a random slope term with a free intercept-slope correlation is not linear in the variance components")
  TRUE
}

#' Denominator degrees of freedom for one or more contrasts
#'
#' Turns a Wald statistic into an F or t test by working out how much the
#' estimated standard error is itself worth trusting. Without this a mixed
#' model reports a chi-square, which treats the variance components as known
#' and is anti-conservative when there are few clusters.
#'
#' @section Which method:
#' `"satterthwaite"` matches a scaled chi-square to the estimated variance of
#' the contrast. It is the default because it costs first derivatives only, it
#' applies to every structure this package fits, and it is what `lmerTest`
#' reports.
#'
#' `"kenward-roger"` also inflates the covariance for the uncertainty in the
#' variance components, which matters most when clusters are few. It costs
#' second derivatives, and it is only defined where the marginal covariance is
#' linear in the variance components -- random intercepts and slopes with a
#' single residual variance. Asked for elsewhere it stops and says why.
#'
#' `"residual"` is `N - p`, which is exact with nothing integrated out and
#' optimistic otherwise. `"asymptotic"` returns `Inf`, recovering the
#' chi-square.
#'
#' @section What it is not for:
#' Both approximations are derived for LINEAR mixed models. For a non-gaussian
#' family the small-sample problem is a different one and neither answers it;
#' [ilm_pb_lrt()] simulates the null rather than approximating its reference,
#' and is the remedy there.
#'
#' @param object A fitted `"ilm_model"` object.
#' @param L A contrast vector of length `p`, or a `q`-by-`p` contrast matrix.
#' @param method One of `"auto"`, `"satterthwaite"`, `"kenward-roger"`,
#'   `"residual"`, `"asymptotic"`. `"auto"` gives the exact residual df when
#'   nothing was integrated out, Satterthwaite for a gaussian mixed model, and
#'   `Inf` otherwise.
#' @param h Step size for the numerical derivatives.
#' @return A list with `df`, the `method` actually used, and for
#'   Kenward-Roger the adjusted covariance `V`.
#' @seealso [ilm_anova()], [ilm_emmeans()], [ilm_pb_lrt()].
#' @examples
#' set.seed(1)
#' d <- ilm_sim(n_id = 25)
#' f <- ilm_model(score ~ income + (1 | id), data = d, family = "gaussian",
#'                verbose = FALSE)
#' L <- c(0, 1)                      # the slope on income
#' ilm_denom_df(f, L)
#' @export
ilm_denom_df <- function(object, L,
                         method = c("auto", "satterthwaite", "kenward-roger",
                                    "residual", "asymptotic"),
                         h = 1e-5) {
  if (!inherits(object, "ilm_model"))
    stop("`object` must be a fitted ilm_model object, not ", class(object)[1],
         call. = FALSE)
  method <- match.arg(method)
  L <- if (is.matrix(L)) L else matrix(L, nrow = 1L)
  p <- length(object$beta)
  if (ncol(L) != p)
    stop("`L` has ", ncol(L), " column(s) but the model has ", p,
         " fixed-effect coefficient(s)", call. = FALSE)

  if (method == "auto") {
    method <- if (isTRUE(object$exact_df)) "residual"
              else if (identical(object$family$name, "gaussian")) "satterthwaite"
              else "asymptotic"
  }
  if (method == "residual")
    return(list(df = if (isTRUE(object$exact_df)) object$resid_df
                     else nrow(object$X) - p,
                method = "residual", V = NULL))
  if (method == "asymptotic")
    return(list(df = Inf, method = "asymptotic", V = NULL))

  ## Neither approximation is derived for a non-gaussian family. Saying so is
  ## better than returning a number that looks like an answer.
  if (!identical(object$family$name, "gaussian"))
    stop(method, " degrees of freedom are derived for LINEAR mixed models, and ",
         "this is a ", object$family$name, " fit. Use `method = \"asymptotic\"` ",
         "for the large-sample reference, or ilm_pb_lrt() to simulate the null ",
         "instead of approximating it.", call. = FALSE)
  if (isTRUE(object$exact_df))
    return(list(df = object$resid_df, method = "residual", V = NULL))

  if (method == "kenward-roger") {
    ok <- ilm_kr_applicable(object)
    if (!isTRUE(ok))
      stop("Kenward-Roger is not defined for this model: ", ok,
           ". Use `method = \"satterthwaite\"`, which does apply here.",
           call. = FALSE)
    return(ilm_df_kr(object, L, h = h))
  }
  ilm_df_satt(object, L, h = h)
}

#' The gradient of a contrast's variance with respect to the variance components
#'
#' Central differences on `theta`, each step costing one beta-block Hessian.
#' `theta` is short -- a handful of entries even for an elaborate model -- so
#' this stays cheap where differencing over `beta` would not.
#'
#' @param object A fitted `"ilm_model"` object.
#' @param L Contrast matrix.
#' @param h Step size.
#' @return A list with `V0` (the covariance at the estimate) and `G`, a list of
#'   derivative matrices, one per variance component.
#' @keywords internal
#' @noRd
ilm_vbeta_grad <- function(object, L, h = 1e-5) {
  ix <- ilm_par_blocks(object)
  p0 <- object$opt$par
  V0 <- ilm_vbeta_at(object, p0, h = h)
  if (is.null(V0)) return(NULL)
  G <- vector("list", length(ix$theta))
  for (j in seq_along(ix$theta)) {
    up <- p0; up[ix$theta[j]] <- up[ix$theta[j]] + h
    dn <- p0; dn[ix$theta[j]] <- dn[ix$theta[j]] - h
    Vu <- ilm_vbeta_at(object, up, h = h)
    Vd <- ilm_vbeta_at(object, dn, h = h)
    if (is.null(Vu) || is.null(Vd)) return(NULL)
    G[[j]] <- (Vu - Vd) / (2 * h)
  }
  list(V0 = V0, G = G, theta_idx = ix$theta)
}

#' The sampling covariance of the variance components
#'
#' sdreport already estimates it: `cov.fixed` covers the whole optimised
#' vector, and the block outside `beta` is exactly what is wanted.
#'
#' @param object A fitted `"ilm_model"` object.
#' @return A covariance matrix for `theta`, or `NULL` when it is unavailable.
#' @keywords internal
#' @noRd
ilm_theta_vcov <- function(object) {
  cf <- object$sdr$cov.fixed
  if (is.null(cf)) return(NULL)
  ix <- ilm_par_blocks(object)
  if (!length(ix$theta)) return(NULL)
  A <- cf[ix$theta, ix$theta, drop = FALSE]
  ## a term held at its boundary is treated as known, which is what holding
  ## it means: its variance parameters contribute no uncertainty to the df
  if (length(object$hessian_held)) A[!is.finite(A)] <- 0
  if (!all(is.finite(A))) return(NULL)
  A
}

#' Satterthwaite degrees of freedom
#'
#' For one contrast this is the scaled chi-square match directly. For several
#' at once the contrast matrix is first rotated to the basis in which
#' `L Vb L'` is diagonal, so that the rows become independent one-degree
#' contrasts, each gets its own `nu`, and the F denominator is assembled from
#' them by the usual expectation argument. That is Fai and Cornelius' method
#' and it is what `lmerTest` reports, so the two are directly comparable.
#'
#' @param object A fitted `"ilm_model"` object.
#' @param L Contrast matrix.
#' @param h Step size.
#' @return A list with `df` and `method`.
#' @keywords internal
#' @noRd
ilm_df_satt <- function(object, L, h = 1e-5) {
  gr <- ilm_vbeta_grad(object, L, h = h)
  A  <- ilm_theta_vcov(object)
  if (is.null(gr) || is.null(A))
    return(list(df = Inf, method = "asymptotic", V = NULL))
  V0 <- gr$V0; G <- gr$G
  q <- nrow(L)

  ## one row: no rotation needed
  nu_of <- function(l) {
    d <- as.numeric(l %*% V0 %*% l)
    if (!is.finite(d) || d <= 0) return(NA_real_)
    g <- vapply(G, function(Gj) as.numeric(l %*% Gj %*% l), numeric(1))
    den <- as.numeric(t(g) %*% A %*% g)
    if (!is.finite(den) || den <= 0) return(NA_real_)
    2 * d^2 / den
  }
  if (q == 1L) {
    nu <- nu_of(as.numeric(L))
    return(list(df = if (is.na(nu)) Inf else nu, method = "satterthwaite",
                V = NULL))
  }

  ## several rows: rotate to the basis where L Vb L' is diagonal
  LVL <- L %*% V0 %*% t(L)
  ev <- eigen((LVL + t(LVL)) / 2, symmetric = TRUE)
  keep <- ev$values > max(ev$values) * 1e-10
  if (!any(keep)) return(list(df = Inf, method = "asymptotic", V = NULL))
  P <- t(ev$vectors[, keep, drop = FALSE])       # rows are the new contrasts
  Lr <- P %*% L
  nus <- vapply(seq_len(nrow(Lr)), function(m) nu_of(Lr[m, ]), numeric(1))
  nus <- nus[is.finite(nus) & nus > 2]
  if (!length(nus)) return(list(df = Inf, method = "asymptotic", V = NULL))
  E <- sum(nus / (nus - 2))
  qq <- length(nus)
  df <- if (E > qq) 2 * E / (E - qq) else Inf
  list(df = df, method = "satterthwaite", V = NULL)
}

#' Kenward-Roger degrees of freedom and adjusted covariance
#'
#' Only reached when [ilm_kr_applicable()] has agreed that the marginal
#' covariance is linear in the variance components, which makes the second
#' derivative of the beta block of the Hessian exactly twice the symmetric part
#' of `Q_ij` and lets the adjustment be assembled without ever forming `V`.
#'
#' The adjustment has two halves. The covariance is inflated, because treating
#' `theta_hat` as though it were `theta` understates the spread of the fixed
#' effects; and the degrees of freedom are then computed from the inflated
#' version, which is what makes the two halves consistent with one another.
#'
#' @param object A fitted `"ilm_model"` object.
#' @param L Contrast matrix.
#' @param h Step size.
#' @return A list with `df`, `method` and the adjusted covariance `V`.
#' @keywords internal
#' @noRd
ilm_df_kr <- function(object, L, h = 1e-5) {
  ix <- ilm_par_blocks(object)
  A  <- ilm_theta_vcov(object)
  p0 <- object$opt$par
  V0 <- ilm_vbeta_at(object, p0, h = h)
  if (is.null(V0) || is.null(A))
    return(list(df = Inf, method = "asymptotic", V = NULL))
  nt <- length(ix$theta)

  ## first derivatives of Vb
  Vp <- vector("list", nt); Vm <- vector("list", nt); G1 <- vector("list", nt)
  for (j in seq_len(nt)) {
    up <- p0; up[ix$theta[j]] <- up[ix$theta[j]] + h
    dn <- p0; dn[ix$theta[j]] <- dn[ix$theta[j]] - h
    Vp[[j]] <- ilm_vbeta_at(object, up, h = h)
    Vm[[j]] <- ilm_vbeta_at(object, dn, h = h)
    if (is.null(Vp[[j]]) || is.null(Vm[[j]]))
      return(list(df = Inf, method = "asymptotic", V = NULL))
    G1[[j]] <- (Vp[[j]] - Vm[[j]]) / (2 * h)
  }

  ## second derivatives, by differencing the first
  G2 <- matrix(list(), nt, nt)
  for (j in seq_len(nt)) for (k in j:nt) {
    if (j == k) {
      G2[[j, k]] <- (Vp[[j]] - 2 * V0 + Vm[[j]]) / (h^2)
    } else {
      pp <- p0; pp[ix$theta[j]] <- pp[ix$theta[j]] + h
                pp[ix$theta[k]] <- pp[ix$theta[k]] + h
      pm <- p0; pm[ix$theta[j]] <- pm[ix$theta[j]] + h
                pm[ix$theta[k]] <- pm[ix$theta[k]] - h
      mp <- p0; mp[ix$theta[j]] <- mp[ix$theta[j]] - h
                mp[ix$theta[k]] <- mp[ix$theta[k]] + h
      mm <- p0; mm[ix$theta[j]] <- mm[ix$theta[j]] - h
                mm[ix$theta[k]] <- mm[ix$theta[k]] - h
      Vs <- lapply(list(pp, pm, mp, mm), function(z) ilm_vbeta_at(object, z, h = h))
      if (any(vapply(Vs, is.null, TRUE)))
        return(list(df = Inf, method = "asymptotic", V = NULL))
      G2[[j, k]] <- (Vs[[1]] - Vs[[2]] - Vs[[3]] + Vs[[4]]) / (4 * h^2)
    }
    G2[[k, j]] <- G2[[j, k]]
  }

  ## Kenward and Roger (1997) eq. 2: the inflation, written with Vb and its
  ## derivatives rather than with V, which the linearity check has licensed.
  Vi <- solve(V0)
  adj <- matrix(0, nrow(V0), ncol(V0))
  for (j in seq_len(nt)) for (k in seq_len(nt)) {
    Pj <- -Vi %*% G1[[j]] %*% Vi
    Pk <- -Vi %*% G1[[k]] %*% Vi
    Qjk <- Vi %*% G2[[j, k]] %*% Vi / 2
    adj <- adj + A[j, k] * (Qjk - Pj %*% V0 %*% Pk)
  }
  Vkr <- V0 + 2 * V0 %*% adj %*% V0
  Vkr <- (Vkr + t(Vkr)) / 2

  ## df from the inflated covariance, by the same scaled chi-square match
  nu_of <- function(l, Vuse) {
    d <- as.numeric(l %*% Vuse %*% l)
    if (!is.finite(d) || d <= 0) return(NA_real_)
    g <- vapply(G1, function(Gj) as.numeric(l %*% Gj %*% l), numeric(1))
    den <- as.numeric(t(g) %*% A %*% g)
    if (!is.finite(den) || den <= 0) return(NA_real_)
    2 * d^2 / den
  }
  q <- nrow(L)
  if (q == 1L) {
    nu <- nu_of(as.numeric(L), Vkr)
    return(list(df = if (is.na(nu)) Inf else nu, method = "kenward-roger",
                V = Vkr))
  }
  LVL <- L %*% Vkr %*% t(L)
  ev <- eigen((LVL + t(LVL)) / 2, symmetric = TRUE)
  keep <- ev$values > max(ev$values) * 1e-10
  if (!any(keep)) return(list(df = Inf, method = "asymptotic", V = Vkr))
  Lr <- t(ev$vectors[, keep, drop = FALSE]) %*% L
  nus <- vapply(seq_len(nrow(Lr)), function(m) nu_of(Lr[m, ], Vkr), numeric(1))
  nus <- nus[is.finite(nus) & nus > 2]
  if (!length(nus)) return(list(df = Inf, method = "asymptotic", V = Vkr))
  E <- sum(nus / (nus - 2)); qq <- length(nus)
  list(df = if (E > qq) 2 * E / (E - qq) else Inf,
       method = "kenward-roger", V = Vkr)
}
