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
## KENWARD-ROGER goes further: it inflates Vb for the uncertainty in theta,
## and gives the F test its own df and a scale factor. It is computed as
## pbkrtest computes it, from the covariance of the observations written as
## Sigma = sum_r s_r G_r -- each grouping term's variances and covariances,
## then the residual variance -- which is linear in s for random intercepts
## and slopes, correlated or not, nested or crossed. See ilm_kr_parts().
##
## It used to be assembled from differences of Vb(theta) on the log-Cholesky
## scale, and the algebra was wrong: the Q term carried the wrong sign and the
## P Phi P term was missing, so for a single variance component the "inflated"
## covariance came out smaller than Phi, by a factor 1 - 4/n. The df were
## Satterthwaite's formula applied to that matrix, not Kenward and Roger's.
## Differences of Vb cannot give Q and R separately in any case, which is why
## this works from Sigma instead.
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

#' Is Kenward-Roger available for this fit?
#'
#' Kenward and Roger derived their adjustment for a LINEAR mixed model fitted
#' by REML, and it is computed here from the covariance of the observations,
#' written as a sum of variance parameters times known matrices. That is
#' linear in the parameters for random intercepts and slopes -- correlated or
#' not, nested or crossed -- with one residual variance. A correlation over
#' time, a modelled dispersion or a smooth's penalty is not in that form, and
#' the covariance of all N observations is formed at once, so it stops at a
#' size.
#'
#' @param object A fitted `"ilm_model"` object.
#' @return `TRUE` when Kenward-Roger is available, otherwise a character
#'   string saying why not.
#' @keywords internal
#' @noRd
ilm_kr_applicable <- function(object) {
  if (identical(object$ar$type, "rw1"))
    return("the model has a random-walk term, which the Kenward-Roger computation here does not cover")
  if (!is.null(object$ar))
    return("the model has an AR(1) or CAR(1) term, whose correlation parameter makes the marginal covariance non-linear in the variance components")
  if (!is.null(object$disp_formula))
    return("the model has a dispersion formula, so the residual variance is not a single variance component")
  if (!identical(object$family$link, "identity"))
    return("its link is not the identity, so it is not a linear mixed model")
  if (!is.null(object$Zzi))
    return("the model has a zero part, so it is not a linear mixed model")
  if (isTRUE(object$n_censored > 0L))
    return("the response is censored, so the model is not a linear mixed model")
  if (!is.null(object$weights))
    return("the model has frequency weights, which the covariance of the observations here does not carry")
  grp <- vapply(object$re, function(e) !identical(e$kind, "basis"), TRUE)
  if (!length(grp))
    return("the model has no random effects")
  if (!all(grp))
    return("the model has a penalised smooth, whose variance is a smoothing parameter rather than a variance component of the design")
  if (!isTRUE(object$reml))
    return("it is derived for REML estimates of the variance components, and this model was fitted by maximum likelihood; refit it with `reml = TRUE`")
  if (nrow(object$X) > ilm_kr_max_n)
    return(sprintf("it forms the covariance of all %d observations at once, and above %d rows that takes more memory than it should",
                   nrow(object$X), ilm_kr_max_n))
  TRUE
}

## Kenward-Roger forms the N-by-N covariance of the observations and its
## inverse, as pbkrtest does. At 4000 rows that is 128 MB for each.
ilm_kr_max_n <- 4000L

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
#' variance components, which matters most when clusters are few, and gives
#' an F test its own denominator df and scale factor (Kenward and Roger 1997).
#' It is computed as `pbkrtest` computes it, and agrees with it. It needs a
#' REML fit (`reml = TRUE`) of a gaussian model whose random terms are
#' intercepts and slopes -- correlated or not, nested or crossed -- with one
#' residual variance and no weights, and it forms the covariance of all the
#' observations, so it stops above 4000 rows. Asked for elsewhere it stops and
#' says why. For one contrast its scale is exactly 1: the t statistic is the
#' estimate over its standard error from the adjusted covariance `V`, on the
#' df returned.
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
#'   Kenward-Roger the adjusted covariance `V` and `scale`, the factor by which
#'   the Wald F formed with `V` is multiplied before it is referred to
#'   F(q, df).
#' @references
#' Kenward, M. G., & Roger, J. H. (1997). Small sample inference for fixed
#' effects from restricted maximum likelihood. *Biometrics*, 53(3), 983--997.
#'
#' Halekoh, U., & Hojsgaard, S. (2014). A Kenward-Roger approximation and
#' parametric bootstrap methods for tests in linear mixed models: the R
#' package pbkrtest. *Journal of Statistical Software*, 59(9), 1--30.
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
      stop("Kenward-Roger is not available for this model: ", ok,
           ". `method = \"satterthwaite\"` does apply here.",
           call. = FALSE)
    return(ilm_df_kr(object, L))
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
#' Only reached when [ilm_kr_applicable()] has agreed. The adjusted
#' covariance comes from [ilm_kr_parts()]; the F test's denominator df and
#' scale from [ilm_kr_ftest()], after `L` is reduced to independent rows --
#' the test is the same for any rows spanning the same hypothesis.
#'
#' @param object A fitted `"ilm_model"` object.
#' @param L Contrast matrix.
#' @param parts Optional output of [ilm_kr_parts()], to reuse.
#' @return A list with `df`, `method`, the adjusted covariance `V` and
#'   `scale`, the factor for the Wald F formed with `V`.
#' @keywords internal
#' @noRd
ilm_df_kr <- function(object, L, parts = NULL) {
  if (is.null(parts)) parts <- ilm_kr_parts(object)
  LVL <- L %*% parts$Phi %*% t(L)
  ev <- eigen((LVL + t(LVL)) / 2, symmetric = TRUE)
  keep <- ev$values > max(ev$values) * 1e-10
  if (!any(keep))
    return(list(df = Inf, method = "asymptotic", V = parts$PhiA, scale = 1))
  Lr <- t(ev$vectors[, keep, drop = FALSE]) %*% L
  ft <- ilm_kr_ftest(parts, Lr)
  if (!is.finite(ft$df) || ft$df <= 0) {
    warning("Kenward-Roger's denominator df could not be computed for this ",
            "contrast, so the large-sample reference is used.", call. = FALSE)
    return(list(df = Inf, method = "asymptotic", V = parts$PhiA, scale = 1))
  }
  list(df = ft$df, method = "kenward-roger", V = parts$PhiA, scale = ft$scale)
}

#' The pieces of Kenward and Roger's adjustment
#'
#' At the fitted variance parameters: the fixed-effect covariance `Phi`, its
#' adjusted version `PhiA`, and the `P_i` and `W` the degrees of freedom need.
#' Computed as `pbkrtest` computes them, from the covariance of the
#' observations written as
#'
#'     Sigma = sum_r s_r G_r ,
#'
#' with `s` each grouping term's variances and covariances and then the
#' residual variance, and each `G_r` a known matrix. `Sigma` is linear in `s`,
#' so the second-derivative terms `R_ij` are zero, and
#'
#'     P_i  = -X' S^-1 G_i S^-1 X
#'     Q_ij =  X' S^-1 G_i S^-1 G_j S^-1 X
#'     W    =  the inverse of the expected REML information for s
#'     PhiA =  Phi + 2 Phi { sum_ij W_ij (Q_ij - P_i Phi P_j) } Phi
#'
#' (Kenward and Roger 1997, section 2; Halekoh and Hojsgaard 2014, section 3).
#' Each `G_r` is a sum of outer products of sparse N-by-levels blocks, so
#' nothing N-by-N is formed but `Sigma` and its inverse.
#'
#' @param object A fitted `"ilm_model"` object that [ilm_kr_applicable()]
#'   accepts.
#' @return A list with `Phi`, `PhiA`, `P` (a list) and `W`.
#' @keywords internal
#' @noRd
ilm_kr_parts <- function(object) {
  X <- object$X; n <- nrow(X)
  s2 <- unname(object$dispersion[[1]])^2
  blk <- list()            # the N-by-levels blocks: one per term and dimension
  prs <- list()            # each parameter's G, as pairs (u, v) meaning u v'
  Sig <- diag(s2, n)
  for (k in seq_along(object$re)) {
    e <- object$re[[k]]; nm <- names(object$re)[k]; d <- e$d
    ## a term's covariance on the variance scale: the category scale times
    ## the within-group structure, whose first SD the fit holds at one
    S <- object$Sigma[[nm]][1, 1] *
      (if (d > 1L) object$Sigma_d[[nm]] else matrix(1, 1L, 1L))
    id <- integer(d)
    for (a in seq_len(d)) {
      blk[[length(blk) + 1L]] <- Matrix::sparseMatrix(
        i = seq_len(n), j = e$group, x = e$Z[, a], dims = c(n, e$nl))
      id[a] <- length(blk)
    }
    for (a in seq_len(d)) for (b in seq_len(d))
      if (S[a, b] != 0)
        Sig <- Sig + S[a, b] *
          as.matrix(Matrix::tcrossprod(blk[[id[a]]], blk[[id[b]]]))
    ## the parameters in pbkrtest's order, the lower triangle by columns, with
    ## a covariance only where the structure estimates one
    free <- d > 1L && isTRUE(object$dcor[k])
    for (b in seq_len(d)) for (a in b:d)
      if (a == b || free)
        prs[[length(prs) + 1L]] <- if (a == b) list(c(id[a], id[a]))
                                   else list(c(id[a], id[b]), c(id[b], id[a]))
  }
  prs[[length(prs) + 1L]] <- list()          # the residual variance: G = I
  R <- length(prs)

  SI <- chol2inv(chol(Sig))
  SX <- SI %*% X
  Phi <- solve(crossprod(X, SX)); Phi <- (Phi + t(Phi)) / 2
  ZX  <- lapply(blk, function(Z) as.matrix(Matrix::crossprod(Z, SX)))
  SIZ <- lapply(blk, function(Z) as.matrix(SI %*% Z))

  ## P_r, and G_r S^-1 X for the Q terms
  P <- vector("list", R); O <- vector("list", R)
  for (r in seq_len(R)) {
    if (!length(prs[[r]])) { P[[r]] <- -crossprod(SX); O[[r]] <- SX; next }
    Pr <- 0; Or <- 0
    for (uv in prs[[r]]) {
      Pr <- Pr - crossprod(ZX[[uv[1]]], ZX[[uv[2]]])
      Or <- Or + as.matrix(blk[[uv[1]]] %*% ZX[[uv[2]]])
    }
    P[[r]] <- (Pr + t(Pr)) / 2; O[[r]] <- Or
  }
  SIO <- lapply(O, function(o) SI %*% o)
  Q <- matrix(list(), R, R)
  for (r in seq_len(R)) for (s in seq_len(R)) Q[[r, s]] <- crossprod(O[[r]], SIO[[s]])

  ## the expected REML information, 1/2 tr(Pi G_r Pi G_s), through
  ## tr(S^-1 G_r S^-1 G_s) - 2 tr(Phi Q_rs) + tr(Phi P_r Phi P_s)
  Kt <- matrix(0, R, R)
  for (r in seq_len(R)) for (s in r:R) {
    a <- prs[[r]]; b <- prs[[s]]; val <- 0
    if (!length(a) && !length(b)) {
      val <- sum(SI * SI)
    } else if (!length(a) || !length(b)) {
      ## tr(S^-1 S^-1 u v') = sum((S^-1 v) * (S^-1 u))
      for (uv in c(a, b)) val <- val + sum(SIZ[[uv[1]]] * SIZ[[uv[2]]])
    } else {
      ## tr(S^-1 u v' S^-1 w z') = tr((v' S^-1 w)(z' S^-1 u))
      for (uv in a) for (wz in b) {
        m1 <- as.matrix(Matrix::crossprod(blk[[uv[2]]], SIZ[[wz[1]]]))
        m2 <- as.matrix(Matrix::crossprod(blk[[wz[2]]], SIZ[[uv[1]]]))
        val <- val + sum(m1 * t(m2))
      }
    }
    Kt[r, s] <- Kt[s, r] <- val
  }
  IE2 <- matrix(0, R, R)
  for (r in seq_len(R)) for (s in r:R)
    IE2[r, s] <- IE2[s, r] <- Kt[r, s] - 2 * sum(Phi * Q[[r, s]]) +
      sum((Phi %*% P[[r]]) * (P[[s]] %*% Phi))
  W <- 2 * ilm_kr_ginv(IE2)

  UU <- 0
  for (r in seq_len(R)) for (s in seq_len(R))
    UU <- UU + W[r, s] * (Q[[r, s]] - P[[r]] %*% Phi %*% P[[s]])
  PhiA <- Phi + 2 * Phi %*% UU %*% Phi
  dimnames(Phi) <- dimnames(PhiA) <- list(colnames(X), colnames(X))
  list(Phi = Phi, PhiA = (PhiA + t(PhiA)) / 2, P = P, W = W)
}

## The information's inverse, or its pseudo-inverse when a variance parameter
## is not identified -- pbkrtest's rule, with MASS::ginv()'s tolerance.
#' @keywords internal
#' @noRd
ilm_kr_ginv <- function(M) {
  ev <- eigen((M + t(M)) / 2, symmetric = TRUE)
  if (min(abs(ev$values)) > 1e-10) return(solve(M))
  keep <- abs(ev$values) > max(abs(ev$values)) * sqrt(.Machine$double.eps)
  V <- ev$vectors[, keep, drop = FALSE]
  V %*% (t(V) / ev$values[keep])
}

#' Kenward and Roger's F test: its denominator df and scale
#'
#' For `H0: L beta = 0` with `L` of full row rank `q`, as `pbkrtest` computes
#' them (Kenward and Roger 1997, section 3). The Wald F is formed with the
#' adjusted covariance, multiplied by `scale` and referred to F(q, df). For one
#' row the scale is exactly 1 and the df is Satterthwaite's formula with the
#' unadjusted covariance and the expected information, so a t test is the
#' estimate over its adjusted standard error, on this df.
#'
#' @param parts Output of [ilm_kr_parts()].
#' @param L Contrast matrix of full row rank.
#' @return A list with `df`, `scale` and `q`.
#' @keywords internal
#' @noRd
ilm_kr_ftest <- function(parts, L) {
  Phi <- parts$Phi; P <- parts$P; W <- parts$W
  LPL <- L %*% Phi %*% t(L)
  if (any(!is.finite(LPL)) || min(eigen((LPL + t(LPL)) / 2, symmetric = TRUE,
                                        only.values = TRUE)$values) <= 0)
    return(list(df = Inf, scale = 1, q = nrow(L)))
  Theta <- crossprod(L, solve(LPL, L))
  TP <- Theta %*% Phi
  U <- lapply(P, function(Pi) TP %*% Pi %*% Phi)
  A1 <- A2 <- 0
  for (i in seq_along(U)) for (j in seq_along(U)) {
    A1 <- A1 + W[i, j] * sum(diag(U[[i]])) * sum(diag(U[[j]]))
    A2 <- A2 + W[i, j] * sum(U[[i]] * t(U[[j]]))
  }
  q <- nrow(L)
  B <- (A1 + 6 * A2) / (2 * q)
  g <- ((q + 1) * A1 - (q + 4) * A2) / ((q + 2) * A2)
  den <- 3 * q + 2 * (1 - g)
  c1 <- g / den; c2 <- (q - g) / den; c3 <- (q + 2 - g) / den
  V0 <- 1 + c1 * B; V1 <- 1 - c2 * B; V2 <- 1 - c3 * B
  if (abs(V0) < 1e-10) V0 <- 0
  rho <- (1 / q) * ((1 - A2 / q) / V1)^2 * V0 / V2
  m <- 4 + (q + 2) / (q * rho - 1)
  scale <- if (abs(m - 2) < 1e-2) 1 else m * (1 - A2 / q) / (m - 2)
  list(df = m, scale = scale, q = q)
}
