## ---- rescaling the fixed-effect columns inside the fit (item 6) -------------
##
## Columns on very different scales -- income in dollars beside an indicator
## -- make the optimiser's steps and the Hessian's differencing badly
## conditioned, though the model is the same model in any units. Measured: a
## zero-inflated negative binomial with income in dollars failed 12 fits of
## 12 (the Hessian not positive definite, three standard errors NaN) where
## income in thousands gave 11 of 12 clean. So ilm_fit() divides each
## column of the fixed design, the zero part's and the dispersion model's
## whose SD is outside ilm_rescale_band by that SD, fits there, and converts
## back EXACTLY:
## eta = X beta = X_s beta_s, so the random effects, the variance parameters
## and everything reported from them are identical in both coordinates, and
## only the coefficients move, by a linear map:
##   beta = A beta_s (+ the intercept's shift when centring), V = A V_s A',
##   a precision P = A^-T P_s A^-1, a Hessian likewise.
## Under REML beta is integrated out, and the change of variable scales the
## restricted likelihood by |det A|: the objective is put back by the log
## of it, so logLik() does not depend on the units either.
##
## Left alone, as the help says: offsets (a fixed coefficient of 1), the
## columns of a flexible baseline (functions of the response), constant
## columns and 0/1 indicators (already well scaled), a smooth's penalised
## part (it is not in X; its unpenalised columns are, and are rescaled like
## any other), and random-slope covariates, whose scale sits in the
## random-effect covariance.

## ONLY COLUMNS IN EXTREME UNITS. Rescaling exists to fix columns like income
## in dollars (SD 9,000) or a rate of SD 1e-4, not to move fits whose
## columns are already well scaled. Dividing every column by its SD moved
## those too, by factors like 0.97: the optimiser took another path and
## stopped elsewhere within its tolerance -- enough to fail tests pinned
## near it, and on a beta whose precision ran towards 1e12, to land in a
## basin 3.4 log-likelihood units worse. A column whose SD is inside the band
## is left exactly as it is, so a fit whose columns all are runs the
## unrescaled path bit for bit. The band is where nlminb is comfortable.
ilm_rescale_band <- c(1e-2, 1e2)

## Which columns to scale and by what. `centre` only where the matrix has an
## intercept to take up the shift, and only on the columns rescaled.
#' @keywords internal
#' @noRd
ilm_col_scales <- function(M, centre = FALSE, skip = integer(0)) {
  k <- ncol(M)
  s <- rep(1, k); m <- rep(0, k)
  if (is.null(M) || !k || !nrow(M)) return(list(s = s, m = m, int = NA_integer_, any = FALSE))
  cst <- vapply(seq_len(k), function(j) {
    v <- M[, j]; all(v == v[1])
  }, TRUE)
  int <- which(cst & M[1, ] != 0)[1]
  for (j in seq_len(k)) {
    if (cst[j] || j %in% skip) next
    v <- M[, j]
    if (all(v %in% c(0, 1))) next                 # an indicator: already fine
    sj <- stats::sd(v)
    if (!is.finite(sj) || sj <= 0) next
    if (sj >= ilm_rescale_band[1] && sj <= ilm_rescale_band[2]) next
    s[j] <- sj
    if (centre && !is.na(int)) m[j] <- mean(v)
  }
  list(s = s, m = m, int = int, any = any(s != 1 | m != 0))
}

## The scaled matrix.
#' @keywords internal
#' @noRd
ilm_apply_scales <- function(M, sc) {
  if (is.null(M) || !sc$any) return(M)
  out <- M                                         # names and attributes kept
  out[] <- sweep(sweep(M, 2L, sc$m, "-"), 2L, sc$s, "/")
  out
}

## The map from scaled coefficients to the user's, for one block of k
## coefficients: A (k x k) and the shift is carried by A's intercept row.
#' @keywords internal
#' @noRd
ilm_scale_map <- function(sc) {
  k <- length(sc$s)
  A <- diag(1 / sc$s, k, k)
  ## centred: the intercept takes back what the centring moved,
  ## beta_1 = beta_s,1 - sum_j m_j beta_s,j / s_j (its own s is 1, m 0)
  if (!is.na(sc$int) && any(sc$m != 0))
    A[sc$int, ] <- A[sc$int, ] - sc$m / sc$s
  A
}

## One map over a whole named parameter vector: the blocks called `beta`
## (C categories of p, category by category), `gzi` and `gamma` transform;
## everything else is the identity. Returns A and its inverse.
#' @keywords internal
#' @noRd
ilm_param_map <- function(nms, maps) {
  n <- length(nms)
  A <- diag(1, n, n)
  for (b in names(maps)) {
    ix <- which(nms == b)
    if (!length(ix)) next
    Ab <- maps[[b]]; k <- nrow(Ab)
    ncat <- length(ix) %/% k
    for (cc in seq_len(ncat)) {
      jj <- ix[(cc - 1L) * k + seq_len(k)]
      A[jj, jj] <- Ab
    }
  }
  list(A = A, Ainv = solve(A))
}
