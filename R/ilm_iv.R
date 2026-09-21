## ---------------------------------------------------------------------------
## Instrumental variables, by two-stage least squares.
##
## A regressor correlated with the error term makes its coefficient biased, and
## no amount of data fixes it. An instrument is a variable that moves the
## regressor without acting on the outcome any other way, and the part of the
## regressor it explains is clean.
##
## Three things about this are routinely got wrong.
##
## THE STANDARD ERRORS. Running lm() twice by hand gives the right point
## estimate and the WRONG standard errors, because the second stage computes
## its residuals from the FITTED regressor rather than from the actual one.
## Writing x = xhat + v, the naive residual is u + v*beta where u is the
## structural one, so the error is of size beta^2 Var(v) + 2 beta Cov(u, v)
## -- and Cov(u, v) is not zero, because the regressor being endogenous is
## exactly the statement that it is not. The direction therefore depends on
## the sign of the coefficient and of that covariance: in the example in
## ?ilm_iv the naive standard errors come out 24% too LARGE, and flipping the
## sign of beta sends them the other way. "Wrong" is the claim, not "small".
##
## The structural residual is y - X beta with the regressor as observed, and
## that is what goes into sigma^2 here.
##
## THE STRENGTH OF THE INSTRUMENT. A weak instrument does not merely widen the
## interval, it biases the estimate towards the one it was meant to fix and
## leaves the interval too narrow at the same time. The usual rule of thumb is
## a first-stage F above 10. Lee, McCrary, Moreira and Porter (2022) show that
## threshold is far too lenient: a conventional 5% t-test needs an F above
## about 104 before it really has 5% size. Both numbers are reported.
##
## WHETHER IT WAS NEEDED. If the regressor is not actually endogenous, least
## squares is unbiased and much more precise, and using an instrument throws
## that away. The Durbin-Wu-Hausman test asks, and is reported unasked.
##
## References:
##   Wooldridge, J. M. (2010). Econometric Analysis of Cross Section and Panel
##     Data, 2nd ed. MIT Press, chapter 5.
##   Lee, D. S., McCrary, J., Moreira, M. J. and Porter, J. (2022). Valid t-ratio
##     inference for IV. American Economic Review 112, 3260-3290.
##   Anderson, T. W. and Rubin, H. (1949). Estimation of the parameters of a
##     single equation in a complete system of stochastic equations. Annals of
##     Mathematical Statistics 20, 46-63.
## ---------------------------------------------------------------------------

#' Split a two-part instrumental-variables formula
#'
#' @keywords internal
#' @noRd
ilm_iv_parts <- function(formula) {
  if (!inherits(formula, "formula") || length(formula) != 3L)
    stop("`formula` must be two-sided, of the form y ~ x + w | z + w",
         call. = FALSE)
  rhs <- formula[[3L]]
  if (!(is.call(rhs) && identical(rhs[[1L]], as.name("|"))))
    stop("`formula` needs a `|` separating the model from the instruments: ",
         "y ~ x + w | z + w, where everything to the right of the bar is the ",
         "FULL list of exogenous variables -- the instruments AND any ",
         "covariate that appears on both sides.", call. = FALSE)
  env <- environment(formula)
  txt <- function(e) paste(deparse(e), collapse = " ")
  list(main = stats::as.formula(paste(txt(formula[[2L]]), "~", txt(rhs[[2L]])),
                                env),
       inst = stats::as.formula(paste("~", txt(rhs[[3L]])), env))
}

#' Instrumental-variables regression by two-stage least squares
#'
#' For a linear model in which one or more regressors are correlated with the
#' error -- a treatment people chose rather than were assigned, a price set in
#' response to demand, an exposure measured with error. An instrument moves the
#' offending regressor without acting on the outcome any other way, and the
#' part of the regressor it explains is used in its place.
#'
#' @section Writing the formula:
#'
#' `y ~ x + w | z + w`. To the left of the bar is the model you want. To the
#' right is the **full** list of exogenous variables: the instruments, plus
#' every covariate that appears on both sides. Leaving `w` out of the right
#' side does not mean "w is fine", it means "w is endogenous too and has no
#' instrument", which is a different and usually unintended model. That case is
#' an error here rather than a silent one.
#'
#' @section What it reports without being asked:
#'
#' * **First-stage F**, one per endogenous regressor, on the excluded
#'   instruments. Above 10 is the familiar rule of thumb; above about 104 is
#'   what a conventional 5% t-test actually needs, after Lee et al. (2022).
#'   Both are shown, because the gap between them is the point.
#' * **Durbin-Wu-Hausman**, testing whether the regressor was endogenous at
#'   all. If it was not, least squares is unbiased and far more precise, and
#'   instrumenting has thrown that away for nothing.
#' * **Sargan's J**, when there are more instruments than endogenous
#'   regressors, testing whether they agree with each other. It cannot test
#'   whether they are all invalid together, which is the assumption that
#'   matters and the one no test reaches.
#'
#' When the instrument is weak the remedy is [ilm_iv_ar()], which builds a
#' confidence set that stays valid however weak it is.
#'
#' @param formula `y ~ x + w | z + w`; see above.
#' @param data A data frame.
#' @param weights Optional frequency weights.
#' @param cluster Optional one-sided formula or vector for cluster-robust
#'   standard errors; see [ilm_vcov_cluster()].
#' @param robust Use heteroskedasticity-robust standard errors. Implied by
#'   `cluster`.
#' @param na.action How to handle missing values.
#' @return An object of class `"ilm_iv"`: the coefficients, their standard
#'   errors and tests, the diagnostics above, and the pieces they were built
#'   from.
#' @seealso [ilm_iv_ar()] for weak-instrument-robust intervals,
#'   [ilm_dag_model()] for choosing an adjustment set when no instrument is
#'   available.
#' @references Lee, D. S., McCrary, J., Moreira, M. J. and Porter, J. (2022).
#'   Valid t-ratio inference for IV. *American Economic Review* 112, 3260-3290.
#' @examples
#' set.seed(1); n <- 500
#' z <- rnorm(n); u <- rnorm(n)
#' x <- 0.8 * z + u + rnorm(n)          # x is correlated with the error u
#' y <- 1 + 0.5 * x + u + rnorm(n)
#' d <- data.frame(y = y, x = x, z = z)
#' coef(lm(y ~ x, data = d))            # biased upwards by the shared u
#' ilm_iv(y ~ x | z, data = d)
#' @export
ilm_iv <- function(formula, data, weights = NULL, cluster = NULL,
                   robust = FALSE, na.action = stats::na.omit) {
  pp <- ilm_iv_parts(formula)
  cl <- match.call()
  ## one model frame for both sides, so the same rows survive na.action
  allv <- unique(c(all.vars(pp$main), all.vars(pp$inst)))
  wexpr <- substitute(weights)
  wnm <- if (is.null(wexpr)) NULL else all.vars(wexpr)
  full <- stats::reformulate(setdiff(c(allv[-1L], wnm), allv[1L]),
                             response = as.name(allv[1L]),
                             env = environment(formula))
  mf <- stats::model.frame(full, data, na.action = na.action)
  y <- as.numeric(stats::model.response(mf))
  X <- stats::model.matrix(stats::terms(pp$main, data = mf), mf)
  Z <- stats::model.matrix(stats::terms(pp$inst, data = mf), mf)
  n <- length(y); k <- ncol(X)
  w <- if (is.null(wnm)) rep(1, n) else as.numeric(mf[[wnm]])
  if (any(w <= 0)) stop("`weights` must be positive", call. = FALSE)

  ## Anything in the model but not in the instrument list is being treated as
  ## endogenous. That is a real model, but it is almost never the one someone
  ## meant to write, and it is invisible in the output.
  endo <- setdiff(colnames(X), colnames(Z))
  excl <- setdiff(colnames(Z), colnames(X))
  if (!length(endo))
    stop("every regressor also appears in the instrument list, so nothing is ",
         "being instrumented and this is ordinary least squares. Put the ",
         "endogenous regressor on the left of the bar only.", call. = FALSE)
  if (length(excl) < length(endo))
    stop("there are ", length(endo), " endogenous regressor",
         if (length(endo) > 1L) "s" else "", " (",
         paste(endo, collapse = ", "), ") but only ", length(excl),
         " excluded instrument", if (length(excl) != 1L) "s" else "",
         ". A coefficient cannot be identified with fewer instruments than ",
         "endogenous regressors. Did a covariate get left out of the right ",
         "side of the bar? Everything exogenous belongs on both sides.",
         call. = FALSE)

  sw <- sqrt(w)
  Xw <- X * sw; Zw <- Z * sw; yw <- y * sw
  qz <- qr(Zw)
  if (qz$rank < ncol(Zw))
    stop("the instrument matrix is rank deficient (rank ", qz$rank, " of ",
         ncol(Zw), "); two of these columns carry the same information.",
         call. = FALSE)
  ## the projection of the regressors onto the instruments -- the part of each
  ## regressor the instruments explain, which is the part that is clean
  Xhat <- qr.fitted(qz, Xw)
  A <- crossprod(Xhat)
  if (rcond(A) < 1e-12)
    stop("the projected design is singular: the instruments do not move the ",
         "endogenous regressor enough to tell its coefficient apart from the ",
         "others.", call. = FALSE)
  Ainv <- solve(A)
  beta <- as.numeric(Ainv %*% crossprod(Xhat, yw))
  names(beta) <- colnames(X)

  ## THE STRUCTURAL RESIDUAL: y - X beta, with X as observed. Using the
  ## second stage's own residuals (y - Xhat beta) is what makes a hand-rolled
  ## 2SLS report the wrong standard errors -- in either direction, depending
  ## on the sign of the coefficient and of the endogeneity.
  res <- y - as.numeric(X %*% beta)
  rw <- res * sw
  df <- n - k
  s2 <- sum(rw^2) / df
  V <- s2 * Ainv

  vtype <- "classical"
  if (!is.null(cluster) || isTRUE(robust)) {
    S <- Xhat * rw
    if (!is.null(cluster)) {
      g <- if (inherits(cluster, "formula")) {
        v <- all.vars(cluster)
        src <- if (all(v %in% names(mf))) mf else data
        gg <- interaction(src[v], drop = TRUE, sep = ":")
        if (length(gg) != n && !is.null(attr(mf, "na.action")))
          gg <- gg[-as.integer(attr(mf, "na.action"))]
        droplevels(gg)
      } else factor(cluster)
      if (length(g) != n)
        stop("`cluster` has ", length(g), " values against ", n, " rows used ",
             "by the fit.", call. = FALSE)
      G <- nlevels(g)
      Sg <- do.call(rbind, lapply(split(seq_len(n), g), function(i)
        colSums(S[i, , drop = FALSE])))
      meat <- crossprod(Sg) * (G / (G - 1)) * ((n - 1) / (n - k))
      vtype <- sprintf("cluster-robust (%d clusters)", G)
    } else {
      meat <- crossprod(S) * (n / (n - k))
      vtype <- "heteroskedasticity-robust"
    }
    V <- Ainv %*% meat %*% Ainv
  }
  dimnames(V) <- list(colnames(X), colnames(X))
  se <- sqrt(pmax(diag(V), 0))

  ## ---- diagnostics ---------------------------------------------------------
  Wx <- setdiff(colnames(X), endo)          # exogenous regressors
  fs <- lapply(endo, function(v) ilm_iv_first_stage(mf, v, Wx, excl, X, Z, w))
  names(fs) <- endo
  dwh <- ilm_iv_dwh(y, X, Z, w, endo, fs)
  sar <- ilm_iv_sargan(rw, Zw, length(excl) - length(endo))

  structure(list(
    coefficients = beta, se = se, vcov = V, vcov_type = vtype,
    residuals = res, fitted = as.numeric(X %*% beta),
    n = n, df = df, sigma = sqrt(s2),
    endogenous = endo, instruments = excl, exogenous = Wx,
    first_stage = fs, dwh = dwh, sargan = sar,
    X = X, Z = Z, Xhat = Xhat, y = y, weights = w, model = mf,
    formula = formula, call = cl), class = "ilm_iv")
}

#' First-stage fit for one endogenous regressor
#'
#' The F statistic is on the EXCLUDED instruments only: whether the exogenous
#' covariates predict the regressor is beside the point, since they are in the
#' structural equation too.
#'
#' @keywords internal
#' @noRd
ilm_iv_first_stage <- function(mf, v, Wx, excl, X, Z, w) {
  xv <- X[, v]
  sw <- sqrt(w)
  Zr <- Z[, unique(c(Wx, excl)), drop = FALSE] * sw
  Zn <- Z[, Wx, drop = FALSE] * sw
  if (!ncol(Zn)) Zn <- matrix(0, nrow(Z), 0L)
  xw <- xv * sw
  rss1 <- sum(stats::lsfit(Zr, xw, intercept = FALSE)$residuals^2)
  rss0 <- if (ncol(Zn))
    sum(stats::lsfit(Zn, xw, intercept = FALSE)$residuals^2) else sum(xw^2)
  q <- length(excl); dfr <- nrow(Z) - ncol(Zr)
  Fv <- ((rss0 - rss1) / q) / (rss1 / dfr)
  list(F = Fv, df1 = q, df2 = dfr,
       p = stats::pf(Fv, q, dfr, lower.tail = FALSE),
       partial_r2 = (rss0 - rss1) / rss0,
       resid = as.numeric(stats::lsfit(Zr, xw, intercept = FALSE)$residuals) / sw)
}

#' Durbin-Wu-Hausman: was the regressor endogenous at all?
#'
#' Add the first-stage residuals to the structural equation and test them
#' jointly. If they carry nothing, the regressor was exogenous and least
#' squares was the better estimator.
#'
#' @keywords internal
#' @noRd
ilm_iv_dwh <- function(y, X, Z, w, endo, fs) {
  sw <- sqrt(w)
  R <- do.call(cbind, lapply(fs, `[[`, "resid"))
  Xa <- cbind(X, R) * sw
  yw <- y * sw
  rss1 <- sum(stats::lsfit(Xa, yw, intercept = FALSE)$residuals^2)
  rss0 <- sum(stats::lsfit(X * sw, yw, intercept = FALSE)$residuals^2)
  q <- ncol(R); dfr <- nrow(X) - ncol(Xa)
  Fv <- ((rss0 - rss1) / q) / (rss1 / dfr)
  list(F = Fv, df1 = q, df2 = dfr,
       p = stats::pf(Fv, q, dfr, lower.tail = FALSE))
}

#' Sargan's J: do the instruments agree with each other?
#'
#' @keywords internal
#' @noRd
ilm_iv_sargan <- function(rw, Zw, df) {
  if (df <= 0L) return(list(statistic = NA_real_, df = 0L, p = NA_real_))
  fit <- stats::lsfit(Zw, rw, intercept = FALSE)
  r2 <- 1 - sum(fit$residuals^2) / sum(rw^2)
  J <- length(rw) * r2
  list(statistic = J, df = df, p = stats::pchisq(J, df, lower.tail = FALSE))
}

#' @export
print.ilm_iv <- function(x, ...) {
  cat("Instrumental-variables regression (2SLS)\n")
  cat("  endogenous:  ", paste(x$endogenous, collapse = ", "), "\n", sep = "")
  cat("  instruments: ", paste(x$instruments, collapse = ", "), "\n", sep = "")
  cat("  ", x$n, " observations; ", x$vcov_type, " standard errors\n\n", sep = "")
  tv <- x$coefficients / x$se
  ct <- cbind(Estimate = x$coefficients, `Std. Error` = x$se,
              `t value` = tv,
              `Pr(>|t|)` = 2 * stats::pt(-abs(tv), x$df))
  stats::printCoefmat(ct, signif.stars = TRUE, has.Pvalue = TRUE,
                      P.values = TRUE)
  cat("\nDiagnostics\n")
  for (v in names(x$first_stage)) {
    f <- x$first_stage[[v]]
    cat(sprintf("  first stage, %-10s F = %8.2f on %d and %d df (p = %.3g)\n",
                v, f$F, f$df1, f$df2, f$p))
  }
  minF <- min(vapply(x$first_stage, `[[`, 0, "F"))
  cat(sprintf("    rule of thumb F > 10: %s; F > 104.7 for a valid 5%% t-test: %s\n",
              if (minF > 10) "met" else "NOT met",
              if (minF > 104.7) "met" else "NOT met"))
  if (minF <= 104.7)
    cat("    A conventional t-test is not 5% here. ilm_iv_ar() gives an\n",
        "    interval that stays valid however weak the instrument is.\n",
        sep = "")
  cat(sprintf("  Durbin-Wu-Hausman    F = %8.2f on %d and %d df (p = %.3g)\n",
              x$dwh$F, x$dwh$df1, x$dwh$df2, x$dwh$p))
  if (x$dwh$p > 0.1)
    cat("    No evidence of endogeneity. If the regressor is exogenous, least\n",
        "    squares is unbiased and much more precise than this.\n", sep = "")
  if (x$sargan$df > 0L) {
    cat(sprintf("  Sargan J             %8.2f on %d df (p = %.3g)\n",
                x$sargan$statistic, x$sargan$df, x$sargan$p))
    cat("    This asks whether the instruments agree with each other. It\n",
        "    cannot ask whether they are all invalid together, which is the\n",
        "    assumption that matters.\n", sep = "")
  }
  invisible(x)
}

#' @export
coef.ilm_iv <- function(object, ...) object$coefficients

#' @export
vcov.ilm_iv <- function(object, ...) object$vcov

#' @export
residuals.ilm_iv <- function(object, ...) object$residuals

#' @export
fitted.ilm_iv <- function(object, ...) object$fitted

#' Anderson-Rubin confidence set for an instrumented coefficient
#'
#' A confidence set that is valid however weak the instrument is, which the
#' usual estimate-plus-or-minus-two-standard-errors is not. It works by testing
#' each candidate value directly: if the coefficient really is `b`, then
#' `y - x b` has the instrument's influence removed from it, so the instruments
#' should predict nothing of what is left. The set is every `b` that survives.
#'
#' @section It can be unbounded, and that is information:
#'
#' When the instrument is weak the set can run to infinity in one or both
#' directions, or be empty. An unbounded set is the honest statement that the
#' data cannot rule out arbitrarily large effects; a Wald interval in the same
#' situation reports a tidy finite range and is simply wrong. An empty set means
#' no value of the coefficient reconciles the instruments with each other, which
#' is evidence against the model rather than against any particular value.
#'
#' @param object An [ilm_iv()] fit with exactly one endogenous regressor.
#' @param level Confidence level.
#' @param range Optional `c(lo, hi)` to search over; by default a wide band
#'   around the 2SLS estimate.
#' @param n_grid Number of grid points.
#' @return A list with `lower`, `upper`, whether the set is `bounded`,
#'   `empty`, and the grid it was built from.
#' @references Anderson, T. W. and Rubin, H. (1949). Estimation of the
#'   parameters of a single equation in a complete system of stochastic
#'   equations. *Annals of Mathematical Statistics* 20, 46-63.
#' @seealso [ilm_iv()].
#' @examples
#' set.seed(2); n <- 400
#' z <- rnorm(n); u <- rnorm(n)
#' x <- 0.15 * z + u + rnorm(n)          # a weak instrument
#' y <- 1 + 0.5 * x + u + rnorm(n)
#' fit <- ilm_iv(y ~ x | z, data = data.frame(y = y, x = x, z = z))
#' ilm_iv_ar(fit)
#' @export
ilm_iv_ar <- function(object, level = 0.95, range = NULL, n_grid = 2001L) {
  if (!inherits(object, "ilm_iv"))
    stop("`object` must be an ilm_iv() fit, not ", class(object)[1],
         call. = FALSE)
  if (length(object$endogenous) != 1L)
    stop("an Anderson-Rubin set is built here for a single endogenous ",
         "regressor; this fit has ", length(object$endogenous), ".",
         call. = FALSE)
  v <- object$endogenous
  X <- object$X; Z <- object$Z; y <- object$y; w <- object$weights
  sw <- sqrt(w)
  xv <- X[, v] * sw
  Wx <- object$exogenous
  Zn <- if (length(Wx)) Z[, Wx, drop = FALSE] * sw else
    matrix(0, nrow(Z), 0L)
  Zf <- Z * sw
  q <- length(object$instruments)
  dfr <- nrow(Z) - ncol(Zf)
  crit <- stats::qf(level, q, dfr)

  b <- unname(object$coefficients[v]); s <- unname(object$se[v])
  if (is.null(range)) range <- c(b - 40 * max(s, 1e-6), b + 40 * max(s, 1e-6))
  grid <- seq(range[1L], range[2L], length.out = n_grid)
  yw <- y * sw
  keep <- vapply(grid, function(b0) {
    e <- yw - b0 * xv
    rss1 <- sum(stats::lsfit(Zf, e, intercept = FALSE)$residuals^2)
    rss0 <- if (ncol(Zn))
      sum(stats::lsfit(Zn, e, intercept = FALSE)$residuals^2) else sum(e^2)
    Fv <- ((rss0 - rss1) / q) / (rss1 / dfr)
    Fv <= crit
  }, TRUE)

  if (!any(keep))
    return(structure(list(lower = NA_real_, upper = NA_real_, bounded = TRUE,
                          empty = TRUE, level = level, grid = grid,
                          keep = keep, term = v), class = "ilm_iv_ar"))
  lo <- grid[which(keep)[1L]]
  hi <- grid[rev(which(keep))[1L]]
  ## a set that reaches either end of the search was not closed by it
  bounded <- !(keep[1L] || keep[length(keep)])
  structure(list(lower = if (keep[1L]) -Inf else lo,
                 upper = if (keep[length(keep)]) Inf else hi,
                 bounded = bounded, empty = FALSE, level = level,
                 grid = grid, keep = keep, term = v,
                 contiguous = all(diff(which(keep)) == 1L)),
            class = "ilm_iv_ar")
}

#' @export
print.ilm_iv_ar <- function(x, ...) {
  cat(sprintf("Anderson-Rubin %.0f%% confidence set for %s\n",
              100 * x$level, x$term))
  if (isTRUE(x$empty)) {
    cat("  EMPTY. No value of the coefficient reconciles the instruments with\n",
        "  each other, which is evidence against the model rather than\n",
        "  against any particular value.\n", sep = "")
    return(invisible(x))
  }
  cat(sprintf("  [%s, %s]\n",
              format(x$lower, digits = 4), format(x$upper, digits = 4)))
  if (!x$bounded)
    cat("  UNBOUNDED. The instrument is too weak to rule out arbitrarily\n",
        "  large effects. A Wald interval would report a tidy finite range\n",
        "  here and would be wrong.\n", sep = "")
  if (isFALSE(x$contiguous))
    cat("  The set is not a single interval; the endpoints above bracket it.\n")
  invisible(x)
}
