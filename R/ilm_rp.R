## ---------------------------------------------------------------------------
## Royston-Parmar flexible parametric survival.
##
## The parametric survival models already here fix the SHAPE of the baseline:
## Weibull says the log cumulative hazard is a straight line in log time,
## log-logistic and lognormal say their own equivalents. When the shape is
## wrong, the fit is wrong everywhere, and ilm_plot_survival() exists precisely
## to catch that.
##
## Royston and Parmar replace the straight line with a restricted cubic spline:
##
##     g(S(t | x))  =  s(log t ; gamma)  +  x' beta
##
## with g one of three links, each of which makes a familiar model the special
## case where the spline is linear:
##
##   hazard   log(-log S)      proportional hazards      linear -> Weibull
##   odds     log((1 - S) / S) proportional odds         linear -> log-logistic
##   normal   -qnorm(S)        probit                    linear -> lognormal
##
## So the coefficients keep an interpretation -- a hazard ratio, an odds ratio
## -- while the baseline is free to bend. That is the whole appeal, and it is
## why the hazard scale is the default: exp(beta) is a hazard ratio, which is
## what a Cox model reports, without a Cox model's inability to say anything
## about the baseline.
##
## HOW IT IS FITTED HERE. The spline basis is a function of the OBSERVED log
## time, so it is data, and its columns simply join the model matrix. The
## spline coefficients are then part of beta and everything downstream --
## standard errors, the anova, the bootstrap -- works unchanged. The one thing
## the likelihood needs beyond eta is d(eta)/d(log t), which is the derivative
## basis times the same beta, so it is passed as a second design matrix with
## zeros in the covariate columns.
## ---------------------------------------------------------------------------

#' Restricted cubic spline basis, and its derivative
#'
#' The Durrleman-Simon natural cubic spline: cubic between the boundary knots
#' and linear outside them, which stops the fitted hazard doing something wild
#' in the tails where there is no data to restrain it.
#'
#' With `m` knots the basis has `m - 1` columns: a linear term and `m - 2`
#' cubic ones. The intercept is not included, since the model matrix already
#' has one.
#'
#' @param x Log time.
#' @param knots Increasing vector of knot positions on the log-time scale, the
#'   first and last being the boundaries.
#' @param deriv Return the derivative with respect to `x` instead of the basis.
#' @return A matrix with `length(knots) - 1` columns.
#' @keywords internal
#' @noRd
ilm_rcs <- function(x, knots, deriv = FALSE) {
  m <- length(knots)
  k1 <- knots[1L]; km <- knots[m]
  out <- matrix(0, length(x), m - 1L)
  out[, 1L] <- if (deriv) 1 else x
  if (m > 2L) {
    pw <- if (deriv) 2L else 3L
    mult <- if (deriv) 3 else 1
    for (j in 2:(m - 1L)) {
      lam <- (km - knots[j]) / (km - k1)
      out[, j] <- mult * (pmax(x - knots[j], 0)^pw -
                          lam * pmax(x - k1, 0)^pw -
                          (1 - lam) * pmax(x - km, 0)^pw)
    }
  }
  colnames(out) <- paste0("rcs", seq_len(m - 1L))
  out
}

#' Knot positions for a flexible parametric survival model
#'
#' Boundary knots at the extremes of the uncensored log event times, interior
#' knots at evenly spaced quantiles of them. Censored times say when someone
#' was last seen, not when the hazard changed, so they do not place knots.
#'
#' @param logt Log times.
#' @param event `1` where the event was observed.
#' @param df Degrees of freedom the spline spends beyond the intercept. `1` is
#'   a straight line and so the corresponding parametric model.
#' @return An increasing vector of `df + 1` knots.
#' @keywords internal
#' @noRd
ilm_rp_knots <- function(logt, event, df) {
  le <- logt[event == 1L]
  if (length(le) < 2L)
    stop("at least two events are needed to place knots for a flexible ",
         "parametric baseline; ", length(le), " observed.", call. = FALSE)
  if (df < 1L) stop("`rp_df` must be at least 1", call. = FALSE)
  nu <- length(unique(le))
  if (df + 1L > nu)
    stop("`rp_df` is ", df, ", which needs ", df + 1L, " knots, but there are ",
         "only ", nu, " distinct event times to place them at.", call. = FALSE)
  ## Knots are placed at quantiles of the events, so each interval between
  ## them has to hold enough events to say anything about the shape there.
  if (df > 1L && nu / df < 10)
    warning("`rp_df` is ", df, " with only ", nu, " distinct event times, ",
            "about ", round(nu / df, 1), " per degree of freedom. The baseline ",
            "will follow the noise between knots; rp_df of 3 to 5 is usual, ",
            "and ilm_rp_lrt() says whether even that is earning its place.",
            call. = FALSE)
  if (df == 1L) return(range(le))
  p <- seq(0, 1, length.out = df + 1L)
  k <- unname(stats::quantile(le, p, names = FALSE))
  k[1L] <- min(le); k[df + 1L] <- max(le)
  if (any(diff(k) <= 0))
    stop("the event times do not have ", df + 1L, " distinct quantiles, so a ",
         "spline with rp_df = ", df, " has repeated knots. Use a smaller ",
         "rp_df, or pass `rp_knots` yourself.", call. = FALSE)
  k
}

#' Is the fitted survival curve monotone?
#'
#' A spline baseline can bend enough to make the fitted survivor function rise,
#' which is not a survivor function. The derivative of the linear predictor
#' with respect to log time has to stay positive, and there is nothing in the
#' likelihood that forces it to.
#'
#' @param object A fitted `"ilm_model"` with a flexible parametric baseline.
#' @return The proportion of observations where it is violated.
#' @keywords internal
#' @noRd
ilm_rp_monotone <- function(object) {
  if (is.null(object$rp)) return(0)
  d <- as.vector(object$rp$D %*% object$beta[, 1L])
  mean(d <= 0)
}

#' Flexible parametric survival baselines
#'
#' The three Royston-Parmar links. Each takes the linear predictor and returns
#' the log density and log survivor contributions; the `log(deta)` Jacobian and
#' the `-log(t)` that turns a density in log time into one in time are added by
#' the likelihood, which is where the derivative lives.
#'
#' @param name `"rp"`, `"rp_odds"` or `"rp_normal"`.
#' @return A family list, as [ilm_family()] returns.
#' @keywords internal
#' @noRd
ilm_rp_family <- function(name) {
  scale <- switch(name, rp = "hazard", rp_odds = "odds", rp_normal = "normal")
  ## The same three shapes the accelerated failure time families use, which is
  ## not a coincidence: a linear spline makes each of them the corresponding
  ## parametric model exactly.
  logf <- switch(scale,
    hazard = function(e) -exp(e),
    odds   = function(e) -logspace_add(0 * e, e),
    normal = function(e) dnorm(e, log = TRUE))
  logS <- switch(scale,
    hazard = function(e) -exp(e),
    odds   = function(e) -logspace_add(0 * e, e),
    normal = function(e) log(pnorm(-e) + 1e-300))
  logF <- switch(scale,
    hazard = function(e) log(1 - exp(-exp(e)) + 1e-300),
    odds   = function(e) e - logspace_add(0 * e, e),
    normal = function(e) log(pnorm(e) + 1e-300))
  ## the density is the survivor function differentiated, and on the hazard and
  ## odds scales that adds a term the survivor function does not have
  extra <- switch(scale,
    hazard = function(e) e,
    odds   = function(e) e - logspace_add(0 * e, e),
    normal = function(e) 0 * e)

  list(name = name, link = "identity", n_disp = 0L,
       disp_names = character(0), C_of = function(J) 1L,
       censorable = TRUE, positive = TRUE, rp = TRUE, rp_scale = scale,
       nll = function(eta, y, w, disp, cens = NULL, etad = NULL, ...) {
         e <- eta[, 1]; lt <- log(y)
         ## d(eta)/d(log t) must be positive for a survivor function to fall,
         ## and nothing in the likelihood forces it to. pmax() is not available
         ## on an automatic-differentiation type, and would give a zero
         ## gradient where it bit even if it were, so the floor is a softplus:
         ## logspace_add(0, k x) / k, which is positive everywhere, smooth, and
         ## within 5e-6 of x itself once x exceeds about 0.3. A typical value
         ## here is 1 / scale, around 1.4, so the approximation costs nothing
         ## in the valid region and the likelihood falls away smoothly outside
         ## it -- which is the penalty, so no separate one is needed.
         ## ilm_rp_monotone() reports whether the answer landed inside.
         dpos <- logspace_add(0 * etad, 20 * etad) / 20
         pen <- 0
         if (is.null(cens) || !any(cens != 0L))
           return(pen - sum(w * (logf(e) + extra(e) + log(dpos) - lt)))
         o <- cens == 0L; l <- cens < 0L; r <- cens > 0L
         out <- pen
         if (any(o))
           out <- out - sum(w[o] * (logf(e[o]) + extra(e[o]) +
                                    log(dpos[o]) - lt[o]))
         if (any(r)) out <- out - sum(w[r] * logS(e[r]))
         if (any(l)) out <- out - sum(w[l] * logF(e[l]))
         out
       },
       linkinv = function(e) e,
       ## Simulating needs the inverse of the fitted survivor function, which
       ## has no closed form once the baseline is a spline. ilm_simulate() and
       ## every envelope built on it therefore go through ilm_rp_sim(), which
       ## inverts it numerically.
       sim = function(eta, w, disp, ...)
         stop("simulating from a flexible parametric baseline needs the ",
              "fitted knots, so it goes through ilm_rp_sim(); this should not ",
              "be reached", call. = FALSE),
       surv = function(e) exp(logS(e)),
       ## The residual machinery asks the family for this rather than keeping
       ## its own switch. eta already carries the spline at the observed time,
       ## so neither y nor a scale is needed here.
       cdf = function(y, eta, sc) 1 - exp(logS(eta)))
}

## Invert the fitted survivor function by bisection on log time. S is monotone
## in log t wherever the fit is valid, so bisection is reliable and needs no
## derivative; 60 steps take the bracket below 1e-15 of its width.
## The linear predictor with the baseline removed: what stays fixed while the
## survivor function is inverted for a time. Built from the conditional eta so
## that random effects come along.
#' @keywords internal
#' @noRd
ilm_rp_offset <- function(object, conditional = TRUE) {
  rp <- object$rp
  eta <- as.numeric(ilm_eta_hat(object, conditional)[, 1L])
  base <- as.vector(ilm_rcs(log(as.numeric(object$y)), rp$knots) %*%
                    object$beta[rp$cols, 1L])
  eta - base
}

#' @keywords internal
#' @noRd
ilm_rp_sim <- function(object, n = NULL, u = NULL, off = NULL) {
  rp <- object$rp
  if (is.null(rp)) stop("not a flexible parametric model", call. = FALSE)
  if (is.null(off)) {
    Xc <- object$X[, -rp$cols, drop = FALSE]
    bc <- object$beta[-rp$cols, 1L, drop = FALSE]
    off <- as.vector(Xc %*% bc)
  }
  gam <- object$beta[rp$cols, 1L]
  N <- length(off)
  if (is.null(u)) u <- stats::runif(N)
  Sfun <- object$family$surv
  lo <- rep(rp$lo, N); hi <- rep(rp$hi, N)
  for (k in seq_len(60L)) {
    mid <- (lo + hi) / 2
    s <- Sfun(as.vector(ilm_rcs(mid, rp$knots) %*% gam) + off)
    up <- s > u                    # still alive at mid, so the event is later
    lo <- ifelse(up, mid, lo)
    hi <- ifelse(up, hi, mid)
  }
  exp((lo + hi) / 2)
}

#' Was the flexible baseline worth it?
#'
#' Compares a flexible parametric fit against the straight-line special case it
#' generalises -- Weibull on the hazard scale, log-logistic on the odds scale,
#' lognormal on the normal scale -- by refitting with `rp_df = 1`.
#'
#' @details
#' The comparison is a likelihood-ratio test on `rp_df - 1` degrees of freedom.
#' It is a test of the baseline SHAPE, not of any covariate: a large statistic
#' says the straight line was wrong, which is the thing the flexible model
#' exists to fix. Because the knots are placed from the data rather than fixed
#' in advance, treat the p-value as indicative.
#'
#' @param object A fitted `"ilm_model"` with a flexible parametric baseline.
#' @param verbose Print the verdict.
#' @return Invisibly, a list with both log-likelihoods, the statistic, its
#'   degrees of freedom, a p-value and a `status`.
#' @seealso [ilm_plot_survival()], which checks the same thing against the
#'   Kaplan-Meier estimate rather than against the special case.
#' @examples
#' set.seed(1)
#' d <- data.frame(x = rnorm(300))
#' d$time <- exp(1.5 + 0.8 * d$x + 0.6 * log(rexp(300)))
#' d$event <- 1L
#' f <- ilm_model(time ~ x, data = d, family = "rp", rp_df = 3,
#'                censor = ilm_surv(d$time, d$event), verbose = FALSE)
#' ilm_rp_lrt(f)
#' @export
ilm_rp_lrt <- function(object, verbose = TRUE) {
  if (!inherits(object, "ilm_model") || is.null(object$rp))
    stop("`object` must be a model fitted with family = \"rp\", \"rp_odds\" ",
         "or \"rp_normal\"", call. = FALSE)
  df <- object$rp$df
  if (df <= 1L)
    stop("this model already has a straight-line baseline (rp_df = 1), so ",
         "there is nothing to compare it against.", call. = FALSE)
  cl <- stats::getCall(object)
  cl$rp_df <- 1L
  cl$verbose <- FALSE
  base <- try(suppressWarnings(eval(cl, parent.frame())), silent = TRUE)
  if (inherits(base, "try-error"))
    stop("the straight-line baseline would not fit, so the comparison cannot ",
         "be made: ", conditionMessage(attr(base, "condition")), call. = FALSE)
  l1 <- as.numeric(stats::logLik(object)); l0 <- as.numeric(stats::logLik(base))
  stat <- 2 * (l1 - l0); k <- df - 1L
  p <- stats::pchisq(stat, k, lower.tail = FALSE)
  simple <- switch(object$family$rp_scale, hazard = "weibull",
                   odds = "loglogistic", normal = "lognormal")
  status <- if (!is.finite(p)) "INCONCLUSIVE" else
            if (p < 0.01) "FLEXIBLE" else if (p < 0.05) "UNCLEAR" else "SIMPLE"
  if (verbose) {
    cat(sprintf("\nflexible baseline against the %s special case\n", simple))
    cat(sprintf("  logLik %.4f vs %.4f, chi-square %.2f on %d df, p = %s\n",
                l1, l0, stat, k, format.pval(p, digits = 3)))
    cat(switch(status,
      FLEXIBLE = sprintf("  the straight line is too rigid; keep the spline\n"),
      UNCLEAR  = sprintf("  borderline; the spline is doing a little work\n"),
      SIMPLE   = sprintf("  a straight line fits as well: family = \"%s\" is simpler and says the same\n",
                         simple),
      "  no verdict\n"))
  }
  invisible(list(logLik = l1, logLik_simple = l0, statistic = stat, df = k,
                 p.value = p, simple = simple, status = status))
}

## The survival curve for a flexible baseline. The linear predictor depends on
## time through the spline, so unlike an accelerated failure time model there
## is no single eta to evaluate: the basis is rebuilt at every time asked for.
#' @keywords internal
#' @noRd
ilm_rp_survival <- function(object, newdata, times, conf) {
  rp <- object$rp
  if (is.null(times)) {
    y <- as.numeric(object$y)
    times <- seq(min(y[y > 0]), max(y), length.out = 100L)
  }
  times <- sort(unique(times[times > 0]))
  ## the covariate part, without the baseline, for each requested pattern
  Xn <- ilm_newX(object, newdata)$X
  keep <- setdiff(seq_len(ncol(object$X)), rp$cols)
  if (ncol(Xn) != length(keep))
    stop("the new data give ", ncol(Xn), " covariate columns where the fit ",
         "has ", length(keep), call. = FALSE)
  off <- as.vector(Xn %*% object$beta[keep, 1L])
  V <- object$sdr$cov.fixed
  gam <- object$beta[rp$cols, 1L]
  Bs <- ilm_rcs(log(times), rp$knots)
  crit <- stats::qnorm(1 - (1 - conf) / 2)
  Sfun <- object$family$surv
  ## The variance of eta at a given time comes from the whole coefficient
  ## vector, spline included, so the interval widens where the baseline is
  ## least well determined as well as where the covariates are.
  ## cov.fixed is indexed by TMB's own parameter names, which are the names of
  ## the parameter list -- "beta", "theta" -- not the display names the package
  ## shows the user. Using the latter matched nothing and gave a zero-width
  ## interval.
  pn <- names(object$sdr$par.fixed)
  bidx <- if (!is.null(V) && all(is.finite(V)) && !is.null(pn))
    which(pn == "beta") else integer(0)
  out <- lapply(seq_along(off), function(i) {
    e <- as.vector(Bs %*% gam) + off[i]
    S <- Sfun(e)
    se <- rep(0, length(times))
    if (length(bidx) == ncol(object$X)) {
      Z <- cbind(Bs, matrix(rep(Xn[i, ], each = length(times)),
                            nrow = length(times)))
      Vb <- V[bidx, bidx, drop = FALSE]
      se <- sqrt(pmax(rowSums((Z %*% Vb) * Z), 0))
    }
    data.frame(row = i, time = times, surv = S,
               lower = Sfun(e + crit * se), upper = Sfun(e - crit * se))
  })
  do.call(rbind, out)
}
