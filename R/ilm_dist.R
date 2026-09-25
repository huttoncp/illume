## ---------------------------------------------------------------------------
## A fitted model's response distribution, as functions.
##
## Everything family-specific about the response -- the inverse link, the
## dispersion's transform, the zero part's mixing or truncation, a censored
## row's contribution, an ordered response's thresholds -- lives here and in
## the family objects, so code scoring or simulating from a fit hands over
## linear predictors and internal parameters and never re-derives illume's
## parameterisation. The density agrees with the fit's own likelihood: summed
## over a fixed-effects fit's rows it is the optimised objective, exactly.
## ---------------------------------------------------------------------------

#' A fitted model's response distribution, as functions
#'
#' The density, distribution function, quantile function, random generator
#' and mean of a fit's response, in illume's own parameterisation, taking
#' linear predictors and internal parameters -- a draw from [ilm_draws()], say
#' -- rather than means.
#'
#' @details
#' **Inputs.** `lp` is a named list of linear predictors on the link scale:
#' `mu`, the mean's (an `n`-by-`C` matrix for a multinomial response); `zi`,
#' the zero part's, for a zero-inflated or hurdle model; and `disp`, the
#' dispersion model's, for a model with `dispformula` -- its data columns
#' only, since the dispersion's power of the mean is added here from
#' `par["mu_pow"]`. `par` holds the family's own parameters, named as the
#' blocks of [ilm_draws()]: `logdisp` (the log dispersion), `zeta_raw` (an
#' ordered response's thresholds, as fitted: a first value and log
#' increments) and `mu_pow`. Each entry of `lp` and `par` may have one value
#' or one per element, so a vector of rows by draws is scored in one call.
#'
#' **Censoring.** `status` is `-1` for a value known only to be at or below
#' `y`, `1` for one at or above it, and `0` for an observed one, as
#' [ilm_censor()] codes them; a censored row's density is the probability of
#' its interval.
#'
#' **A binomial response** is a proportion, with `size` its number of trials.
#' The density then includes the binomial coefficient, which the fit's
#' likelihood leaves out (see [logLik.ilm_model()]).
#'
#' **Categories.** For an ordered or multinomial response `y` is a category:
#' its position among `categories`, or its label. `mean()` gives each
#' category's probability. A multinomial response has no ordering, so no
#' distribution or quantile function.
#'
#' The flexible parametric survival families are not covered: their linear
#' predictor is a function of time itself; see [ilm_survival()].
#'
#' @param object A fitted `"ilm_model"`.
#' @return An object of class `"ilm_dist"`: a list with the functions
#'   `d(y, lp, par, log = FALSE, status = NULL, size = 1)`,
#'   `p(q, lp, par, lower.tail = TRUE, size = 1)` (`NULL` for a multinomial
#'   response), `q(prob, lp, par, size = 1)` (likewise), `r(lp, par, size = 1)`
#'   and `mean(lp, par, size = 1)`, and `family`, `discrete`, `zero` (the kind
#'   of zero part, or `NULL`), `categories`, `lp_names` and `par_names`.
#' @seealso [ilm_draws()], [ilm_matrices()].
#' @examples
#' set.seed(1)
#' d <- data.frame(x = rnorm(200))
#' d$y <- rpois(200, exp(0.5 + 0.4 * d$x))
#' fit <- ilm_model(y ~ x, data = d, family = "poisson", verbose = FALSE)
#' dd <- ilm_dist(fit)
#' lp <- list(mu = as.numeric(fit$X %*% fit$beta))
#' ## the fit's log-likelihood, row by row
#' sum(dd$d(d$y, lp, par = NULL, log = TRUE))
#' logLik(fit)
#' @export
ilm_dist <- function(object) {
  if (!inherits(object, "ilm_model"))
    stop("`object` must be a fitted ilm_model, not ", class(object)[1],
         call. = FALSE)
  fam <- object$family
  name <- if (is.null(fam)) "multinomial" else fam$name
  if (name %in% c("rp", "rp_odds", "rp_normal"))
    stop("ilm_dist() does not cover the flexible parametric survival ",
         "families: their linear predictor is a function of the time itself. ",
         "Use ilm_survival() for the survival curve.", call. = FALSE)
  zero <- if (!is.null(object$Zzi)) (if (is.null(object$zi_type)) "inflated"
                                     else object$zi_type) else NULL
  has_dm <- !is.null(object$disp_formula) || isTRUE(object$disp_mu)
  disp_mu <- isTRUE(object$disp_mu)
  cats <- if (name == "multinomial" || isTRUE(fam$ordinal)) object$ylevels
          else NULL
  J <- length(cats)
  cont_zero <- isTRUE(fam$continuous)
  aft <- isTRUE(fam$aft)

  ## ---- parameters, and the per-element log dispersion -------------------
  getp <- function(par, nm) {
    if (is.list(par)) return(par[[nm]])
    if (is.null(par)) return(NULL)
    v <- par[names(par) == nm]
    if (length(v)) unname(v) else NULL
  }
  linkinv <- if (!is.null(fam)) fam$linkinv else function(e) e
  ldisp <- function(lp, par) {
    if (is.null(fam) || fam$n_disp == 0L) return(NULL)
    if (has_dm) {
      ld <- if (!is.null(lp$disp)) lp$disp else 0
      if (disp_mu)
        ld <- ld + getp(par, "mu_pow") * log(abs(linkinv(lp$mu)) + 1e-8)
      return(ld)
    }
    ld <- getp(par, "logdisp")
    if (is.null(ld))
      stop("`par` needs \"logdisp\", the ", name, " family's log dispersion",
           call. = FALSE)
    ld
  }
  zeta_of <- function(par) {
    zr <- getp(par, "zeta_raw")
    if (is.null(zr))
      stop("`par` needs \"zeta_raw\", the thresholds as fitted", call. = FALSE)
    zr <- if (is.matrix(zr)) zr else matrix(zr, nrow = 1L)
    ## a first threshold and log increments: zeta_j = zr_1 + sum of the first
    ## j - 1 increments, as the fit builds them
    z <- zr
    if (ncol(zr) > 1L)
      for (j in 2:ncol(zr)) z[, j] <- z[, j - 1L] + exp(zr[, j])
    z
  }
  zprob <- function(lp) {
    z <- if (!is.null(lp$zi)) lp$zi else lp$hu
    if (is.null(z))
      stop("this model has a zero part: give its linear predictor as lp$zi",
           call. = FALSE)
    stats::plogis(z)
  }
  lvec <- function(y) {
    if (is.null(cats)) return(y)
    if (is.factor(y) || is.character(y)) match(as.character(y), cats)
    else as.integer(y)
  }

  ## ---- the families without a zero part ----------------------------------
  ## each takes the mean's linear predictor and the log dispersion (NULL when
  ## there is none), vectorised over elements
  base <- switch(name,
    gaussian = list(
      d = function(y, e, ld, log, size) stats::dnorm(y, e, exp(ld), log = log),
      p = function(q, e, ld, lt, size) stats::pnorm(q, e, exp(ld), lower.tail = lt),
      q = function(p, e, ld, size) stats::qnorm(p, e, exp(ld)),
      r = function(e, ld, size) stats::rnorm(length(e), e, exp(ld)),
      m = function(e, ld, size) e),
    binomial = list(
      d = function(y, e, ld, log, size)
        stats::dbinom(round(y * size), size, stats::plogis(e), log = log),
      p = function(q, e, ld, lt, size)
        stats::pbinom(floor(q * size + 1e-9), size, stats::plogis(e), lower.tail = lt),
      q = function(p, e, ld, size) stats::qbinom(p, size, stats::plogis(e)) / size,
      r = function(e, ld, size)
        stats::rbinom(length(e), size, stats::plogis(e)) / size,
      m = function(e, ld, size) stats::plogis(e)),
    poisson = list(
      d = function(y, e, ld, log, size) stats::dpois(y, exp(e), log = log),
      p = function(q, e, ld, lt, size) stats::ppois(q, exp(e), lower.tail = lt),
      q = function(p, e, ld, size) stats::qpois(p, exp(e)),
      r = function(e, ld, size) stats::rpois(length(e), exp(e)),
      m = function(e, ld, size) exp(e)),
    nbinom = list(
      d = function(y, e, ld, log, size)
        stats::dnbinom(y, size = exp(ld), mu = exp(e), log = log),
      p = function(q, e, ld, lt, size)
        stats::pnbinom(q, size = exp(ld), mu = exp(e), lower.tail = lt),
      q = function(p, e, ld, size) stats::qnbinom(p, size = exp(ld), mu = exp(e)),
      r = function(e, ld, size)
        stats::rnbinom(length(e), size = exp(ld), mu = exp(e)),
      m = function(e, ld, size) exp(e)),
    beta = list(
      d = function(y, e, ld, log, size) {
        mu <- stats::plogis(e); ph <- exp(ld)
        stats::dbeta(y, mu * ph, (1 - mu) * ph, log = log)
      },
      p = function(q, e, ld, lt, size) {
        mu <- stats::plogis(e); ph <- exp(ld)
        stats::pbeta(q, mu * ph, (1 - mu) * ph, lower.tail = lt)
      },
      q = function(p, e, ld, size) {
        mu <- stats::plogis(e); ph <- exp(ld)
        stats::qbeta(p, mu * ph, (1 - mu) * ph)
      },
      r = function(e, ld, size) {
        mu <- stats::plogis(e); ph <- exp(ld)
        stats::rbeta(length(e), mu * ph, (1 - mu) * ph)
      },
      m = function(e, ld, size) stats::plogis(e)),
    lognormal = list(
      d = function(y, e, ld, log, size) stats::dlnorm(y, e, exp(ld), log = log),
      p = function(q, e, ld, lt, size) stats::plnorm(q, e, exp(ld), lower.tail = lt),
      q = function(p, e, ld, size) stats::qlnorm(p, e, exp(ld)),
      r = function(e, ld, size) stats::rlnorm(length(e), e, exp(ld)),
      m = function(e, ld, size) exp(e + exp(2 * ld) / 2)),
    ## log T = eta + scale * log(E), E standard exponential: a Weibull with
    ## shape 1 / scale and scale exp(eta)
    weibull = list(
      d = function(y, e, ld, log, size)
        stats::dweibull(y, shape = exp(-ld), scale = exp(e), log = log),
      p = function(q, e, ld, lt, size)
        stats::pweibull(q, shape = exp(-ld), scale = exp(e), lower.tail = lt),
      q = function(p, e, ld, size) stats::qweibull(p, shape = exp(-ld), scale = exp(e)),
      r = function(e, ld, size)
        stats::rweibull(length(e), shape = exp(-ld), scale = exp(e)),
      m = function(e, ld, size) exp(e) * gamma(1 + exp(ld))),
    loglogistic = list(
      d = function(y, e, ld, log, size) {
        s <- exp(ld); z <- (log(y) - e) / s
        v <- stats::dlogis(z, log = TRUE) - log(s) - log(y)
        if (log) v else exp(v)
      },
      p = function(q, e, ld, lt, size)
        stats::plogis((log(q) - e) / exp(ld), lower.tail = lt),
      q = function(p, e, ld, size) exp(e + exp(ld) * stats::qlogis(p)),
      r = function(e, ld, size) exp(e + exp(ld) * stats::rlogis(length(e))),
      ## the mean exists only for a scale below one
      m = function(e, ld, size) {
        s <- exp(ld)
        ifelse(s < 1, exp(e) * (pi * s) / sin(pi * s), Inf)
      }),
    NULL)

  ## ---- an ordered response ------------------------------------------------
  if (isTRUE(fam$ordinal)) {
    pfun <- fam$pfun; qfun <- fam$qfun
    cum <- function(e, par) {                    # P(Y <= j), j = 1..J-1
      z <- zeta_of(par)
      z[rep_len(seq_len(nrow(z)), length(e)), , drop = FALSE] - e
    }
    probs <- function(e, par) {
      cm <- cbind(0, pfun(cum(e, par)), 1)
      pmax(cm[, -1L, drop = FALSE] - cm[, -(J + 1L), drop = FALSE], 0)
    }
  }

  ## ---- the public functions ---------------------------------------------
  d <- function(y, lp, par, log = FALSE, status = NULL, size = 1) {
    e <- lp$mu
    if (isTRUE(fam$ordinal)) e <- rep_len(e, max(length(e), length(y)))
    if (name == "multinomial" || isTRUE(fam$ordinal)) {
      P <- if (name == "multinomial") ilm_dist_softmax(e, J) else probs(e, par)
      v <- P[cbind(seq_len(nrow(P)), rep_len(lvec(y), nrow(P)))]
      return(if (log) log(v) else v)
    }
    ld <- ldisp(lp, par)
    n <- max(length(e), length(y))
    y <- rep_len(y, n); e <- rep_len(e, n)
    if (!is.null(ld)) ld <- rep_len(ld, n)
    size <- rep_len(size, n)
    if (is.null(zero)) {
      v <- base$d(y, e, ld, TRUE, size)
      if (!is.null(status)) {
        st <- rep_len(status, n)
        lo <- st < 0L; hi <- st > 0L
        ## the probability of the interval, floored as the fit's likelihood
        ## floors it
        if (any(lo))
          v[lo] <- log(base$p(y[lo], e[lo], ld[lo], TRUE, size[lo]) + 1e-300)
        if (any(hi))
          v[hi] <- log(base$p(y[hi], e[hi], ld[hi], FALSE, size[hi]) + 1e-300)
      }
      return(if (log) v else exp(v))
    }
    pz <- rep_len(zprob(lp), n)
    l0 <- if (cont_zero) rep(-Inf, n) else base$d(rep(0, n), e, ld, TRUE, size)
    v <- numeric(n)
    z0 <- y == 0
    if (zero == "hurdle") {
      v[z0] <- log(pz[z0])
      v[!z0] <- log1p(-pz[!z0]) + base$d(y[!z0], e[!z0], ld[!z0], TRUE, size[!z0]) -
        (if (cont_zero) 0 else log1p(-exp(l0[!z0])))
    } else {
      v[z0] <- log(pz[z0] + (1 - pz[z0]) * exp(l0[z0]))
      v[!z0] <- log1p(-pz[!z0]) + base$d(y[!z0], e[!z0], ld[!z0], TRUE, size[!z0])
    }
    if (log) v else exp(v)
  }

  p <- if (name == "multinomial") NULL else function(q, lp, par, lower.tail = TRUE,
                                                    size = 1) {
    e <- lp$mu
    if (isTRUE(fam$ordinal)) {
      qq <- lvec(q)
      e <- rep_len(e, max(length(e), length(qq)))
      cm <- cbind(pfun(cum(e, par)), 1)
      v <- cm[cbind(seq_len(nrow(cm)), pmin(pmax(rep_len(qq, nrow(cm)), 1L), J))]
      v[rep_len(qq, nrow(cm)) < 1] <- 0
      return(if (lower.tail) v else 1 - v)
    }
    ld <- ldisp(lp, par)
    n <- max(length(e), length(q))
    q <- rep_len(q, n); e <- rep_len(e, n)
    if (!is.null(ld)) ld <- rep_len(ld, n)
    size <- rep_len(size, n)
    if (is.null(zero)) return(base$p(q, e, ld, lower.tail, size))
    pz <- rep_len(zprob(lp), n)
    Fc <- base$p(q, e, ld, TRUE, size)
    v <- if (zero == "hurdle" && !cont_zero) {
      f0 <- exp(base$d(rep(0, n), e, ld, TRUE, size))
      pz + (1 - pz) * pmax(Fc - f0, 0) / (1 - f0)
    } else pz + (1 - pz) * Fc
    v[q < 0] <- 0
    if (lower.tail) v else 1 - v
  }

  q <- if (name == "multinomial") NULL else function(prob, lp, par, size = 1) {
    e <- lp$mu
    if (isTRUE(fam$ordinal)) {
      e <- rep_len(e, max(length(e), length(prob)))
      cm <- cbind(pfun(cum(e, par)), 1)
      pr <- rep_len(prob, nrow(cm))
      return(as.integer(rowSums(cm < pr)) + 1L)
    }
    ld <- ldisp(lp, par)
    n <- max(length(e), length(prob))
    prob <- rep_len(prob, n); e <- rep_len(e, n)
    if (!is.null(ld)) ld <- rep_len(ld, n)
    size <- rep_len(size, n)
    if (is.null(zero)) return(base$q(prob, e, ld, size))
    pz <- rep_len(zprob(lp), n)
    out <- numeric(n)
    up <- prob > pz
    a <- (prob[up] - pz[up]) / (1 - pz[up])
    if (zero == "hurdle" && !cont_zero) {
      f0 <- exp(base$d(rep(0, sum(up)), e[up], ld[up], TRUE, size[up]))
      a <- f0 + a * (1 - f0)
    }
    out[up] <- base$q(a, e[up], ld[up], size[up])
    out
  }

  r <- function(lp, par, size = 1) {
    e <- lp$mu
    if (name == "multinomial") {
      P <- ilm_dist_softmax(e, J)
      return(as.integer(rowSums(stats::runif(nrow(P)) > t(apply(P, 1L, cumsum)))) + 1L)
    }
    if (isTRUE(fam$ordinal)) {
      z <- e + qfun(stats::runif(length(e)))
      zt <- zeta_of(par)
      zt <- zt[rep_len(seq_len(nrow(zt)), length(e)), , drop = FALSE]
      return(as.integer(rowSums(zt < z)) + 1L)
    }
    n <- length(e)
    ld <- ldisp(lp, par)
    if (!is.null(ld)) ld <- rep_len(ld, n)
    size <- rep_len(size, n)
    if (is.null(zero)) return(base$r(e, ld, size))
    pz <- rep_len(zprob(lp), n)
    out <- numeric(n)
    pos <- stats::runif(n) >= pz
    if (zero == "hurdle" && !cont_zero) {
      f0 <- exp(base$d(rep(0, sum(pos)), e[pos], ld[pos], TRUE, size[pos]))
      out[pos] <- base$q(f0 + stats::runif(sum(pos)) * (1 - f0), e[pos], ld[pos],
                         size[pos])
    } else out[pos] <- base$r(e[pos], ld[pos], size[pos])
    out
  }

  mean <- function(lp, par, size = 1) {
    e <- lp$mu
    if (name == "multinomial") return(ilm_dist_softmax(e, J))
    if (isTRUE(fam$ordinal)) return(probs(e, par))
    n <- length(e)
    ld <- ldisp(lp, par)
    if (!is.null(ld)) ld <- rep_len(ld, n)
    m <- base$m(e, ld, rep_len(size, n))
    if (is.null(zero)) return(m)
    pz <- rep_len(zprob(lp), n)
    if (zero == "hurdle" && !cont_zero) {
      f0 <- exp(base$d(rep(0, n), e, ld, TRUE, rep_len(size, n)))
      (1 - pz) * m / (1 - f0)
    } else (1 - pz) * m
  }

  structure(list(
    d = d, p = p, q = q, r = r, mean = mean,
    family = name,
    discrete = name %in% c("binomial", "poisson", "nbinom", "multinomial") ||
      isTRUE(fam$ordinal),
    zero = zero, categories = cats,
    lp_names = c("mu", if (!is.null(zero)) "zi",
                 if (!is.null(object$disp_formula)) "disp"),
    par_names = c(if (!is.null(fam) && fam$n_disp > 0L && !has_dm) "logdisp",
                  if (isTRUE(fam$ordinal)) "zeta_raw",
                  if (disp_mu) "mu_pow")),
    class = "ilm_dist")
}

## Category probabilities from a multinomial fit's C sum-to-zero linear
## predictors, as the fit maps them: onto all J categories through the
## sum-to-zero contrasts, then softmax.
#' @keywords internal
#' @noRd
ilm_dist_softmax <- function(e, J) {
  e <- if (is.matrix(e)) e else matrix(e, ncol = J - 1L)
  eJ <- e %*% t(stats::contr.sum(J))
  eJ <- eJ - apply(eJ, 1L, max)
  P <- exp(eJ)
  P / rowSums(P)
}

#' @export
print.ilm_dist <- function(x, ...) {
  cat("<ilm_dist> the response distribution of a", x$family, "fit",
      if (!is.null(x$zero)) paste0("with a zero part (", x$zero, ")") else "",
      "\n  linear predictors:", paste(x$lp_names, collapse = ", "),
      "\n  parameters:", if (length(x$par_names))
        paste(x$par_names, collapse = ", ") else "none",
      "\n  functions: d", if (!is.null(x$p)) ", p, q", ", r, mean\n", sep = " ")
  invisible(x)
}
