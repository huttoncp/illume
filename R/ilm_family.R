#' Response distributions supported by illume
#'

## ---------------------------------------------------------------------------
## Accelerated failure time.
##
## All three are the same model with a different error distribution:
##
##     log T  =  X beta  +  scale * W
##
## and everything that differs between them is the distribution of W. Writing
## z = (log t - eta) / scale, the density and survivor function on the LOG scale
## are what the likelihood needs, and the -log(t) that turns a density in log t
## into one in t is the Jacobian.
##
## A coefficient is a log time ratio: +0.5 means survival times are exp(0.5)
## times longer, whatever the baseline hazard does. That is the appeal of the
## accelerated failure time reading over a proportional hazards one, and it is
## the same interpretation in all three.
##
## Censoring is not a special case here, it is the point: a subject still alive
## at the end of follow-up contributes the probability of surviving that long.
## The mechanism is the one ilm_censor() already provides.
## ---------------------------------------------------------------------------

#' @keywords internal
#' @noRd
ilm_aft_family <- function(name) {
  ## log density and log survivor of the standardised error W.
  ## logspace_add(0, z) is log(1 + exp(z)) computed without overflowing, which
  ## plain log1p(exp(z)) does for z beyond about 700.
  logf <- switch(name,
    lognormal   = function(z) dnorm(z, log = TRUE),
    loglogistic = function(z) z - 2 * logspace_add(0 * z, z),
    weibull     = function(z) z - exp(z))
  logS <- switch(name,
    lognormal   = function(z) log(pnorm(-z) + 1e-300),
    loglogistic = function(z) -logspace_add(0 * z, z),
    weibull     = function(z) -exp(z))
  logF <- switch(name,
    lognormal   = function(z) log(pnorm(z) + 1e-300),
    loglogistic = function(z) z - logspace_add(0 * z, z),
    weibull     = function(z) log(1 - exp(-exp(z)) + 1e-300))
  rw <- switch(name,
    lognormal   = function(n) stats::rnorm(n),
    loglogistic = function(n) stats::rlogis(n),
    ## log of a standard exponential is the smallest extreme value
    ## distribution, which is what makes T Weibull
    weibull     = function(n) log(stats::rexp(n)))

  list(name = name, link = "log", n_disp = 1L, disp_names = "log_scale",
       C_of = function(J) 1L, censorable = TRUE, positive = TRUE, aft = TRUE,
       nll = function(eta, y, w, disp, cens = NULL, logsig = NULL, ...) {
         sc <- if (is.null(logsig)) exp(disp[1]) else exp(logsig)
         lt <- log(y); z <- (lt - eta[, 1]) / sc
         if (is.null(cens) || !any(cens != 0L))
           return(-sum(w * (logf(z) - log(sc) - lt)))
         o <- cens == 0L; l <- cens < 0L; r <- cens > 0L
         out <- 0
         if (any(o)) out <- out - sum(w[o] * (logf(z[o]) - log(sc) - lt[o]))
         ## right-censored: known only to have survived past t
         if (any(r)) out <- out - sum(w[r] * logS(z[r]))
         ## left-censored: known only to have failed before t
         if (any(l)) out <- out - sum(w[l] * logF(z[l]))
         out
       },
       linkinv = function(e) exp(e),
       sim = function(eta, w, disp, logsig = NULL)
         exp(eta[, 1] +
             (if (is.null(logsig)) exp(disp[1]) else exp(logsig)) * rw(nrow(eta))),
       ## The residual machinery asks the family for these rather than carrying
       ## its own switch. Two separate switches over families is how ilm_rqr()
       ## and ilm_pearson_ovr() both came to assume a multinomial response.
       cdf = function(y, eta, sc) {
         z <- (log(y) - eta) / sc
         switch(name,
           lognormal   = stats::pnorm(z),
           loglogistic = stats::plogis(z),
           weibull     = -expm1(-exp(z)))
       })
}
#' Describes how a response distribution enters the model. The random-effect
#' machinery -- covariance structures, random slopes, smooths, AR(1) -- is
#' identical for every family and never sees the response; a family supplies only
#' the piece that turns a linear predictor into a log-likelihood.
#'
#' @details
#' Supported families:
#'
#' \describe{
#'   \item{`"gaussian"`}{continuous response, identity link, one dispersion
#'     parameter (the residual standard deviation).}
#'   \item{`"binomial"`}{two-category response, logit link. The response may be
#'     0/1, a two-level factor, `TRUE`/`FALSE`, or a character column -- all
#'     four are treated identically. For a factor or character response the
#'     **second** level is the one being modelled, as in [stats::glm()], so
#'     `factor(c("no", "yes"))` models the probability of `"yes"`. A proportion
#'     with `weights` giving the number of trials also works.}
#'   \item{`"poisson"`}{non-negative counts, log link, no dispersion parameter.}
#'   \item{`"nbinom"`}{counts with more variability than Poisson allows, log
#'     link, one dispersion parameter. Variance is `mu + mu^2/k`, so smaller `k`
#'     means more overdispersion; as `k` grows it approaches the Poisson.}
#'   \item{`"multinomial"`}{three or more unordered categories. The only family
#'     needing more than one linear predictor: with `J` categories it uses
#'     `J - 1` dimensions with sum-to-zero coding.}
#' }
#'
#' Zero-inflation, hurdle models, Tweedie and other specialised families are
#' deliberately out of scope. `glmmTMB` covers those well and there is nothing to
#' gain from a weaker reimplementation.
#'
#' @section On numerical stability:
#' The binomial log-likelihood uses `logspace_add()` rather than
#' `log(1 + exp(eta))`. They are mathematically identical, but the naive form
#' overflows to infinity once `eta` exceeds about 709, which an optimiser can
#' easily reach while exploring. The stable form is exact throughout.
#'
#' @param family Character: one of `"gaussian"`, `"binomial"`, `"poisson"`,
#'   `"nbinom"`, `"multinomial"`.
#' @return A list describing the family, with elements `name`, `link`,
#'   `n_disp` (number of dispersion parameters), `disp_names`, `C_of()` (number
#'   of linear predictor dimensions), `nll()`, `linkinv()` and `sim()`.
#' @references
#' McCullagh, P., & Nelder, J. A. (1989). *Generalized Linear Models*, 2nd ed.
#' Chapman & Hall.
#'
#' Hilbe, J. M. (2011). *Negative Binomial Regression*, 2nd ed. Cambridge
#' University Press.
#' @examples
#' ilm_family("poisson")$n_disp      # Poisson has no dispersion parameter
#' ilm_family("nbinom")$disp_names
#' @export
ilm_family <- function(family = c("gaussian", "binomial", "poisson",
                                  "nbinom", "beta", "multinomial",
                                  "ordinal", "ordinal_probit",
                                  "ordinal_cloglog",
                                  "weibull", "lognormal", "loglogistic",
                                  "rp", "rp_odds", "rp_normal")) {
  family <- match.arg(family)
  if (family %in% c("weibull", "lognormal", "loglogistic"))
    return(ilm_aft_family(family))
  if (startsWith(family, "ordinal")) return(ilm_ordinal_family(family))
  if (family %in% c("rp", "rp_odds", "rp_normal"))
    return(ilm_rp_family(family))

  ## Each nll() receives eta (N x C), the response, the frequency weights and
  ## the dispersion parameters on the log scale, and returns the negative
  ## log-likelihood.  Constants that do not involve parameters may be dropped.
  switch(family,

    gaussian = list(
      name = "gaussian", link = "identity", n_disp = 1L,
      disp_names = "log_sigma", C_of = function(J) 1L, censorable = TRUE,
      nll = function(eta, y, w, disp, cens = NULL, logsig = NULL, ...) {
        ## logsig is a dispersion MODEL: one log standard deviation per row,
        ## in place of the single parameter. dnorm vectorises over it, so
        ## nothing else in the likelihood changes.
        mu <- eta[, 1]; s <- if (is.null(logsig)) exp(disp[1]) else exp(logsig)
        if (is.null(cens) || !any(cens != 0L))
          return(-sum(w * dnorm(y, mu, s, log = TRUE)))
        ## A censored row contributes the probability of the interval it is
        ## known to lie in, not a density at a value it was never observed at.
        ## The upper tail is written as pnorm(mu - y) rather than
        ## 1 - pnorm(y - mu), which loses all its precision out in the tail.
        o <- cens == 0L; l <- cens < 0L; r <- cens > 0L
        out <- 0
        if (any(o)) out <- out - sum(w[o] * dnorm(y[o], mu[o], s, log = TRUE))
        if (any(l)) out <- out - sum(w[l] * log(pnorm((y[l] - mu[l]) / s) + 1e-300))
        if (any(r)) out <- out - sum(w[r] * log(pnorm((mu[r] - y[r]) / s) + 1e-300))
        out
      },
      linkinv = function(e) e,
      sim = function(eta, w, disp, logsig = NULL)
        stats::rnorm(nrow(eta), eta[, 1],
                     if (is.null(logsig)) exp(disp[1]) else exp(logsig))),

    binomial = list(
      name = "binomial", link = "logit", n_disp = 0L,
      disp_names = character(0), C_of = function(J) 1L,
      nll = function(eta, y, w, disp, ...) {
        e <- eta[, 1]
        ## y may be 0/1, or a proportion with w = number of trials
        -sum(w * (y * e - logspace_add(0 * e, e)))
      },
      linkinv = function(e) 1 / (1 + exp(-e)),
      sim = function(eta, w, disp, ...)
        stats::rbinom(nrow(eta), size = pmax(1, round(w)),
                      prob = 1 / (1 + exp(-eta[, 1]))) / pmax(1, round(w))),

    poisson = list(
      name = "poisson", link = "log", n_disp = 0L, zi_ok = TRUE,
      disp_names = character(0), C_of = function(J) 1L,
      nll = function(eta, y, w, disp, ...) {
        -sum(w * dpois(y, exp(eta[, 1]), log = TRUE))
      },
      ## Per row rather than summed, which is what a zero-inflated or hurdle
      ## likelihood needs: it has to weigh the density at the observed count
      ## against the density at zero, one row at a time.
      logden = function(eta, y, disp, logsig = NULL)
        dpois(y, exp(eta[, 1]), log = TRUE),
      p0 = function(mu, disp) exp(-mu),
      pcdf = function(q, mu, disp) stats::ppois(q, mu),
      ## zero-truncated draw by inverse CDF -- exact, and unlike rejection
      ## sampling it does not slow to a crawl when the mean is near zero
      rtrunc = function(n, mu, disp) {
        f0 <- exp(-mu)
        stats::qpois(f0 + stats::runif(n) * (1 - f0), mu)
      },
      linkinv = function(e) exp(e),
      sim = function(eta, w, disp, ...) stats::rpois(nrow(eta), exp(eta[, 1]))),

    nbinom = list(
      name = "nbinom", link = "log", n_disp = 1L, zi_ok = TRUE,
      disp_names = "log_k", C_of = function(J) 1L,
      nll = function(eta, y, w, disp, logsig = NULL, ...) {
        mu <- exp(eta[, 1]); k <- if (is.null(logsig)) exp(disp[1]) else exp(logsig)
        -sum(w * dnbinom2(y, mu = mu, var = mu + mu * mu / k, log = TRUE))
      },
      logden = function(eta, y, disp, logsig = NULL) {
        mu <- exp(eta[, 1]); k <- if (is.null(logsig)) exp(disp[1]) else exp(logsig)
        dnbinom2(y, mu = mu, var = mu + mu * mu / k, log = TRUE)
      },
      p0 = function(mu, disp) stats::dnbinom(0, size = disp, mu = mu),
      pcdf = function(q, mu, disp) stats::pnbinom(q, size = disp, mu = mu),
      rtrunc = function(n, mu, disp) {
        f0 <- stats::dnbinom(0, size = disp, mu = mu)
        stats::qnbinom(f0 + stats::runif(n) * (1 - f0), size = disp, mu = mu)
      },
      linkinv = function(e) exp(e),
      sim = function(eta, w, disp, logsig = NULL)
        stats::rnbinom(nrow(eta),
                       size = if (is.null(logsig)) exp(disp[1]) else exp(logsig),
                       mu = exp(eta[, 1]))),

    ## A proportion that is not a count of anything. Percent cover, share of
    ## time, a score already scaled onto the unit interval: a binomial fit
    ## needs a denominator these have none of, and a gaussian one puts mass
    ## outside [0, 1] and assumes a spread that does not shrink at the ends.
    ##
    ## Parameterised by the mean and a PRECISION, after Ferrari and
    ## Cribari-Neto: y ~ Beta(mu * phi, (1 - mu) * phi), so the mean is mu and
    ## the variance is mu (1 - mu) / (1 + phi). Larger phi means LESS spread,
    ## which is the opposite of every other dispersion parameter here and is
    ## worth saying twice.
    beta = list(
      name = "beta", link = "logit", n_disp = 1L,
      disp_names = "phi", C_of = function(J) 1L, unit = TRUE,
      ## A zero part is available, but only as a hurdle. A CONTINUOUS density
      ## has no mass at zero, so there is no "the count produced a zero on its
      ## own" route for a mixture to add -- every zero came from the zero
      ## process, which is what a hurdle says. The truncation factor is then 1
      ## rather than 1 - f(0), and `continuous` is the flag that says so.
      zi_ok = TRUE, continuous = TRUE,
      nll = function(eta, y, w, disp, logsig = NULL, ...) {
        mu <- 1 / (1 + exp(-eta[, 1]))
        phi <- if (is.null(logsig)) exp(disp[1]) else exp(logsig)
        -sum(w * dbeta(y, mu * phi, (1 - mu) * phi, log = TRUE))
      },
      logden = function(eta, y, disp, logsig = NULL) {
        mu <- 1 / (1 + exp(-eta[, 1]))
        phi <- if (is.null(logsig)) exp(disp[1]) else exp(logsig)
        dbeta(y, mu * phi, (1 - mu) * phi, log = TRUE)
      },
      p0 = function(mu, disp) rep(0, length(mu)),
      pcdf = function(q, mu, disp)
        stats::pbeta(pmin(pmax(q, 0), 1), mu * disp, (1 - mu) * disp),
      ## continuous, so the quantile residual is the distribution function at
      ## the observation with nothing to randomise across
      cdf = function(y, eta, sc) {
        mu <- stats::plogis(eta)
        stats::pbeta(y, mu * sc, (1 - mu) * sc)
      },
      linkinv = function(e) 1 / (1 + exp(-e)),
      sim = function(eta, w, disp, logsig = NULL, ...) {
        mu <- stats::plogis(eta[, 1])
        phi <- if (is.null(logsig)) exp(disp[1]) else exp(logsig)
        stats::rbeta(nrow(eta), mu * phi, (1 - mu) * phi)
      }),

    multinomial = list(
      name = "multinomial", link = "logit", n_disp = 0L,
      disp_names = character(0), C_of = function(J) J - 1L,
      ## y arrives as a count matrix (N x J) and Tct projects the C free
      ## dimensions onto all J categories under sum-to-zero coding
      nll = function(eta, y, w, disp, Tct, ...) {
        etaJ <- eta %*% Tct
        -(sum(y * etaJ) - sum(w * log(rowSums(exp(etaJ)))))
      },
      linkinv = function(e) e,
      sim = function(eta, w, disp, ..., Tct) {
        P <- exp(eta %*% Tct); P <- P / rowSums(P)
        cp <- t(apply(P, 1, cumsum))
        as.integer(rowSums(stats::runif(nrow(P)) > cp)) + 1L
      })
  )
}

## ---------------------------------------------------------------------------
## Cumulative link models for an ordered response.
##
## One linear predictor and a set of thresholds:
##
##     P(Y <= j | x)  =  F(theta_j - x'beta),    theta_1 < ... < theta_{J-1}
##
## so a single beta shifts every cut at once. That is the proportional-odds
## assumption under a logit link, and it is what buys an ordered outcome its
## economy: J - 1 categories described by one coefficient per predictor rather
## than J - 1 of them. It is also an assumption, which ilm_check_proportional()
## is here to test.
##
## The sign convention is the usual one for this model and it is a trap worth
## naming: the linear predictor is SUBTRACTED from the threshold, so a positive
## coefficient shifts probability towards the HIGHER categories. Reading it the
## other way reverses every conclusion.
##
## There is no intercept in x -- the thresholds are the intercepts, and a
## column of ones would be perfectly confounded with all of them at once.
##
## The thresholds are fitted as a first value and a set of log increments, so
## the ordering holds by construction rather than by a constraint the optimiser
## has to respect. A boundary is then unreachable rather than something to
## detect: an increment can approach zero only by its logarithm running to
## minus infinity, which the convergence checks see.
## ---------------------------------------------------------------------------

#' @keywords internal
#' @noRd
ilm_ordinal_family <- function(name) {
  link <- switch(name, ordinal = "logit", ordinal_probit = "probit",
                 ordinal_cloglog = "cloglog")
  ## the distribution function of the latent error, and its numeric twin for
  ## everything outside the likelihood
  pf <- switch(link,
    logit   = function(z) 1 / (1 + exp(-z)),
    probit  = function(z) pnorm(z),
    cloglog = function(z) 1 - exp(-exp(z)))
  pfn <- switch(link,
    logit   = stats::plogis, probit = stats::pnorm,
    cloglog = function(z) 1 - exp(-exp(z)))
  qfn <- switch(link,
    logit   = stats::qlogis, probit = stats::qnorm,
    cloglog = function(p) log(-log(1 - p)))
  list(
    name = name, link = link, ordinal = TRUE, n_disp = 0L,
    disp_names = character(0), C_of = function(J) 1L,
    ## `zeta` arrives as the thresholds themselves, already cumulated and
    ## ordered. y is the category index, 1 to J.
    nll = function(eta, y, w, disp, zeta = NULL, ...) {
      ## The outermost cuts are at plus and minus infinity, and there is no
      ## way to write that as a number in a vector the tape can hold: pasting
      ## sentinels onto `zeta` with c() drops its class and the tape breaks.
      ## Instead every row indexes a real threshold and a MASK, built from the
      ## response and so constant, selects whether that value or the limit is
      ## used. Nothing infinite is ever evaluated.
      e1 <- eta[, 1]
      idx <- attr(y, "ord_idx")
      Fu <- idx$mu * pf(zeta[idx$iu] - e1) + (1 - idx$mu)
      Fl <- idx$ml * pf(zeta[idx$il] - e1)
      -sum(w * log(Fu - Fl + 1e-300))
    },
    ## the cumulative probabilities a fit implies, for prediction and residuals
    pcum = function(eta, zeta) {
      outer(-as.numeric(eta), c(-Inf, as.numeric(zeta), Inf), `+`) * -1
    },
    pfun = pfn, qfun = qfn,
    linkinv = function(e) e,
    sim = function(eta, w, disp, zeta = NULL, ...) {
      ## draw the latent variable and see which band it falls in, which is the
      ## model's own account of where a category comes from
      n <- nrow(eta)
      z <- as.numeric(eta[, 1]) + qfn(stats::runif(n))
      as.integer(rowSums(outer(z, as.numeric(zeta), `>`))) + 1L
    })
}

#' Category probabilities implied by a cumulative link fit
#'
#' @param eta Linear predictor, one value per row.
#' @param zeta Thresholds, `J - 1` of them, increasing.
#' @param pfun The link's distribution function.
#' @return An `N x J` matrix of probabilities whose rows sum to one.
#' @keywords internal
#' @noRd
ilm_ord_probs <- function(eta, zeta, pfun) {
  eta <- as.numeric(eta); zeta <- as.numeric(zeta)
  N <- length(eta); J <- length(zeta) + 1L
  cm <- cbind(0, pfun(outer(-eta, zeta, `+`)), 1)
  p <- cm[, -1L, drop = FALSE] - cm[, -(J + 1L), drop = FALSE]
  ## a cut can round to the same double as the one below it in the far tail,
  ## which shows up as a zero or a whisker below it
  pmax(p, 0)
}

#' Check that a response is valid for its family
#'
#' Catches the common mistakes early and explains them, rather than letting the
#' optimiser fail obscurely later: counts that are negative or fractional,
#' binomial responses outside 0 to 1, and so on.
#'
#' @param y The response, after any conversion.
#' @param family A family object from [ilm_family()].
#' @return `TRUE`, invisibly; called for its error messages.
#' @keywords internal
#' @noRd
ilm_check_response <- function(y, family, has_zero_part = FALSE) {
  nm <- family$name
  if (nm == "gaussian" && !is.numeric(y))
    stop("gaussian family needs a numeric response", call. = FALSE)
  if (nm %in% c("poisson", "nbinom")) {
    if (!is.numeric(y) || any(y < 0))
      stop(nm, " family needs non-negative counts", call. = FALSE)
    if (any(abs(y - round(y)) > 1e-8))
      stop(nm, " family needs whole-number counts; the response has fractional ",
           "values. If these are rates, use an offset instead.", call. = FALSE)
  }
  if (isTRUE(family$ordinal)) {
    if (!is.numeric(y) || any(y < 1) || any(abs(y - round(y)) > 1e-8))
      stop("an ordinal response must be a factor or whole numbers giving the ",
           "category, from 1 upwards", call. = FALSE)
  }
  if (nm == "beta") {
    if (!is.numeric(y))
      stop("the beta family needs a numeric proportion", call. = FALSE)
    ## With a zero part in the model the zeros are accounted for, so only
    ## the ones at the other end are still without a likelihood.
    bad <- if (isTRUE(has_zero_part)) y >= 1 else y <= 0 | y >= 1
    if (any(bad))
      stop("the beta family is defined on the OPEN interval (0, 1), and ",
           sum(bad), " value", if (sum(bad) > 1L) "s are" else " is",
           " exactly ", if (isTRUE(has_zero_part)) "1" else "0 or 1",
           ". The beta density has no mass at either endpoint, so those rows ",
           "have no likelihood. Either the boundary values are a separate ",
           "process, in which case model them as one -- ",
           if (isTRUE(has_zero_part))
             "this model already does that for the zeros, and a value at 1 needs the same treatment, which is done by modelling 1 - y so the boundary is at zero"
           else
             "ilm_model(ziformula = ~ 1, zi_type = \"hurdle\") puts the zeros in their own part, and a boundary at 1 is the same problem on 1 - y",
           " -- or they are a rounding of interior values, in which case ",
           "ilm_squeeze() shifts them inside by the Smithson-Verkuilen ",
           "amount and says how far it moved them.", call. = FALSE)
  }
  if (nm == "binomial" && (any(y < 0) || any(y > 1)))
    stop("binomial family needs a 0/1 response, or a proportion between 0 and 1 ",
         "with weights giving the number of trials", call. = FALSE)
  invisible(TRUE)
}
