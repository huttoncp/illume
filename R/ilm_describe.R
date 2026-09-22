## ---------------------------------------------------------------------------
## Class-aware description of variables, with a gaussian agreement index.
##
## Built on collapse primitives and base stats, so the only cost beyond base R
## is collapse itself, which depends on nothing but Rcpp.
## ---------------------------------------------------------------------------

## ---- distribution fitting --------------------------------------------------
##
## Candidate families are fitted by maximum likelihood and ranked by AIC.
##
## Discrete and continuous candidates are NEVER compared against each other: a
## probability mass and a probability density are not on the same scale, so
## their AICs are not comparable. The data decide which branch is used.

#' @keywords internal
#' @noRd
ilm_fit_families <- function(x) {
  x <- x[is.finite(x)]
  n <- length(x)
  if (n < 5L || fsd(x) == 0) return(NULL)
  out <- list()
  add <- function(family, ll, k) {
    if (is.finite(ll)) out[[length(out) + 1L]] <<-
      data.frame(family = family, logLik = ll, k = k, stringsAsFactors = FALSE)
  }
  nll_opt <- function(par0, nll, lower) {
    fit <- try(stats::nlminb(par0, nll, lower = lower), silent = TRUE)
    if (inherits(fit, "try-error") || fit$convergence != 0) return(NULL)
    -fit$objective
  }

  is_count <- all(abs(x - round(x)) < 1e-8) && all(x >= 0)

  if (is_count) {
    mu <- fmean(x)
    if (mu > 0) add("poisson", sum(stats::dpois(x, mu, log = TRUE)), 1L)
    if (all(x %in% c(0, 1)))
      add("bernoulli", sum(stats::dbinom(x, 1, mu, log = TRUE)), 1L)
    if (mu > 0) {
      p <- 1 / (1 + mu)
      add("geometric", sum(stats::dgeom(x, p, log = TRUE)), 1L)
      ## negative binomial: mu fixed at the mean, size profiled
      ll <- nll_opt(1, function(s)
        -sum(stats::dnbinom(x, size = s, mu = mu, log = TRUE)), 1e-8)
      if (!is.null(ll)) add("negbinomial", ll, 2L)
    }
  } else {
    m <- fmean(x); s <- sqrt(fmean((x - m)^2))
    add("normal", sum(stats::dnorm(x, m, s, log = TRUE)), 2L)
    rng <- diff(range(x))
    if (rng > 0) add("uniform", -n * log(rng), 2L)

    if (all(x > 0)) {
      lx <- log(x); lm_ <- fmean(lx); ls <- sqrt(fmean((lx - lm_)^2))
      if (ls > 0) add("lognormal", sum(stats::dlnorm(x, lm_, ls, log = TRUE)), 2L)
      add("exponential", sum(stats::dexp(x, 1 / m, log = TRUE)), 1L)
      ## gamma, started from method of moments
      v <- fvar(x); sh0 <- m^2 / v; rt0 <- m / v
      ll <- nll_opt(c(sh0, rt0), function(p)
        -sum(stats::dgamma(x, shape = p[1], rate = p[2], log = TRUE)), c(1e-8, 1e-8))
      if (!is.null(ll)) add("gamma", ll, 2L)
      ll <- nll_opt(c(1, m), function(p)
        -sum(stats::dweibull(x, shape = p[1], scale = p[2], log = TRUE)), c(1e-8, 1e-8))
      if (!is.null(ll)) add("weibull", ll, 2L)
    }
    if (all(x > 0 & x < 1)) {
      ll <- nll_opt(c(1, 1), function(p)
        -sum(stats::dbeta(x, p[1], p[2], log = TRUE)), c(1e-8, 1e-8))
      if (!is.null(ll)) add("beta", ll, 2L)
    }
  }
  if (!length(out)) return(NULL)
  r <- do.call(rbind, out)
  r$AIC <- -2 * r$logLik + 2 * r$k
  r <- r[order(r$AIC), , drop = FALSE]

  ## Some candidates NEST others: negbinomial contains poisson as size -> Inf,
  ## and gamma and weibull both contain exponential at shape = 1. A nested
  ## model's likelihood can never be worse than its special case, so it wins on
  ## AIC whenever noise exceeds the 2-point penalty -- which makes genuinely
  ## Poisson data look overdispersed. Decide these pairs with a likelihood-ratio
  ## test instead, and keep the simpler model unless the extra parameter earns
  ## its place.
  nested <- list(
    list(simple = "poisson",     complex = "negbinomial", boundary = TRUE),
    list(simple = "exponential", complex = "gamma",       boundary = FALSE),
    list(simple = "exponential", complex = "weibull",     boundary = FALSE))
  for (nz in nested) {
    if (!identical(r$family[1], nz$complex) || !(nz$simple %in% r$family)) next
    lr <- 2 * (r$logLik[r$family == nz$complex] - r$logLik[r$family == nz$simple])
    if (!is.finite(lr) || lr < 0) lr <- 0
    ## the Poisson/negbinomial comparison sits on a parameter boundary, where
    ## the null distribution is a 50:50 mixture of chi-square(0) and chi-square(1)
    p <- stats::pchisq(lr, 1, lower.tail = FALSE) * if (nz$boundary) 0.5 else 1
    if (p > 0.05) {
      i <- which(r$family == nz$simple)
      r <- rbind(r[i, , drop = FALSE], r[-i, , drop = FALSE])
    }
  }
  r
}

## A guess, plus an honest statement of how much to trust it.
##
## Reporting one family and a probability would invite over-reading. Two
## families separated by less than 2 AIC are not distinguishable by these data,
## and below ~30 observations almost nothing is.
#' @keywords internal
#' @noRd
ilm_dist_guess <- function(x, min_n = 30L) {
  f <- ilm_fit_families(x)
  if (is.null(f)) return(list(dist = NA_character_, dist_delta = NA_real_,
                              dist_note = "not assessable"))
  n <- sum(is.finite(x))
  best <- f$family[1]
  delta <- if (nrow(f) > 1L) f$AIC[2] - f$AIC[1] else Inf
  note <- if (n < min_n) {
    sprintf("n < %d, not discriminable", min_n)
  } else if (delta < 0) {
    ## a nested rival had lower AIC but did not survive its likelihood-ratio
    ## test, so the simpler model was kept on parsimony
    sprintf("simpler than %s, whose extra parameter is not justified", f$family[2])
  } else if (delta < 2) {
    sprintf("or %s; not discriminable", f$family[2])
  } else if (delta < 10) {
    sprintf("weakly preferred over %s", f$family[2])
  } else "clearly preferred"
  delta <- abs(delta)
  list(dist = best, dist_delta = round(delta, 1), dist_note = note)
}

## ---- gaussian plausibility, and why not ------------------------------------
##
## This deliberately does NOT report a normality test p-value. Any test of
## "exactly normal" rejects everything once n is large, so the p-value answers
## a question nobody asked. What a researcher wants to know is how FAR from
## normal the variable is, which is an effect size.
##
## The measure is the Kolmogorov distance between the empirical distribution
## and the best-fitting normal -- the largest amount, on the cumulative
## probability scale, by which a normal model misstates the data. It converges
## to a population quantity instead of degenerating with n.
##
## Sampling noise alone produces a distance of about 0.6/sqrt(n) (the Lilliefors
## scale), so the score charges only the excess beyond what normal data of the
## same size would plausibly give. `gauss` is therefore an AGREEMENT INDEX on
## 0-1, not a probability, and is documented as such.

#' @keywords internal
#' @noRd
ilm_gauss_d <- function(x) {
  x <- sort(x[is.finite(x)]); n <- length(x)
  s <- fsd(x); if (!is.finite(s) || s == 0) return(NA_real_)
  p <- stats::pnorm(x, fmean(x), s)
  max(pmax(abs(p - seq_len(n) / n), abs(p - (seq_len(n) - 1) / n)))
}

## how many modes does a kernel density estimate show?  This is what catches
## mixtures -- the failure mode a single-family fit gets confidently wrong.
#' @keywords internal
#' @noRd
ilm_n_modes <- function(x, rel = 0.10, drop = 0.5) {
  x <- x[is.finite(x)]
  if (length(x) < 20 || fsd(x) == 0) return(1L)
  d <- try(stats::density(x, n = 512), silent = TRUE)
  if (inherits(d, "try-error")) return(1L)
  pk <- which(diff(sign(diff(d$y))) == -2L) + 1L
  pk <- pk[d$y[pk] > rel * max(d$y)]
  if (length(pk) <= 1L) return(max(1L, length(pk)))
  ## Height alone counts ripple as peaks, which made uniform, exponential and
  ## Poisson data all look "multimodal". Require PROMINENCE: two peaks are
  ## distinct only if the valley between them drops well below both.
  keep <- pk[1L]
  for (i in 2:length(pk)) {
    last <- keep[length(keep)]
    valley <- min(d$y[last:pk[i]])
    if (valley < drop * min(d$y[last], d$y[pk[i]])) {
      keep <- c(keep, pk[i])
    } else if (d$y[pk[i]] > d$y[last]) {
      keep[length(keep)] <- pk[i]
    }
  }
  length(keep)
}

## `cap` is the excess departure at which alignment reaches 0. It is not
## plucked from the air: canonical clear departures cluster just above 0.11
## (lognormal 0.110, t(3) 0.105, bimodal 0.108, Poisson 0.148), so 0.12 puts
## "obviously not normal" at the bottom of the scale instead of compressing
## everything interesting into the top of it. It is exposed as an argument.
#' How far a variable is from gaussian, and why
#'
#' Reports an agreement index on 0-1 and, when it is low, the reason.
#'
#' This deliberately does not report a normality test p-value. Any test of
#' exact normality rejects everything once `n` is large, so the p-value answers
#' a question nobody asked. What matters is how far from normal a variable is,
#' which is an effect size.
#'
#' The measure is the Kolmogorov distance between the data and the best-fitting
#' normal: the largest amount, on the cumulative probability scale, by which a
#' normal model misstates the data. Sampling noise alone produces a distance of
#' about `0.6/sqrt(n)`, so the index charges only the excess beyond what normal
#' data of the same size would give. It is an agreement index, not a
#' probability.
#'
#' @param x A numeric vector.
#' @param min_n Below this many observations the index is not computed.
#' @param cap The excess distance at which agreement reaches 0. The default
#'   0.12 is set so clear departures (lognormal 0.110, t(3) 0.105, bimodal
#'   0.108, Poisson 0.148) land at the bottom of the scale.
#' @return A list with `gauss` (0-1), `ks_d` (the raw distance) and
#'   `gauss_note` (empty when agreement is high).
#' @references
#' Lilliefors, H. W. (1967). On the Kolmogorov-Smirnov test for normality with
#' mean and variance unknown. Journal of the American Statistical Association,
#' 62(318), 399-402.
#' @examples
#' set.seed(1)
#' ilm_gauss_check(rnorm(500))
#' ilm_gauss_check(rlnorm(500))
#' ilm_gauss_check(c(rnorm(250), rnorm(250, 5)))
#' @export
ilm_gauss_check <- function(x, min_n = 20L, cap = 0.12) {
  x <- x[is.finite(x)]; n <- length(x)
  if (n < min_n || fsd(x) == 0)
    return(list(gauss = NA_real_, ks_d = NA_real_,
                gauss_note = if (n >= min_n) "constant (zero variance)"
                             else "n too small to assess"))
  D <- ilm_gauss_d(x)
  ## 95th percentile of D for genuinely normal data of this size, so the score
  ## charges only the departure that exceeds sampling noise
  Dref <- 0.9 / sqrt(n)
  score <- 1 - min(1, max(0, (D - Dref) / cap))

  if (score >= 0.95)
    return(list(gauss = round(score, 3), ks_d = round(D, 4), gauss_note = ""))

  m <- fmean(x); s <- fsd(x); z <- (x - m) / s
  ku <- fmean(z^4) - 3
  ## Moment skewness is unstable under heavy tails -- it called a symmetric
  ## t(3) sample "right-skewed". Bowley's quartile skewness is bounded and
  ## robust, so it reports direction rather than tail noise.
  q <- stats::quantile(x, c(.25, .5, .75), names = FALSE)
  bow <- if (q[3] > q[1]) (q[3] + q[1] - 2 * q[2]) / (q[3] - q[1]) else 0
  nu <- fndistinct(x)
  discrete <- all(abs(x - round(x)) < 1e-8) && nu <= max(20L, n / 50)
  reasons <- character(0)
  ## ordered by how much the reason should change what the user does next
  ## (a kernel density over discrete values is spiky, so skip modes there)
  if (!discrete && ilm_n_modes(x) >= 2L)
    reasons <- c(reasons, "multimodal (check for subgroups)")
  if (discrete) reasons <- c(reasons, sprintf("discrete (%d distinct values)", nu))
  ## A stack of identical values at one end of a continuous variable is not a
  ## shape a distribution produces; it is a limit. A detection floor, an
  ## instrument ceiling, a capped scale. It is worth naming before skewness is,
  ## because the remedy is different: ilm_censor(), not a transformation.
  if (!discrete && nu > 20L) {
    pmin_ <- fmean(x == min(x)); pmax_ <- fmean(x == max(x))
    if (pmin_ >= 0.02 && n * pmin_ >= 5)
      reasons <- c(reasons, sprintf("%.0f%% of values sit exactly at the minimum (%s): a floor, see ilm_censor()",
                                    100 * pmin_, format(min(x), digits = 4)))
    if (pmax_ >= 0.02 && n * pmax_ >= 5)
      reasons <- c(reasons, sprintf("%.0f%% of values sit exactly at the maximum (%s): a ceiling, see ilm_censor()",
                                    100 * pmax_, format(max(x), digits = 4)))
  }
  if (all(x >= 0) && (min(x) - 0) < 0.05 * s) reasons <- c(reasons, "bounded at zero")
  if (bow > 0.1) reasons <- c(reasons, "right-skewed")
  if (bow < -0.1) reasons <- c(reasons, "left-skewed")
  if (ku > 1) reasons <- c(reasons, "heavy-tailed")
  if (ku < -1) reasons <- c(reasons, "light-tailed / flat")
  if (!length(reasons)) reasons <- "departs from normal, no single dominant cause"
  list(gauss = round(score, 3), ks_d = round(D, 4),
       gauss_note = paste(utils::head(reasons, 2L), collapse = "; "))
}

## ---- describe --------------------------------------------------------------
##
## Every column here answers "what about this variable will break a model?",
## not merely "what does this variable look like". The default output is kept
## narrow on the principle that a summary nobody reads is worse than a short
## one: quartiles, skew and kurtosis are available on request, but `gauss` and
## `gauss_note` already say whether they are worth asking for.

## percentile column names follow whatever `probs` the user asked for:
## c(0, .5, 1) -> p0, p50, p100;  c(.025, .975) -> p2.5, p97.5
#' @keywords internal
#' @noRd
ilm_prob_names <- function(probs)
  ## formatC(digits = 6) pads to that width, which produced names like
  ## "p      0"; as.character on a rounded value gives p0 / p50 / p2.5
  paste0("p", trimws(as.character(round(probs * 100, 6))))

## ---- numeric ---------------------------------------------------------------

#' @keywords internal
#' @noRd
ilm_describe_num <- function(x, digits = 3,
                             gauss = c("index", "ks_d", "both", "none"),
                             probs = c(0, 0.5, 1), skew = FALSE, kurt = FALSE,
                             dispersion = TRUE, cap = 0.12, min_n = 20L) {
  gauss <- match.arg(gauss)
  if (!is.numeric(probs) || any(probs < 0 | probs > 1))
    stop("`probs` must be numeric values between 0 and 1", call. = FALSE)
  obs <- length(x); nn <- fnobs(x); na <- obs - nn
  xv <- x[is.finite(x)]
  s <- if (nn > 1) fsd(xv) else NA_real_
  m <- if (nn > 0) fmean(xv) else NA_real_
  z <- if (nn > 1 && is.finite(s) && s > 0) (xv - m) / s else numeric(0)

  out <- data.frame(
    obs = obs, n = nn, na = na,
    sum = round(fsum(xv), digits), mean = round(m, digits),
    sd = round(s, digits), se = round(s / sqrt(nn), digits),
    stringsAsFactors = FALSE)

  ## quantiles honour `digits` like every other numeric column -- they did not,
  ## which is why a large-scale variable printed unrounded percentiles
  if (length(probs)) {
    q <- if (nn > 0) round(fquantile(xv, probs), digits)
         else rep(NA_real_, length(probs))
    qd <- as.data.frame(as.list(q), stringsAsFactors = FALSE)
    names(qd) <- ilm_prob_names(probs)
    out <- cbind(out, qd)
  }
  if (skew) out$skew <- round(if (length(z) > 2) fmean(z^3) else NA_real_, digits)
  if (kurt) out$kurt <- round(if (length(z) > 3) fmean(z^4) - 3 else NA_real_, digits)

  out$p_zero <- round(if (nn > 0) fmean(xv == 0) else NA_real_, digits)
  ## Variance-to-mean ratio: 1 is Poisson, above 1 points at family="nbinom".
  ## Reported for any non-negative integer variable, with no guessing about
  ## which of those are "really" counts. A rule keyed on the support starting
  ## near zero was tried and produced an inconsistency users would notice: x1
  ## (1-100) got a ratio while x2 (101-200) did not, though they are the same
  ## kind of variable. The ratio is only meaningful if you intend to model the
  ## variable as a count, which the user knows and the function cannot, so
  ## `dispersion = FALSE` turns the column off.
  if (dispersion) {
    is_count <- nn > 0 && all(abs(xv - round(xv)) < 1e-8) && all(xv >= 0)
    out$dispersion <- if (is_count && is.finite(m) && m > 0)
      round(fvar(xv) / m, digits) else NA_real_
  }

  if (gauss != "none") {
    g <- ilm_gauss_check(xv, min_n = min_n, cap = cap)
    if (gauss %in% c("index", "both")) out$gauss <- g$gauss
    if (gauss %in% c("ks_d", "both"))  out$ks_d  <- g$ks_d
    out$gauss_note <- g$gauss_note
  }
  out
}

## ---- factor / character ----------------------------------------------------

#' @keywords internal
#' @noRd
ilm_describe_cat <- function(x, digits = 3, sep = "_", top = 3L, rare_n = 5L) {
  obs <- length(x); nn <- fnobs(x); na <- obs - nn
  chr <- as.character(x); vals <- chr[!is.na(chr)]
  tb <- sort(table(vals), decreasing = TRUE); nu <- length(tb)

  out <- data.frame(
    obs = obs, n = nn, na = na,
    ## "" is not NA to R but is almost always missing to the analyst
    n_empty = sum(trimws(vals) == ""),
    n_unique = nu,
    ordered = is.ordered(x),
    p_max = round(if (nu) max(tb) / sum(tb) else NA_real_, digits),
    ## the usual reason a factor model will not fit
    n_rare = sum(tb < rare_n),
    n_unused = if (is.factor(x)) sum(!(levels(x) %in% vals)) else 0L,
    case_variants = { u <- unique(vals); sum(table(tolower(trimws(u))) > 1L) },
    counts_tb = paste(sprintf("%s%s%d", names(tb)[seq_len(min(top, nu))],
                              sep, as.integer(tb)[seq_len(min(top, nu))]),
                      collapse = ", "),
    stringsAsFactors = FALSE)
  z <- character(0)
  if (out$n_empty > 0) z <- c(z, sprintf("%d empty strings (not NA)", out$n_empty))
  if (out$n_unused > 0) z <- c(z, sprintf("%d unused levels", out$n_unused))
  if (out$case_variants > 0) z <- c(z, sprintf("%d case/space variants", out$case_variants))
  if (out$n_rare > 0) z <- c(z, sprintf("%d rare levels (separation risk)", out$n_rare))
  if (is.finite(out$p_max) && out$p_max > 0.99) z <- c(z, "near-constant")
  out$note <- paste(utils::head(z, 2L), collapse = "; ")
  out
}

## ---- logical ---------------------------------------------------------------

#' @keywords internal
#' @noRd
ilm_describe_lgl <- function(x, digits = 3) {
  obs <- length(x); nn <- fnobs(x); na <- obs - nn
  p <- if (nn > 0) fmean(x) else NA_real_
  data.frame(obs = obs, n = nn, na = na,
             n_TRUE = sum(x, na.rm = TRUE), n_FALSE = sum(!x, na.rm = TRUE),
             p_TRUE = round(p, digits),
             note = if (is.finite(p) && (p < 0.01 || p > 0.99))
               "near-constant (separation risk)" else "",
             stringsAsFactors = FALSE)
}

## ---- Date / POSIXct --------------------------------------------------------

#' @keywords internal
#' @noRd
ilm_label_spacing <- function(md, secs) {
  if (!is.finite(md)) return(NA_character_)
  if (secs) {
    if (md == 1) return("1 sec"); if (md == 60) return("minutely")
    if (md == 3600) return("hourly"); if (md == 86400) return("daily")
    if (md == 604800) return("weekly")
    return(sprintf("%g secs", md))
  }
  if (isTRUE(all.equal(md, 1))) return("daily")
  if (isTRUE(all.equal(md, 7))) return("weekly")
  if (md >= 28 && md <= 31) return("monthly")
  if (md >= 90 && md <= 92) return("quarterly")
  if (md >= 365 && md <= 366) return("yearly")
  sprintf("%g days", md)
}

#' @keywords internal
#' @noRd
ilm_describe_time <- function(x, digits = 3) {
  obs <- length(x); nn <- fnobs(x); na <- obs - nn
  xv <- x[!is.na(x)]; u <- sort(unique(xv)); secs <- inherits(x, "POSIXt")
  d <- if (length(u) > 1L) as.numeric(diff(u), units = if (secs) "secs" else "days")
       else numeric(0)
  md <- if (length(d)) as.numeric(names(sort(table(round(d, 6)), decreasing = TRUE))[1]) else NA_real_
  n_gaps <- if (length(d) && is.finite(md) && md > 0) sum(d > 1.5 * md) else 0L
  ## calendar steps are not constant in days, so allow variation proportional
  ## to the usual step; real gaps are still caught
  regular <- length(d) > 0 && is.finite(md) &&
    all(abs(d - md) <= pmax(1e-6, 0.15 * md))

  out <- data.frame(
    obs = obs, n = nn, na = na, n_unique = length(u),
    start = if (length(u)) as.character(min(u)) else NA_character_,
    end   = if (length(u)) as.character(max(u)) else NA_character_,
    span_days = if (length(u) > 1L)
      round(as.numeric(difftime(max(u), min(u), units = "days")), digits) else 0,
    spacing = ilm_label_spacing(md, secs), regular = regular, n_gaps = n_gaps,
    n_dup_times = nn - length(u),
    tz = if (secs) { t0 <- attr(x, "tzone")
      if (is.null(t0) || !nzchar(t0[1])) "unset" else t0[1] } else NA_character_,
    stringsAsFactors = FALSE)
  z <- character(0)
  if (!regular && n_gaps > 0)
    z <- c(z, sprintf("irregular, %d gaps (AR(1) needs regular spacing)", n_gaps))
  else if (!regular) z <- c(z, "irregular spacing")
  ## Repeated dates are the norm in long-format panel data -- many units share
  ## a date -- so this is described, not warned about. Duplicates only break
  ## AR(1) indexing WITHIN a group, which is visible when `by` is used.
  if (out$n_dup_times > 0)
    z <- c(z, sprintf("repeated dates: %d unique across %d rows (long format?)",
                      length(u), nn))
  if (secs && identical(out$tz, "unset")) z <- c(z, "timezone unset")
  out$note <- paste(utils::head(z, 2L), collapse = "; ")
  out
}

## ---- dispatch --------------------------------------------------------------

#' @keywords internal
#' @noRd
ilm_is_time <- function(v) inherits(v, "Date") || inherits(v, "POSIXt")

#' @keywords internal
#' @noRd
ilm_class_of <- function(v) {
  if (is.logical(v)) "logical"
  else if (ilm_is_time(v)) "time"
  else if (is.numeric(v)) "numeric"
  else "categorical"
}

ILM_CLASSES <- c("numeric", "categorical", "logical", "time")

#' Class-aware description of one variable
#'
#' Summarises a vector, or one column of a data frame, with statistics chosen
#' for its class. Numeric variables also get a gaussian agreement index and,
#' where that is low, a plain-language reason.
#'
#' The default output is deliberately narrow. Quartiles, skewness and kurtosis
#' are available on request, but `gauss` and `gauss_note` already say whether
#' they are worth asking for.
#'
#' @param data A data frame, or a vector when `y` is `NULL`.
#' @param y Name of the column to summarise, or several names. With none,
#'   `data` is taken as the thing to describe: a bare vector, or a data frame
#'   whose columns are all of one kind. For a frame of mixed kinds use
#'   [ilm_describe_all()], which returns one table per kind.
#' @param by Optional character vector of grouping columns.
#' @param digits Rounding for numeric columns, quantiles included.
#' @param gauss One of `"index"` (the 0-1 agreement index), `"ks_d"` (the raw
#'   Kolmogorov distance behind it), `"both"`, or `"none"`.
#' @param probs Quantiles to report, as proportions. `c(0, 0.5, 1)` gives
#'   `p0`, `p50` and `p100`; any vector works and the column names follow it.
#' @param skew,kurt Add moment-based skewness and excess kurtosis.
#' @param dispersion Add the variance-to-mean ratio for non-negative integer
#'   variables. Meaningful only if you intend to model the variable as a count.
#' @param rare_n Levels with fewer observations than this count as rare.
#' @param cap,min_n Control the `gauss` index; see [ilm_gauss_check()].
#' @return A one-row data frame, or one row per group when `by` is used.
#' @seealso [ilm_describe_all()] for every column, [ilm_gauss_check()] for the
#'   index, [ilm_describe_na()] for missingness.
#' @examples
#' d <- ilm_sim()
#' ilm_describe(d, "score")
#' ilm_describe(d, "income", probs = c(0.1, 0.5, 0.9))
#' ilm_describe(d, "score", by = "grp")
#' @export
ilm_describe <- function(data, y = NULL, by = NULL, digits = 3,
                          gauss = c("index", "ks_d", "both", "none"),
                          probs = c(0, 0.5, 1), skew = FALSE, kurt = FALSE,
                          dispersion = TRUE, rare_n = 5L,
                          cap = 0.12, min_n = 20L) {
  gauss <- match.arg(gauss)
  ## An ilm_anomaly() result is accepted directly: the flagged rows are what
  ## the user wants to look at next, and rebuilding that subset by hand from
  ## `row` is both a papercut and a chance to line the wrong rows up.
  if (inherits(data, "ilm_anomaly"))
    data <- ilm_from_anomaly(data, "ilm_describe", score = TRUE)
  one <- function(v) {
    if (is.logical(v)) ilm_describe_lgl(v, digits)
    else if (ilm_is_time(v)) ilm_describe_time(v, digits)
    else if (is.numeric(v)) ilm_describe_num(v, digits, gauss, probs, skew, kurt, dispersion,
                                             cap, min_n)
    else ilm_describe_cat(v, digits, rare_n = rare_n)
  }

  ## Resolving what is actually being described.
  ##
  ## Passing a whole data frame and no `y` used to hand the frame itself to
  ## one(), which is not logical, not a time and not numeric, so it fell
  ## through to the categorical branch and failed with "the condition has
  ## length > 1" for any frame of more than one column -- and returned nonsense
  ## rather than failing for a frame of exactly one.
  ##
  ## Columns of different kinds do not share a set of statistics. That is
  ## precisely why ilm_describe_all() returns one table per kind, so a mixed
  ## frame is an error that names that function rather than a mangled table.
  cols <- NULL
  if (!is.null(y)) {
    miss <- setdiff(y, names(data))
    if (length(miss))
      stop("column", if (length(miss) > 1L) "s" else "", " not found in the ",
           "data: ", paste(miss, collapse = ", "), ". Available: ",
           paste(utils::head(names(data), 12), collapse = ", "),
           if (length(names(data)) > 12L) ", ..." else "", call. = FALSE)
    if (length(y) > 1L) cols <- y else x <- data[[y]]
  } else if (is.data.frame(data)) {
    kinds <- unique(vapply(data, ilm_class_of, ""))
    if (length(kinds) > 1L)
      stop("`data` holds columns of more than one kind (",
           paste(sort(kinds), collapse = ", "), "), which do not share a set ",
           "of statistics. Name a column with `y`, or use ilm_describe_all(), ",
           "which returns one table per kind.", call. = FALSE)
    if (!ncol(data))
      stop("`data` has no columns to describe", call. = FALSE)
    cols <- names(data)
  } else x <- data

  if (!is.null(cols)) {
    r <- do.call(rbind, lapply(cols, function(cn)
      ilm_describe(data, cn, by = by, digits = digits, gauss = gauss,
                   probs = probs, skew = skew, kurt = kurt,
                   dispersion = dispersion, rare_n = rare_n,
                   cap = cap, min_n = min_n)))
    nrep <- nrow(r) / length(cols)
    return(cbind(variable = rep(cols, each = nrep), r,
                 stringsAsFactors = FALSE))
  }

  if (is.null(by)) return(one(x))
  miss <- setdiff(by, names(data))
  if (length(miss))
    stop("`by` variable(s) not found in the data: ", paste(miss, collapse = ", "),
         call. = FALSE)
  g <- interaction(data[by], drop = TRUE)
  parts <- lapply(split(x, g), one)
  cbind(setNames(data.frame(names(parts), stringsAsFactors = FALSE),
                 paste(by, collapse = ".")),
        do.call(rbind, parts), stringsAsFactors = FALSE)
}

## A variable with one value tells you nothing a table of statistics can show,
## and it will break a model matrix. It is reported separately rather than
## padding every other section with a row of NAs.
#' @keywords internal
#' @noRd
ilm_constant_tbl <- function(data, cols) {
  if (!length(cols)) return(NULL)
  do.call(rbind, lapply(cols, function(cn) {
    v <- data[[cn]]; nn <- fnobs(v); u <- unique(v[!is.na(v)])
    data.frame(variable = cn, class = ilm_class_of(v), obs = length(v),
               n = nn, na = length(v) - nn,
               value = if (length(u)) as.character(u[1]) else NA_character_,
               stringsAsFactors = FALSE)
  }))
}

#' Describe every column of a data frame
#'
#' Applies [ilm_describe()] to each column and returns one table per class, so
#' columns with different summaries are not forced into a common shape.
#'
#' Constant columns are split into their own `constant` section rather than
#' padding every other table with rows of `NA`: a variable with one value tells
#' you nothing a statistic can show, and it will break a model matrix.
#'
#' @inheritParams ilm_describe
#' @param class `"all"`, or one or more of `"numeric"`, `"categorical"`,
#'   `"logical"`, `"time"`.
#' @return A named list of data frames, one per class present, plus `constant`
#'   when any column has a single value; a single data frame if only one
#'   section results.
#' @seealso [ilm_frame_issues()] for problems belonging to pairs of columns.
#' @examples
#' d <- ilm_sim()
#' ilm_describe_all(d)
#' ilm_describe_all(d, class = "numeric", skew = TRUE)
#' @export
ilm_describe_all <- function(data, by = NULL, digits = 3,
                              gauss = c("index", "ks_d", "both", "none"),
                              probs = c(0, 0.5, 1), skew = FALSE, kurt = FALSE,
                              dispersion = TRUE, class = "all", rare_n = 5L,
                              cap = 0.12, min_n = 20L) {
  gauss <- match.arg(gauss)
  ## An ilm_anomaly() result is accepted directly: the flagged rows are what
  ## the user wants to look at next, and rebuilding that subset by hand from
  ## `row` is both a papercut and a chance to line the wrong rows up.
  if (inherits(data, "ilm_anomaly"))
    data <- ilm_from_anomaly(data, "ilm_describe_all", score = TRUE)
  bad <- setdiff(class, c("all", ILM_CLASSES))
  if (length(bad))
    stop("unknown `class`: ", paste(sQuote(bad), collapse = ", "),
         ". Options are ", paste(sQuote(c("all", ILM_CLASSES)), collapse = ", "),
         ".", call. = FALSE)
  miss <- setdiff(by, names(data))
  if (length(miss))
    stop("`by` variable(s) not found in the data: ", paste(miss, collapse = ", "),
         call. = FALSE)

  ## grouping variables describe the split, so they are not also described
  cand <- setdiff(names(data), by)
  cl <- vapply(data[cand], ilm_class_of, "")
  ## constancy is judged on the whole column, so the same variables are set
  ## aside whether or not `by` is used
  is_const <- vapply(data[cand], function(v) fndistinct(v) <= 1L, TRUE)
  const <- ilm_constant_tbl(data, cand[is_const])
  cand <- cand[!is_const]; cl <- cl[!is_const]

  want <- if (identical(class, "all")) unique(cl) else intersect(class, unique(cl))
  res <- lapply(want, function(k) {
    cols <- cand[cl == k]
    if (!length(cols)) return(NULL)
    r <- do.call(rbind, lapply(cols, function(cn)
      ilm_describe(data, cn, by = by, digits = digits, gauss = gauss,
                    probs = probs, skew = skew, kurt = kurt,
                    dispersion = dispersion, rare_n = rare_n,
                    cap = cap, min_n = min_n)))
    nrep <- nrow(r) / length(cols)
    cbind(variable = rep(cols, each = nrep), r, stringsAsFactors = FALSE)
  })
  names(res) <- want
  res <- res[!vapply(res, is.null, TRUE)]
  ## The constant section is attached only for class = "all". Asking for one
  ## class should return that data frame, not a list of it plus something else:
  ## a return type that changes depending on whether the data happen to contain
  ## a constant column is a trap.
  if (identical(class, "all")) {
    if (!is.null(const)) res$constant <- const
  } else if (!is.null(const)) {
    message("ilm_describe_all: ", nrow(const),
            " constant column(s) set aside (", paste(const$variable, collapse = ", "),
            "); use class = \"all\" to see them")
  }
  if (length(res) == 1L) res[[1L]] else res
}

## ---- whole-frame issues ----------------------------------------------------

#' Problems that belong to pairs of columns
#'
#' Constant columns, identifier-like columns, duplicated columns and
#' near-perfect collinearity. These cannot live in a per-variable table because
#' they are properties of the data frame as a whole.
#'
#' A continuous variable is unique per row by construction, so uniqueness is
#' reported as identifier-like only for discrete-valued columns.
#'
#' @param data A data frame.
#' @param cor_cut Absolute correlation at or above which a numeric pair is
#'   reported as collinear.
#' @return A data frame with `issue`, `columns` and `detail`; zero rows when
#'   nothing is found.
#' @examples
#' ilm_frame_issues(ilm_sim())
#' @export
ilm_frame_issues <- function(data, cor_cut = 0.999) {
  out <- list()
  add <- function(issue, cols, detail)
    out[[length(out) + 1L]] <<- data.frame(issue = issue, columns = cols,
                                           detail = detail, stringsAsFactors = FALSE)
  n <- nrow(data)
  for (cn in names(data)) {
    v <- data[[cn]]; nu <- fndistinct(v)
    ## a continuous variable is unique per row by construction, so uniqueness
    ## only signals an identifier for discrete-valued columns
    discrete_like <- is.character(v) || is.factor(v) ||
      (is.numeric(v) && !is.logical(v) &&
         all(abs(v[!is.na(v)] - round(v[!is.na(v)])) < 1e-8))
    if (nu <= 1L) add("constant", cn, "zero variance; breaks the model matrix")
    else if (discrete_like && nu == sum(!is.na(v)) && nu == n)
      add("id_like", cn, "unique per row; a key, not a predictor")
  }
  nm <- names(data)
  if (length(nm) > 1L) for (i in seq_len(length(nm) - 1L)) for (j in (i + 1L):length(nm))
    if (identical(unname(data[[i]]), unname(data[[j]])))
      add("duplicate_columns", paste(nm[i], nm[j], sep = " = "), "identical values")
  num <- names(data)[vapply(data, function(v) is.numeric(v) && !is.logical(v), TRUE)]
  if (length(num) > 1L) {
    cm <- suppressWarnings(stats::cor(data[num], use = "pairwise.complete.obs"))
    for (i in seq_len(length(num) - 1L)) for (j in (i + 1L):length(num)) {
      r <- cm[i, j]
      if (is.finite(r) && abs(r) >= cor_cut)
        add("collinear", paste(num[i], num[j], sep = " ~ "),
            sprintf("r = %.4f; coefficients will be unstable", r))
    }
  }
  if (!length(out))
    return(data.frame(issue = character(0), columns = character(0),
                      detail = character(0), stringsAsFactors = FALSE))
  do.call(rbind, out)
}
