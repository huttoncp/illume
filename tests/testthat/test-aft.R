# Accelerated failure time: weibull, lognormal and log-logistic.
#
# All three are the same model with a different error distribution,
#
#     log T = X beta + scale * W,
#
# so a coefficient is a log time ratio: +0.5 means survival times are exp(0.5)
# times longer. Censoring is not a special case, it is the point -- a subject
# still event-free at the end of follow-up contributes the probability of
# having survived that long.
#
# The fixed-effects case is pinned to survival::survreg(), which fits exactly
# these three models.

surv_gen <- function(seed, dist, n = 900, scale = 0.7, cens_rate = 0.35) {
  set.seed(seed)
  d <- data.frame(x = stats::rnorm(n),
                  grp = factor(sample(c("a", "b"), n, TRUE)))
  eta <- 1.5 + 0.8 * d$x - 0.6 * (d$grp == "b")
  w <- switch(dist,
    weibull     = log(stats::rexp(n)),
    lognormal   = stats::rnorm(n),
    loglogistic = stats::rlogis(n))
  tt <- exp(eta + scale * w)
  cmax <- stats::quantile(tt, 1 - cens_rate)
  ct <- pmin(cmax, stats::rexp(n, rate = 1 / (3 * cmax)))
  d$time <- pmin(tt, ct)
  d$event <- as.integer(tt <= ct)
  d
}

## ---- against survreg -------------------------------------------------------

test_that("every AFT family matches survreg", {
  skip_if_not_installed("survival")
  for (dist in c("weibull", "lognormal", "loglogistic")) {
    d <- surv_gen(7, dist)
    f <- ilm_model(time ~ x + grp, data = d, family = dist,
                   censor = ilm_surv(d$time, d$event), verbose = FALSE)
    sv <- survival::survreg(survival::Surv(time, event) ~ x + grp, data = d,
                            dist = dist)
    expect_equal(unname(coef(f)), unname(coef(sv)), tolerance = 1e-4,
                 label = dist)
    expect_equal(unname(f$dispersion[[1]]), sv$scale, tolerance = 1e-4,
                 label = dist)
    expect_equal(as.numeric(logLik(f)), as.numeric(logLik(sv)),
                 tolerance = 1e-6, label = dist)
    expect_equal(unname(sqrt(diag(vcov(f)))[["x"]]),
                 unname(sqrt(diag(vcov(sv)))[["x"]]), tolerance = 1e-3,
                 label = dist)
  }
})

test_that("the three families are genuinely different models", {
  skip_if_not_installed("survival")
  d <- surv_gen(7, "weibull")
  ll <- vapply(c("weibull", "lognormal", "loglogistic"), function(dist)
    as.numeric(logLik(ilm_model(time ~ x + grp, data = d, family = dist,
                                censor = ilm_surv(d$time, d$event),
                                verbose = FALSE))), 0)
  # the data were generated Weibull, so that one should win
  expect_equal(names(which.max(ll)), "weibull")
  expect_gt(diff(range(ll)), 5)
})

test_that("a frailty model recovers what survreg has no answer for", {
  set.seed(21)
  ng <- 80; nt <- 8
  d <- data.frame(id = factor(rep(seq_len(ng), each = nt)),
                  x = stats::rnorm(ng * nt))
  b <- stats::rnorm(ng, 0, 0.5)[as.integer(d$id)]
  tt <- exp(1.2 + 0.7 * d$x + b + 0.6 * log(stats::rexp(ng * nt)))
  ct <- stats::rexp(ng * nt, rate = 1 / (2 * stats::median(tt)))
  d$time <- pmin(tt, ct); d$event <- as.integer(tt <= ct)
  f <- ilm_model(time ~ x + (1 | id), data = d, family = "weibull",
                 censor = ilm_surv(d$time, d$event), verbose = FALSE)
  expect_equal(unname(coef(f)[["x"]]), 0.7, tolerance = 0.12)
  expect_equal(unname(f$dispersion[[1]]), 0.6, tolerance = 0.12)
  expect_equal(sqrt(f$Sigma[["id"]][1, 1]), 0.5, tolerance = 0.2)
  # ignoring the unit pushes the between-subject spread into the scale
  g <- ilm_model(time ~ x, data = d, family = "weibull",
                 censor = ilm_surv(d$time, d$event), verbose = FALSE)
  expect_gt(g$dispersion[[1]], f$dispersion[[1]])
})

## ---- residuals -------------------------------------------------------------

test_that("quantile residuals are uniform, censored or not", {
  for (dist in c("weibull", "lognormal", "loglogistic")) {
    for (cz in c(FALSE, TRUE)) {
      d <- surv_gen(7, dist, n = 2000, cens_rate = if (cz) 0.4 else 0)
      cs <- if (cz) ilm_surv(d$time, d$event) else NULL
      f <- ilm_model(time ~ x, data = d, family = dist, censor = cs,
                     verbose = FALSE)
      p <- suppressWarnings(stats::ks.test(ilm_rqr(f, seed = 1L), "punif")$p.value)
      expect_gt(p, 0.01, label = paste(dist, "censored =", cz))
    }
  }
})

test_that("the AFT residual is centred and unit-scaled even when censored", {
  # the obvious alternative, (log t - eta) / scale, is evaluated at the
  # censoring time rather than the event for a censored row, so it sits low:
  # measured at 40% censoring it pulled the mean to -0.64, -0.84 and -0.56
  for (dist in c("weibull", "lognormal", "loglogistic")) {
    d <- surv_gen(7, dist, n = 2000, cens_rate = 0.4)
    f <- ilm_model(time ~ x, data = d, family = dist,
                   censor = ilm_surv(d$time, d$event), verbose = FALSE)
    r <- illume:::ilm_pearson_ovr(f)
    expect_lt(abs(mean(r)), 0.1, label = paste(dist, "mean"))
    expect_lt(abs(stats::sd(r) - 1), 0.15, label = paste(dist, "sd"))
  }
})

test_that("the residual machinery reaches an AFT fit", {
  d <- surv_gen(7, "weibull", n = 600)
  f <- ilm_model(time ~ x, data = d, family = "weibull",
                 censor = ilm_surv(d$time, d$event), verbose = FALSE)
  expect_silent(ilm_coef_table(f))
  expect_true(is.finite(AIC(f)))
  expect_true(all(is.finite(ilm_simulate(f, 3, seed = 1))))
  ff <- tempfile(fileext = ".png")
  grDevices::png(ff, width = 1000, height = 650)
  on.exit({grDevices::dev.off(); unlink(ff)}, add = TRUE)
  suppressWarnings(suppressMessages(ilm_appraise(f, B = 60L)))
  grDevices::dev.off()
  on.exit(unlink(ff), add = FALSE)
  expect_gt(file.size(ff), 5000)
})

## ---- the survival wrapper --------------------------------------------------

test_that("ilm_surv flips the convention so the user does not have to", {
  s <- ilm_surv(c(5, 9, 12, 3, 20), c(1, 0, 1, 1, 0))
  # survival::Surv uses 1 for the event; ilm_censor stores 1 for censored
  expect_equal(as.integer(s), c(0L, 1L, 0L, 0L, 1L))
  expect_true(is.na(attr(s, "lower")))
  expect_equal(as.integer(ilm_surv(1:4, c(TRUE, FALSE, TRUE, TRUE))),
               c(0L, 1L, 0L, 0L))
})

test_that("misspecified follow-up is named", {
  expect_error(ilm_surv(1:5), "`event` is required")
  expect_error(ilm_surv(1:5, c(1, 0)), "same length")
  expect_error(ilm_surv(1:5, rep(2, 5)), "must be 1")
  expect_error(ilm_surv(c(0, 1, 2, 3, 4), c(1, 1, 1, 0, 1)), "strictly positive")
  expect_error(ilm_surv(1:5, rep(0, 5)), "no subject had the event")
  expect_warning(ilm_surv(1:100, c(1, rep(0, 99))), "not the other way round")
})

test_that("a non-positive time is refused with the reason", {
  d <- surv_gen(7, "weibull", n = 200)
  d$time[1] <- 0
  expect_error(ilm_model(time ~ x, data = d, family = "weibull",
                         verbose = FALSE),
               "must be strictly positive")
  expect_error(ilm_model(time ~ x, data = d, family = "weibull",
                         verbose = FALSE),
               "before the first assessment")
  # the gaussian family has no such restriction
  expect_silent(ilm_model(x ~ time, data = d, family = "gaussian",
                          verbose = FALSE))
})

test_that("the coefficient reads as a log time ratio", {
  # a doubling of survival time is a coefficient of log(2)
  set.seed(31)
  n <- 4000
  d <- data.frame(x = rep(0:1, each = n / 2))
  d$time <- exp(1 + log(2) * d$x + 0.5 * log(stats::rexp(n)))
  d$event <- 1L
  f <- ilm_model(time ~ x, data = d, family = "weibull", verbose = FALSE)
  expect_equal(unname(coef(f)[["x"]]), log(2), tolerance = 0.06)
  # and the ratio of median times matches
  md <- tapply(d$time, d$x, stats::median)
  expect_equal(unname(md[2] / md[1]), 2, tolerance = 0.15)
})
