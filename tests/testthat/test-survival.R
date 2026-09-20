# Survival curves, and the check that the family was the right choice.
#
# A parametric survival model buys a smooth curve and interpretable
# coefficients. What it costs is an assumption about the SHAPE of the baseline.
# The Kaplan-Meier makes no such assumption, so laying the fitted curve over it
# is the check.
#
# The load-bearing test here is that the verdict PASSES a correctly specified
# model. An earlier version failed all four cases it was shown, including the
# right one, because simulated replicates were not censored the way the study
# was: 141 of 190 censored subjects drew a time beyond their own censoring
# time, and the simulated Kaplan-Meier sat above the observed one everywhere.

surv_d <- function(seed, dist = "weibull", n = 500, mult = 2.5) {
  set.seed(seed)
  d <- data.frame(x = stats::rnorm(n))
  w <- switch(dist, weibull = log(stats::rexp(n)),
              lognormal = stats::rnorm(n), loglogistic = stats::rlogis(n))
  tt <- exp(1.5 + 0.8 * d$x + 0.7 * w)
  ct <- stats::rexp(n, rate = 1 / (mult * stats::median(tt)))
  d$time <- pmin(tt, ct); d$event <- as.integer(tt <= ct)
  d
}

fit_surv <- function(d, dist = "weibull")
  ilm_model(time ~ x, data = d, family = dist,
            censor = ilm_surv(d$time, d$event), verbose = FALSE)

## ---- the Kaplan-Meier ------------------------------------------------------

test_that("the product-limit estimator matches survival::survfit", {
  skip_if_not_installed("survival")
  set.seed(5)
  n <- 300
  tt <- stats::rexp(n, 0.1); ct <- stats::rexp(n, 0.05)
  tm <- pmin(tt, ct); ev <- as.integer(tt <= ct)
  mine <- illume:::ilm_km(tm, ev)
  theirs <- survival::survfit(survival::Surv(tm, ev) ~ 1)
  ix <- match(mine$time, theirs$time)
  expect_equal(mine$surv, theirs$surv[ix], tolerance = 1e-12)
  expect_gt(length(mine$time), 100L)
})

test_that("the estimator handles the edges", {
  # nobody had the event
  e <- illume:::ilm_km(1:5, rep(0L, 5))
  expect_length(e$time, 0L)
  expect_equal(illume:::ilm_km_at(e, c(1, 3)), c(NA_real_, NA_real_))
  # everybody did, so the curve reaches zero and Greenwood's last term is
  # undefined rather than infinite
  k <- illume:::ilm_km(1:5, rep(1L, 5))
  expect_equal(k$surv[5], 0)
  expect_true(all(is.finite(k$se)))
  # before the first event the estimate is 1
  expect_equal(illume:::ilm_km_at(k, 0.5), 1)
})

## ---- the fitted curve ------------------------------------------------------

test_that("predicted survival agrees with survreg's quantiles", {
  skip_if_not_installed("survival")
  d <- surv_d(9)
  f <- fit_surv(d)
  sv <- survival::survreg(survival::Surv(time, event) ~ x, d, dist = "weibull")
  med_sv <- stats::predict(sv, newdata = data.frame(x = c(-1, 0, 1)),
                           type = "quantile", p = 0.5)
  med_me <- vapply(c(-1, 0, 1), function(xx) {
    s <- ilm_survival(f, data.frame(x = xx), times = seq(0.05, 200, by = 0.05))
    s$time[which.min(abs(s$surv - 0.5))]
  }, 0)
  expect_equal(unname(med_me), unname(med_sv), tolerance = 0.02)
})

test_that("the curve decreases, stays in range and widens with uncertainty", {
  d <- surv_d(9)
  f <- fit_surv(d)
  s <- ilm_survival(f, newdata = data.frame(x = c(-1, 1)),
                    times = c(1, 2, 5, 10, 20))
  expect_equal(nrow(s), 10L)
  expect_true(all(s$surv >= 0 & s$surv <= 1))
  expect_true(all(s$lower <= s$surv + 1e-9))
  expect_true(all(s$upper >= s$surv - 1e-9))
  # survival falls with time within a covariate pattern
  for (r in unique(s$row)) expect_true(all(diff(s$surv[s$row == r]) <= 0))
  # a larger x means longer survival, since the coefficient is positive
  expect_true(all(s$surv[s$row == 2] > s$surv[s$row == 1]))
})

test_that("a survival curve needs a survival family", {
  set.seed(2)
  d <- data.frame(x = stats::rnorm(200))
  d$y <- stats::rnorm(200)
  g <- ilm_model(y ~ x, data = d, family = "gaussian", verbose = FALSE)
  # the message now names the flexible baseline too, since that is also an
  # answer to "how do I get a survival curve"
  expect_error(ilm_survival(g), "needs a time-to-event family")
  expect_error(ilm_survival(g), "rp", fixed = TRUE)
  expect_error(ilm_plot_survival(g), "needs a time-to-event family")
  expect_error(ilm_survival(d), "must be a fitted ilm_model")
})

## ---- the verdict -----------------------------------------------------------

draws_to <- function(expr, w = 700, h = 450) {
  ff <- tempfile(fileext = ".png")
  grDevices::png(ff, width = w, height = h)
  r <- force(expr)
  grDevices::dev.off()
  list(size = file.size(ff), value = r)
}

test_that("a correctly specified model passes", {
  d <- surv_d(9, "weibull")
  f <- fit_surv(d, "weibull")
  r <- draws_to(ilm_plot_survival(f, d$time, d$event, B = 50, verbose = FALSE))
  expect_equal(r$value$status, "OK")
  expect_lt(r$value$outside, 0.05)
  expect_gt(r$size, 4000)
  # the raw gap alone would have condemned it: the Kaplan-Meier of 500
  # subjects wobbles by more than a tenth in the tail, which is exactly why
  # the envelope exists rather than a threshold on the gap
  expect_gt(r$value$gap, 0.05)
})

test_that("the wrong family is caught", {
  set.seed(3)
  n <- 500
  tt <- exp(2 + 0.9 * stats::rlogis(n) * 1.4)
  ct <- stats::rexp(n, rate = 1 / (3 * stats::median(tt)))
  d <- data.frame(time = pmin(tt, ct), event = as.integer(tt <= ct),
                  x = stats::rnorm(n))
  fw <- fit_surv(d, "weibull")
  fl <- fit_surv(d, "loglogistic")
  rw <- draws_to(ilm_plot_survival(fw, d$time, d$event, B = 50, verbose = FALSE))
  rl <- draws_to(ilm_plot_survival(fl, d$time, d$event, B = 50, verbose = FALSE))
  expect_equal(rw$value$status, "FAIL")
  expect_equal(rl$value$status, "OK")
  # and the likelihood agrees with the picture
  expect_lt(AIC(fl), AIC(fw))
})

test_that("replicates are censored the way the study was", {
  # this is what makes the verdict mean anything
  d <- surv_d(9)
  cs <- ilm_surv(d$time, d$event)
  f <- ilm_model(time ~ x, data = d, family = "weibull", censor = cs,
                 verbose = FALSE)
  # the observed codes are recovered from the observed response
  expect_identical(as.integer(illume:::ilm_censor_for(cs, d$time)),
                   as.integer(1L - d$event))
  ys <- illume:::ilm_sim_cond(f, 20, seed = 1)
  rate <- apply(ys, 2, function(v) mean(illume:::ilm_censor_for(cs, v) != 0L))
  # near the observed rate, and varying, because whether a draw outlives its
  # censoring time is itself random
  expect_equal(mean(rate), mean(d$event == 0), tolerance = 0.08)
  expect_gt(stats::sd(rate), 0)
  # no replicate outlives the follow-up by a wide margin
  expect_lt(max(ys), max(d$time) * 1.5)
})

test_that("grouped curves are drawn and validated per level", {
  d <- surv_d(11)
  d$grp <- factor(rep(c("a", "b"), length.out = nrow(d)))
  d$time <- d$time * ifelse(d$grp == "b", 2, 1)
  f <- ilm_model(time ~ x + grp, data = d, family = "weibull",
                 censor = ilm_surv(d$time, d$event), verbose = FALSE)
  r <- draws_to(ilm_plot_survival(f, d$time, d$event, by = d$grp, B = 30,
                                  verbose = FALSE))
  expect_equal(r$value$levels, c("a", "b"))
  expect_equal(ncol(r$value$km), 2L)
  expect_equal(ncol(r$value$fitted), 2L)
  expect_gt(r$size, 4000)
})

test_that("the event indicator is recovered from the model when not given", {
  d <- surv_d(9)
  f <- fit_surv(d)
  a <- draws_to(ilm_plot_survival(f, B = 0, verbose = FALSE))
  b <- draws_to(ilm_plot_survival(f, d$time, d$event, B = 0, verbose = FALSE))
  expect_equal(a$value$km, b$value$km)
  # with no envelope there is no verdict to give
  expect_equal(a$value$status, "INCONCLUSIVE")
  expect_true(is.na(a$value$outside))
})

test_that("misspecification is named", {
  d <- surv_d(9, n = 200)
  f <- fit_surv(d)
  expect_error(ilm_plot_survival(f, d$time[-1], d$event), "but the model has")
  expect_error(ilm_plot_survival(f, d$time, d$event[-1]), "but the model has")
  expect_error(ilm_plot_survival(f, d$time, d$event, by = d$x[-1]),
               "but the model has")
  g <- ilm_model(time ~ x, data = d, family = "weibull", verbose = FALSE)
  expect_error(ilm_plot_survival(g), "`event` is required")
})

test_that("the printed report names the remedy when the shape is wrong", {
  set.seed(3)
  n <- 400
  tt <- exp(2 + 0.9 * stats::rlogis(n) * 1.4)
  ct <- stats::rexp(n, rate = 1 / (3 * stats::median(tt)))
  d <- data.frame(time = pmin(tt, ct), event = as.integer(tt <= ct),
                  x = stats::rnorm(n))
  f <- fit_surv(d, "weibull")
  ff <- tempfile(fileext = ".png")
  grDevices::png(ff, width = 700, height = 450)
  on.exit({grDevices::dev.off(); unlink(ff)}, add = TRUE)
  out <- capture.output(ilm_plot_survival(f, d$time, d$event, B = 40,
                                          verbose = TRUE))
  expect_true(any(grepl("largest vertical gap", out)))
  expect_true(any(grepl("outside the envelope", out)))
  expect_true(any(grepl("another family", out)))
})
