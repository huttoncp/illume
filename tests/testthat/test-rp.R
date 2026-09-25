# Royston-Parmar flexible parametric survival.
#
# The baseline is a restricted cubic spline in log time instead of a straight
# line, so the shape is free while the coefficients keep their interpretation:
# a hazard ratio on the hazard scale, an odds ratio on the odds scale.
#
# There is no flexsurv or rstpm2 here to check against, so the load-bearing
# test is the identity: with rp_df = 1 the spline IS a straight line, and each
# scale must reproduce the parametric model it generalises. Those are pinned to
# survival::survreg() in test-aft.R, so the validation carries through.

rp_gen <- function(seed, dist = "weibull", n = 900, scale = 0.7, cr = 0.3) {
  set.seed(seed)
  d <- data.frame(x = stats::rnorm(n),
                  grp = factor(sample(c("a", "b"), n, TRUE)))
  eta <- 1.5 + 0.8 * d$x - 0.6 * (d$grp == "b")
  w <- switch(dist, weibull = log(stats::rexp(n)),
              lognormal = stats::rnorm(n), loglogistic = stats::rlogis(n))
  tt <- exp(eta + scale * w)
  q <- stats::quantile(tt, 1 - cr)
  ct <- pmin(q, stats::rexp(n, rate = 1 / (3 * q)))
  d$time <- pmin(tt, ct); d$event <- as.integer(tt <= ct)
  d
}

## a piecewise-constant hazard: 0.15 before t = 3, then 0.9. No parametric
## family in the package can bend like that, which is the point.
rp_piecewise <- function(seed, n = 800, b = 0.7) {
  set.seed(seed)
  d <- data.frame(x = stats::rnorm(n))
  lp <- b * d$x
  e <- -log(stats::runif(n))
  t1 <- e / (0.15 * exp(lp))
  tt <- pmax(ifelse(t1 <= 3, t1, 3 + (e - 0.45 * exp(lp)) / (0.9 * exp(lp))),
             1e-4)
  ct <- stats::rexp(n, rate = 1 / (2 * stats::median(tt)))
  d$time <- pmin(tt, ct); d$event <- as.integer(tt <= ct)
  d
}

rp_true_S <- function(t, x, b = 0.7) {
  H <- ifelse(t <= 3, 0.15 * t, 0.45 + 0.9 * (t - 3))
  exp(-H * exp(b * x))
}

## ---- the spline basis ------------------------------------------------------

test_that("the basis derivative is the derivative of the basis", {
  kn <- c(0, 1, 2, 3.5)
  x <- seq(-0.5, 4, length.out = 200)
  D <- illume:::ilm_rcs(x, kn, deriv = TRUE)
  h <- 1e-6
  num <- (illume:::ilm_rcs(x + h, kn) - illume:::ilm_rcs(x - h, kn)) / (2 * h)
  expect_equal(D, num, tolerance = 1e-6, ignore_attr = TRUE)
  expect_equal(ncol(illume:::ilm_rcs(x, kn)), length(kn) - 1L)
})

test_that("the basis is linear outside the boundary knots", {
  # which is what stops the fitted hazard doing something wild in the tails
  kn <- c(0, 1, 2, 3.5)
  x <- c(seq(-3, -0.1, length.out = 50), seq(3.6, 7, length.out = 50))
  B <- illume:::ilm_rcs(x, kn)
  for (part in list(1:50, 51:100)) {
    dd <- apply(B[part, , drop = FALSE], 2, function(z) diff(diff(z)))
    expect_lt(max(abs(dd)), 1e-8)
  }
  # with two knots there is only the linear term
  expect_equal(ncol(illume:::ilm_rcs(x, c(0, 1))), 1L)
  expect_equal(as.vector(illume:::ilm_rcs(x, c(0, 1))), x)
})

test_that("knots come from the events, not from the censored times", {
  lt <- log(c(1, 2, 3, 4, 100, 200))
  ev <- c(1L, 1L, 1L, 1L, 0L, 0L)
  # only four events here, so the thin-spline warning is expected
  k <- suppressWarnings(illume:::ilm_rp_knots(lt, ev, 3L))
  expect_length(k, 4L)
  # the censored times at 100 and 200 say when someone was last seen, not when
  # the hazard changed, so they do not stretch the boundary knots out to them
  expect_equal(range(k), range(lt[ev == 1L]))
  expect_true(all(diff(k) > 0))
  # df = 1 is the two boundaries and nothing between
  expect_equal(illume:::ilm_rp_knots(lt, ev, 1L), range(lt[ev == 1L]))
  expect_error(illume:::ilm_rp_knots(lt, c(1L, rep(0L, 5)), 3L),
               "at least two events")
  # three of the four events are tied, so there are two distinct times and a
  # four-knot spline has nowhere to put them
  expect_error(illume:::ilm_rp_knots(log(c(1, 1, 1, 2)), rep(1L, 4), 3L),
               "distinct event times")
})

## ---- the identity that carries the validation ------------------------------

test_that("rp_df = 1 is exactly the parametric model it generalises", {
  for (pr in list(c("rp", "weibull"), c("rp_odds", "loglogistic"),
                  c("rp_normal", "lognormal"))) {
    d <- rp_gen(7, pr[2])
    cs <- ilm_surv(d$time, d$event)
    fr <- ilm_model(time ~ x + grp, data = d, family = pr[1], rp_df = 1,
                    censor = cs, verbose = FALSE)
    fa <- ilm_model(time ~ x + grp, data = d, family = pr[2], censor = cs,
                    verbose = FALSE)
    expect_equal(as.numeric(logLik(fr)), as.numeric(logLik(fa)),
                 tolerance = 1e-7, label = pr[1])
    sc <- fa$dispersion[[1]]
    # the flexible parameterisation is the accelerated one divided by -scale
    expect_equal(unname(coef(fr)[["rcs1"]]), 1 / sc, tolerance = 1e-5,
                 label = pr[1])
    expect_equal(unname(coef(fr)[["x"]]), -unname(coef(fa)[["x"]]) / sc,
                 tolerance = 1e-5, label = pr[1])
    # and so the two models have the same number of parameters
    expect_equal(attr(logLik(fr), "df"), attr(logLik(fa), "df"))
  }
})

## ---- the shape it exists for -----------------------------------------------

test_that("a flexible baseline fits what no parametric family can", {
  d <- rp_piecewise(12)
  cs <- ilm_surv(d$time, d$event)
  par_aic <- vapply(c("weibull", "lognormal", "loglogistic"), function(fam)
    AIC(ilm_model(time ~ x, data = d, family = fam, censor = cs,
                  verbose = FALSE)), 0)
  f5 <- ilm_model(time ~ x, data = d, family = "rp", rp_df = 5, censor = cs,
                  verbose = FALSE)
  expect_lt(AIC(f5), min(par_aic) - 50)
  # rp_df = 1 is the Weibull, so it must land on the Weibull's AIC
  f1 <- ilm_model(time ~ x, data = d, family = "rp", rp_df = 1, censor = cs,
                  verbose = FALSE)
  expect_equal(AIC(f1), unname(par_aic[["weibull"]]), tolerance = 1e-4)
})

test_that("the hazard ratio survives a baseline the model cannot draw", {
  # the appeal of the model: the coefficient is right even where the shape is
  # hard. Measured across rp_df of 3, 5 and 7 it came back 0.711, 0.699, 0.701
  # against a truth of 0.700
  d <- rp_piecewise(12)
  cs <- ilm_surv(d$time, d$event)
  for (k in c(3L, 5L, 7L)) {
    f <- ilm_model(time ~ x, data = d, family = "rp", rp_df = k, censor = cs,
                   verbose = FALSE)
    expect_equal(unname(coef(f)[["x"]]), 0.7, tolerance = 0.06,
                 label = paste("rp_df", k))
    expect_equal(illume:::ilm_rp_monotone(f), 0)
  }
})

test_that("more degrees of freedom get closer to the truth", {
  d <- rp_piecewise(12)
  cs <- ilm_surv(d$time, d$event)
  tg <- c(1, 3, 5, 8)
  err <- vapply(c(3L, 5L, 7L), function(k) {
    f <- ilm_model(time ~ x, data = d, family = "rp", rp_df = k, censor = cs,
                   verbose = FALSE)
    s <- ilm_survival(f, newdata = data.frame(x = 0), times = tg)
    max(abs(s$surv - rp_true_S(tg, 0)))
  }, 0)
  expect_lt(err[2], err[1])          # 5 beats 3
  expect_lt(err[3], err[1])          # 7 beats 3
  expect_lt(err[3], 0.06)
})

test_that("the survival-curve check catches a baseline that is too rigid", {
  # this is the diagnostic earning its place: rp_df = 3 on a discontinuous
  # hazard came back FAIL with 82% of the curve outside the envelope, and
  # rp_df = 5 came back OK with 2%
  d <- rp_piecewise(12)
  cs <- ilm_surv(d$time, d$event)
  draw <- function(k) {
    f <- ilm_model(time ~ x, data = d, family = "rp", rp_df = k, censor = cs,
                   verbose = FALSE)
    ff <- tempfile(fileext = ".png")
    grDevices::png(ff, width = 600, height = 400)
    r <- ilm_plot_survival(f, d$time, d$event, B = 40, verbose = FALSE)
    grDevices::dev.off(); unlink(ff)
    r
  }
  expect_equal(draw(3L)$status, "FAIL")
  expect_equal(draw(5L)$status, "OK")
})

## ---- is the spline earning its place? --------------------------------------

test_that("the flexible baseline is tested against its own special case", {
  d <- rp_piecewise(12)
  f <- ilm_model(time ~ x, data = d, family = "rp", rp_df = 3,
                 censor = ilm_surv(d$time, d$event), verbose = FALSE)
  r <- ilm_rp_lrt(f, verbose = FALSE)
  expect_equal(r$status, "FLEXIBLE")
  expect_equal(r$df, 2L)
  expect_lt(r$p.value, 1e-6)
  expect_equal(r$simple, "weibull")

  # and on data that really is Weibull it says the simple model suffices
  dw <- rp_gen(7, "weibull")
  fw <- ilm_model(time ~ x + grp, data = dw, family = "rp", rp_df = 3,
                  censor = ilm_surv(dw$time, dw$event), verbose = FALSE)
  rw <- ilm_rp_lrt(fw, verbose = FALSE)
  expect_equal(rw$status, "SIMPLE")
  expect_gt(rw$p.value, 0.05)

  out <- capture.output(ilm_rp_lrt(fw, verbose = TRUE))
  expect_true(any(grepl("weibull", out)))
  expect_true(any(grepl("simpler", out)))
  expect_error(ilm_rp_lrt(
    ilm_model(time ~ x, data = dw, family = "rp", rp_df = 1,
              censor = ilm_surv(dw$time, dw$event), verbose = FALSE)),
    "nothing to compare")
})

## ---- the rest of the package copes -----------------------------------------

test_that("the residual and simulation machinery reaches a flexible fit", {
  d <- rp_piecewise(12, n = 1000)
  f <- ilm_model(time ~ x, data = d, family = "rp", rp_df = 5,
                 censor = ilm_surv(d$time, d$event), verbose = FALSE)
  u <- ilm_rqr(f, seed = 1L)
  expect_gt(suppressWarnings(stats::ks.test(u, "punif")$p.value), 0.01)
  r <- illume:::ilm_pearson_ovr(f)
  expect_lt(abs(mean(r)), 0.1)
  expect_lt(abs(stats::sd(r) - 1), 0.15)
  # simulating inverts the fitted survivor function numerically, since a
  # spline baseline has no closed-form draw
  ys <- ilm_simulate(f, 20, seed = 1)
  expect_true(all(is.finite(ys)) && all(ys > 0))
  expect_equal(stats::median(ys), stats::median(d$time), tolerance = 0.25)
  ff <- tempfile(fileext = ".png")
  grDevices::png(ff, width = 1000, height = 650)
  on.exit({grDevices::dev.off(); unlink(ff)}, add = TRUE)
  suppressWarnings(suppressMessages(ilm_appraise(f, B = 40L)))
  grDevices::dev.off()
  on.exit(unlink(ff), add = FALSE)
  expect_gt(file.size(ff), 5000)
})

test_that("the survival curve carries an interval that is not degenerate", {
  d <- rp_piecewise(12)
  f <- ilm_model(time ~ x, data = d, family = "rp", rp_df = 5,
                 censor = ilm_surv(d$time, d$event), verbose = FALSE)
  s <- ilm_survival(f, newdata = data.frame(x = c(-1, 1)),
                    times = c(1, 3, 5, 8))
  expect_equal(nrow(s), 8L)
  expect_true(all(s$lower < s$surv & s$surv < s$upper))
  expect_true(all(s$surv >= 0 & s$surv <= 1))
  for (r in unique(s$row)) expect_true(all(diff(s$surv[s$row == r]) <= 0))
  # a larger x means a larger hazard, so shorter survival
  expect_true(all(s$surv[s$row == 1] > s$surv[s$row == 2]))
})

test_that("predict says what to use instead of a fitted value", {
  d <- rp_piecewise(12, n = 300)
  f <- ilm_model(time ~ x, data = d, family = "rp", rp_df = 3,
                 censor = ilm_surv(d$time, d$event), verbose = FALSE)
  expect_error(stats::predict(f, newdata = data.frame(x = 0)),
               "no fitted value for a covariate pattern")
  expect_error(stats::predict(f, newdata = data.frame(x = 0)),
               "ilm_survival", fixed = TRUE)
})

test_that("misspecification is named", {
  d <- rp_piecewise(12, n = 300)
  d0 <- d; d0$time[1] <- 0
  expect_error(ilm_model(time ~ x, data = d0, family = "rp", verbose = FALSE),
               "strictly positive")
  expect_error(ilm_model(time ~ x, data = d, family = "rp", rp_df = 4000,
                         censor = ilm_surv(d$time, d$event), verbose = FALSE),
               "distinct event times to place them at")
  # and a spline with too little behind each knot is warned about rather than
  # silently fitted to the noise between them
  expect_warning(ilm_model(time ~ x, data = d, family = "rp", rp_df = 40,
                           censor = ilm_surv(d$time, d$event), verbose = FALSE),
                 "per degree of freedom")
  expect_error(ilm_model(time ~ x, data = d, family = "rp", rp_knots = 1,
                         censor = ilm_surv(d$time, d$event), verbose = FALSE),
               "at least two knots")
  g <- ilm_model(time ~ x, data = d, family = "weibull",
                 censor = ilm_surv(d$time, d$event), verbose = FALSE)
  expect_error(ilm_rp_lrt(g), "family = ")
})

test_that("explicit knots are honoured", {
  d <- rp_piecewise(12, n = 400)
  kn <- quantile(log(d$time[d$event == 1L]), c(0, 0.5, 1))
  f <- ilm_model(time ~ x, data = d, family = "rp", rp_knots = kn,
                 censor = ilm_surv(d$time, d$event), verbose = FALSE)
  expect_equal(unname(f$rp$knots), unname(sort(kn)))
  expect_equal(f$rp$df, 2L)
  expect_true("rcs2" %in% names(coef(f)))
})
