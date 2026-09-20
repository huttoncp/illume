make_panel <- function(nunit = 40L, nper = 8L, tstart = 5L, att = 0.8,
                       pre_trend = 0, seed = 1L) {
  set.seed(seed)
  d <- expand.grid(unit = seq_len(nunit), time = seq_len(nper))
  d$treated <- as.integer(d$unit <= nunit / 2)
  d$post <- as.integer(d$time >= tstart)
  d$y <- 1 + 0.3 * d$time + pre_trend * d$treated * d$time +
    rnorm(nunit)[d$unit] + att * d$treated * d$post + rnorm(nrow(d))
  d
}

test_that("the ATT is the interaction and recovers the truth", {
  d <- make_panel()
  f <- ilm_did(d, "y", "unit", "time", treated = "treated", post = "post",
               verbose = FALSE)
  expect_s3_class(f, "ilm_did")
  expect_equal(f$n_treated, 20L); expect_equal(f$n_control, 20L)
  expect_equal(f$treat_time, 5)
  expect_equal(f$family, "gaussian")
  expect_lt(f$att$lower, 0.8); expect_gt(f$att$upper, 0.8)
  ## the fit is an ordinary ilm_model, so the rest of the package works on it
  expect_s3_class(f$fit, "ilm_model")
  expect_true(is.finite(stats::AIC(f$fit)))

  ## treat_time is an alternative to a post column
  f2 <- ilm_did(d, "y", "unit", "time", treated = "treated", treat_time = 5,
                verbose = FALSE)
  expect_equal(f2$att$estimate, f$att$estimate)
})

test_that("parallel trends passes when true and fails when violated", {
  ok <- ilm_did(make_panel(pre_trend = 0), "y", "unit", "time",
                treated = "treated", post = "post", verbose = FALSE)
  expect_equal(ok$parallel$status, "OK")
  expect_equal(ok$parallel$n_pre, 4L)

  ## the treated group already climbing before anything happened
  bad <- ilm_did(make_panel(pre_trend = 0.35, seed = 2L), "y", "unit", "time",
                 treated = "treated", post = "post", verbose = FALSE)
  expect_equal(bad$parallel$status, "FAIL")
  expect_gt(abs(bad$parallel$diff_slope), 0.1)
  expect_output(print(bad), "already diverging")

  ## with one pre-period there is nothing to test, and that is said
  one <- ilm_did(make_panel(nper = 5L, tstart = 2L), "y", "unit", "time",
                 treated = "treated", post = "post", verbose = FALSE)
  expect_equal(one$parallel$status, "UNTESTED")
})

test_that("the event study brackets treatment and uses a pre-period reference", {
  f <- ilm_did(make_panel(), "y", "unit", "time", treated = "treated",
               post = "post", verbose = FALSE)
  ev <- f$event
  expect_false(is.null(ev))
  expect_true(any(ev$rel_time < 0)); expect_true(any(ev$rel_time >= 0))
  ## exactly one reference period, pinned at zero with no uncertainty
  ref <- ev[ev$se == 0, ]
  expect_equal(nrow(ref), 1L)
  expect_equal(ref$estimate, 0)
  expect_lt(ref$rel_time, 0)
  expect_true(all(diff(ev$rel_time) > 0))
})

test_that("staggered adoption is refused unless asked for", {
  set.seed(3)
  d <- expand.grid(unit = 1:40, time = 1:8)
  start <- ifelse(d$unit <= 10, 4, ifelse(d$unit <= 20, 6, Inf))
  d$treatment <- as.integer(d$time >= start)
  d$y <- 1 + 0.3 * d$time + rnorm(40)[d$unit] + 0.8 * d$treatment + rnorm(nrow(d))
  expect_error(ilm_did(d, "y", "unit", "time", treatment = "treatment",
                       verbose = FALSE), "staggered")
  f <- ilm_did(d, "y", "unit", "time", treatment = "treatment",
               allow_staggered = TRUE, verbose = FALSE)
  expect_true(f$staggered)
  expect_output(print(f), "CAVEAT")

  ## common timing given the same way is NOT flagged
  d2 <- d; d2$treatment <- as.integer(d2$unit <= 20 & d2$time >= 5)
  expect_false(ilm_did(d2, "y", "unit", "time", treatment = "treatment",
                       verbose = FALSE)$staggered)
})

test_that("did takes other families and other treatment codings", {
  set.seed(6)
  d <- expand.grid(unit = 1:40, time = 1:8)
  d$treated <- d$unit <= 20            # logical
  d$post <- factor(ifelse(d$time >= 5, "after", "before"),
                   levels = c("before", "after"))
  eta <- -0.2 + 0.1 * d$time + rnorm(40)[d$unit] * 0.4 +
    0.9 * as.integer(d$treated) * as.integer(d$post == "after")
  d$y <- rbinom(nrow(d), 1, plogis(eta))
  f <- ilm_did(d, "y", "unit", "time", treated = "treated", post = "post",
               verbose = FALSE)
  expect_equal(f$family, "binomial")
  expect_true(is.finite(f$att$estimate))
})

test_that("did refuses what it cannot do", {
  d <- make_panel()
  expect_error(ilm_did(d, "y", "unit", "time", verbose = FALSE),
               "give either")
  expect_error(ilm_did(d, "y", "unit", "time", treated = "treated",
                       verbose = FALSE), "give `post` or `treat_time`")
  expect_error(ilm_did(d, "nope", "unit", "time", treated = "treated",
                       post = "post", verbose = FALSE), "not found")
  expect_error(ilm_did("nope", "y", "unit", "time", verbose = FALSE),
               "must be a data frame")
  ## a group marker with only one value cannot mark two groups
  d2 <- d; d2$treated <- 1L
  expect_error(ilm_did(d2, "y", "unit", "time", treated = "treated",
                       post = "post", verbose = FALSE), "1 distinct value")
})

test_that("the event-study plot draws and needs an event study", {
  f <- ilm_did(make_panel(), "y", "unit", "time", treated = "treated",
               post = "post", verbose = FALSE)
  pf <- file.path(tempdir(), "did.png")
  grDevices::png(pf, width = 600, height = 400); on.exit(unlink(pf))
  expect_silent(ilm_plot_did(f))
  grDevices::dev.off()
  expect_true(file.exists(pf))
  expect_error(ilm_plot_did(structure(list(), class = "ilm_did")),
               "no event study")
  expect_error(ilm_plot_did(1), "must be an ilm_did")
})
