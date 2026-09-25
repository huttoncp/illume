test_that("a single core returns no pool and runs sequentially", {
  expect_null(ilm_pool(1L))
  expect_equal(ilm_lapply(NULL, 1:3, function(i) i * 2), list(2, 4, 6))
})

test_that("consistency checking runs and reports a refit rate", {
  fit <- fit_basic()
  cc <- ilm_consistency(fit, B = 6L, verbose = FALSE)
  expect_true(all(c("table", "n_ok", "conv_rate") %in% names(cc)))
  expect_true(cc$conv_rate >= 0 && cc$conv_rate <= 1)
  expect_true(all(cc$table$status %in%
                    c("OK", "WARN", "FAIL", "BOUNDARY", "INCONCLUSIVE")))
})

test_that("a refit held at a covariance boundary counts as a refit", {
  ## pdHess is FALSE by design when a direction is held, and the fixed
  ## effects are usable; ilm_consistency() used to drop such refits
  set.seed(3); n <- 240
  d <- data.frame(g = factor(rep(1:24, each = 10)), x = rnorm(n))
  u <- rnorm(24, sd = 0.8)[d$g]
  eta <- cbind(0, 0.5 * d$x + u, u)
  d$k <- factor(apply(exp(eta) / rowSums(exp(eta)), 1, function(p)
    sample(c("a", "b", "c"), 1, prob = p)))
  f <- suppressMessages(ilm_model(k ~ x + (1 | g), data = d,
                                  family = "multinomial", verbose = FALSE))
  expect_identical(f$hessian_how, "boundary")
  expect_false(isTRUE(f$sdr$pdHess))
  ## refitted to its own response it lands on the same boundary, and counts
  r <- ilm_refit_many(f, matrix(f$y))
  expect_false(is.null(r[[1]]))
  b <- grep("^beta", names(r[[1]]))
  expect_equal(r[[1]][b], ilm_summarise_fit(f)[b], tolerance = 1e-3)
})

test_that("a lone variance fitted at zero is a boundary, not a bias", {
  set.seed(11); n <- 300
  d <- data.frame(g = factor(rep(1:30, each = 10)), x = rnorm(n))
  d$y <- 1 + 0.5 * d$x + rnorm(30, sd = 0.1)[d$g] + rnorm(n)
  f <- suppressMessages(ilm_model(y ~ x + (1 | g), data = d,
                                  family = "gaussian", verbose = FALSE))
  expect_lt(ilm_summarise_fit(f)[["sd_g1"]], 1e-3)
  cc <- ilm_consistency(f, B = 12L, verbose = FALSE)
  expect_identical(cc$table$status[cc$table$parameter == "sd_g1"], "BOUNDARY")
})

test_that("get_coef and set_coef round-trip, and theta actually propagates", {
  # if set_coef only stored the vector, nudging a covariance parameter would
  # change nothing and marginaleffects would report zero uncertainty from it
  fit <- fit_basic()
  b0 <- get_coef.ilm_model(fit)
  expect_identical(b0, coef(fit, full = TRUE))
  p0 <- get_predict.ilm_model(fit, groups = "typical")$estimate
  expect_equal(get_predict.ilm_model(set_coef.ilm_model(fit, b0), groups = "typical")$estimate, p0)

  bp <- b0; i <- grep("^subj:", names(bp))[1]; bp[i] <- bp[i] + 0.3
  pp <- get_predict.ilm_model(set_coef.ilm_model(fit, bp), groups = "population", ndraw = 40)$estimate
  p0m <- get_predict.ilm_model(fit, groups = "population", ndraw = 40)$estimate
  expect_gt(max(abs(pp - p0m)), 1e-6)
})

test_that("get_predict returns long format with rowid and group", {
  fit <- fit_basic()
  g <- get_predict.ilm_model(fit, groups = "typical")
  expect_true(all(c("rowid", "group", "estimate") %in% names(g)))
  expect_equal(nrow(g), nrow(fit$X) * fit$J)
  expect_setequal(unique(g$group), fit$ylevels)
})

test_that("model_performance returns the indices it promises", {
  fit <- fit_basic()
  mp <- model_performance.ilm_model(fit, verbose = FALSE)
  expect_true(all(c("AIC", "BIC", "logLik", "df", "Log_score", "Brier",
                    "Accuracy") %in% names(mp)))
  expect_false("R2_conditional" %in% names(mp))   # deliberately absent
})

test_that("parametric bootstrap runs and respects its p-value floor", {
  fit <- fit_basic()
  r <- ilm_pb_lrt(fit, "x1", B = 8L, verbose = FALSE)
  expect_true(r$p_boot >= 1 / (r$n_ok + 1) - 1e-12)
  expect_true(r$p_boot <= 1)
  expect_equal(r$df, fit$C)     # x1 is one column
})
