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

test_that("get_coef and set_coef round-trip, and theta actually propagates", {
  # if set_coef only stored the vector, nudging a covariance parameter would
  # change nothing and marginaleffects would report zero uncertainty from it
  fit <- fit_basic()
  b0 <- get_coef.ilm_model(fit)
  expect_identical(b0, coef(fit, full = TRUE))
  p0 <- get_predict.ilm_model(fit, marginal = FALSE)$estimate
  expect_equal(get_predict.ilm_model(set_coef.ilm_model(fit, b0), marginal = FALSE)$estimate, p0)

  bp <- b0; i <- grep("^subj:", names(bp))[1]; bp[i] <- bp[i] + 0.3
  pp <- get_predict.ilm_model(set_coef.ilm_model(fit, bp), marginal = TRUE, ndraw = 40)$estimate
  p0m <- get_predict.ilm_model(fit, marginal = TRUE, ndraw = 40)$estimate
  expect_gt(max(abs(pp - p0m)), 1e-6)
})

test_that("get_predict returns long format with rowid and group", {
  fit <- fit_basic()
  g <- get_predict.ilm_model(fit, marginal = FALSE)
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
