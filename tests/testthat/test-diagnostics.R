test_that("quantile residuals lie on the unit interval", {
  fit <- fit_basic()
  u <- ilm_rqr(fit)
  expect_length(u, nrow(fit$X))
  expect_true(all(u >= 0 & u <= 1))
  expect_equal(u, ilm_rqr(fit, seed = 1))     # reproducible given the seed
})

test_that("the residual construction is uniform when probabilities are known", {
  # the oracle case: no estimation anywhere, so the transform must be exact
  set.seed(99)
  N <- 4000; J <- 3
  P <- matrix(runif(N * J, 0.2, 1), N, J); P <- P / rowSums(P)
  cp <- t(apply(P, 1, cumsum))
  y <- as.integer(rowSums(runif(N) > cp)) + 1L
  S <- -log(P); obs <- S[cbind(seq_len(N), y)]
  lo <- rowSums(P * (S < obs - 1e-10)); hi <- rowSums(P * (S <= obs + 1e-10))
  u <- lo + runif(N) * (hi - lo)
  expect_equal(mean(u), 0.5, tolerance = 0.03)
  expect_gt(suppressWarnings(ks.test(u, "punif")$p.value), 0.01)
})

test_that("scoring rules are on their expected scales", {
  fit <- fit_basic()
  s <- ilm_scores(fit)
  expect_named(s, c("log_score", "brier", "accuracy"))
  expect_gt(s[["log_score"]], 0)
  expect_true(s[["brier"]] >= 0 && s[["brier"]] <= 2)
  expect_true(s[["accuracy"]] >= 0 && s[["accuracy"]] <= 1)
})

test_that("calibration returns one entry per category", {
  fit <- fit_basic()
  cal <- ilm_calibration(fit, nbins = 4L, B = 20L)
  expect_length(cal, fit$J)
  ok <- Filter(Negate(is.null), cal)
  expect_true(all(vapply(ok, function(z) all(z$lo <= z$hi, na.rm = TRUE), TRUE)))
})

test_that("random-effect distances are non-negative with the right df", {
  fit <- fit_basic()
  rm_ <- ilm_re_mahalanobis(fit)
  expect_true(all(rm_$d >= 0))
  expect_equal(rm_$df, fit$wk[[1]])
  expect_equal(length(rm_$d), fit$nlk[[1]])
})

test_that("simulation produces valid category codes and is reproducible", {
  fit <- fit_basic()
  ys <- ilm_simulate(fit, 3L, seed = 7L)
  expect_equal(dim(ys), c(nrow(fit$X), 3L))
  expect_true(all(ys >= 1L & ys <= fit$J))
  expect_equal(ys, ilm_simulate(fit, 3L, seed = 7L))
})

test_that("appraise draws without error", {
  fit <- fit_basic()
  pdf(NULL); on.exit(dev.off())
  expect_silent(invisible(ilm_appraise(fit, nbins = 4L, B = 20L)))
})

test_that("a targeted covariate check runs and reports a status", {
  fit <- fit_basic()
  r <- ilm_check_covariate(fit, fit$model$grp, name = "grp", B = 12L,
                            verbose = FALSE)
  expect_true(r$status %in% c("OK", "WARN", "FAIL", "INCONCLUSIVE"))
})

test_that("autocorrelation helpers behave on known input", {
  z <- c(1, 2, 3, 4, 1, 2, 3, 4)
  g <- c(1, 1, 1, 1, 2, 2, 2, 2)
  tt <- c(1, 2, 3, 4, 1, 2, 3, 4)
  a <- ilm_resid_acf(z, g, tt, maxlag = 2L, min_pairs = 2L)
  expect_length(a, 2L)
  expect_gt(a[1], 0.9)     # perfectly increasing within group
})

test_that("one-vs-rest Pearson residuals have the right shape", {
  fit <- fit_basic()
  R <- ilm_pearson_ovr(fit)
  expect_equal(dim(R), c(nrow(fit$X), fit$J))
  expect_true(all(is.finite(R)))
})
