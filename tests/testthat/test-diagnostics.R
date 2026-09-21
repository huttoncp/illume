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
  # B is kept small for speed; the envelope is then unstable and says so, which
  # is not what this test is about
  cal <- suppressWarnings(ilm_calibration(fit, nbins = 4L, B = 20L))
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
  # as above: a deliberately cheap B draws a legitimate warning from the
  # calibration panel, so the silence being checked is at the default
  expect_silent(suppressWarnings(invisible(ilm_appraise(fit, nbins = 4L, B = 20L))))
})

test_that("appraise passes a small B through to the calibration warning", {
  fit <- fit_basic()
  pdf(NULL); on.exit(dev.off())
  expect_warning(invisible(ilm_appraise(fit, nbins = 4L, B = 20L)),
                 "unstable envelope")
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

## --- missing data and the row alignment it forces ---------------------------
##
## A fit drops rows with missing values, so its residuals are shorter than the
## data frame a covariate is passed from. Both of these used to fail with
## "arguments must have same length", which names no remedy and does not say
## which two lengths disagreed. Found by calling illume from outside it -- not
## by this suite, which is the usual story.

na_fit <- function(seed = 4, n = 400, n_na = 9) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n), z = rnorm(n), w = rnorm(n))
  d$y <- 0.5 * d$x + 1.2 * d$w + rnorm(n, 0, 0.6)
  d$fac <- factor(sample(c("p", "q", "r"), n, TRUE))
  d$x[seq(3, by = 7, length.out = n_na)] <- NA
  list(d = d, fit = ilm_model(y ~ x, data = d, family = "gaussian",
                              verbose = FALSE))
}

test_that("a covariate check works when the fit dropped rows", {
  z <- na_fit()
  expect_s3_class(z$fit, "ilm_model")
  expect_lt(nrow(z$fit$model), nrow(z$d))          # rows really were dropped
  expect_silent(r <- ilm_check_covariate(z$fit, z$d$w, name = "w", B = 20L,
                                         verbose = FALSE))
  expect_true(r$status %in% c("OK", "WARN", "FAIL"))
})

test_that("aligning automatically gives the SAME answer as subsetting by hand", {
  ## The point of the test: a wrong alignment also runs without erroring, so
  ## "it no longer crashes" is not evidence. Dropping the last k rows instead
  ## of the right k gives z = 0.41 where the truth is z = 23.5.
  z <- na_fit()
  om <- as.integer(attr(z$fit$model, "na.action"))
  auto <- ilm_check_covariate(z$fit, z$d$w, name = "w", B = 40L, seed = 7L,
                              verbose = FALSE)
  hand <- ilm_check_covariate(z$fit, z$d$w[-om], name = "w", B = 40L, seed = 7L,
                              verbose = FALSE)
  expect_equal(auto$observed, hand$observed)
  expect_equal(auto$z, hand$z)
  wrong <- ilm_check_covariate(z$fit,
             z$d$w[seq_len(nrow(z$d) - length(om))],
             name = "w", B = 40L, seed = 7L, verbose = FALSE)
  expect_false(isTRUE(all.equal(auto$observed, wrong$observed)))
})

test_that("the check keeps its power when rows have been dropped", {
  z <- na_fit()
  hit <- ilm_check_covariate(z$fit, z$d$w, name = "w", B = 40L, seed = 7L,
                             verbose = FALSE)
  nul <- ilm_check_covariate(z$fit, z$d$z, name = "z", B = 40L, seed = 7L,
                             verbose = FALSE)
  expect_identical(hit$status, "FAIL")             # w really is omitted
  expect_gt(abs(hit$z), 5)
  expect_identical(nul$status, "OK")               # z really is not
})

test_that("ilm_check_omitted takes a frame longer than the model frame", {
  z <- na_fit()
  expect_silent(o <- ilm_check_omitted(z$fit, z$d[c("w", "z", "fac")],
                                       B = 20L, verbose = FALSE))
  expect_named(o, c("w", "z", "fac"))
  expect_identical(o$w$status, "FAIL")
})

test_that("a covariate with missing values of its own is handled", {
  ## it is a variable the model did NOT use, so nothing has filtered it
  z <- na_fit()
  z$d$wna <- z$d$w
  z$d$wna[c(1, 2, 5, 400)] <- NA
  expect_silent(r <- ilm_check_covariate(z$fit, z$d$wna, name = "wna",
                                         B = 20L, verbose = FALSE))
  expect_true(r$status %in% c("OK", "WARN", "FAIL"))
})

test_that("a length that matches nothing is refused, with the remedy named", {
  z <- na_fit()
  expect_error(ilm_check_covariate(z$fit, rnorm(13), name = "junk",
                                   verbose = FALSE), "13 value")
  expect_error(ilm_check_covariate(z$fit, rnorm(13), name = "junk",
                                   verbose = FALSE), "same data frame")
})
