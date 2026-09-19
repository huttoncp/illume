# Models with no random or smooth terms.  A gaussian model of this kind is an
# ordinary linear model, and the whole point of the special case is that illume
# should then agree with lm() exactly rather than approximately.

sim_lm <- function(seed = 1, n = 200) {
  set.seed(seed)
  dd <- data.frame(x = stats::rnorm(n),
                   z = factor(sample(c("a", "b", "c"), n, TRUE)))
  dd$y <- 2 + 1.5 * dd$x - 0.8 * (dd$z == "b") + stats::rnorm(n, 0, 1.2)
  dd
}

test_that("a fixed-effects-only model fits at all", {
  dd <- sim_lm()
  f <- lum_model(y ~ x + z, data = dd, family = "gaussian", verbose = FALSE)
  expect_s3_class(f, "lum_model")
  expect_equal(f$opt$convergence, 0L)
  expect_equal(length(f$re), 0L)
  expect_true(isTRUE(f$exact_df))
  expect_equal(f$resid_df, nrow(dd) - ncol(f$X))
})

test_that("an intercept-only model fits", {
  dd <- sim_lm(2)
  f <- lum_model(y ~ 1, data = dd, family = "gaussian", verbose = FALSE)
  expect_equal(length(coef(f)), 1L)
  expect_equal(unname(coef(f)), unname(coef(stats::lm(y ~ 1, dd))), tolerance = 1e-5)
})

test_that("anova on an intercept-only model is empty, not an error", {
  # there are no non-intercept terms to test, which is a valid answer rather
  # than a failure; rbind of nothing gives NULL, which used to error
  dd <- sim_lm(2)
  f <- lum_model(y ~ 1, data = dd, family = "gaussian", verbose = FALSE)
  a <- lum_anova(f, type = 3)
  expect_s3_class(a, "anova")
  expect_equal(nrow(a), 0L)
  expect_true("F value" %in% names(a))
  expect_silent(invisible(capture.output(print(a))))
  a2 <- lum_anova(f, type = 2)
  expect_equal(nrow(a2), 0L)
})

test_that("an intercept-only model summarises and predicts", {
  dd <- sim_lm(2)
  f <- lum_model(y ~ 1, data = dd, family = "gaussian", verbose = FALSE)
  out <- capture.output(print(summary(f)))
  expect_true(any(grepl("Linear model fit", out)))
  expect_equal(unname(predict(f)[1, 1]), unname(coef(f)[1]), tolerance = 1e-6)
})

test_that("coefficients and standard errors match lm() exactly", {
  dd <- sim_lm(3)
  f <- lum_model(y ~ x + z, data = dd, family = "gaussian", verbose = FALSE)
  m <- stats::lm(y ~ x + z, data = dd)
  cl <- summary(m)$coefficients
  expect_equal(unname(coef(f)), unname(cl[, 1]), tolerance = 1e-5)
  expect_equal(unname(sqrt(diag(vcov(f)))), unname(cl[, 2]), tolerance = 1e-5)
})

test_that("the residual SD is the unbiased one, as lm reports", {
  # maximum likelihood divides by n; lm divides by n - p.  Reporting the ML
  # value alongside rescaled standard errors would be internally inconsistent.
  dd <- sim_lm(4)
  f <- lum_model(y ~ x + z, data = dd, family = "gaussian", verbose = FALSE)
  m <- stats::lm(y ~ x + z, data = dd)
  expect_equal(unname(f$dispersion), summary(m)$sigma, tolerance = 1e-5)
})

test_that("tests are t, not z, when they can be exact", {
  dd <- sim_lm(5)
  f <- lum_model(y ~ x + z, data = dd, family = "gaussian", verbose = FALSE)
  ct <- lum_coef_table(f)
  expect_true("t value" %in% names(ct))
  expect_false("z value" %in% names(ct))
  m <- stats::lm(y ~ x + z, data = dd)
  expect_equal(unname(ct[["t value"]]), unname(summary(m)$coefficients[, 3]),
               tolerance = 1e-4)
  expect_equal(unname(ct[["Pr(>|t|)"]]), unname(summary(m)$coefficients[, 4]),
               tolerance = 1e-4)
})

test_that("anova uses F, and matches car::Anova on the equivalent lm", {
  skip_if_not_installed("car")
  dd <- sim_lm(6)
  f <- lum_model(y ~ x + z, data = dd, family = "gaussian", verbose = FALSE)
  m <- stats::lm(y ~ x + z, data = dd)
  a <- lum_anova(f, type = 3)
  ca <- car::Anova(m, type = 3)
  expect_true("F value" %in% names(a))
  expect_equal(a[["F value"]], ca[rownames(a), "F value"], tolerance = 1e-4)
  expect_equal(a[["Pr(>F)"]], ca[rownames(a), "Pr(>F)"], tolerance = 1e-4)
})

test_that("a mixed model still gets z and chi-square, not t and F", {
  # the exact reference exists only when nothing is integrated out
  dd <- sim_lm(7); dd$g <- factor(sample(20, nrow(dd), TRUE))
  f <- lum_model(y ~ x + (1 | g), data = dd, family = "gaussian", verbose = FALSE)
  expect_false(isTRUE(f$exact_df))
  expect_true("z value" %in% names(lum_coef_table(f)))
  expect_true("Chisq" %in% names(lum_anova(f, type = 3)))
})

test_that("non-gaussian fixed-effects models do not claim exact inference", {
  # a Poisson GLM has no exact small-sample t or F reference
  dd <- sim_lm(8)
  dd$cnt <- stats::rpois(nrow(dd), exp(0.4 + 0.3 * dd$x))
  f <- lum_model(cnt ~ x, data = dd, family = "poisson", verbose = FALSE)
  expect_false(isTRUE(f$exact_df))
  expect_true("z value" %in% names(lum_coef_table(f)))
})

test_that("the latent budget is not reported when nothing is latent", {
  dd <- sim_lm(9)
  f <- lum_model(y ~ x + z, data = dd, family = "gaussian", verbose = FALSE)
  expect_false("latent_budget" %in% f$checks$check)
  expect_true(all(f$checks$status %in% c("OK", "WARN", "FAIL", "BOUNDARY",
                                         "INCONCLUSIVE")))
})

test_that("summary describes a linear model rather than a mixed one", {
  dd <- sim_lm(10)
  f <- lum_model(y ~ x + z, data = dd, family = "gaussian", verbose = FALSE)
  out <- capture.output(print(summary(f)))
  expect_true(any(grepl("Linear model fit", out)))
  expect_true(any(grepl("Residual degrees of freedom", out)))
  expect_false(any(grepl("Random effects", out)))
  expect_false(any(grepl("groups:", out)))
})

test_that("predict works without random effects", {
  dd <- sim_lm(11)
  f <- lum_model(y ~ x + z, data = dd, family = "gaussian", verbose = FALSE)
  p <- predict(f)
  expect_equal(dim(p), c(nrow(dd), 1L))
  m <- stats::lm(y ~ x + z, data = dd)
  expect_equal(unname(p[, 1]), unname(stats::fitted(m)), tolerance = 1e-5)
  nd <- data.frame(x = c(-1, 0, 1), z = factor("a", levels = levels(dd$z)))
  expect_equal(nrow(predict(f, newdata = nd)), 3L)
})

test_that("simulate works without random effects", {
  dd <- sim_lm(12)
  f <- lum_model(y ~ x + z, data = dd, family = "gaussian", verbose = FALSE)
  ys <- lum_simulate(f, 3L, seed = 1L)
  expect_equal(dim(ys), c(nrow(dd), 3L))
  expect_true(all(is.finite(ys)))
  # a gaussian draw is a measurement, not a category index.  Checking only
  # finiteness would pass on a matrix of 1s and 2s, which is what a
  # multinomial-only simulator returns for a gaussian fit.
  expect_false(is.integer(ys))
  expect_gt(length(unique(as.vector(ys))), nrow(dd))
  expect_gt(sd(as.vector(ys)), 0.5 * sd(dd$y))
})
