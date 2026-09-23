beta_data <- function(n = 800L, seed = 1L, phi = 10) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n), g = factor(sample(c("a", "b"), n, TRUE)),
                  z = rnorm(n))
  mu <- stats::plogis(0.4 + 0.8 * d$x - 0.5 * (d$g == "b"))
  d$y <- stats::rbeta(n, mu * phi, (1 - mu) * phi)
  d
}
zbeta_data <- function(n = 1200L, seed = 7L, phi = 10) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n), z = rnorm(n))
  mu <- stats::plogis(0.4 + 0.8 * d$x)
  d$y <- ifelse(runif(n) < stats::plogis(-0.6 + 0.9 * d$z), 0,
                stats::rbeta(n, mu * phi, (1 - mu) * phi))
  d
}

test_that("beta regression recovers the mean and the precision", {
  d <- beta_data()
  f <- ilm_model(y ~ x + g, data = d, family = "beta", verbose = FALSE)
  expect_true(f$ok)
  expect_equal(f$family$name, "beta")
  expect_equal(unname(coef(f)), c(0.4, 0.8, -0.5), tolerance = 0.12)
  expect_equal(unname(f$dispersion), 10, tolerance = 0.2)
  ## the mean is on the logit scale and the fitted values stay inside (0, 1)
  p <- as.numeric(predict(f))
  expect_true(all(p > 0 & p < 1))
  expect_equal(mean(p), mean(d$y), tolerance = 0.02)
  ## larger phi means LESS spread, which is the opposite of every other
  ## dispersion parameter here
  f2 <- ilm_model(y ~ x + g, data = beta_data(seed = 2L, phi = 60),
                  family = "beta", verbose = FALSE)
  expect_gt(unname(f2$dispersion), unname(f$dispersion))
  expect_lt(stats::sd(f2$y), stats::sd(f$y))
})

test_that("beta agrees with glmmTMB, with and without a precision model", {
  skip_if_not_installed("glmmTMB")
  d <- beta_data()
  f <- ilm_model(y ~ x + g, data = d, family = "beta", verbose = FALSE)
  m <- glmmTMB::glmmTMB(y ~ x + g, data = d, family = glmmTMB::beta_family())
  expect_equal(unname(coef(f)), unname(glmmTMB::fixef(m)$cond),
               tolerance = 1e-4)
  expect_equal(unname(sqrt(diag(vcov(f)))), unname(sqrt(diag(vcov(m)$cond))),
               tolerance = 1e-5)
  expect_equal(unname(f$dispersion), glmmTMB::sigma(m), tolerance = 1e-4)
  expect_equal(as.numeric(logLik(f)), as.numeric(logLik(m)), tolerance = 1e-6)

  ## variable precision: the spread changes with a covariate
  set.seed(2); n <- 800L
  d2 <- data.frame(x = rnorm(n), z = rnorm(n))
  mu2 <- stats::plogis(0.3 + 0.7 * d2$x); phi2 <- exp(2.5 + 0.8 * d2$z)
  d2$y <- stats::rbeta(n, mu2 * phi2, (1 - mu2) * phi2)
  f2 <- ilm_model(y ~ x, data = d2, family = "beta", dispformula = ~ z,
                  verbose = FALSE)
  m2 <- glmmTMB::glmmTMB(y ~ x, dispformula = ~ z, data = d2,
                         family = glmmTMB::beta_family())
  expect_equal(unname(coef(f2)), unname(glmmTMB::fixef(m2)$cond),
               tolerance = 1e-3)
  expect_equal(unname(f2$disp_coef), unname(glmmTMB::fixef(m2)$disp),
               tolerance = 1e-3)
  expect_equal(as.numeric(logLik(f2)), as.numeric(logLik(m2)), tolerance = 1e-5)
  expect_equal(unname(f2$disp_coef), c(2.5, 0.8), tolerance = 0.2)
})

test_that("beta works with a random effect and gives uniform residuals", {
  set.seed(4); ng <- 50L; ni <- 12L; n <- ng * ni
  d <- data.frame(id = factor(rep(seq_len(ng), each = ni)), x = rnorm(n))
  b <- rnorm(ng, 0, 0.7)
  mu <- stats::plogis(0.3 + 0.6 * d$x + b[as.integer(d$id)])
  d$y <- stats::rbeta(n, mu * 10, (1 - mu) * 10)
  f <- ilm_model(y ~ x + (1 | id), data = d, family = "beta", verbose = FALSE)
  expect_true(f$ok)
  expect_equal(unname(coef(f)), c(0.3, 0.6), tolerance = 0.15)
  expect_equal(sqrt(f$Sigma$id[1, 1]), 0.7, tolerance = 0.2)
  expect_equal(unname(f$dispersion), 10, tolerance = 2)

  fb <- ilm_model(y ~ x + g, data = beta_data(1000L), family = "beta",
                  verbose = FALSE)
  u <- ilm_rqr(fb, seed = 1L)
  expect_true(all(u > 0 & u < 1))
  expect_gt(suppressWarnings(stats::ks.test(u, "punif")$p.value), 0.01)
  ys <- ilm_simulate(fb, 50L, seed = 2L)
  expect_true(all(ys > 0 & ys < 1))
  expect_equal(mean(ys), mean(fb$y), tolerance = 0.02)
})

test_that("a value at a boundary is refused, with both remedies named", {
  d <- beta_data(300L)
  d0 <- d; d0$y[1:3] <- 0
  e <- tryCatch(ilm_model(y ~ x, data = d0, family = "beta", verbose = FALSE),
                error = conditionMessage)
  expect_match(e, "OPEN interval")
  expect_match(e, "ziformula")          # separate process
  expect_match(e, "ilm_squeeze")        # rounding
  d1 <- d; d1$y[1:2] <- 1
  expect_error(ilm_model(y ~ x, data = d1, family = "beta", verbose = FALSE),
               "OPEN interval")
  ## and the squeeze is a real remedy: it moves everything and says how far
  expect_message(ilm_squeeze(c(0, 0.5, 1)), "moved towards")
  s <- ilm_squeeze(d0$y, quiet = TRUE)
  expect_true(all(s > 0 & s < 1))
  expect_silent(ilm_model(s ~ x, data = transform(d0, s = s), family = "beta",
                          verbose = FALSE))
  expect_error(ilm_squeeze(c(-1, 0.5)), "outside")
  expect_error(ilm_squeeze(c(0, 1), n = 1), "at least 2")
})

test_that("a zero part on a beta response is a hurdle and nothing else", {
  d <- zbeta_data()
  expect_gt(mean(d$y == 0), 0.2)
  f <- ilm_model(y ~ x, data = d, family = "beta", ziformula = ~ z,
                 zi_type = "hurdle", verbose = FALSE)
  expect_true(f$ok)
  expect_equal(unname(coef(f)), c(0.4, 0.8), tolerance = 0.1)
  expect_equal(ilm_zi_coef(f)$estimate, c(-0.6, 0.9), tolerance = 0.15)
  ## a MIXTURE is not defined for a continuous response and says so, rather
  ## than fitting a worse approximation of something that does not exist
  expect_error(ilm_model(y ~ x, data = d, family = "beta", ziformula = ~ z,
                         verbose = FALSE),
               "continuous")
  ## the marginal mean carries the zeros
  expect_equal(mean(as.numeric(predict(f))), mean(d$y), tolerance = 0.02)
  ## and it is NOT the count part's mean, which is the mean of the positives
  expect_gt(mean(stats::plogis(as.numeric(predict(f, type = "link")))),
            mean(as.numeric(predict(f))) * 1.2)
})

test_that("a continuous zero part simulates and residualises correctly", {
  d <- zbeta_data(800L, seed = 9L)
  f <- ilm_model(y ~ x, data = d, family = "beta", ziformula = ~ z,
                 zi_type = "hurdle", verbose = FALSE)
  ## the only jump is the point mass at zero; above it the response is smooth,
  ## so F(y) - F(y - 1) is meaningless and using it silently destroys the
  ## uniformity while leaving the residuals looking like residuals
  u <- ilm_rqr(f, seed = 1L)
  expect_true(all(u > 0 & u < 1))
  expect_gt(max(u), 0.97)
  expect_gt(suppressWarnings(stats::ks.test(u, "punif")$p.value), 0.01)

  ys <- ilm_simulate(f, 100L, seed = 3L)
  expect_true(all(ys >= 0 & ys < 1))
  expect_equal(mean(ys == 0), mean(d$y == 0), tolerance = 0.03)
  expect_equal(mean(ys), mean(d$y), tolerance = 0.03)
  ## the zero check applies here, because the model has a part for them
  z <- suppressMessages(ilm_check_zeros(f, B = 200L))
  expect_equal(z$status, "OK")
  ## but not to a plain beta fit, which could not hold a zero anyway
  fb <- ilm_model(y ~ x + g, data = beta_data(300L), family = "beta",
                  verbose = FALSE)
  expect_error(ilm_check_zeros(fb), "cannot contain a zero")
})

test_that("zero-inflated beta agrees with glmmTMB", {
  skip_if_not_installed("glmmTMB")
  d <- zbeta_data()
  f <- ilm_model(y ~ x, data = d, family = "beta", ziformula = ~ z,
                 zi_type = "hurdle", verbose = FALSE)
  m <- glmmTMB::glmmTMB(y ~ x, ziformula = ~ z, data = d,
                        family = glmmTMB::beta_family())
  expect_equal(unname(coef(f)), unname(glmmTMB::fixef(m)$cond), tolerance = 1e-4)
  expect_equal(ilm_zi_coef(f)$estimate, unname(glmmTMB::fixef(m)$zi),
               tolerance = 1e-4)
  expect_equal(ilm_zi_coef(f)$se, unname(sqrt(diag(vcov(m)$zi))),
               tolerance = 1e-5)
  expect_equal(unname(f$dispersion), glmmTMB::sigma(m), tolerance = 1e-4)
  expect_equal(as.numeric(logLik(f)), as.numeric(logLik(m)), tolerance = 1e-6)
})

test_that("cluster-robust standard errors cover the beta family", {
  set.seed(11); G <- 40L; ni <- 10L; n <- G * ni
  d <- data.frame(g = factor(rep(seq_len(G), each = ni)), x = rnorm(n))
  mu <- stats::plogis(0.3 + 0.6 * d$x + rep(rnorm(G, 0, 0.8), each = ni))
  d$y <- stats::rbeta(n, mu * 12, (1 - mu) * 12)
  f <- ilm_model(y ~ x, data = d, family = "beta", verbose = FALSE)
  ## the score sums to zero at the maximum, which is what being there means
  expect_equal(unname(colSums(ilm_estfun(f))), c(0, 0), tolerance = 1e-3)
  r <- ilm_robust(f, ~ g)
  expect_equal(nrow(r), 2L)
  expect_true(all(r$se > 0))
  ## ignoring a cluster effect this size understates the intercept's error
  expect_gt(r$se[r$term == "(Intercept)"] / r$se_model[r$term == "(Intercept)"],
            1.5)
})
