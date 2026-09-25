# The dispersion formula: the remedy for what ilm_check_variance() diagnoses.
#
# The check reports two patterns -- spread that trends with the fitted value,
# and spread that differs between groups -- and both are the same fix: model
# the logarithm of the dispersion. `~ group` for the second, `~ mu` for the
# first, `~ x` for anything else. `mu` is a reserved name meaning the fitted
# mean.
#
# Pinned to glmmTMB, whose dispformula is the same idea, and to nlme's
# varIdent and varPower, which are the same idea in different clothes.

disp_group <- function(seed = 31, n = 1200) {
  set.seed(seed)
  d <- data.frame(x = stats::rnorm(n),
                  g = factor(sample(c("a", "b", "c"), n, TRUE)))
  sdv <- c(a = 0.5, b = 1.5, c = 3.0)[as.character(d$g)]
  d$y <- 1 + 0.8 * d$x + stats::rnorm(n, 0, sdv)
  d
}

## ---- against the references ------------------------------------------------

test_that("a separate variance per group matches glmmTMB and gls", {
  skip_if_not_installed("glmmTMB")
  skip_if_not_installed("nlme")
  d <- disp_group()
  f <- ilm_model(y ~ x, data = d, family = "gaussian", dispformula = ~ g,
                 verbose = FALSE)
  g1 <- glmmTMB::glmmTMB(y ~ x, data = d, family = stats::gaussian(),
                         dispformula = ~ g)
  n1 <- nlme::gls(y ~ x, data = d, weights = nlme::varIdent(form = ~ 1 | g),
                  method = "ML")
  expect_equal(as.numeric(logLik(f)), as.numeric(logLik(g1)), tolerance = 1e-5)
  expect_equal(as.numeric(logLik(f)), as.numeric(logLik(n1)), tolerance = 1e-5)
  expect_equal(unname(coef(f)), unname(glmmTMB::fixef(g1)$cond), tolerance = 1e-4)
  expect_equal(unname(sqrt(diag(vcov(f)))[["x"]]),
               unname(sqrt(diag(vcov(g1)$cond))[["x"]]), tolerance = 1e-3)
  # and the fitted spreads are the ones the data were made with
  dv <- illume:::ilm_disp_vec(f)
  by_g <- sort(vapply(split(dv, d$g), function(z) z[1], 0))
  expect_equal(unname(by_g), c(0.5, 1.5, 3.0), tolerance = 0.12)
})

test_that("a dispersion that changes with a covariate matches glmmTMB", {
  skip_if_not_installed("glmmTMB")
  set.seed(32)
  n <- 1200
  d <- data.frame(x = stats::rnorm(n), z = stats::rnorm(n))
  d$y <- 1 + 0.8 * d$x + stats::rnorm(n, 0, exp(-0.2 + 0.5 * d$z))
  f <- ilm_model(y ~ x, data = d, family = "gaussian", dispformula = ~ z,
                 verbose = FALSE)
  g <- glmmTMB::glmmTMB(y ~ x, data = d, family = stats::gaussian(),
                        dispformula = ~ z)
  expect_equal(as.numeric(logLik(f)), as.numeric(logLik(g)), tolerance = 1e-5)
  expect_equal(unname(f$disp_coef), unname(glmmTMB::fixef(g)$disp),
               tolerance = 1e-3)
  # the truth, recovered
  expect_equal(unname(f$disp_coef[[1]]), -0.2, tolerance = 0.06)
  expect_equal(unname(f$disp_coef[[2]]), 0.5, tolerance = 0.08)
})

test_that("a power of the fitted mean matches nlme's varPower", {
  skip_if_not_installed("nlme")
  set.seed(33)
  n <- 1200
  d <- data.frame(x = stats::runif(n, 1, 5))
  mu <- 2 + 1.5 * d$x
  d$y <- stats::rnorm(n, mu, 0.25 * mu)
  f <- ilm_model(y ~ x, data = d, family = "gaussian", dispformula = ~ mu,
                 verbose = FALSE)
  g <- nlme::gls(y ~ x, data = d,
                 weights = nlme::varPower(form = ~ fitted(.)), method = "ML")
  expect_equal(unname(f$disp_coef[["disp:mu_power"]]),
               unname(coef(g$modelStruct$varStruct, unconstrained = FALSE)[[1]]),
               tolerance = 0.02)
  expect_equal(as.numeric(logLik(f)), as.numeric(logLik(g)), tolerance = 1e-2)
  expect_equal(f$opt$convergence, 0L)
})

test_that("the mean-power term is warm started, not left to run away", {
  # from beta = 0 the fitted mean is zero for every row, log|mu| is the log of
  # the numerical floor, and the power multiplying it can go anywhere. A cold
  # start converged falsely with the power at -115818 and the slope at zero.
  set.seed(33)
  n <- 1200
  d <- data.frame(x = stats::runif(n, 1, 5))
  mu <- 2 + 1.5 * d$x
  d$y <- stats::rnorm(n, mu, 0.25 * mu)
  f <- ilm_model(y ~ x, data = d, family = "gaussian", dispformula = ~ mu,
                 verbose = FALSE)
  expect_equal(f$opt$convergence, 0L)
  expect_lt(abs(f$disp_coef[["disp:mu_power"]] - 1), 0.3)
  expect_equal(unname(coef(f)[["x"]]), 1.5, tolerance = 0.1)
})

test_that("a mean that crosses zero is warned about", {
  set.seed(35)
  n <- 400
  d <- data.frame(x = stats::rnorm(n))
  d$y <- stats::rnorm(n, 2 * d$x, 0.5)     # the mean spans zero
  expect_warning(ilm_model(y ~ x, data = d, family = "gaussian",
                           dispformula = ~ mu, verbose = FALSE),
                 "changes sign")
})

## ---- it fixes what the check complained about ------------------------------

test_that("the diagnosis clears once the remedy is applied", {
  set.seed(34)
  ng <- 30; nt <- 30
  d <- data.frame(id = factor(rep(seq_len(ng), each = nt)),
                  x = stats::rnorm(ng * nt),
                  s = factor(sample(c("a", "b"), ng * nt, TRUE)))
  b <- stats::rnorm(ng, 0, 0.5)[as.integer(d$id)]
  d$y <- 1 + 0.8 * d$x + b +
    stats::rnorm(ng * nt, 0, ifelse(d$s == "a", 0.5, 2.0))
  m0 <- ilm_model(y ~ x + s + (1 | id), data = d, family = "gaussian",
                  verbose = FALSE)
  m1 <- ilm_model(y ~ x + s + (1 | id), data = d, family = "gaussian",
                  dispformula = ~ s, verbose = FALSE)
  r0 <- suppressWarnings(ilm_check_variance(m0, by = "s", B = 60, plot = FALSE,
                                            verbose = FALSE))
  r1 <- suppressWarnings(ilm_check_variance(m1, by = "s", B = 60, plot = FALSE,
                                            verbose = FALSE))
  expect_true(r0$status %in% c("WARN", "FAIL"))
  expect_equal(r1$status, "OK")
  expect_gt(r0$ratio, 4)
  expect_lt(r1$ratio, 2)
  expect_lt(AIC(m1), AIC(m0) - 100)
})

test_that("the check names the remedy it is asking for", {
  set.seed(34)
  n <- 600
  d <- data.frame(x = stats::rnorm(n), s = factor(sample(c("a", "b"), n, TRUE)))
  d$y <- 1 + 0.8 * d$x + stats::rnorm(n, 0, ifelse(d$s == "a", 0.5, 2.5))
  m <- ilm_model(y ~ x + s, data = d, family = "gaussian", verbose = FALSE)
  r <- suppressWarnings(ilm_check_variance(m, by = "s", B = 60, plot = FALSE,
                                           verbose = FALSE))
  expect_match(r$suggestion, "dispformula", fixed = TRUE)
})

## ---- the rest of the package copes -----------------------------------------

test_that("the per-row dispersion reaches the residuals and the simulators", {
  d <- disp_group(n = 800)
  f <- ilm_model(y ~ x, data = d, family = "gaussian", dispformula = ~ g,
                 verbose = FALSE)
  g0 <- ilm_model(y ~ x, data = d, family = "gaussian", verbose = FALSE)
  # uniform with the model, hopeless without: the residual is being divided by
  # each row's own standard deviation
  expect_gt(suppressWarnings(stats::ks.test(ilm_rqr(f, seed = 1L), "punif")$p.value),
            0.01)
  expect_lt(suppressWarnings(stats::ks.test(ilm_rqr(g0, seed = 1L), "punif")$p.value),
            1e-6)
  # Pearson residuals are unit-scaled for the same reason
  expect_lt(abs(stats::sd(illume:::ilm_pearson_ovr(f)) - 1), 0.1)
  # simulated data inherits the varying spread
  ys <- ilm_simulate(f, 20, seed = 1)
  sds <- vapply(split(as.vector(ys), rep(d$g, 20)), stats::sd, 0)
  expect_equal(unname(sort(sds)), c(0.5, 1.5, 3.0), tolerance = 0.25)
})

test_that("a dispersion model rules out exact t and F inference", {
  # exact inference rests on a constant variance. With the variance itself
  # estimated as a function of covariates that independence goes, and the t
  # distribution becomes an approximation; gls() reports t here by convention.
  d <- disp_group(n = 600)
  f <- ilm_model(y ~ x, data = d, family = "gaussian", dispformula = ~ g,
                 verbose = FALSE)
  g0 <- ilm_model(y ~ x, data = d, family = "gaussian", verbose = FALSE)
  expect_false(isTRUE(f$exact_df))
  expect_true(isTRUE(g0$exact_df))
  expect_true("z value" %in% names(ilm_coef_table(f)))
  expect_true("t value" %in% names(ilm_coef_table(g0)))
})

test_that("summary says the spread is not constant", {
  d <- disp_group(n = 600)
  f <- ilm_model(y ~ x, data = d, family = "gaussian", dispformula = ~ g,
                 verbose = FALSE)
  out <- capture.output(summary(f))
  expect_true(any(grepl("Dispersion model", out)))
  expect_true(any(grepl("median row", out)))
})

## ---- informative failure ---------------------------------------------------

test_that("misspecification is named", {
  d <- disp_group(n = 300)
  expect_error(ilm_model(y ~ x, data = d, family = "gaussian",
                         dispformula = y ~ g, verbose = FALSE),
               "one-sided formula")
  expect_error(ilm_model(y ~ x, data = d, family = "gaussian",
                         dispformula = ~ nope, verbose = FALSE),
               "not in the data")
  expect_error(ilm_model(y ~ x, data = d, family = "gaussian",
                         dispformula = ~ nope, verbose = FALSE),
               "reserved name")
  # a family with no dispersion has nothing to model
  d$cnt <- stats::rpois(nrow(d), 3)
  expect_error(ilm_model(cnt ~ x, data = d, family = "poisson",
                         dispformula = ~ g, verbose = FALSE),
               "no dispersion parameter to model")
  expect_error(ilm_model(cnt ~ x, data = d, family = "poisson",
                         dispformula = ~ g, verbose = FALSE),
               "nbinom", fixed = TRUE)
})

test_that("it works for families other than gaussian", {
  set.seed(36)
  n <- 1500
  d <- data.frame(x = stats::rnorm(n),
                  g = factor(sample(c("a", "b"), n, TRUE)))
  k <- ifelse(d$g == "a", 8, 1)
  d$cnt <- stats::rnbinom(n, size = k, mu = exp(1.2 + 0.4 * d$x))
  f <- ilm_model(cnt ~ x, data = d, family = "nbinom", dispformula = ~ g,
                 verbose = FALSE)
  expect_equal(f$opt$convergence, 0L)
  dv <- illume:::ilm_disp_vec(f)
  by_g <- vapply(split(dv, d$g), function(z) z[1], 0)
  # the group with less overdispersion has the larger k
  expect_gt(by_g[["a"]], by_g[["b"]])
  expect_lt(AIC(f),
            AIC(ilm_model(cnt ~ x, data = d, family = "nbinom", verbose = FALSE)))
})
