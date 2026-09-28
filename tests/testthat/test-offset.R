## Offsets: offset(log(exposure)) in the formula, as glm() and glmmTMB read it.
## Each fit is compared with the reference implementation on the same data.

q <- function(e) suppressMessages(suppressWarnings(e))
rates <- function(n = 300, seed = 1) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n), e = runif(n, 0.5, 20),
                  g = factor(rep(seq_len(30), length.out = n)))
  d$y <- rpois(n, d$e * exp(-1 + 0.4 * d$x + rnorm(30, 0, 0.5)[d$g]))
  d
}

test_that("a Poisson rate matches glm() with the same offset", {
  d <- rates()
  f <- q(ilm_model(y ~ x + offset(log(e)), data = d, family = "poisson", verbose = FALSE))
  g <- stats::glm(y ~ x + offset(log(e)), family = poisson, data = d)
  expect_equal(unname(coef(f)), unname(coef(g)), tolerance = 1e-5)
  expect_equal(as.numeric(logLik(f)), as.numeric(logLik(g)), tolerance = 1e-6)
  expect_equal(unname(sqrt(diag(vcov(f)))), unname(sqrt(diag(vcov(g)))), tolerance = 1e-4)
  ## the fitted counts carry each row's exposure
  expect_equal(as.numeric(predict(f, type = "response")), unname(fitted(g)),
               tolerance = 1e-5)
  ## and leaving the offset out is a different model
  f0 <- q(ilm_model(y ~ x, data = d, family = "poisson", verbose = FALSE))
  expect_gt(abs(coef(f0)[[1]] - coef(f)[[1]]), 1)
})

test_that("a negative binomial rate matches MASS::glm.nb()", {
  skip_if_not_installed("MASS")
  set.seed(2); n <- 400
  d <- data.frame(x = rnorm(n), e = runif(n, 1, 10))
  d$y <- MASS::rnegbin(n, mu = d$e * exp(0.5 + 0.3 * d$x), theta = 3)
  f <- q(ilm_model(y ~ x + offset(log(e)), data = d, family = "nbinom", verbose = FALSE))
  g <- MASS::glm.nb(y ~ x + offset(log(e)), data = d)
  expect_equal(unname(coef(f)), unname(coef(g)), tolerance = 1e-4)
  expect_equal(as.numeric(logLik(f)), as.numeric(logLik(g)), tolerance = 1e-5)
})

test_that("a binomial offset matches glm()", {
  ## a logit offset -- a known log-odds shift, as a case-control sampling
  ## correction is -- since illume's binomial family has the logit link only
  set.seed(3); n <- 500
  d <- data.frame(x = rnorm(n), s = rnorm(n, 0, 0.5))
  d$y <- rbinom(n, 1, plogis(-0.5 + 0.5 * d$x + d$s))
  f <- q(ilm_model(y ~ x + offset(s), data = d, family = "binomial", verbose = FALSE))
  g <- stats::glm(y ~ x + offset(s), family = binomial, data = d)
  expect_equal(unname(coef(f)), unname(coef(g)), tolerance = 1e-5)
})

test_that("a gaussian offset matches lm()", {
  set.seed(4); n <- 100
  d <- data.frame(x = rnorm(n), b = rnorm(n))
  d$y <- 1 + 0.5 * d$x + d$b + rnorm(n)
  f <- q(ilm_model(y ~ x + offset(b), data = d, family = "gaussian", verbose = FALSE))
  g <- stats::lm(y ~ x + offset(b), data = d)
  expect_equal(unname(coef(f)), unname(coef(g)), tolerance = 1e-5)
  expect_equal(unname(sqrt(diag(vcov(f)))), unname(sqrt(diag(vcov(g)))), tolerance = 1e-5)
})

test_that("a mixed Poisson rate matches glmmTMB", {
  skip_if_not_installed("glmmTMB")
  d <- rates()
  f <- q(ilm_model(y ~ x + offset(log(e)) + (1 | g), data = d, family = "poisson",
                   verbose = FALSE))
  g <- glmmTMB::glmmTMB(y ~ x + offset(log(e)) + (1 | g), family = poisson, data = d)
  expect_equal(unname(coef(f)), unname(glmmTMB::fixef(g)$cond), tolerance = 1e-4)
  expect_equal(as.numeric(logLik(f)), as.numeric(logLik(g)), tolerance = 1e-5)
  expect_equal(sqrt(unname(f$Sigma[[1]][1, 1])),
               unname(attr(glmmTMB::VarCorr(g)$cond$g, "stddev")), tolerance = 1e-3)
})

test_that("new data predicts at its own exposure", {
  d <- rates()
  f <- q(ilm_model(y ~ x + offset(log(e)), data = d, family = "poisson", verbose = FALSE))
  nd <- data.frame(x = c(0, 0), e = c(1, 2))
  p <- as.numeric(predict(f, newdata = nd, type = "response"))
  expect_equal(p[2] / p[1], 2, tolerance = 1e-10)
  expect_equal(p[1], exp(coef(f)[[1]]), tolerance = 1e-10)
  ## new data without the exposure has nothing to predict at, and says which
  expect_error(predict(f, newdata = data.frame(x = 0)), "has no column .e.")
})

test_that("a refit and a simulation keep the offset", {
  d <- rates()
  f <- q(ilm_model(y ~ x + offset(log(e)), data = d, family = "poisson", verbose = FALSE))
  d2 <- rates(seed = 9)
  f2 <- q(ilm_refit(f, data = d2))
  g2 <- stats::glm(y ~ x + offset(log(e)), family = poisson, data = d2)
  expect_equal(unname(coef(f2)), unname(coef(g2)), tolerance = 1e-5)
  ## simulated counts are centred on the fitted counts, exposure and all
  s <- q(ilm_simulate(f, nsim = 200, seed = 1))
  expect_equal(mean(as.matrix(s)), mean(fitted(stats::glm(y ~ x + offset(log(e)),
                                                           family = poisson, data = d))),
               tolerance = 0.02)
})

test_that("an offset's missing values drop rows like any other column", {
  d <- rates(); d$e[1:4] <- NA
  f <- q(ilm_model(y ~ x + offset(log(e)), data = d, family = "poisson", verbose = FALSE))
  r <- ilm_rows_used(f)
  expect_identical(c(r$n_used, r$n_dropped), c(296L, 4L))
  expect_identical(r$dropped_by, c(e = 4L))
  expect_length(f$offset, 296L)
})

test_that("an offset is refused where there is no single linear predictor", {
  set.seed(5); n <- 200
  d <- data.frame(x = rnorm(n), e = runif(n, 1, 5))
  d$k <- factor(sample(c("a", "b", "c"), n, replace = TRUE))
  expect_error(ilm_model(k ~ x + offset(log(e)), data = d, family = "multinomial",
                         verbose = FALSE), "one linear predictor")
  d$o <- factor(sample(1:3, n, replace = TRUE), ordered = TRUE)
  expect_error(ilm_model(o ~ x + offset(log(e)), data = d, family = "ordinal",
                         verbose = FALSE), "one linear predictor")
})

test_that("an offset of zero is no offset", {
  d <- rates(); d$z0 <- 0
  f <- q(ilm_model(y ~ x + offset(z0), data = d, family = "poisson", verbose = FALSE))
  expect_null(f$offset)
})

test_that("fractional counts are pointed at an offset, which then fits", {
  d <- rates()
  d$rate <- d$y / d$e
  expect_error(ilm_model(rate ~ x, data = d, family = "poisson", verbose = FALSE),
               "offset(log(exposure))", fixed = TRUE)
  f <- q(ilm_model(y ~ x + offset(log(e)), data = d, family = "poisson", verbose = FALSE))
  expect_true(f$ok)
})
