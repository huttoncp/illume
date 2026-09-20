# Randomised quantile residuals and the homoscedasticity check.
#
# The residual test here is uniformity, not finiteness. ilm_rqr() previously
# errored for gaussian and poisson and returned plausible-looking but
# NON-UNIFORM values for binomial, and a test that only checked is.finite()
# passed on all of it.

sim_fam <- function(seed = 1, n = 1200) {
  set.seed(seed)
  d <- data.frame(g = factor(rep(seq_len(40), each = n / 40)),
                  x = stats::rnorm(n))
  b <- stats::rnorm(40, 0, 0.5)[as.integer(d$g)]
  e <- 0.3 + 0.8 * d$x + b
  d$gaus <- e + stats::rnorm(n)
  d$pois <- stats::rpois(n, exp(pmin(e, 4)))
  d$nb   <- stats::rnbinom(n, size = 2, mu = exp(pmin(e, 4)))
  d$bin  <- stats::rbinom(n, 1, 1 / (1 + exp(-e)))
  P <- exp(cbind(e, 0.4 * e) %*% t(stats::contr.sum(3))); P <- P / rowSums(P)
  d$cat <- factor(c("a", "b", "c")[apply(P, 1, function(p)
    sample.int(3, 1, prob = p))])
  d
}

fit_fam <- function(d, fam) {
  y <- c(gaussian = "gaus", poisson = "pois", nbinom = "nb",
         binomial = "bin", multinomial = "cat")[[fam]]
  ilm_model(stats::as.formula(paste(y, "~ x + (1 | g)")), data = d,
            family = fam, verbose = FALSE)
}

test_that("quantile residuals are uniform for every family", {
  d <- sim_fam()
  for (fam in c("gaussian", "poisson", "nbinom", "binomial", "multinomial")) {
    r <- ilm_rqr(fit_fam(d, fam), TRUE, 1L)
    expect_length(r, nrow(d))
    expect_true(all(r >= 0 & r <= 1), info = fam)
    p <- suppressWarnings(stats::ks.test(r, "punif")$p.value)
    expect_gt(p, 0.01, label = paste(fam, "uniformity p"))
  }
})

test_that("quantile residuals still detect a wrong family", {
  d <- sim_fam()
  # overdispersed counts fitted as Poisson
  f <- ilm_model(nb ~ x + (1 | g), data = d, family = "poisson",
                 verbose = FALSE)
  p <- suppressWarnings(stats::ks.test(ilm_rqr(f, TRUE, 1L), "punif")$p.value)
  expect_lt(p, 1e-5)
})

## ---- homoscedasticity ------------------------------------------------------

sim_var <- function(kind, seed, n = 900) {
  set.seed(seed)
  d <- data.frame(g = factor(rep(seq_len(30), each = n / 30)),
                  x = stats::rnorm(n),
                  s = factor(sample(c("a", "b"), n, TRUE)))
  b <- stats::rnorm(30, 0, 0.5)[as.integer(d$g)]
  mu <- 1 + 0.8 * d$x + b
  sdv <- switch(kind,
    constant = rep(1, n),
    trend    = 0.3 + 0.7 * exp(0.5 * d$x),
    bygroup  = ifelse(d$s == "a", 0.5, 2.0))
  d$y <- stats::rnorm(n, mu, sdv)
  d
}

test_that("constant variance is not flagged", {
  st <- vapply(1:3, function(s) {
    f <- ilm_model(y ~ x + (1 | g), data = sim_var("constant", s),
                   family = "gaussian", verbose = FALSE)
    suppressWarnings(ilm_check_variance(f, B = 100, plot = FALSE,
                                        verbose = FALSE))$status
  }, "")
  expect_false(any(st == "FAIL"))
})

test_that("spread trending with the fitted value is detected", {
  st <- vapply(1:3, function(s) {
    f <- ilm_model(y ~ x + (1 | g), data = sim_var("trend", s),
                   family = "gaussian", verbose = FALSE)
    suppressWarnings(ilm_check_variance(f, B = 100, plot = FALSE,
                                        verbose = FALSE))$status
  }, "")
  expect_true(all(st %in% c("WARN", "FAIL")))
})

test_that("variance differing between strata is detected, with the remedy named", {
  f <- ilm_model(y ~ x + s + (1 | g), data = sim_var("bygroup", 1),
                 family = "gaussian", verbose = FALSE)
  r <- suppressWarnings(ilm_check_variance(f, by = "s", B = 100, plot = FALSE,
                                           verbose = FALSE))
  expect_true(r$status %in% c("WARN", "FAIL"))
  expect_gt(r$ratio, 4)
  expect_match(r$suggestion, "separate variance per level")
})

test_that("the by argument is validated", {
  f <- ilm_model(y ~ x + (1 | g), data = sim_var("constant", 1),
                 family = "gaussian", verbose = FALSE)
  expect_error(suppressWarnings(ilm_check_variance(f, by = "nope", B = 20)),
               "not a column of the model frame")
  expect_error(suppressWarnings(ilm_check_variance(f, by = rep(1, 3), B = 20)),
               "values but the model has")
  expect_error(ilm_check_variance(data.frame(x = 1)), "must be a fitted ilm_model")
})

test_that("too few replicates is warned about", {
  f <- ilm_model(y ~ x + (1 | g), data = sim_var("constant", 1),
                 family = "gaussian", verbose = FALSE)
  expect_warning(ilm_check_variance(f, B = 20, plot = FALSE, verbose = FALSE),
                 "FAIL verdict is unreachable")
})
