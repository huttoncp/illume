pow_fit <- function(n = 200L, seed = 1L, b = 0.2) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n))
  d$y <- b * d$x + rnorm(n)
  ilm_model(y ~ x, data = d, family = "gaussian", verbose = FALSE)
}

test_that("simulated power matches the closed form for a linear model", {
  skip_on_cran()
  f <- pow_fit()
  ## the closed form uses the REALISED spread of x and the FITTED residual
  ## scale, not the population values they were drawn from: the simulation
  ## resamples this dataset, so those are the quantities it inherits
  sdx <- stats::sd(f$model$x); sig <- unname(f$dispersion)
  nn <- c(100L, 200L, 400L, 800L)
  cf <- vapply(nn, function(k) {
    ncp <- 0.2 * sqrt(k) * sdx / sig
    stats::pnorm(-1.959964 + ncp) + stats::pnorm(-1.959964 - ncp)
  }, 0)
  p <- ilm_power(f, n = nn, effect = 0.2, sims = 600L, seed = 1L,
                 progress = FALSE)
  expect_s3_class(p, "ilm_power")
  expect_equal(nrow(p), 4L)
  expect_lt(max(abs(p$power - cf)), 0.04)
  ## power rises with n, and every replicate converged for a model this simple
  expect_false(is.unsorted(p$power))
  expect_true(all(p$converged == 1))
  expect_equal(p$power, p$power_converged)
})

test_that("the rejection rate under the null is alpha", {
  skip_on_cran()
  f <- pow_fit()
  p <- ilm_power(f, n = c(200L, 800L), effect = 0, sims = 800L, seed = 2L,
                 progress = FALSE)
  expect_lt(max(abs(p$power - 0.05)), 0.025)
  ## and the nominal level is inside the Monte Carlo interval
  expect_true(all(p$mc_lower < 0.05 & 0.05 < p$mc_upper))
  ## a tighter level rejects less
  p1 <- ilm_power(f, n = 800L, effect = 0, sims = 800L, alpha = 0.01,
                  seed = 2L, progress = FALSE)
  expect_lt(p1$power, p$power[2])
})

test_that("the required sample size is a range, not a number", {
  skip_on_cran()
  f <- pow_fit()
  p <- ilm_power(f, n = c(100L, 200L, 400L, 800L), effect = 0.2, sims = 600L,
                 seed = 1L, progress = FALSE)
  r <- ilm_power_n(p, target = 0.8)
  expect_equal(nrow(r), 1L)
  ## against the closed-form n with the realised quantities
  sdx <- stats::sd(f$model$x); sig <- unname(f$dispersion)
  n_cf <- (1.959964 + 0.8416212)^2 * sig^2 / (0.2^2 * sdx^2)
  expect_lt(abs(r$n - n_cf) / n_cf, 0.15)
  ## the Monte Carlo error makes it an interval, and the point sits inside
  expect_lt(r$n_lower, r$n)
  expect_gt(r$n_upper, r$n)
  expect_error(ilm_power_n(p[1, , drop = FALSE]), "at least two sample sizes")
  expect_error(ilm_power_n(mtcars), "must be an ilm_power")
})

test_that("the Monte Carlo interval is reported and behaves", {
  skip_on_cran()
  f <- pow_fit()
  a <- ilm_power(f, n = 200L, effect = 0.2, sims = 100L, seed = 3L,
                 progress = FALSE)
  b <- ilm_power(f, n = 200L, effect = 0.2, sims = 1000L, seed = 3L,
                 progress = FALSE)
  ## more replicates narrow the interval; nothing else does
  expect_lt(b$mc_upper - b$mc_lower, a$mc_upper - a$mc_lower)
  expect_true(a$mc_lower < a$power && a$power < a$mc_upper)
  ## a Wilson interval stays inside [0, 1] even at the extremes, where a Wald
  ## one does not
  p1 <- ilm_power(f, n = 4000L, effect = 0.5, sims = 100L, seed = 4L,
                  progress = FALSE)
  expect_lte(p1$mc_upper, 1)
  expect_gte(p1$mc_lower, 0)
  expect_output(print(a), "Wilson interval")
})

test_that("it covers other families and effect grids", {
  skip_on_cran()
  set.seed(3); n <- 300L
  d <- data.frame(x = rnorm(n))
  d$y <- rbinom(n, 1L, stats::plogis(-0.3 + 0.5 * d$x))
  f <- ilm_model(y ~ x, data = d, family = "binomial", verbose = FALSE)
  p <- ilm_power(f, n = c(100L, 400L), effect = c(0.3, 0.5), sims = 200L,
                 seed = 3L, progress = FALSE)
  expect_equal(nrow(p), 4L)
  ## bigger effect and bigger n both raise power
  expect_gt(p$power[p$n == 400 & p$effect == 0.5],
            p$power[p$n == 100 & p$effect == 0.5])
  expect_gt(p$power[p$n == 400 & p$effect == 0.5],
            p$power[p$n == 400 & p$effect == 0.3])
  expect_equal(attr(p, "family"), "binomial")

  ## a count model
  set.seed(9); dp <- data.frame(x = rnorm(400L))
  dp$y <- rpois(400L, exp(0.5 + 0.3 * dp$x))
  fp <- ilm_model(y ~ x, data = dp, family = "poisson", verbose = FALSE)
  pp <- ilm_power(fp, n = c(100L, 400L), sims = 150L, seed = 9L,
                  progress = FALSE)
  expect_true(all(pp$power >= 0 & pp$power <= 1))
  expect_gt(pp$power[2], pp$power[1])
})

test_that("a mixed design resamples clusters, not rows", {
  skip_on_cran()
  set.seed(4); ng <- 30L; ni <- 8L; n <- ng * ni
  d <- data.frame(id = factor(rep(seq_len(ng), each = ni)), x = rnorm(n))
  d$y <- 0.4 * d$x + rep(rnorm(ng, 0, 0.8), each = ni) + rnorm(n)
  f <- ilm_model(y ~ x + (1 | id), data = d, family = "gaussian",
                 verbose = FALSE)
  p <- ilm_power(f, n = c(120L, 480L), sims = 100L, seed = 4L,
                 progress = FALSE)
  expect_true(attr(p, "grouped"))
  expect_output(print(p), "clusters are resampled whole")
  expect_gte(p$power[2], p$power[1])
  ## the resampler gives back roughly the size asked for, and relabels a
  ## cluster drawn twice as two clusters -- otherwise the design has fewer
  ## independent groups than it looks like it has
  set.seed(1)
  r <- ilm_power_resample(d, 480L, "id")
  expect_gt(nrow(r), 400L)
  expect_gt(nlevels(factor(r$id)), 40L)
  expect_true(all(table(r$id) == ni))
})

test_that("failures to converge are counted as non-detections", {
  skip_on_cran()
  ## a small, badly identified design, so some replicates fail
  set.seed(6); ns <- 60L
  d <- data.frame(x = rnorm(ns), g = factor(rep(1:6, each = 10)))
  d$y <- rbinom(ns, 1L, stats::plogis(-2 + 1.5 * d$x +
                                        rep(rnorm(6, 0, 2), each = 10)))
  f <- ilm_model(y ~ x + (1 | g), data = d, family = "binomial",
                 verbose = FALSE)
  p <- ilm_power(f, n = c(40L, 80L), sims = 80L, seed = 6L, progress = FALSE)
  expect_true(all(p$converged <= 1))
  ## power counting failures as misses is never above power conditional on
  ## convergence, and they differ exactly when something failed
  expect_true(all(p$power <= p$power_converged + 1e-12))
  if (any(p$converged < 1)) {
    expect_true(any(p$power < p$power_converged))
    expect_output(print(p), "CONDITIONAL ON CONVERGENCE")
  }
})

test_that("it refuses what it cannot simulate", {
  f <- pow_fit(100L)
  expect_error(ilm_power(f, term = "zz", sims = 5L), "is not a coefficient")
  expect_error(ilm_power(f, n = 5L, sims = 5L), "not a study")
  expect_error(ilm_power(mtcars), "must be a fitted ilm_model")
  set.seed(5); d <- data.frame(y = rnorm(100L))
  fi <- ilm_model(y ~ 1, data = d, family = "gaussian", verbose = FALSE)
  expect_error(ilm_power(fi, sims = 5L), "intercept only")
})
