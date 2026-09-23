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

test_that("the test counted is the one the analysis reports", {
  skip_on_cran()
  ## a gaussian fit with nothing integrated out reports t, not z; at 20 rows
  ## the difference is the difference between 0.395 and 0.43
  set.seed(21); d <- data.frame(g = factor(rep(c("a", "b"), 10)))
  d$y <- 0.8 * (d$g == "b") + rnorm(20)
  f <- ilm_model(y ~ g, data = d, family = "gaussian", verbose = FALSE)
  p <- ilm_power(f, n = 20L, effect = 0.8, sims = 20L, progress = FALSE)
  expect_identical(attr(p, "test"), "t")
  expect_output(print(p), "counting the t test")
  fb <- ilm_model(yb ~ x, family = "binomial", verbose = FALSE,
                  data = transform(data.frame(x = rnorm(200)),
                                   yb = rbinom(200, 1, 0.5)))
  expect_identical(attr(ilm_power(fb, n = 100L, sims = 5L, progress = FALSE),
                        "test"), "Wald z")
})

test_that("a term with several coefficients is tested jointly", {
  skip_on_cran()
  set.seed(22); d <- data.frame(g = factor(rep(c("a", "b", "c"), 40)))
  d$y <- c(0, 0.4, 0.8)[as.integer(d$g)] + rnorm(120)
  f <- ilm_model(y ~ g, data = d, family = "gaussian", verbose = FALSE)
  p <- ilm_power(f, n = c(60L, 120L), term = "g", sims = 40L,
                 progress = FALSE)
  expect_identical(attr(p, "test"), "F")
  expect_identical(attr(p, "coefs"), c("gb", "gc"))
  expect_identical(attr(p, "effect_scale"), "multiple")
  expect_equal(p$effect, c(1, 1))
  expect_output(print(p), "MULTIPLE of the assumed")
  ## one of its coefficients can still be followed on its own
  p1 <- ilm_power(f, n = 60L, term = "gc", sims = 10L, progress = FALSE)
  expect_identical(attr(p1, "coefs"), "gc")
  expect_identical(attr(p1, "effect_scale"), "link")
})

test_that("contrasts survive the refit", {
  skip_on_cran()
  ## sum-to-zero coding names the coefficient g1, and a refit that dropped the
  ## contrasts came back with gb instead: every replicate failed and the
  ## power was exactly zero
  set.seed(23); d <- data.frame(g = factor(rep(c("a", "b"), 50)))
  d$y <- 0.5 * (d$g == "b") + rnorm(100)
  f1 <- ilm_model(y ~ g, data = d, family = "gaussian", verbose = FALSE)
  f2 <- ilm_model(y ~ g, data = d, family = "gaussian", verbose = FALSE,
                  contrasts = list(g = "contr.sum"))
  p1 <- ilm_power(f1, n = c(50L, 100L), sims = 60L, seed = 4L,
                  progress = FALSE)
  p2 <- ilm_power(f2, n = c(50L, 100L), sims = 60L, seed = 4L,
                  progress = FALSE)
  expect_true(all(p2$converged == 1))
  ## the same study under either coding
  expect_equal(p2$power, p1$power)
})

test_that("a transformed response or predictor is simulated on the model's scale", {
  skip_on_cran()
  ## the model frame holds log(yp) and log(x), not yp and x; a refit through
  ## the formula found neither and failed on every replicate
  set.seed(24); d <- data.frame(x = stats::rexp(150) + 0.1)
  d$yp <- exp(0.2 + 0.3 * log(d$x) + rnorm(150, 0, 0.5))
  f <- ilm_model(log(yp) ~ log(x), data = d, family = "gaussian",
                 verbose = FALSE)
  p <- ilm_power(f, n = c(50L, 150L), sims = 40L, progress = FALSE)
  expect_true(all(p$converged == 1))
  expect_gt(p$power[2], p$power[1])
})

test_that("a multinomial model is tested across its categories", {
  skip_on_cran()
  set.seed(25); n <- 300
  d <- data.frame(x = rnorm(n))
  eta <- cbind(0, 0.2 + 0.6 * d$x, -0.1 - 0.4 * d$x)
  pr <- exp(eta) / rowSums(exp(eta))
  d$y <- factor(apply(pr, 1, function(p) sample(c("a", "b", "c"), 1, prob = p)))
  f <- ilm_model(y ~ x, data = d, family = "multinomial", verbose = FALSE)
  ## by default the first TERM, across both of its categories -- the first
  ## coefficient would have been a category's intercept
  p <- ilm_power(f, n = c(60L, 150L), sims = 40L, progress = FALSE)
  expect_identical(attr(p, "term"), "x")
  expect_identical(attr(p, "coefs"), c("a:x", "b:x"))
  expect_identical(attr(p, "test"), "Wald chi-square")
  expect_true(all(p$converged > 0.9))
  expect_gt(p$power[2], p$power[1])
  ## and one category on its own
  p1 <- ilm_power(f, n = 150L, term = "b:x", sims = 20L, progress = FALSE)
  expect_identical(attr(p1, "test"), "Wald z")
})

test_that("frequency weights are drawn as the observations they stand for", {
  skip_on_cran()
  set.seed(26); d <- data.frame(x = rep(c(-1, 0, 1), each = 40))
  d$y <- rpois(nrow(d), exp(0.5 + 0.25 * d$x))
  agg <- stats::aggregate(list(w = rep(1, nrow(d))), by = list(x = d$x, y = d$y),
                          FUN = sum)
  fw <- ilm_model(y ~ x, data = agg, family = "poisson", weights = w,
                  verbose = FALSE)
  fe <- ilm_model(y ~ x, data = d, family = "poisson", verbose = FALSE)
  pw <- ilm_power(fw, n = c(60L, 120L), effect = 0.25, sims = 300L,
                  progress = FALSE)
  pe <- ilm_power(fe, n = c(60L, 120L), effect = 0.25, sims = 300L,
                  progress = FALSE)
  ## one weighted row is several observations, and a study of 60 draws 60 of
  ## them, not 60 rows of the aggregate
  expect_lt(max(abs(pw$power - pe$power)), 0.1)
})

test_that("zero parts, censoring, serial correlation and survival come along", {
  skip_on_cran()
  set.seed(27); dz <- data.frame(x = rnorm(300))
  dz$y <- ifelse(runif(300) < 0.25, 0, rpois(300, exp(0.6 + 0.3 * dz$x)))
  fz <- ilm_model(y ~ x, data = dz, family = "poisson", ziformula = ~ 1,
                  verbose = FALSE)
  pz <- ilm_power(fz, n = 300L, sims = 20L, progress = FALSE)
  expect_gt(pz$converged, 0.9)

  dc <- data.frame(x = rnorm(300))
  dc$y <- pmax(0.3 * dc$x + rnorm(300), -0.5)
  fc <- ilm_model(y ~ x, data = dc, family = "gaussian", verbose = FALSE,
                  censor = ilm_censor(dc$y, lower = -0.5))
  pc <- ilm_power(fc, n = 300L, sims = 20L, progress = FALSE)
  expect_gt(pc$converged, 0.9)

  ## eight times per unit and a clear random intercept: with five times the
  ## AR(1) absorbed the intercept, whose fitted SD went to 6e-5, and studies
  ## drawn from a fit at that boundary fail as often as the analysis would --
  ng <- 30L; nt <- 8L
  da <- expand.grid(t = seq_len(nt), id = factor(seq_len(ng)))
  da$x <- rep(rnorm(ng), each = nt)
  ## and a residual of its own, which the model has alongside the AR(1): data
  ## without one put the fitted residual SD at zero, where the likelihood is
  ## flat and a sixth of the refits could not be inverted
  da$y <- 0.3 * da$x + rep(rnorm(ng, 0, 0.8), each = nt) +
    unlist(lapply(seq_len(ng), function(i) stats::arima.sim(list(ar = 0.6), nt))) +
    rnorm(nrow(da), 0, 0.7)
  fa <- ilm_model(y ~ x + (1 | id), data = da, family = "gaussian",
                  ar = ilm_ar1(da$t, da$id, verbose = FALSE), verbose = FALSE)
  pa <- ilm_power(fa, n = 240L, sims = 10L, progress = FALSE)
  ## most of the refits must work, not nine in ten: this model is hard enough
  ## that the count moves with the platform's arithmetic -- Windows on R 4.6.1
  ## converged 8 of these 10, which "> 0.8" refused. At a true rate of 0.9,
  ## ten replicates come in at 8 or fewer 26% of the time and at 6 or fewer
  ## 1.3%, so seven is the bar
  expect_gte(pa$converged, 0.7)
  ## the replicate's AR(1) specification is rebuilt for the rows it drew
  r <- ilm_power_rows(240L, ilm_ar_group(fa$ar), seq_len(nrow(da)))
  a2 <- ilm_ar_rows(fa$ar, r$rows, r$copy)
  expect_identical(a2$Tt, fa$ar$Tt)
  expect_identical(a2$n_group, max(r$copy))

  ds <- data.frame(x = rnorm(200))
  tt <- stats::rweibull(200, 1.4, exp(1 - 0.4 * ds$x)); cs <- runif(200, 1, 6)
  ds$time <- pmin(tt, cs); ds$event <- as.integer(tt <= cs)
  fr <- ilm_model(time ~ x, data = ds, family = "rp", verbose = FALSE,
                  censor = ilm_surv(ds$time, ds$event))
  pr <- ilm_power(fr, n = 200L, sims = 10L, progress = FALSE)
  expect_gt(pr$converged, 0.8)
})

test_that("a model fitted to a complex sample is refused, with the way round it", {
  set.seed(28); d <- data.frame(x = rnorm(200), w = runif(200, 1, 3))
  d$y <- 0.5 * d$x + rnorm(200)
  f <- ilm_model(y ~ x, data = d, family = "gaussian", verbose = FALSE,
                 design = ilm_design(d, weights = ~ w))
  expect_error(ilm_power(f, sims = 5L), "design effect")
})
