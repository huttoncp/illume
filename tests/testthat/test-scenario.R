scen_data <- function(n = 600L, seed = 1L, binary = TRUE) {
  set.seed(seed)
  d <- data.frame(dose = runif(n, 0, 20), age = rnorm(n, 50, 15))
  d$y <- if (binary)
    rbinom(n, 1L, stats::plogis(-3 + 0.12 * d$dose + 0.05 * d$age))
  else 2 + 0.15 * d$dose + 0.05 * d$age + rnorm(n)
  d
}

test_that("under an identity link the two averagings agree exactly", {
  d <- scen_data(binary = FALSE)
  f <- ilm_model(y ~ dose + age, data = d, family = "gaussian",
                 verbose = FALSE)
  a <- ilm_scenario(f, dose = c(0, 10, 20), sims = 300L, progress = FALSE)
  b <- ilm_scenario(f, dose = c(0, 10, 20), over = "reference", sims = 300L,
                    progress = FALSE)
  expect_equal(a$estimate, b$estimate, tolerance = 1e-10)
})

test_that("under a non-linear link they differ, and each is what it claims", {
  d <- scen_data()
  f <- ilm_model(y ~ dose + age, data = d, family = "binomial",
                 verbose = FALSE)
  a <- ilm_scenario(f, dose = 10, sims = 600L, progress = FALSE)
  b <- ilm_scenario(f, dose = 10, over = "reference", sims = 600L,
                    progress = FALSE)
  expect_false(isTRUE(all.equal(a$estimate, b$estimate)))
  ## the sample version is standardisation: set dose for EVERY unit, keep each
  ## unit's own age, predict, average
  nd <- d; nd$dose <- 10
  hand <- mean(stats::plogis(as.numeric(
    stats::model.matrix(~ dose + age, nd) %*% stats::coef(f))))
  expect_equal(a$estimate, hand, tolerance = 0.01)
  ## the reference version is the prediction at the average covariate
  nr <- data.frame(dose = 10, age = mean(d$age))
  hand2 <- stats::plogis(as.numeric(
    stats::model.matrix(~ dose + age, nr) %*% stats::coef(f)))
  expect_equal(b$estimate, hand2, tolerance = 0.01)
  ## and they are not the same number
  expect_gt(abs(hand - hand2), 0.01)
  expect_equal(attr(a, "over"), "sample")
  expect_output(print(b), "NOT the")
})

test_that("the grid, the intervals and the contrasts behave", {
  d <- scen_data()
  f <- ilm_model(y ~ dose + age, data = d, family = "binomial",
                 verbose = FALSE)
  s <- ilm_scenario(f, dose = c(0, 10, 20), age = c(40, 60), sims = 400L,
                    contrast = "first", progress = FALSE)
  expect_s3_class(s, "ilm_scenario")
  expect_equal(nrow(s), 6L)                     # every combination
  expect_true(all(c("dose", "age", "estimate", "lower", "upper") %in% names(s)))
  expect_true(all(s$lower < s$estimate & s$estimate < s$upper))
  expect_true(all(s$estimate > 0 & s$estimate < 1))
  ## a probability rises with dose at fixed age
  sub <- s[s$age == 40, ]
  expect_false(is.unsorted(sub$estimate))
  ct <- attr(s, "contrasts")
  expect_equal(nrow(ct), 5L)                    # each against the first
  expect_true(all(ct$lower < ct$estimate & ct$estimate < ct$upper))
  ## the estimate is the difference at the fit; the interval is computed
  ## WITHIN each draw, so it is not built from two independently-summarised
  ## numbers
  dr <- attr(s, "draws")
  expect_equal(ct$estimate[1], s$estimate[2] - s$estimate[1], tolerance = 1e-12)
  expect_equal(ct$lower[1], unname(stats::quantile(dr[, 2] - dr[, 1], 0.025)),
               tolerance = 1e-12)
  ## pairwise gives every pair
  s2 <- ilm_scenario(f, dose = c(0, 10, 20), sims = 200L,
                     contrast = "pairwise", progress = FALSE)
  expect_equal(nrow(attr(s2, "contrasts")), 3L)
  expect_null(attr(ilm_scenario(f, dose = 10, sims = 100L,
                                progress = FALSE), "contrasts"))
})

test_that("a value outside the observed range is reported", {
  d <- scen_data()
  f <- ilm_model(y ~ dose + age, data = d, family = "binomial",
                 verbose = FALSE)
  s <- ilm_scenario(f, dose = c(10, 60), sims = 200L, progress = FALSE)
  expect_false(s$outside[1])
  expect_true(s$outside[2])
  expect_output(print(s), "outside the observed")
  ## a round number just past the edge is NOT flagged: the smallest observed
  ## dose is about 0.01, and calling dose = 0 an extrapolation every time is
  ## how a warning gets ignored
  s0 <- ilm_scenario(f, dose = c(0, 20), sims = 200L, progress = FALSE)
  expect_false(any(s0$outside))
  ## a factor level that was never seen is refused by model.matrix, and a
  ## level that was is fine
  dd <- d; dd$g <- factor(sample(c("a", "b"), nrow(d), TRUE))
  fg <- ilm_model(y ~ dose + g, data = dd, family = "binomial",
                  verbose = FALSE)
  expect_silent(ilm_scenario(fg, g = "b", sims = 100L, progress = FALSE))
})

test_that("an impossible COMBINATION of possible values is reported", {
  ## each value inside its own range, the pair in a corner nobody occupies
  set.seed(5); n <- 500L
  age <- rnorm(n, 45, 8)
  d <- data.frame(age = age,
                  service = pmax(0, age - 22 + rnorm(n, 0, 2)))
  d$y <- 1 + 0.02 * d$age + 0.05 * d$service + rnorm(n)
  f <- ilm_model(y ~ age + service, data = d, family = "gaussian",
                 verbose = FALSE)
  s <- ilm_scenario(f, age = c(45, 25), service = c(23, 35), sims = 200L,
                    progress = FALSE)
  ## a 25-year-old with 35 years of service: both values ordinary
  expect_gt(25, min(d$age)); expect_lt(35, max(d$service))
  i <- which(s$age == 25 & s$service == 35)
  expect_true(s$outside[i])
  ## while the sensible combination is not flagged
  j <- which(s$age == 45 & s$service == 23)
  expect_false(s$outside[j])
  sup <- attr(s, "support")
  expect_gt(sup[[i]]$distance, sup[[i]]$typical)
  expect_gt(sup[[i]]$distance, sup[[j]]$distance)
  expect_output(print(s), "COMBINATION")
})

test_that("a mixed model is averaged over the random effects", {
  set.seed(6); ng <- 50L; ni <- 12L; n <- ng * ni
  d <- data.frame(id = factor(rep(seq_len(ng), each = ni)),
                  dose = runif(n, 0, 20))
  d$y <- rbinom(n, 1L, stats::plogis(-2 + 0.12 * d$dose +
                                       rep(rnorm(ng, 0, 1.5), each = ni)))
  f <- ilm_model(y ~ dose + (1 | id), data = d, family = "binomial",
                 verbose = FALSE)
  s <- ilm_scenario(f, dose = c(0, 20), sims = 200L, progress = FALSE)
  expect_true(attr(s, "mixed"))
  expect_output(print(s), "population mean")
  ## the projection brackets what was actually seen at low and high dose
  lo <- mean(d$y[d$dose < 7]); hi <- mean(d$y[d$dose > 13])
  expect_lt(s$estimate[1], lo)
  expect_gt(s$estimate[2], hi)
  ## a population mean is pulled towards 0.5 relative to the value for a
  ## cluster whose random effect is zero
  cond <- stats::plogis(as.numeric(c(1, 20) %*% stats::coef(f)))
  expect_lt(s$estimate[2], cond)
})

test_that("it says what it is holding fixed, and refuses what it cannot do", {
  d <- scen_data(400L)
  f <- ilm_model(y ~ dose + age, data = d, family = "binomial",
                 verbose = FALSE)
  s <- ilm_scenario(f, dose = 10, sims = 100L, progress = FALSE)
  expect_equal(attr(s, "held"), "age")
  expect_output(print(s), "held as observed: age")
  expect_error(ilm_scenario(f, sims = 10L), "name the predictors")
  expect_error(ilm_scenario(f, nope = 1, sims = 10L), "not in the model")
  expect_error(ilm_scenario(mtcars, dose = 1), "must be a fitted ilm_model")
  set.seed(7); n <- 300L
  dm <- data.frame(x = rnorm(n))
  dm$y <- factor(sample(c("a", "b", "c"), n, TRUE))
  fm <- ilm_model(y ~ x, data = dm, family = "multinomial", verbose = FALSE)
  expect_error(ilm_scenario(fm, x = 0, sims = 10L), "vector of probabilities")
})

test_that("a zero-inflated fit projects the response, zeros included", {
  set.seed(8); n <- 600L
  d <- data.frame(x = runif(n, 0, 5), z = rnorm(n))
  d$y <- ifelse(runif(n) < stats::plogis(-0.5 + 0.8 * d$z), 0,
                rpois(n, exp(0.3 + 0.4 * d$x)))
  f <- ilm_model(y ~ x, data = d, family = "poisson", ziformula = ~ z,
                 verbose = FALSE)
  s <- ilm_scenario(f, x = c(0, 5), sims = 200L, progress = FALSE)
  expect_true(all(s$estimate > 0))
  expect_gt(s$estimate[2], s$estimate[1])
  ## the projection is the mean of the RESPONSE, so it sits below the count
  ## part's own mean -- the zeros are in it
  cnt <- exp(as.numeric(c(1, 5) %*% stats::coef(f)))
  expect_lt(s$estimate[2], cnt)
})
