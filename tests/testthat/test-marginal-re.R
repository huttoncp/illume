## Averaging over the random effects
##
## `marginal = TRUE` computes E_u[g^-1(eta + z'u)]. Two things make that easy to
## get wrong in a way nothing complains about:
##
##   1. A random effect is a dk x C MATRIX per group, not one number. Drawing
##      only the intercept variance and adding it as a shift is right when
##      dk == 1 and silently wrong otherwise -- it drops the slope variance AND
##      holds the shift constant across rows when the true spread grows with
##      distance from where the slope is centred.
##   2. Every answer is a probability in a plausible range either way.
##
## So these tests check the INTEGRAL against a brute-force average computed
## independently, rather than checking that the function returns something.

skip_slow <- function() skip_on_cran()

test_that("a random intercept integrates to the right number", {
  ## this case was already correct; the test is that it stayed correct
  skip_slow()
  set.seed(3); n <- 600
  d <- data.frame(x = rnorm(n), g = factor(rep(1:60, each = 10)))
  d$y <- rbinom(n, 1, plogis(-0.4 + 0.8 * d$x +
                             rnorm(60, 0, 1.2)[as.integer(d$g)]))
  f <- ilm_model(y ~ x + (1 | g), data = d, family = "binomial",
                 verbose = FALSE)
  gr <- data.frame(x = seq(-2, 2, length.out = 5))
  pm <- predict(f, newdata = gr, type = "response", marginal = TRUE,
                ndraw = 20000)
  set.seed(99)
  u <- rnorm(200000, 0, sqrt(f$Sigma[[1]][1, 1]))
  bf <- vapply(gr$x, function(z)
    mean(plogis(sum(c(1, z) * stats::coef(f)) + u)), 0)
  expect_equal(as.numeric(pm), bf, tolerance = 0.005)
  ## and it is a genuinely different curve from the conditional one
  pc <- predict(f, newdata = gr, type = "response", marginal = FALSE)
  expect_lt(diff(range(pm)), diff(range(pc)))
})

test_that("a random SLOPE is integrated, not flattened to an intercept", {
  skip_slow()
  set.seed(21); ni <- 120; nt <- 5
  d <- expand.grid(time = 0:(nt - 1), id = factor(seq_len(ni)))
  d$y <- rbinom(nrow(d), 1, plogis(-0.3 + 0.35 * d$time +
    rnorm(ni, 0, 0.8)[as.integer(d$id)] +
    rnorm(ni, 0, 0.9)[as.integer(d$id)] * d$time))
  f <- ilm_model(y ~ time + (1 + time | id), data = d, family = "binomial",
                 verbose = FALSE)
  V <- f$Sigma[[1]][1, 1] * f$Sigma_d[["id"]]
  expect_gt(sqrt(V[2, 2]), 0.3)              # there really is slope variance

  g <- data.frame(time = 0:4)
  pm <- predict(f, newdata = g, type = "response", marginal = TRUE,
                ndraw = 4000)
  set.seed(7)
  U <- matrix(rnorm(4e5 * 2), 4e5, 2) %*% ilm_msqrt(V)
  bf <- vapply(g$time, function(tt) mean(plogis(
    sum(c(1, tt) * stats::coef(f)) + U[, 1] + tt * U[, 2])), 0)
  expect_equal(as.numeric(pm), bf, tolerance = 0.008)

  ## the intercept-only average is a DIFFERENT number, so passing the check
  ## above is not something the old behaviour could also have done
  set.seed(7)
  u1 <- rnorm(4e5, 0, sqrt(V[1, 1]))
  old <- vapply(g$time, function(tt)
    mean(plogis(sum(c(1, tt) * stats::coef(f)) + u1)), 0)
  expect_gt(max(abs(bf - old)), 0.05)
  expect_gt(abs(pm[5] - old[5]), 0.05)
})

test_that("an identity link is averaged exactly, not simulated", {
  ## E[eta + z'u] = eta when u has mean zero, so there is an exact answer and
  ## simulating it would return a noisy estimate of a number already known --
  ## with 200 draws and a random slope that noise reached 0.14
  set.seed(5)
  d <- expand.grid(time = 0:3, id = factor(1:80))
  d$y <- 2 + 0.5 * d$time + rnorm(80, 0, 1.5)[as.integer(d$id)] +
    rnorm(80, 0, 0.7)[as.integer(d$id)] * d$time + rnorm(nrow(d), 0, 1)
  f <- ilm_model(y ~ time + (1 + time | id), data = d, family = "gaussian",
                 verbose = FALSE)
  g <- data.frame(time = 0:3)
  expect_identical(
    predict(f, newdata = g, type = "response", marginal = TRUE),
    predict(f, newdata = g, type = "response", marginal = FALSE))
  ## and it does not depend on ndraw, because nothing is drawn
  expect_identical(
    predict(f, newdata = g, type = "response", marginal = TRUE, ndraw = 5),
    predict(f, newdata = g, type = "response", marginal = TRUE, ndraw = 5000))
})

test_that("a bar over a column the prediction data lacks is said out loud", {
  ## averaging over part of a term is a different quantity from averaging over
  ## the term, so it cannot be done quietly
  skip_slow()
  set.seed(31)
  d <- expand.grid(time = 0:3, id = factor(1:60))
  d$x <- rnorm(nrow(d))
  d$y <- rbinom(nrow(d), 1, plogis(-0.2 + 0.5 * d$x +
    rnorm(60, 0, 0.8)[as.integer(d$id)] +
    rnorm(60, 0, 0.6)[as.integer(d$id)] * d$time))
  f <- ilm_model(y ~ x + (1 + time | id), data = d, family = "binomial",
                 verbose = FALSE)
  expect_warning(predict(f, newdata = data.frame(x = c(-1, 0, 1)),
                         type = "response", marginal = TRUE),
                 "integrates its intercept")
})

test_that("the draw factorisation is shared with the power simulator", {
  ## ilm_power_draw() and predict(marginal = TRUE) integrate the same
  ## distribution for different purposes. They went out of step once already;
  ## this pins that they read it from one place.
  set.seed(12)
  d <- expand.grid(time = 0:3, id = factor(1:50))
  d$y <- 1 + 0.4 * d$time + rnorm(50, 0, 1.1)[as.integer(d$id)] +
    rnorm(50, 0, 0.5)[as.integer(d$id)] * d$time + rnorm(nrow(d))
  f <- ilm_model(y ~ time + (1 + time | id), data = d, family = "gaussian",
                 verbose = FALSE)
  fac <- ilm_re_factors(f, 1L)
  expect_identical(fac$d, 2L)
  ## A A' = Sigma_d and B'B = Sigma, so (Zb %*% A Z B) has the covariance the
  ## objective gave the term
  expect_equal(fac$A %*% t(fac$A), f$Sigma_d[["id"]], tolerance = 1e-10)
  expect_equal(t(fac$B) %*% fac$B, f$Sigma[[1]], tolerance = 1e-10)
  Zb <- ilm_re_design(f, 1L, d, fac$d)
  expect_identical(dim(Zb), c(nrow(d), 2L))
  expect_true(all(Zb[, 1] == 1))
  expect_identical(as.numeric(Zb[, 2]), as.numeric(d$time))
})

test_that("a random intercept still gives a design of ones", {
  set.seed(2)
  d <- data.frame(x = rnorm(200), g = factor(rep(1:20, each = 10)))
  d$y <- d$x + rnorm(20)[as.integer(d$g)] + rnorm(200)
  f <- ilm_model(y ~ x + (1 | g), data = d, family = "gaussian",
                 verbose = FALSE)
  fac <- ilm_re_factors(f, 1L)
  expect_identical(fac$d, 1L)
  expect_identical(fac$A, matrix(1, 1L, 1L))
  expect_identical(dim(ilm_re_design(f, 1L, d, 1L)), c(200L, 1L))
})
