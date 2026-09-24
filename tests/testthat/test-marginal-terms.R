## Averaging over EVERY random term, and saying which effect is reported.
##
## predict(marginal = TRUE) left an AR or CAR latent out of the average;
## ilm_scenario() averaged over each term's intercept alone, dropped a smooth's
## penalised part, and drew the fixed effects alone for its interval; ilm_ame()
## reported the effect for a group at zero while its help implied the
## population's. Each is checked here against a number computed independently.

pois_ar_data <- function(seed = 4, ng = 60, nt = 12) {
  set.seed(seed)
  d <- expand.grid(t = seq_len(nt), id = factor(seq_len(ng)))
  u <- rnorm(ng, 0, sqrt(0.2))[d$id]
  a <- unlist(lapply(seq_len(ng), function(i)
    as.numeric(stats::arima.sim(list(ar = 0.6), nt, sd = sqrt(0.7 * (1 - 0.36))))))
  d$x <- rnorm(nrow(d))
  d$y <- rpois(nrow(d), exp(0.2 + 0.3 * d$x + u + a))
  d
}

test_that("an AR latent is averaged over, not left out", {
  skip_on_cran()
  d <- pois_ar_data()
  f <- ilm_model(y ~ x + (1 | id), data = d, family = "poisson",
                 ar = ilm_ar1(d$t, d$id, verbose = FALSE), verbose = FALSE)
  expect_gt(f$Sigma[["ar"]][1, 1], 0.2)          # there really is AR variance
  pm <- as.numeric(predict(f, type = "response", marginal = TRUE))
  eta <- as.numeric(f$X %*% f$beta)
  ## a log link has the average in closed form, over both terms
  cf <- exp(eta + (f$Sigma[["id"]][1, 1] + f$Sigma[["ar"]][1, 1]) / 2)
  expect_equal(pm, cf, tolerance = 1e-6)
  ## the grouping term alone is a different, lower number -- the old answer
  expect_lt(mean(exp(eta + f$Sigma[["id"]][1, 1] / 2)), 0.9 * mean(pm))
  ## and the model's own simulations agree with the average
  ys <- ilm_simulate(f, 200L, 1L)
  expect_equal(mean(pm), mean(ys), tolerance = 0.03)
})

test_that("a CAR(1) latent is averaged over too", {
  skip_on_cran()
  d <- pois_ar_data(seed = 5, ng = 40, nt = 10)
  d$t <- d$t + stats::runif(nrow(d), 0, 0.5)     # irregular times
  f <- ilm_model(y ~ x + (1 | id), data = d, family = "poisson",
                 ar = ilm_car1(d$t, d$id, verbose = FALSE), verbose = FALSE)
  pm <- as.numeric(predict(f, type = "response", marginal = TRUE))
  eta <- as.numeric(f$X %*% f$beta)
  cf <- exp(eta + (f$Sigma[["id"]][1, 1] + f$Sigma[["ar"]][1, 1]) / 2)
  expect_equal(pm, cf, tolerance = 1e-6)
})

test_that("one linear predictor is averaged exactly, whatever ndraw or the seed", {
  set.seed(3); n <- 400
  d <- data.frame(x = rnorm(n), g = factor(rep(1:40, each = 10)))
  d$y <- rbinom(n, 1, plogis(-0.4 + 0.8 * d$x + rnorm(40, 0, 1.5)[d$g]))
  f <- ilm_model(y ~ x + (1 | g), data = d, family = "binomial", verbose = FALSE)
  gr <- data.frame(x = c(-2, 0, 2))
  pm <- as.numeric(predict(f, newdata = gr, type = "response", marginal = TRUE))
  s <- sqrt(f$Sigma[["g"]][1, 1])
  ex <- vapply(gr$x, function(z) {
    e <- sum(c(1, z) * stats::coef(f))
    stats::integrate(function(u) stats::plogis(e + u) * stats::dnorm(u, 0, s),
                     -Inf, Inf, rel.tol = 1e-10)$value
  }, 0)
  expect_equal(pm, ex, tolerance = 1e-7)
  expect_identical(pm, as.numeric(predict(f, newdata = gr, type = "response",
                                          marginal = TRUE, ndraw = 7, seed = 99)))
})

test_that("a multinomial outcome with an AR term averages over it by draws", {
  skip_on_cran()
  ## five observations per cell: with one, a multinomial AR latent is barely
  ## identified and its fitted variance is too small to test anything with
  set.seed(8); ng <- 30; nt <- 8
  cells <- expand.grid(t = seq_len(nt), id = factor(seq_len(ng)))
  cells$a <- unlist(lapply(seq_len(ng), function(i)
    as.numeric(stats::arima.sim(list(ar = 0.5), nt, sd = 1.5))))
  d <- cells[rep(seq_len(nrow(cells)), each = 5L), ]
  d$x <- rnorm(nrow(d))
  eta <- cbind(0, 0.5 * d$x + d$a, -0.3 * d$x + 0.5 * d$a)
  d$k <- factor(apply(exp(eta) / rowSums(exp(eta)), 1, function(p)
    sample(c("a", "b", "c"), 1, prob = p)))
  f <- suppressMessages(ilm_model(k ~ x, data = d, family = "multinomial",
                                  ar = ilm_ar1(d$t, d$id, verbose = FALSE),
                                  verbose = FALSE))
  gr <- data.frame(x = c(-1, 1))
  pm <- predict(f, newdata = gr, type = "response", marginal = TRUE,
                ndraw = 20000)
  ## brute force: the AR latent at a row is N(0, Sigma$ar), in the category
  ## dimensions, mapped through the sum-to-zero coding
  set.seed(1)
  A <- matrix(rnorm(2e5 * f$C), ncol = f$C) %*% ilm_msqrt(f$Sigma[["ar"]])
  Tc <- contr.sum(f$J)
  bf <- t(vapply(gr$x, function(z) {
    e <- as.numeric(c(1, z) %*% f$beta)
    E <- sweep(A, 2L, e, "+")
    colMeans(ilm_softmax_J(E, Tc))
  }, numeric(f$J)))
  expect_lt(max(abs(pm - bf)), 0.006)
  ## and leaving the AR term out -- the old answer -- is a different one
  pc <- predict(f, newdata = gr, type = "response", marginal = FALSE)
  expect_gt(max(abs(pm - pc)), 0.02)
})

test_that("a scenario averages over a random slope, not its intercept alone", {
  skip_on_cran()
  set.seed(21); ni <- 120; nt <- 5
  d <- expand.grid(time = 0:(nt - 1), id = factor(seq_len(ni)))
  d$y <- rbinom(nrow(d), 1, plogis(-0.3 + 0.35 * d$time +
    rnorm(ni, 0, 0.8)[as.integer(d$id)] +
    rnorm(ni, 0, 0.9)[as.integer(d$id)] * d$time))
  f <- ilm_model(y ~ time + (1 + time | id), data = d, family = "binomial",
                 verbose = FALSE)
  sc <- ilm_scenario(f, time = 4, sims = 200, progress = FALSE)
  nd <- f$model; nd$time <- 4
  pa <- mean(predict(f, newdata = nd, type = "response", marginal = TRUE))
  expect_equal(sc$estimate, pa, tolerance = 0.01)
  ## the intercept-only average is a different number
  V <- f$Sigma[[1]][1, 1] * f$Sigma_d[["id"]]
  e4 <- sum(c(1, 4) * stats::coef(f))
  io <- stats::integrate(function(u) stats::plogis(e4 + u) *
                           stats::dnorm(u, 0, sqrt(V[1, 1])), -Inf, Inf)$value
  expect_gt(abs(sc$estimate - io), 0.03)
  ## and the interval draws the variance components as well as the fixed effects
  set.seed(1)
  P <- ilm_scen_par_draws(f, 50L)
  expect_identical(ncol(P), length(stats::coef(f, full = TRUE)))
  expect_gt(stats::sd(P[, ncol(P)]), 0)
})

test_that("a scenario keeps a smooth's penalised part", {
  set.seed(2); n <- 400
  d <- data.frame(x = stats::runif(n), z = rnorm(n))
  d$y <- sin(2 * pi * d$x) + 0.3 * d$z + rnorm(n, sd = 0.4)
  f <- ilm_model(y ~ s(x) + z, data = d, family = "gaussian", verbose = FALSE)
  sc <- suppressWarnings(ilm_scenario(f, x = c(0.25, 0.75), sims = 100,
                                      progress = FALSE))
  pr <- vapply(c(0.25, 0.75), function(v) {
    nd <- f$model; nd$x <- v
    mean(predict(f, newdata = nd, type = "response"))
  }, 0)
  expect_equal(sc$estimate, pr, tolerance = 0.02)
  ## the curve's peak and trough, where it used to come out -1.31 and +1.72
  expect_equal(sc$estimate, c(1, -1), tolerance = 0.1)
})

test_that("ilm_ame() reports a typical group's effect, or the population's when asked", {
  set.seed(12); ng <- 80
  d <- data.frame(x = rnorm(ng * 10), g = factor(rep(seq_len(ng), each = 10)))
  d$y <- rbinom(nrow(d), 1, plogis(0.2 + 0.8 * d$x + rnorm(ng, 0, 1.5)[d$g]))
  f <- ilm_model(y ~ x + (1 | g), data = d, family = "binomial", verbose = FALSE)
  b <- stats::coef(f); eta <- b[1] + b[2] * d$x
  s <- sqrt(f$Sigma[["g"]][1, 1])
  dens <- function(e) stats::plogis(e) * (1 - stats::plogis(e))
  ## at a random effect of zero: the slope of the inverse link, averaged over rows
  cond <- mean(b[2] * dens(eta))
  ## over the random effects: the same slope, integrated
  pop <- mean(vapply(eta, function(e) b[2] * stats::integrate(function(u)
    dens(e + u) * stats::dnorm(u, 0, s), -Inf, Inf)$value, 0))
  a0 <- ilm_ame(f, "x")
  a1 <- ilm_ame(f, "x", marginal = TRUE)
  expect_equal(a0$estimate, unname(cond), tolerance = 1e-4)
  expect_equal(a1$estimate, unname(pop), tolerance = 1e-4)
  ## through a logit the population's effect is the flatter one
  expect_lt(a1$estimate, a0$estimate)
  expect_true(is.finite(a1$se) && a1$se > 0)
})
