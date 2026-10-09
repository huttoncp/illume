## A smooth is part of the fitted mean, so simulated data follow the fitted
## curve (review finding 3.2). Its coefficients used to be drawn afresh from
## the penalty's prior, which left only the straight-line part: fitting
## 2 sin(x), the mean of 500 simulations ran opposite to the fit.

test_that("simulated data from a smooth follow the fitted curve", {
  set.seed(6)
  n <- 300
  d <- data.frame(x = stats::runif(n, 0, 2 * pi))
  d$y <- 2 * sin(d$x) + stats::rnorm(n, 0, 0.5)
  f <- ilm_model(y ~ s(x), data = d, family = "gaussian", verbose = FALSE)
  fit <- as.numeric(predict(f))
  ys <- ilm_simulate(f, nsim = 200, seed = 1)
  expect_gt(stats::cor(rowMeans(ys), fit), 0.99)
  ## the mean of 200 draws is within a few Monte Carlo SEs of the fit
  expect_lt(max(abs(rowMeans(ys) - fit)), 5 * unname(f$dispersion) / sqrt(200))
  ## and a refit to one simulated dataset recovers the curve
  r <- ilm_refit(f, y = ys[, 1])
  expect_gt(stats::cor(as.numeric(predict(r)), fit), 0.95)
})

test_that("a Poisson smooth simulates counts around its fitted means", {
  set.seed(6)
  n <- 300
  d <- data.frame(x = stats::runif(n, 0, 2 * pi))
  d$yc <- stats::rpois(n, exp(1 + sin(d$x)))
  f <- ilm_model(yc ~ s(x), data = d, family = "poisson", verbose = FALSE)
  mu <- as.numeric(predict(f, type = "response"))
  ys <- ilm_simulate(f, nsim = 200, seed = 1)
  expect_equal(mean(ys), mean(mu), tolerance = 0.05)
  expect_gt(stats::cor(rowMeans(ys), mu), 0.98)
})
