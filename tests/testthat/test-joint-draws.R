## Joint draws of the fixed and random parameters, behind predict()'s
## intervals for a model with a smooth. They were centred on the wrong
## values whenever the family has a parameter declared after the random
## effects, and every interval missed its own fit. mgcv is the outside check.

smooth_data <- function(n, seed) {
  set.seed(seed)
  d <- data.frame(x = runif(n))
  d$f <- sin(2 * pi * d$x)
  d
}
at <- data.frame(x = c(0.05, 0.25, 0.5, 0.75, 0.95))

test_that("a gaussian smooth's intervals bracket its fit and its SEs match mgcv", {
  skip_if_not_installed("mgcv")
  d <- smooth_data(500, 1); d$y <- d$f + rnorm(nrow(d), sd = 0.5)
  f <- ilm_model(y ~ s(x), data = d, family = "gaussian", verbose = FALSE)
  p <- predict(f, newdata = at, se.fit = TRUE, interval = "confidence")
  expect_true(p$joint)
  expect_true(all(p$lower < p$fit & p$fit < p$upper))
  m <- mgcv::gam(y ~ s(x), data = d, method = "REML")
  pm <- predict(m, newdata = at, se.fit = TRUE)
  expect_equal(as.numeric(p$fit), as.numeric(pm$fit), tolerance = 0.05)
  expect_equal(as.numeric(p$se.fit), as.numeric(pm$se.fit), tolerance = 0.15)
})

test_that("a negative binomial smooth's intervals bracket its fit on the response scale", {
  skip_if_not_installed("mgcv")
  d <- smooth_data(600, 2)
  d$y <- MASS::rnegbin(nrow(d), mu = exp(1 + d$f), theta = 3)
  f <- ilm_model(y ~ s(x), data = d, family = "nbinom", verbose = FALSE)
  p <- predict(f, newdata = at, type = "response", se.fit = TRUE,
               interval = "confidence")
  expect_true(all(p$lower < p$fit & p$fit < p$upper))
  m <- mgcv::gam(y ~ s(x), data = d, family = mgcv::nb(), method = "REML")
  pm <- predict(m, newdata = at, type = "response", se.fit = TRUE)
  expect_equal(as.numeric(p$fit), as.numeric(pm$fit), tolerance = 0.05)
  ## the response-scale SEs were 0.02 to 16 times mgcv's
  expect_equal(as.numeric(p$se.fit), as.numeric(pm$se.fit), tolerance = 0.2)
})

test_that("each block of the draws is centred on its own estimate", {
  d <- smooth_data(300, 3); d$y <- d$f + rnorm(nrow(d), sd = 0.5)
  for (reml in c(FALSE, TRUE)) {
    f <- ilm_model(y ~ s(x), data = d, family = "gaussian", reml = reml,
                   verbose = FALSE)
    jd <- ilm_joint_draws(f, 4000L, 1L)
    est <- c(f$opt$par, f$sdr$par.random)
    for (nm in unique(jd$which)) {
      m <- rowMeans(jd$draws[jd$which == nm, , drop = FALSE])
      s <- apply(jd$draws[jd$which == nm, , drop = FALSE], 1L, stats::sd)
      ## within five Monte Carlo standard errors of the estimate, by name
      expect_true(all(abs(m - est[names(est) == nm]) < 5 * s / sqrt(4000)),
                  info = paste(nm, if (reml) "(REML)" else "(ML)"))
    }
  }
})
