## REML on a model with nothing random and no dispersion parameter. Under
## REML the fixed effects are integrated out with the random effects, and
## such a model then has no parameter left to estimate. It used to stop with
## the optimiser's own words, "'d' must be a nonempty numeric vector", for a
## Poisson or a flexible-parametric fit; it is now refused with a reason. A
## model that keeps a scale or a variance parameter -- gaussian, Weibull --
## still fits by REML, and so does one with a random effect.

rnr_data <- function(seed = 31, n = 400) {
  set.seed(seed)
  s <- data.frame(x = stats::runif(n), g = factor(sample(1:20, n, TRUE)))
  tt <- exp(1 + 1.2 * s$x + 0.5 * log(stats::rexp(n)))
  ct <- stats::rexp(n, 1 / (2 * stats::median(tt)))
  s$time <- pmin(tt, ct); s$event <- as.integer(tt <= ct)
  s$yg <- 1 + s$x + stats::rnorm(n)
  s$yp <- stats::rpois(n, exp(0.2 + 0.5 * s$x + stats::rnorm(20, 0, 0.3)[as.integer(s$g)]))
  s
}
q <- function(e) suppressMessages(suppressWarnings(e))

test_that("REML with nothing random and no dispersion parameter is refused, with a reason", {
  s <- rnr_data()
  expect_error(q(ilm_model(yp ~ x, data = s, family = "poisson", reml = TRUE, verbose = FALSE)),
               "nothing for REML to estimate.*reml = FALSE")
  expect_error(q(ilm_model(time ~ x, data = s, family = "rp", reml = TRUE,
                           censor = ilm_surv(s$time, s$event), verbose = FALSE)),
               "nothing for REML to estimate.*reml = FALSE")
})

test_that("REML still fits a model that keeps a scale or variance parameter, or a random effect", {
  s <- rnr_data()
  expect_s3_class(q(ilm_model(yg ~ x, data = s, family = "gaussian", reml = TRUE, verbose = FALSE)),
                  "ilm_model")
  expect_s3_class(q(ilm_model(time ~ x, data = s, family = "weibull", reml = TRUE,
                              censor = ilm_surv(s$time, s$event), verbose = FALSE)), "ilm_model")
  f <- q(ilm_model(yp ~ x + (1 | g), data = s, family = "poisson", reml = TRUE, verbose = FALSE))
  expect_s3_class(f, "ilm_model")
  expect_true(isTRUE(f$reml))
  ## and by ML the refused models fit
  expect_s3_class(q(ilm_model(yp ~ x, data = s, family = "poisson", reml = FALSE, verbose = FALSE)),
                  "ilm_model")
})
