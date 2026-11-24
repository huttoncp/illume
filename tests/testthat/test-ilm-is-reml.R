## ilm_is_reml(): whether a fit was made by restricted maximum likelihood. A
## package that must not compare likelihoods across fits needs to ask, and
## fit$reml, which it read until now, is not documented.

ir_data <- function(seed = 1) {
  set.seed(seed)
  d <- data.frame(g = factor(rep(1:10, each = 6)), x = stats::rnorm(60))
  d$y <- 1 + 0.5 * d$x + stats::rnorm(10)[d$g] + stats::rnorm(60)
  d
}
ir_fit <- function(...) suppressWarnings(suppressMessages(ilm_model(...)))

test_that("it is TRUE for a REML fit and FALSE for a maximum-likelihood one", {
  d <- ir_data()
  f_ml <- ir_fit(y ~ x + (1 | g), data = d, family = "gaussian", verbose = FALSE)
  f_re <- ir_fit(y ~ x + (1 | g), data = d, family = "gaussian", reml = TRUE,
                 verbose = FALSE)
  expect_identical(ilm_is_reml(f_ml), FALSE)
  expect_identical(ilm_is_reml(f_re), TRUE)
  ## and it says what the fit holds
  expect_identical(ilm_is_reml(f_ml), isTRUE(f_ml$reml))
  expect_identical(ilm_is_reml(f_re), isTRUE(f_re$reml))
})

test_that("it follows the default: REML for a gaussian dispersion model with nothing random", {
  d <- ir_data()
  f <- ir_fit(y ~ x, data = d, family = "gaussian", dispformula = ~ x,
              verbose = FALSE)
  expect_identical(ilm_is_reml(f), TRUE)
  f0 <- ir_fit(y ~ x, data = d, family = "gaussian", dispformula = ~ x,
               reml = FALSE, verbose = FALSE)
  expect_identical(ilm_is_reml(f0), FALSE)
  ## a model with a random effect is by maximum likelihood unless asked
  fr <- ir_fit(y ~ x + (1 | g), data = d, family = "gaussian", dispformula = ~ x,
               verbose = FALSE)
  expect_identical(ilm_is_reml(fr), FALSE)
})

test_that("a refit by maximum likelihood is FALSE, and the REML fit is unchanged", {
  d <- ir_data()
  f <- ir_fit(y ~ x + (1 | g), data = d, family = "gaussian", reml = TRUE,
              verbose = FALSE)
  m <- suppressMessages(illume:::ilm_as_ml(f, "a test"))
  expect_identical(ilm_is_reml(m), FALSE)
  expect_identical(ilm_is_reml(f), TRUE)
})

test_that("a fit that never recorded it gives FALSE, and anything else is refused", {
  d <- ir_data()
  f <- ir_fit(y ~ x + (1 | g), data = d, family = "gaussian", reml = TRUE,
              verbose = FALSE)
  f$reml <- NULL
  expect_identical(ilm_is_reml(f), FALSE)
  expect_error(ilm_is_reml(1), "ilm_model")
  expect_error(ilm_is_reml(d), "data.frame")
  expect_error(ilm_is_reml(NULL), "NULL")
})

test_that("ilm_data() refuses anything but a fit", {
  expect_error(ilm_data(1), "ilm_model")
  expect_error(ilm_data(data.frame(x = 1)), "data.frame")
})
