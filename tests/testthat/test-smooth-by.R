## A smooth with a numeric `by` is not centred, so its unpenalised part spans
## the by variable: s(x, by = z) already carries z's main effect.

sb_data <- function(seed = 11) {
  set.seed(seed); n <- 300
  d <- data.frame(x = stats::runif(n), z = stats::rnorm(n, 1, 0.5))
  d$y <- 0.5 + d$z * sin(2 * pi * d$x) + stats::rnorm(n, 0, 0.5)
  d
}

test_that("a numeric by variable beside its own smooth is refused, saying why", {
  ## it used to fit with a failed Hessian and every standard error NaN
  d <- sb_data()
  expect_error(ilm_model(y ~ z + s(x, by = z), data = d, family = "gaussian",
                         verbose = FALSE),
               "already contains z's main effect", fixed = TRUE)
  expect_error(ilm_model(y ~ z + s(x, by = z), data = d, family = "gaussian",
                         verbose = FALSE),
               "Drop `z` from the formula; s(x, by = z) carries its effect",
               fixed = TRUE)
  ## without it the model is identified, and every standard error is finite
  f <- ilm_model(y ~ s(x, by = z), data = d, family = "gaussian",
                 verbose = FALSE)
  expect_true(all(is.finite(sqrt(diag(vcov(f))))))
  ## a term that is not the by variable is untouched
  d$w <- stats::rnorm(nrow(d))
  expect_s3_class(ilm_model(y ~ w + s(x, by = z), data = d,
                            family = "gaussian", verbose = FALSE), "ilm_model")
})
