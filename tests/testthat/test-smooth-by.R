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

test_that("a smooth with a numeric by is named as mgcv names it, apart from s(x)", {
  ## both used to be called "s(x)": the second overwrote the first's
  ## penalised part, the columns shared names, and predict() on new rows
  ## stopped on a column count
  set.seed(11); n <- 600
  d <- data.frame(x = stats::runif(n), z = stats::rnorm(n, 1, 0.5),
                  g = factor(sample(1:30, n, TRUE)))
  d$y <- 0.5 + d$z * sin(2 * pi * d$x) +
    stats::rnorm(30, 0, 0.4)[as.integer(d$g)] + stats::rnorm(n, 0, 0.5)
  f <- ilm_model(y ~ s(x) + s(x, by = z) + (1 | g), data = d,
                 family = "gaussian", verbose = FALSE)
  expect_identical(names(stats::coef(f)),
                   c("(Intercept)", "s(x).f1", "s(x):z.f1", "s(x):z.f2"))
  expect_identical(names(f$smooths), c("s(x)", "s(x):z"))
  expect_setequal(names(f$Sigma), c("s(x)", "s(x):z", "g"))
  ## the same model with the by-smooth on a copy of x, where nothing collided
  d$x_ <- d$x
  f2 <- ilm_model(y ~ s(x) + s(x_, by = z) + (1 | g), data = d,
                  family = "gaussian", verbose = FALSE)
  expect_equal(as.numeric(stats::logLik(f)), as.numeric(stats::logLik(f2)),
               tolerance = 1e-7)
  nd <- d[1:8, ]
  expect_equal(unname(stats::predict(f, nd, type = "link")),
               unname(stats::predict(f2, nd, type = "link")), tolerance = 1e-5)
  ## two smooths with one name stop rather than overwrite
  expect_error(ilm_model(y ~ s(x) + s(x, k = 5), data = d,
                         family = "gaussian", verbose = FALSE),
               "two smooths are both called 's(x)'", fixed = TRUE)
})
