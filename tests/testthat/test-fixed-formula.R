## formula(terms(fit)) is the fixed-effects formula of a fit: its random-effect
## terms and smooths left out, its ordinary terms, offset and removed intercept
## kept. It is documented under ?"ilm_model-accessors" in place of a function of
## its own, so these tests hold it to what the help says.

ff_data <- function(seed = 1) {
  set.seed(seed)
  n <- 120
  d <- data.frame(g = factor(sample(c("a", "b", "c"), n, TRUE)),
                  t = stats::runif(n), z = stats::runif(n),
                  id = factor(rep(1:12, each = 10)),
                  e = stats::runif(n, 1, 3))
  d$y <- 1 + 0.5 * d$t + sin(2 * pi * d$z) + 0.2 * as.integer(d$g) +
    stats::rnorm(12)[d$id] + stats::rnorm(n)
  d$k <- stats::rpois(n, exp(0.3 + 0.4 * d$t) * d$e)
  d
}
ff_fit <- function(...) suppressWarnings(suppressMessages(ilm_model(...)))
ff_txt <- function(f) paste(deparse(f), collapse = " ")

test_that("random terms and smooths are left out, ordinary terms kept", {
  skip_if_not_installed("mgcv")
  d <- ff_data()
  f <- ff_fit(y ~ g + poly(t, 2) + s(z) + (1 | id), data = d, family = "gaussian",
              verbose = FALSE)
  expect_identical(ff_txt(formula(f)), "y ~ g + poly(t, 2) + s(z) + (1 | id)")
  fx <- formula(terms(f))
  expect_identical(ff_txt(fx), "y ~ g + poly(t, 2)")
  ## the same formula the fit keeps as its fixed part, with the same environment
  expect_identical(fx, f$fixed_formula)
  expect_identical(environment(fx), environment(f$fixed_formula))
  ## the smooth's unpenalised column is a coefficient, not a term
  expect_true("s(z).f1" %in% names(coef(f)))
  expect_false(any(grepl("s(z)", attr(terms(f), "term.labels"), fixed = TRUE)))
  expect_identical(attr(terms(f), "term.labels"), c("g", "poly(t, 2)"))
})

test_that("an offset is kept", {
  d <- ff_data()
  f <- ff_fit(k ~ g + offset(log(e)) + (1 | id), data = d, family = "poisson",
              verbose = FALSE)
  fx <- formula(terms(f))
  expect_identical(ff_txt(fx), "k ~ g + offset(log(e))")
  expect_identical(fx, f$fixed_formula)
})

test_that("a removed intercept and an intercept-only fit are kept as written", {
  d <- ff_data()
  f0 <- ff_fit(y ~ g + t - 1, data = d, family = "gaussian", verbose = FALSE)
  expect_identical(ff_txt(formula(terms(f0))), "y ~ g + t - 1")
  expect_identical(formula(terms(f0)), f0$fixed_formula)
  f1 <- ff_fit(y ~ 1 + (1 | id), data = d, family = "gaussian", verbose = FALSE)
  expect_identical(ff_txt(formula(terms(f1))), "y ~ 1")
  expect_identical(formula(terms(f1)), f1$fixed_formula)
})

test_that("a dispersion formula is not part of it", {
  d <- ff_data()
  f <- ff_fit(y ~ g, data = d, family = "gaussian", dispformula = ~ t,
              verbose = FALSE)
  expect_identical(ff_txt(formula(terms(f))), "y ~ g")
  expect_identical(formula(terms(f)), f$fixed_formula)
})
