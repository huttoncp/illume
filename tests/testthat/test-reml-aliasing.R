## Under REML the coefficients are integrated out, and the parameter_aliasing
## check read OK from the variance parameters alone, never looking at the
## coefficients. It now takes their covariance from the joint precision.

q <- function(e) suppressMessages(suppressWarnings(e))
al <- function(f) f$checks[f$checks$check == "parameter_aliasing", ]

near_data <- function(seed = 4, fam = "binomial") {
  set.seed(seed); A <- 30; n <- 600
  d <- data.frame(area = factor(sample(A, n, TRUE)), x1 = rnorm(n))
  d$x2 <- d$x1 + rnorm(n, 0, 0.02)                  # nearly the same column
  u <- rnorm(A, 0, 0.7)[d$area]
  eta <- -0.2 + 0.4 * d$x1 + u
  d$y <- if (fam == "binomial") rbinom(n, 1, plogis(eta)) else eta + rnorm(n)
  d
}

test_that("a near-aliased pair is caught under REML as under ML", {
  for (fam in c("binomial", "gaussian")) {
    d <- near_data(fam = fam)
    for (reml in c(FALSE, TRUE)) {
      f <- q(ilm_model(y ~ x1 + x2 + (1 | area), data = d, family = fam, reml = reml, verbose = FALSE))
      a <- al(f)
      expect_identical(a$status, "FAIL", info = paste(fam, reml))
      expect_match(a$detail, "x1 <-> x2", fixed = TRUE, info = paste(fam, reml))
    }
  }
})

test_that("an ordinary REML fit's check reads its coefficients", {
  set.seed(5); A <- 30; n <- 600
  d <- data.frame(area = factor(sample(A, n, TRUE)), x = rnorm(n), z = rnorm(n))
  d$y <- rbinom(n, 1, plogis(-0.3 + 0.5 * d$x + 0.2 * d$z + rnorm(A, 0, 0.8)[d$area]))
  f <- q(ilm_model(y ~ x + z + (1 | area), data = d, family = "binomial", reml = TRUE, verbose = FALSE))
  a <- al(f)
  expect_identical(a$status, "OK")
  expect_false(grepl("theta <-> theta", a$detail, fixed = TRUE))
  ## the pair named is a real pair of parameters, among the coefficients or
  ## the variance
  expect_match(a$detail, "[(](\\(Intercept\\)|x|z|area:L\\[1,1\\]) <-> ")
  ## the line says what REML leaves out
  expect_match(a$detail, "under REML, correlations between coefficients and variance parameters are not assessed",
               fixed = TRUE)
  ## maximum likelihood, unchanged: the joint covariance of everything, and
  ## no such caveat
  f0 <- q(ilm_model(y ~ x + z + (1 | area), data = d, family = "binomial", verbose = FALSE))
  expect_identical(al(f0)$status, "OK")
  expect_match(al(f0)$detail, "largest |parameter correlation|", fixed = TRUE)
  expect_false(grepl("under REML", al(f0)$detail, fixed = TRUE))
})

test_that("a fit with one parameter to correlate says so", {
  ## REML, intercept only, the area variance at zero: the fit reduces to the
  ## plain model and its one coefficient is all there is
  set.seed(2); A <- 30; n <- 600
  d <- data.frame(area = factor(sample(A, n, TRUE)))
  d$y <- rbinom(n, 1, 0.4)
  f <- q(ilm_model(y ~ 1 + (1 | area), data = d, family = "binomial", reml = TRUE, verbose = FALSE))
  expect_identical(f$hessian_how, "reduced")
  a <- al(f)
  expect_identical(a$status, "OK")
  expect_match(a$detail, "one parameter only ((Intercept)): nothing to correlate", fixed = TRUE)
  ## and where the variance is not at zero, the pair is a real pair, never a
  ## parameter against itself
  set.seed(4)
  d <- data.frame(area = factor(sample(A, n, TRUE))); d$y <- rbinom(n, 1, 0.4)
  f <- q(ilm_model(y ~ 1 + (1 | area), data = d, family = "binomial", reml = TRUE, verbose = FALSE))
  expect_match(al(f)$detail, "largest |parameter correlation| = 0.000 ((Intercept) <-> area:L[1,1])",
               fixed = TRUE)
})
