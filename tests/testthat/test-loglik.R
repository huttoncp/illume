## The constants in logLik().

test_that("a binomial with trials leaves out only the binomial coefficients", {
  ## glm() includes sum(lchoose(trials, successes)); illume does not. It is a
  ## constant of the data, and the help says so -- this holds it to exactly
  ## that, so a difference of any other size is noticed.
  set.seed(3)
  d <- data.frame(n = rep(c(5, 10, 20), 20), x = stats::rnorm(60))
  d$s <- stats::rbinom(60, d$n, stats::plogis(-0.5 + 0.8 * d$x))
  d$p <- d$s / d$n
  f <- ilm_model(p ~ x, data = d, family = "binomial", weights = n,
                 verbose = FALSE)
  g <- stats::glm(p ~ x, data = d, family = stats::binomial, weights = n)
  expect_equal(as.numeric(logLik(g)) - as.numeric(logLik(f)),
               sum(lchoose(d$n, d$s)), tolerance = 1e-6)
})
