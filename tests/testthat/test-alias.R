## The whole fixed design, the smooths' unpenalised columns included, is
## checked for rank before the fit (item 4, as Craig ruled). Each of these
## used to fit and come back with a failed Hessian and NaN standard errors,
## or a FAIL after the fit, with nothing naming the overlap.

alias_data <- function(n = 300) {
  set.seed(1)
  d <- data.frame(x = stats::runif(n), z = stats::rnorm(n), w = stats::rnorm(n),
                  a = factor(sample(c("p", "q", "r"), n, TRUE)),
                  b = factor(sample(c("u", "v"), n, TRUE)))
  d$x1 <- stats::rnorm(n); d$x2 <- 2 * d$x1
  d$y <- sin(3 * d$x) + 0.5 * d$z + stats::rnorm(n, 0, 0.3)
  d
}
fit_or_msg <- function(f, d) tryCatch({
  suppressWarnings(suppressMessages(ilm_model(f, data = d, family = "gaussian", verbose = FALSE)))
  "fitted"
}, error = function(e) conditionMessage(e))

test_that("two smooths with the same numeric by stop, naming the overlap", {
  m <- fit_or_msg(y ~ s(x, by = z) + s(w, by = z), alias_data())
  expect_match(m, "not all separable", fixed = TRUE)
  expect_match(m, "the unpenalised part of s(w):z", fixed = TRUE)
  expect_match(m, "same numeric `by` variable", fixed = TRUE)
})

test_that("a parametric term inside a smooth's null space stops", {
  m <- fit_or_msg(y ~ x + s(x), alias_data())
  expect_match(m, "the unpenalised part of s(x) is a linear combination of the intercept and `x`",
               fixed = TRUE)
  m2 <- fit_or_msg(y ~ poly(z, 2) + s(x, by = z), alias_data())
  expect_match(m2, "not all separable", fixed = TRUE)
})

test_that("collinear ordinary columns and an empty interaction cell are named", {
  m <- fit_or_msg(y ~ x1 + x2, alias_data())
  expect_match(m, "`x2` is a linear combination of `x1`", fixed = TRUE)
  d <- alias_data(); d <- d[!(d$a == "r" & d$b == "v"), ]
  m2 <- fit_or_msg(y ~ a * b, d)
  expect_match(m2, "is zero in every row", fixed = TRUE)
  expect_match(m2, "combinations with no rows (r:v)", fixed = TRUE)
})

test_that("a design with full rank is untouched, and PR 21's case keeps its own words", {
  expect_identical(fit_or_msg(y ~ w + z:w + s(x, by = z), alias_data()), "fitted")
  expect_identical(fit_or_msg(y ~ x1 + z + s(x), alias_data()), "fitted")
  ## a column in dollars is not mistaken for a dependent one
  d <- alias_data(); d$big <- stats::rnorm(nrow(d)) * 1e6
  expect_identical(fit_or_msg(y ~ big + z, d), "fitted")
  m <- fit_or_msg(y ~ z + s(x, by = z), alias_data())
  expect_match(m, "is a term of the model and also the `by` variable", fixed = TRUE)
})
