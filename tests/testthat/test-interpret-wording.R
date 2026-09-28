# Small wording and type fixes in the interpretation, found by another agent.

wording_fit <- function() {
  set.seed(8); n <- 240
  d <- data.frame(x = stats::rnorm(n),
                  g = factor(sample(c("north", "central", "south"), n, TRUE)))
  d$y <- 1 + 0.3 * d$x + stats::rnorm(n)            # g has no effect
  ilm_model(y ~ x + g, data = d, family = "gaussian", verbose = FALSE)
}

test_that("a numeric predictor's values stay numbers", {
  f <- wording_fit()
  ap <- illume:::ilm_avg_pred(f, "x", c(-1, 1), intervals = FALSE)
  expect_type(ap$pred$value, "double")
  expect_identical(ap$pred$value, c(-1, 1))
  ## a factor's levels stay labels
  ag <- illume:::ilm_avg_pred(f, "g", c("north", "south"), intervals = FALSE)
  expect_type(ag$pred$value, "character")
})

test_that("three levels with no interval shown speak of no interval", {
  f <- wording_fit()
  lines <- unlist(ilm_interpret(f))
  gline <- lines[startsWith(lines, "g: ")]
  expect_length(gline, 1L)
  expect_false(grepl("The interval includes zero", gline, fixed = TRUE))
  expect_match(gline, "The data are consistent with no effect", fixed = TRUE)
})
