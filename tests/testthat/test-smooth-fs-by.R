## A factor-smooth interaction, s(x, id, bs = "fs"), with a numeric `by` used
## to fit the model without the `by` -- the same fitted values to the last
## digit -- under the name "s(x,id):z". It is refused until it is handled.

fsby_data <- function(seed = 10) {
  set.seed(seed)
  id <- factor(rep(1:20, each = 10))
  x <- stats::runif(200); z <- stats::runif(200, 0.5, 1.5)
  y <- sin(2 * pi * x) * z + stats::rnorm(20, 0, 0.6)[id] + stats::rnorm(200, 0, 0.4)
  data.frame(y = y, x = x, z = z, id = id)
}

test_that("an fs smooth with a numeric by is refused, saying why", {
  d <- fsby_data()
  expect_error(ilm_model(y ~ s(x, id, bs = "fs", k = 5, by = z), data = d,
                         family = "gaussian", verbose = FALSE),
               "illume does not handle a `by` variable on a factor-smooth interaction",
               fixed = TRUE)
  expect_error(ilm_model(y ~ s(x, id, bs = "fs", k = 5, by = z), data = d,
                         family = "gaussian", verbose = FALSE),
               "`z` was dropped without a word", fixed = TRUE)
})

test_that("an fs smooth with a factor by is still refused as one smooth per level", {
  d <- fsby_data()
  d$h <- factor(ifelse(as.integer(d$id) %% 2 == 0, "a", "b"))
  expect_error(ilm_model(y ~ h + s(x, id, bs = "fs", k = 5, by = h), data = d,
                         family = "gaussian", verbose = FALSE),
               "makes 2 smooths, one for each level of the factor `h`", fixed = TRUE)
})

test_that("an fs smooth without by, and an ordinary smooth with a numeric by, still fit", {
  d <- fsby_data()
  expect_s3_class(suppressMessages(ilm_model(y ~ s(x, id, bs = "fs", k = 5), data = d,
                                             family = "gaussian", verbose = FALSE)),
                  "ilm_model")
  expect_s3_class(suppressMessages(ilm_model(y ~ s(x, by = z), data = d,
                                             family = "gaussian", verbose = FALSE)),
                  "ilm_model")
})
