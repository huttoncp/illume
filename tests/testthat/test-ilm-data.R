## ilm_data(): the rows a fit was made from, with the variables its
## transformed terms are built from added back. A package that places new
## rows on a fit -- limn, telling which rows are past the data -- needs each
## row's raw time even when the time is only inside a function, as in
## poly(t, 2).

idata <- function(seed = 1) {
  set.seed(seed)
  d <- data.frame(t = rep(1:12, 10), g = factor(rep(1:10, each = 12)))
  d$y <- 0.2 * d$t - 0.01 * d$t^2 + stats::rnorm(10)[d$g] + stats::rnorm(120)
  d
}
qfit <- function(...) suppressWarnings(suppressMessages(ilm_model(...)))

test_that("a time inside poly() comes back raw, on the fit's own rows", {
  d <- idata(); d$y[c(5, 40)] <- NA            # dropped rows: the rows must line up
  f <- qfit(y ~ poly(t, 2) + (1 | g), data = d, family = "gaussian", verbose = FALSE)
  out <- ilm_data(f)
  expect_s3_class(out, "data.frame")
  expect_identical(nrow(out), nrow(f$model))
  expect_true("t" %in% names(out))
  expect_false("t" %in% names(f$model))         # the frame alone does not have it
  kept <- d[!is.na(d$y), ]
  expect_identical(out$t, kept$t)
  expect_identical(rownames(out), rownames(f$model))
  ## the frame's own columns are unchanged
  for (nm in names(f$model)) expect_identical(out[[nm]], f$model[[nm]])
})

test_that("with no extra columns kept, it returns the frame alone", {
  ## data given as a list, not a data frame: the fit keeps no extra columns
  d <- idata()
  f <- qfit(y ~ poly(t, 2), data = as.list(d), family = "gaussian", verbose = FALSE)
  expect_null(f$data_extra)
  expect_identical(ilm_data(f), f$model)
  ## a term with nothing transformed needs nothing added
  f2 <- qfit(y ~ t + (1 | g), data = d, family = "gaussian", verbose = FALSE)
  expect_null(f2$data_extra)
  expect_identical(ilm_data(f2), f2$model)
})

test_that("a correlation term named by its columns brings its time and group", {
  d <- idata()
  f <- qfit(y ~ poly(t, 2), data = d, family = "gaussian",
            ar = ilm_ar1(~ t | g), verbose = FALSE)
  out <- ilm_data(f)
  expect_true(all(c("t", "g") %in% names(out)))
  expect_identical(out$t, d$t)
  expect_identical(out$g, d$g)
})
