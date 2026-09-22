## ilm_check_ar() and ilm_variogram() must line their columns up with the
## rows the fit kept, the way ilm_check_covariate() and ilm_check_omitted()
## do.
##
## Found from outside the package: a fit on data with any missing value drops
## those rows from the model frame, and a column handed in straight from the
## data frame then has more values than the model has rows. Both functions
## stopped, telling the caller to "subset `time` the same way" -- work they
## cannot easily do and the fit already knows how to.
##
## Most real data has a gap somewhere, so this was every call on real data.

make_gappy <- function(n_id = 25L, n_t = 6L, na_frac = 0.08, seed = 4L) {
  set.seed(seed)
  d <- expand.grid(t = seq_len(n_t), id = seq_len(n_id))
  d$id <- factor(sprintf("u%02d", d$id))
  u <- rnorm(n_id, 0, 1.2)[as.integer(d$id)]
  d$y <- 2 + 0.3 * d$t + u + rnorm(nrow(d))
  d$y[sample.int(nrow(d), max(1L, round(nrow(d) * na_frac)))] <- NA
  d
}

test_that("ilm_check_ar accepts columns straight from the data frame", {
  d <- make_gappy()
  fit <- ilm_model(y ~ t + (1 | id), data = d, family = "gaussian",
                   verbose = FALSE)
  expect_lt(nrow(fit$X), nrow(d))          # rows really were dropped

  expect_no_error(
    suppressMessages(ilm_check_ar(fit, time = d$t, group = d$id,
                                  maxlag = 3L, B = 3L, verbose = FALSE)))
})

test_that("ilm_variogram accepts columns straight from the data frame", {
  d <- make_gappy()
  fit <- ilm_model(y ~ t + (1 | id), data = d, family = "gaussian",
                   verbose = FALSE)
  ## B is tiny because this is testing that the call is accepted, not the
  ## verdict; illume rightly warns that no verdict is reachable at B = 3.
  expect_no_error(
    suppressWarnings(suppressMessages(
      ilm_variogram(fit, time = d$t, group = d$id,
                    breaks = 3L, B = 3L, plot = FALSE, verbose = FALSE))))
})

test_that("an already-aligned column still works", {
  ## The pre-existing calling convention must not break.
  d <- make_gappy()
  fit <- ilm_model(y ~ t + (1 | id), data = d, family = "gaussian",
                   verbose = FALSE)
  expect_no_error(
    suppressMessages(ilm_check_ar(fit, time = fit$model$t,
                                  group = fit$model$id,
                                  maxlag = 3L, B = 3L, verbose = FALSE)))
})

test_that("aligning gives the same answer as subsetting by hand", {
  d <- make_gappy()
  fit <- ilm_model(y ~ t + (1 | id), data = d, family = "gaussian",
                   verbose = FALSE)
  keep <- !is.na(d$y)
  a <- suppressMessages(ilm_check_ar(fit, time = d$t, group = d$id,
                                     maxlag = 3L, B = 5L, seed = 2L,
                                     verbose = FALSE))
  b <- suppressMessages(ilm_check_ar(fit, time = d$t[keep], group = d$id[keep],
                                     maxlag = 3L, B = 5L, seed = 2L,
                                     verbose = FALSE))
  expect_equal(a$table, b$table)
})

test_that("a column of the wrong length is still refused, with a reason", {
  d <- make_gappy()
  fit <- ilm_model(y ~ t + (1 | id), data = d, family = "gaussian",
                   verbose = FALSE)
  expect_error(
    suppressMessages(ilm_check_ar(fit, time = d$t[1:5], group = d$id[1:5],
                                  maxlag = 3L, B = 3L, verbose = FALSE)),
    "value")
})
