## Refitting a model fitted inside someone's own function.
##
## Called through a wrapper's `...`, ilm_model()'s recorded call held ..1,
## ..2 -- and a refit to other data, which rebuilds the model from that call,
## stopped with "..4 used in an incorrect context" as soon as the model had a
## correlation over time. Reported by another agent whose cross-validation
## refits per fold. The record is repaired, and the correlation over time is
## rebuilt from the fit itself.

dots_data <- function() {
  set.seed(1)
  d <- data.frame(g = factor(rep(1:8, each = 10)), t = rep(1:10, 8), x = rnorm(80))
  d$y <- rpois(80, exp(0.4 + 0.3 * d$x))
  d
}

## what a refit should equal: the same model fitted directly to the rows
direct <- function(rows, ar) {
  d <- dots_data()[rows, ]
  suppressMessages(suppressWarnings(
    ilm_model(y ~ x + (1 | g), data = d, family = "poisson", ar = ar, verbose = FALSE)))
}

test_that("a model fitted through a wrapper's ... refits to other data", {
  d <- dots_data()
  wrap <- function(...) ilm_model(..., verbose = FALSE)
  for (ctor in list(ilm_ar1, ilm_car1, ilm_rw1)) {
    f <- suppressMessages(suppressWarnings(
      wrap(y ~ x + (1 | g), data = d, family = "poisson", ar = ctor(~ t | g))))
    ## the record names what the user wrote, not the wrapper's ..N
    expect_false(any(grepl("^[.][.][0-9]+$", vapply(as.list(f$call)[-1L],
                                                     function(a) deparse(a)[1L], ""))))
    r <- suppressMessages(suppressWarnings(ilm_refit(f, data = d[1:60, ])))
    expect_s3_class(r, "ilm_model")
    expect_identical(r$ar$type, f$ar$type)
    expect_equal(coef(r), coef(direct(1:60, ctor(~ t | g))), tolerance = 1e-6)
  }
})

test_that("without a correlation over time it refits as it did", {
  d <- dots_data()
  wrap <- function(...) ilm_model(..., verbose = FALSE)
  f <- suppressMessages(wrap(y ~ x + (1 | g), data = d, family = "poisson"))
  r <- suppressMessages(ilm_refit(f, data = d[1:60, ]))
  expect_equal(coef(r), coef(suppressMessages(
    ilm_model(y ~ x + (1 | g), data = d[1:60, ], family = "poisson", verbose = FALSE))),
    tolerance = 1e-6)
})

test_that("purrr-style calls and a wrapper of a wrapper refit too", {
  d <- dots_data()
  ## a function over a list of data sets, as lapply() or purrr::map() call it
  fs <- lapply(list(d), function(z) suppressMessages(suppressWarnings(
    ilm_model(y ~ x + (1 | g), data = z, family = "poisson",
              ar = ilm_ar1(~ t | g), verbose = FALSE))))
  r <- suppressMessages(suppressWarnings(ilm_refit(fs[[1L]], data = d[1:60, ])))
  expect_equal(coef(r), coef(direct(1:60, ilm_ar1(~ t | g))), tolerance = 1e-6)
  ## a wrapper whose own arguments carry the formula and data: those names
  ## mean nothing outside it, and the refit does not need them
  wrap <- function(...) ilm_model(..., verbose = FALSE)
  wrap2 <- function(fml, dat, ...) wrap(fml, data = dat, family = "poisson", ...)
  f2 <- suppressMessages(suppressWarnings(wrap2(y ~ x + (1 | g), d, ar = ilm_car1(~ t | g))))
  r2 <- suppressMessages(suppressWarnings(ilm_refit(f2, data = d[1:60, ])))
  expect_equal(coef(r2), coef(direct(1:60, ilm_car1(~ t | g))), tolerance = 1e-6)
})
