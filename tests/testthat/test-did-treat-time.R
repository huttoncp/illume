## `treat_time` is one of the time column's own values (review finding 5.5).
## A non-numeric time is worked on as the codes of its ordered values, and
## treat_time used to be compared with those codes: a character year "2005"
## sorted after code 3 and started the treatment at the wrong period, silently
## ("treat_time used: 3"), while 2005 or a Date stopped with a misleading
## "not separable" error.

did_years <- function() {
  set.seed(1)
  d <- expand.grid(unit = 1:40, time = 1:8)
  d$treated <- as.integer(d$unit <= 20)
  d$y <- 1 + 0.3 * d$time + stats::rnorm(40)[d$unit] +
         0.8 * d$treated * (d$time >= 5) + stats::rnorm(nrow(d))
  d$year <- as.character(2000 + d$time)
  d$date <- as.Date(sprintf("%d-01-01", 1999 + d$time))
  d
}

test_that("a character or Date time takes treat_time in its own values", {
  d <- did_years()
  ref <- ilm_did(d, "y", "unit", "time", treated = "treated", treat_time = 5,
                 verbose = FALSE)
  for (tt in list(2005, "2005")) {
    f <- ilm_did(d, "y", "unit", "year", treated = "treated", treat_time = tt,
                 verbose = FALSE)
    expect_equal(f$att$estimate, ref$att$estimate, tolerance = 1e-8)
    expect_identical(f$treat_time, "2005")
    expect_equal(f$n_pre, 4L)
  }
  for (tt in list(as.Date("2004-01-01"), "2004-01-01")) {
    f <- ilm_did(d, "y", "unit", "date", treated = "treated", treat_time = tt,
                 verbose = FALSE)
    expect_equal(f$att$estimate, ref$att$estimate, tolerance = 1e-8)
    expect_identical(f$treat_time, as.Date("2004-01-01"))
  }
  expect_message(ilm_did(d, "y", "unit", "year", treated = "treated",
                         treat_time = "2005"),
                 "treatment starts at year = 2005", fixed = TRUE)
})

test_that("a factor time is ordered by its levels, not alphabetically", {
  d <- did_years()
  mo <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug")
  d$month <- factor(mo[d$time], levels = mo)
  ref <- ilm_did(d, "y", "unit", "time", treated = "treated", treat_time = 5,
                 verbose = FALSE)
  f <- ilm_did(d, "y", "unit", "month", treated = "treated", treat_time = "May",
               verbose = FALSE)
  expect_equal(f$att$estimate, ref$att$estimate, tolerance = 1e-8)
  expect_identical(f$treat_time, "May")
})

test_that("a treat_time that is not a period, or leaves no change, is refused", {
  d <- did_years()
  expect_error(ilm_did(d, "y", "unit", "year", treated = "treated",
                       treat_time = 2010, verbose = FALSE),
               "`treat_time` = 2010 is not one of `year`'s values (2001, 2002",
               fixed = TRUE)
  expect_error(ilm_did(d, "y", "unit", "time", treated = "treated",
                       treat_time = "5", verbose = FALSE),
               "`treat_time` must be a number, as `time` is", fixed = TRUE)
  expect_error(ilm_did(d, "y", "unit", "time", treated = "treated",
                       treat_time = 9, verbose = FALSE),
               "leaves no period at or after it", fixed = TRUE)
  expect_error(ilm_did(d, "y", "unit", "year", treated = "treated",
                       treat_time = "2001", verbose = FALSE),
               "leaves no period before it", fixed = TRUE)
})

test_that("with a per-row treatment, the start and first periods are the time's own", {
  d <- did_years()
  d$tr <- as.integer(d$treated == 1L & d$time >= 5)
  f <- ilm_did(d, "y", "unit", "year", treatment = "tr", verbose = FALSE)
  expect_identical(f$treat_time, "2005")
  expect_setequal(unique(stats::na.omit(f$first_treat)), "2005")
})
