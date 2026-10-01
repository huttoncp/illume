## New rows are placed among a group's fitted times with a tolerance of a
## millionth of the smallest gap between them -- never one scaled by the
## times' size, which merged date-times a second apart (see ilm_time_tol).

tm_walk <- function(t, id, seed) {
  set.seed(seed)
  d <- data.frame(id = id, t = t)
  lv <- unique(as.character(id))
  d$y <- stats::rnorm(nrow(d), 0, 0.3) +
    stats::rnorm(length(lv))[match(as.character(id), lv)] +
    stats::ave(stats::rnorm(nrow(d)), id, FUN = cumsum)
  d
}

test_that("date-times one second apart are each their own fitted time", {
  t0 <- as.POSIXct("2026-01-01 00:00:00", tz = "UTC")
  d <- tm_walk(t0 + rep(0:9, 3), factor(rep(c("a", "b", "c"), each = 10)), 3)
  f <- ilm_model(y ~ 1, data = d, family = "gaussian",
                 ar = ilm_rw1(~ t | id), verbose = FALSE)
  cl <- ilm_cells(f); a <- which(cl$group == "a")
  nd <- data.frame(t = t0 + c(1, 2, 0.5), id = "a")
  m <- ilm_matrices(f, nd)$ar
  expect_identical(m$cell, c(a[2], a[3], NA))
  ## half a second in falls between the first two times, not on either
  expect_identical(c(m$prev_cell[3], m$next_cell[3]), c(a[1], a[2]))
  expect_equal(c(m$dt_prev[3], m$dt_next[3]), c(0.5, 0.5))
  fv <- as.vector(ilm_fitted(f))
  expect_equal(as.vector(predict(f, nd[1:2, ], groups = "fitted")), fv[2:3],
               tolerance = 1e-12)
  expect_error(predict(f, nd[3, ], groups = "fitted"),
               "between two of its group's fitted times", fixed = TRUE)
})

test_that("an irregular series keeps a close pair of fitted times apart", {
  ## times near 1e7, a minute apart but for one pair 0.005 apart: a tolerance
  ## of 1e-9 of the times' size (0.01) merged the pair
  tt <- 1e7 + c(0, 60, 120, 120.005, 180, 240, 300, 360)
  d <- tm_walk(rep(tt, 3), factor(rep(c("a", "b", "c"), each = 8)), 11)
  f <- ilm_model(y ~ 1, data = d, family = "gaussian",
                 ar = ilm_car1(~ t | id), verbose = FALSE)
  cl <- ilm_cells(f); b <- which(cl$group == "b")
  nd <- data.frame(t = 1e7 + c(120, 120.005, 120.0025), id = "b")
  m <- ilm_matrices(f, nd)$ar
  expect_identical(m$cell, c(b[3], b[4], NA))
  expect_identical(c(m$prev_cell[3], m$next_cell[3]), c(b[3], b[4]))
  fv <- as.vector(ilm_fitted(f))
  expect_equal(as.vector(predict(f, nd[1:2, ], groups = "fitted")),
               fv[8 + 3:4], tolerance = 1e-12)
  ## the tolerance follows the smallest gap, not a typical one: a median gap
  ## of 1e4 would give 1e-2 and merge a pair 1e-3 apart
  tol <- illume:::ilm_time_tol(c(0, 1e4, 1e4 + 1e-3, 2e4))
  expect_lt(tol, 1e-3 / 2)
  expect_equal(tol, 1e-9, tolerance = 1e-6)
})

test_that("Dates are placed by day, and a group with one fitted time by it", {
  dd <- as.Date("2026-01-01") + c(0, 1, 2, 5, 6, 9)
  d <- tm_walk(c(rep(dd, 3), dd[3]),
               factor(c(rep(c("a", "b", "c"), each = 6), "s")), 5)
  f <- ilm_model(y ~ 1, data = d, family = "gaussian",
                 ar = ilm_rw1(~ t | id), verbose = FALSE)
  cl <- ilm_cells(f); a <- which(cl$group == "a"); s <- which(cl$group == "s")
  expect_length(s, 1L)
  nd <- data.frame(t = as.Date("2026-01-01") + c(5, 3, 10, 2, 3, 1),
                   id = c("a", "a", "a", "s", "s", "s"))
  m <- ilm_matrices(f, nd)$ar
  expect_identical(m$cell, c(a[4], NA, NA, s, NA, NA))
  expect_identical(c(m$prev_cell[2], m$next_cell[2]), c(a[3], a[4]))
  expect_identical(c(m$prev_cell[3], m$next_cell[3]), c(a[6], NA))
  ## one fitted time: a day after it, then a day before it
  expect_identical(c(m$prev_cell[5], m$next_cell[5]), c(s, NA))
  expect_identical(c(m$prev_cell[6], m$next_cell[6]), c(NA, s))
  expect_equal(c(m$dt_prev[5], m$dt_next[6]), c(1, 1))
  expect_identical(illume:::ilm_time_tol(cl$time[s]), 1e-9)
  fv <- as.vector(ilm_fitted(f))
  expect_equal(as.vector(predict(f, nd[c(1, 4), ], groups = "fitted")),
               fv[c(4, 19)], tolerance = 1e-12)
})
