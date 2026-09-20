# counts, copies, wash_df, translate and recode_errors.

test_that("counts are right and ordering works", {
  y <- factor(c(rep("a", 5), rep("b", 3), "c"))
  r <- ilm_counts(y)
  expect_equal(as.character(r$value), c("a", "b", "c"))
  expect_equal(r$n, c(5L, 3L, 1L))
  expect_equal(as.character(ilm_counts(y, order = "a")$value), c("c", "b", "a"))
  expect_equal(as.character(ilm_counts(y, order = "i")$value), c("a", "b", "c"))
})

test_that("the n argument limits the rows", {
  # documented but ignored in the previous implementation
  y <- factor(sample(letters[1:5], 200, TRUE))
  expect_equal(nrow(ilm_counts(y, n = 2)), 2L)
  expect_equal(nrow(ilm_counts(y, n = "all")), 5L)
  expect_error(ilm_counts(y, n = 0), "positive whole number")
})

test_that("na.rm controls whether missing values are counted", {
  y <- c("a", "a", NA)
  expect_equal(sum(ilm_counts(y)$n), 2L)
  expect_equal(sum(ilm_counts(y, na.rm = FALSE)$n), 3L)
})

test_that("counts_tb shows both ends", {
  y <- factor(c(rep("a", 9), rep("b", 5), rep("c", 2), "d"))
  r <- ilm_counts_tb(y, n = 2)
  expect_equal(as.character(r$top_value), c("a", "b"))
  expect_equal(as.character(r$bot_value), c("d", "c"))
})

test_that("counts_all stacks variables and respects by", {
  d <- ilm_sim(n_id = 20)
  r <- ilm_counts_all(d[, c("grp", "site")])
  expect_true(all(c("variable", "value", "n") %in% names(r)))
  expect_setequal(unique(r$variable), c("grp", "site"))
  rb <- ilm_counts_all(d[, c("grp", "flag")], by = "flag")
  expect_true("flag" %in% names(rb))
  expect_error(ilm_counts_all(d, by = "nope"), "not found")
  expect_error(ilm_counts(d), "use `ilm_counts_all\\(\\)`")
})

test_that("copies reports copy number and group size", {
  d <- data.frame(a = c(1, 1, 2, 2, 2, 3), b = c("x", "x", "y", "y", "z", "w"),
                  stringsAsFactors = FALSE)
  r <- ilm_copies(d)
  expect_equal(r$n_copies, c(2L, 2L, 2L, 2L, 1L, 1L))
  expect_equal(r$copy_number, c(1L, 2L, 1L, 2L, 1L, 1L))
  expect_equal(nrow(ilm_copies(d, filter = "dupes")), 4L)
  expect_equal(nrow(ilm_copies(d, filter = "first")), 4L)
  expect_equal(nrow(ilm_copies(d, filter = "unique")), 2L)
})

test_that("copies can key on a subset of columns", {
  d <- data.frame(a = c(1, 1, 2, 2, 2, 3), b = c("x", "x", "y", "y", "z", "w"),
                  stringsAsFactors = FALSE)
  r <- ilm_copies(d, "a")
  expect_equal(max(r$n_copies), 3L)
  expect_error(ilm_copies(d, "nope"), "not found")
})

test_that("wash_df cleans names, drops empties and retypes", {
  m <- data.frame("Col One" = c("1", "2", "", ""),
                  someFlag = c("TRUE", "FALSE", "", ""),
                  txt = c("a", "b", "", ""),
                  allEmpty = c("", "", "", ""),
                  check.names = FALSE, stringsAsFactors = FALSE)
  w <- ilm_wash_df(m)
  expect_equal(names(w), c("col_one", "some_flag", "txt"))
  expect_equal(nrow(w), 2L)
  expect_type(w$col_one, "integer")
  expect_type(w$some_flag, "logical")
  expect_type(w$txt, "character")
})

test_that("a single unconvertible value blocks coercion", {
  # otherwise one stray "n/a" silently turns a column into NAs
  d <- data.frame(x = c("1", "2", "n/a"), stringsAsFactors = FALSE)
  expect_type(ilm_wash_df(d)$x, "character")
  expect_false(anyNA(ilm_wash_df(d)$x))
})

test_that("translate maps through the dictionary and flags mismatches", {
  expect_equal(ilm_translate(c(1, 2, 3, 99), old = 1:3,
                             new = c("a", "b", "c")),
               c("a", "b", "c", NA))
  expect_equal(ilm_translate(c(1, 99), old = 1, new = 10, default = -1),
               c(10, -1))
  expect_error(ilm_translate(1:3, old = 1:2, new = 1:3), "same length")
  expect_warning(ilm_translate(1:3, old = c(1, 1), new = c("a", "b")),
                 "duplicate")
})

test_that("recode_errors handles vectors, factors and restricted columns", {
  expect_equal(ilm_recode_errors(c(1, 2, 999), errors = 999), c(1, 2, NA))
  expect_equal(ilm_recode_errors(c(1, 999), errors = 999, replacement = 0),
               c(1, 0))
  f <- ilm_recode_errors(factor(c("a", "b", "unknown")), errors = "unknown")
  expect_true(is.na(f[3]))
  expect_false("unknown" %in% levels(f))
  d <- data.frame(x = c(1, 999), y = c(999, 2))
  r <- ilm_recode_errors(d, errors = 999, cols = "x")
  expect_true(is.na(r$x[2]))
  expect_equal(r$y[1], 999)             # untouched
})

test_that("recode_errors errors name the right argument for the input type", {
  d <- data.frame(x = c(1, 999))
  expect_error(ilm_recode_errors(d, errors = 999, cols = "z"), "Available")
  expect_error(ilm_recode_errors(d, errors = 999, ind = 1), "`ind` applies to vector")
  expect_error(ilm_recode_errors(c(1, 2), errors = 1, rows = 1), "apply to data frame")
  expect_error(ilm_recode_errors(d, errors = 999, rows = 99), "out of range")
  expect_error(ilm_recode_errors(d, errors = character(0)), "at least one value")
})
