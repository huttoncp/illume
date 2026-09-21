sel_df <- function()
  data.frame(id = 1:20, score_a = rnorm(20), score_b = rnorm(20),
             other = rnorm(20), label = sample(letters[1:3], 20, TRUE),
             stringsAsFactors = FALSE)

test_that("NULL takes everything eligible", {
  d <- sel_df()
  expect_setequal(ilm_resolve_cols(d, NULL), names(d))
  expect_setequal(ilm_resolve_cols(d, NULL, exclude = "id"),
                  setdiff(names(d), "id"))
  expect_setequal(ilm_resolve_cols(d, NULL, eligible = c("id", "other")),
                  c("id", "other"))
})

test_that("a character vector is taken as names", {
  d <- sel_df()
  expect_equal(ilm_resolve_cols(d, c("score_a", "other")),
               c("score_a", "other"))
  ## order follows what was asked for, not the frame
  expect_equal(ilm_resolve_cols(d, c("other", "id")), c("other", "id"))
  ## an excluded column drops out silently; asking for ONLY excluded ones errors
  expect_equal(ilm_resolve_cols(d, c("id", "other"), exclude = "id"), "other")
  expect_error(ilm_resolve_cols(d, "id", exclude = "id"), "only columns that")
})

test_that("a single unknown string is read as a pattern", {
  d <- sel_df()
  expect_equal(ilm_resolve_cols(d, "^score_"), c("score_a", "score_b"))
  expect_setequal(ilm_resolve_cols(d, "_"), c("score_a", "score_b"))
  expect_equal(ilm_resolve_cols(d, "^l"), "label")
  ## a name beats a pattern when both would work: `other` is a column, so it is
  ## that column and not a regex matching it
  expect_equal(ilm_resolve_cols(d, "other"), "other")
  ## matching nothing is an error naming both readings, since a typo and a
  ## pattern that misses look identical from here
  expect_error(ilm_resolve_cols(d, "^nope"), "not found")
  expect_error(ilm_resolve_cols(d, "^nope"), "as a pattern")
})

test_that("a predicate function selects by column content", {
  d <- sel_df()
  expect_setequal(ilm_resolve_cols(d, is.numeric),
                  c("id", "score_a", "score_b", "other"))
  expect_equal(ilm_resolve_cols(d, is.character), "label")
  expect_setequal(ilm_resolve_cols(d, function(v) is.numeric(v) && max(v) < 10),
                  c("score_a", "score_b", "other"))
  ## a predicate that throws on a column counts as not matching it, rather
  ## than taking the whole call down with it
  expect_error(ilm_resolve_cols(d, function(v) stop("boom")),
               "matched no column")
  expect_equal(ilm_resolve_cols(d, function(v) is.numeric(v) || stop("boom")),
               c("id", "score_a", "score_b", "other"))
  expect_error(ilm_resolve_cols(d, function(v) FALSE), "matched no column")
})

test_that("bad input is refused with a useful message", {
  d <- sel_df()
  expect_error(ilm_resolve_cols(d, 1:3), "must be column names")
  expect_error(ilm_resolve_cols(d, list("a")), "must be column names")
  ## several names, some unknown: say which, and what exists
  expect_error(ilm_resolve_cols(d, c("score_a", "nope", "nah")),
               "nope, nah")
  expect_error(ilm_resolve_cols(d, c("score_a", "nope")), "Available")
  ## the argument name in the message is the caller's, not this function's
  expect_error(ilm_resolve_cols(d, function(v) FALSE, arg = "by"), "`by`")
})

test_that("selection is a value, so it can be held in a variable", {
  ## the whole point of not using non-standard evaluation
  d <- sel_df()
  pick <- "^score_"
  expect_equal(ilm_resolve_cols(d, pick), c("score_a", "score_b"))
  pick2 <- is.numeric
  expect_true("id" %in% ilm_resolve_cols(d, pick2))
  picks <- list("^score_", is.character, c("id", "other"))
  got <- lapply(picks, function(p) ilm_resolve_cols(d, p))
  expect_equal(lengths(got), c(2L, 1L, 2L))
})

test_that("the functions that take cols accept all four forms", {
  d <- sel_df()
  nm <- function(x) sort(unique(as.character(x$variable)))
  expect_equal(nm(ilm_outliers_all(d, cols = "^score_", flagged_only = FALSE)),
               c("score_a", "score_b"))
  expect_equal(nm(ilm_outliers_all(d, cols = is.numeric, flagged_only = FALSE)),
               c("id", "other", "score_a", "score_b"))
  expect_equal(nm(ilm_outliers_all(d, cols = c("other"), flagged_only = FALSE)),
               "other")
  expect_gt(length(nm(ilm_outliers_all(d, flagged_only = FALSE))), 1L)
  ## and `by` goes through the same resolver
  expect_true("label" %in% names(ilm_outliers_all(d, by = "label",
                                                  flagged_only = FALSE)))

  skip_if_not_installed("PCAmixdata")
  expect_setequal(ilm_reduce(d, cols = "^score_")$cols,
                  c("score_a", "score_b"))
  expect_setequal(ilm_reduce(d, cols = is.numeric)$cols,
                  c("id", "score_a", "score_b", "other"))
})
