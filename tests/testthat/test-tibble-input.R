## A tibble is a data frame, and has to behave like one
##
## `ilm_reduce()` handed `X.quanti` straight to PCAmix. PCAmix checks its
## columns with `is.numeric(X.quanti[, j])`, and `[` on a TIBBLE does not drop
## to a vector -- so every numeric column looked non-numeric and the call died
## on "All variables in X.quanti must be numeric". Nothing in that message
## mentions tibbles, and a tibble is what anyone gets from readr, dplyr or
## gapminder, so this broke the most ordinary route into the function.
##
## `ilm_profile()` calls `ilm_reduce()`, so it went the same way.
##
## The test is that the two inputs give the SAME ANSWER, not merely that
## neither errors: a coercion that silently reordered or dropped a column
## would also stop the error.

## A REAL tibble, via tibble::as_tibble(). Pinning the class by hand --
## structure(d, class = c("tbl_df", "tbl", "data.frame")) -- does not
## reproduce this: `[.tbl_df` is only dispatched to when tibble's namespace is
## loaded, so without it `[` falls through to `[.data.frame`, which DOES drop
## to a vector, and the bug disappears. A test built that way passes against
## the broken code, which is the only thing worse than no test.
mixed_tbl <- function(n = 120, seed = 4) {
  skip_if_not_installed("tibble")
  set.seed(seed)
  tibble::as_tibble(data.frame(
    num1 = rnorm(n), num2 = runif(n, 0, 10),
    int1 = sample(1:50, n, TRUE),                    # integer, not double
    fac1 = factor(sample(c("a", "b", "c"), n, TRUE)),
    fac2 = factor(sample(c("x", "y"), n, TRUE)),
    stringsAsFactors = FALSE))
}

test_that("ilm_reduce() gives the same answer for a tibble as a data frame", {
  skip_if_not_installed("PCAmixdata")
  tb <- mixed_tbl(); df <- as.data.frame(tb)
  expect_s3_class(tb, "tbl_df")                      # the input really is one
  a <- ilm_reduce(tb)
  b <- ilm_reduce(df)
  expect_identical(a$method, "famd")                 # mixed data took the mixed route
  expect_identical(a$method, b$method)
  expect_equal(a$ind_coord, b$ind_coord)
  expect_equal(a$eig, b$eig)
})

test_that("each of the three column mixes survives a tibble", {
  skip_if_not_installed("PCAmixdata")
  tb <- mixed_tbl()
  expect_identical(ilm_reduce(tb, cols = c("num1", "num2", "int1"))$method, "pca")
  expect_identical(ilm_reduce(tb, cols = c("fac1", "fac2"))$method, "mca")
  expect_identical(ilm_reduce(tb, cols = c("num1", "fac1"))$method, "famd")
})

test_that("ilm_profile() gives the same answer for a tibble as a data frame", {
  skip_if_not_installed("PCAmixdata")
  tb <- mixed_tbl(); df <- as.data.frame(tb)
  a <- suppressWarnings(ilm_profile(tb, k = 3, seed = 1))
  b <- suppressWarnings(ilm_profile(df, k = 3, seed = 1))
  expect_identical(a$cluster$k, b$cluster$k)
  expect_equal(a$reduce$ind_coord, b$reduce$ind_coord)
  expect_equal(a$cluster$ind_cluster, b$cluster$ind_cluster)
  expect_equal(a$characterization, b$characterization)
})

test_that("an integer column counts as numeric, not as a category", {
  ## is.numeric(1L) is TRUE, and if it were treated as categorical a 50-level
  ## factor would appear and MCA would run instead of PCA
  skip_if_not_installed("PCAmixdata")
  r <- ilm_reduce(mixed_tbl(), cols = c("num1", "int1"))
  expect_identical(r$method, "pca")
})

test_that("ilm_cluster() reduces a mixed tibble rather than refusing it", {
  ## It clusters COORDINATES -- a centroid is not defined on a factor -- but
  ## arriving with raw mixed data is the ordinary way to reach this function,
  ## so it goes through ilm_reduce() and says so instead of stopping.
  skip_if_not_installed("PCAmixdata")
  tb <- mixed_tbl()
  expect_message(cl <- ilm_cluster(tb, k = 2L, seed = 1), "reduced to")
  expect_s3_class(cl, "ilm_cluster")
  ## and it is the same answer as doing the two steps by hand from a tibble,
  ## which is what the routing replaces
  by_hand <- ilm_cluster(ilm_reduce(tb), k = 2L, seed = 1)
  expect_identical(as.integer(cl$ind_cluster$cluster),
                   as.integer(by_hand$ind_cluster$cluster))
})

test_that("an all-numeric tibble still clusters directly", {
  tb <- mixed_tbl()[c("num1", "num2")]
  expect_silent(cl <- ilm_cluster(tb, k = 2L, seed = 1))
  expect_identical(cl$k, 2L)
})
