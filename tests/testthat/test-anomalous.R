## Composing an anomaly scan with the exploration tools.
##
## ilm_anomaly() ends where the interesting question starts: are the odd rows
## alike? Answering it meant rebuilding the subset by hand from `row`, which
## invites lining the wrong rows up. Reported from use.

planted <- function(seed = 1) {
  set.seed(seed)
  d <- mtcars
  d$gear <- factor(d$gear)
  d[c(3, 9), "hp"] <- c(600, 700)
  d
}

test_that("the flagged rows come back, in order, with their position", {
  a <- ilm_anomaly(planted(), alpha = 0.05)
  odd <- ilm_anomalous(a)
  expect_s3_class(odd, "data.frame")
  expect_identical(nrow(odd), sum(a$flag))
  expect_identical(odd$.row, sort(a$row[a$flag]))
  ## the values really are the planted ones, i.e. the rows line up
  expect_true(all(odd$hp > 500))
  expect_true(all(odd$.driver == "hp"))
})

test_that("score columns are optional and the data columns are untouched", {
  a <- ilm_anomaly(planted(), alpha = 0.05)
  bare <- ilm_anomalous(a, score = FALSE)
  expect_false(any(c(".score", ".p_adj", ".driver") %in% names(bare)))
  expect_identical(setdiff(names(bare), ".row"), names(planted()))
})

test_that("flagged = FALSE returns the whole scan", {
  a <- ilm_anomaly(planted(), alpha = 0.05)
  expect_identical(nrow(ilm_anomalous(a, flagged = FALSE)), nrow(mtcars))
})

test_that("a scan that kept nothing says so, and takes `data`", {
  a <- ilm_anomaly(planted(), alpha = 0.05, keep_data = FALSE)
  expect_null(attr(a, "data"))
  expect_error(ilm_anomalous(a), "did not keep the data")
  expect_error(ilm_anomalous(a), "keep_data")        # names the remedy
  expect_identical(nrow(ilm_anomalous(a, data = planted())), sum(a$flag))
})

test_that("the wrong frame is caught rather than silently mismatched", {
  a <- ilm_anomaly(planted(), alpha = 0.05, keep_data = FALSE)
  expect_error(ilm_anomalous(a, data = mtcars[1:5, ]), "but `data` has")
})

test_that("an empty scan keeps its shape", {
  ## nothing flagged is an answer, not a failure -- but code that selects
  ## .driver must not break on the data sets that behaved themselves
  e <- suppressWarnings(ilm_anomaly(mtcars, alpha = 1e-9))
  odd <- ilm_anomalous(e)
  expect_identical(nrow(odd), 0L)
  expect_true(all(c(".row", ".score", ".p_adj", ".driver") %in% names(odd)))
})

test_that("the describe functions take a scan directly", {
  a <- ilm_anomaly(planted(), alpha = 0.05)
  expect_message(x <- ilm_describe_all(a), "row\\(s\\) ilm_anomaly\\(\\) flagged")
  expect_true(is.list(x))
  n_flag <- sum(a$flag)
  expect_true(all(x$numeric$n == n_flag))
  ## .driver rides along for describing: a frequency table of it is often the
  ## whole answer, and it cannot distort a summary the way it would a distance
  expect_true(".driver" %in% x$num$variable || ".driver" %in%
                unlist(lapply(x, function(z) z$variable)))
})

test_that("ilm_describe() takes one too", {
  a <- ilm_anomaly(planted(), alpha = 0.05)
  expect_message(d <- ilm_describe(a, y = "hp"), "flagged")
  expect_s3_class(d, "data.frame")
})

test_that("reduce/profile/cluster take one, and refuse too few rows", {
  a <- ilm_anomaly(planted(), alpha = 0.05)
  ## two rows is not a clustering problem; the warning has to arrive before
  ## whatever error the two rows cause downstream
  expect_warning(try(ilm_reduce(a), silent = TRUE), "only 2 row")
  expect_warning(try(ilm_reduce(a), silent = TRUE), "no structure to find")
})

test_that("a scan with enough flagged rows profiles", {
  set.seed(4)
  d <- data.frame(a = rnorm(200), b = rnorm(200), c = rnorm(200))
  d$b <- d$a + rnorm(200, 0, 0.2)                  # structure to depart FROM
  d$c <- d$a - rnorm(200, 0, 0.2)
  d[1:25, "b"] <- d[1:25, "b"] + 8                 # a real subpopulation
  s <- ilm_anomaly(d, alpha = 0.05, B = 199L)
  skip_if(sum(s$flag) < 10)
  expect_message(p <- ilm_profile(s, k = 2, B = 20, var_contrib = FALSE),
                 "row\\(s\\) ilm_anomaly\\(\\) flagged")
  expect_s3_class(p, "ilm_profile")
})

test_that("the iforest path keeps its data too", {
  skip_if_not_installed("isotree")
  a <- ilm_anomaly(planted(), method = "iforest", alpha = 0.1, seed = 1)
  expect_false(is.null(attr(a, "data")))
  expect_identical(nrow(attr(a, "data")), nrow(mtcars))
  ## and the categorical column survives, which is the point of that method
  expect_true("gear" %in% names(ilm_anomalous(a)))
})

test_that("keep_data = FALSE really drops it, on both methods", {
  expect_null(attr(ilm_anomaly(planted(), keep_data = FALSE), "data"))
  skip_if_not_installed("isotree")
  expect_null(attr(ilm_anomaly(planted(), method = "iforest", seed = 1,
                               keep_data = FALSE), "data"))
})

test_that("a non-anomaly object is refused by name", {
  expect_error(ilm_anomalous(mtcars), "must be an ilm_anomaly")
})
