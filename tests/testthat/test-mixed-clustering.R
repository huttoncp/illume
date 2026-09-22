## Mixed data through the clustering and anomaly paths
##
## The measurements behind these choices, all against known structure:
##
##   clustering, ARI    FAMD+kmeans  RF prox+PAM  VarSelLCM  Gower+PAM
##     balanced             0.433       0.257       0.440      0.447
##     within-cluster r=.75 0.265         --        0.047        --
##
## FAMD + k-means is the default because it is the only one that holds up when
## variables are correlated within a cluster, which real ones are.
##
##   anomalies, AUC     reconstruction  iForest  MCD    LOF+Gower
##     marginal              1.000       0.975   1.000   0.998
##     combination           0.999       0.950   1.000   0.997
##     cross-type            0.771       0.961   0.937   0.990
##     rare category combo   0.505       0.886   0.487   0.363
##
## The reconstruction stays the default: unbeaten at what it is for, and the
## only one with a calibrated null. iForest is there for the bottom two rows.

mixed_df <- function(n = 240, seed = 4) {
  set.seed(seed)
  g <- sample(1:3, n, TRUE)
  data.frame(n1 = rnorm(n, 1.4 * (g - 2)), n2 = rnorm(n, 1.4 * (g == 3)),
             n3 = rnorm(n), n4 = rnorm(n),
             f1 = factor(sample(c("a", "b", "c"), n, TRUE)),
             stringsAsFactors = FALSE)
}

test_that("ilm_cluster() takes raw mixed data and says what it did", {
  skip_if_not_installed("PCAmixdata")
  d <- mixed_df()
  expect_message(cl <- ilm_cluster(d, k = 3L, seed = 1), "reduced to")
  expect_message(ilm_cluster(d, k = 3L, seed = 1), "ilm_profile\\(\\)")
  expect_s3_class(cl, "ilm_cluster")
  expect_identical(cl$k, 3L)
  ## same answer as doing it by hand, which is what the routing replaces
  by_hand <- ilm_cluster(ilm_reduce(d), k = 3L, seed = 1)
  expect_identical(as.integer(cl$ind_cluster$cluster),
                   as.integer(by_hand$ind_cluster$cluster))
})

test_that("an all-numeric frame is untouched and silent", {
  d <- mixed_df()[c("n1", "n2", "n3")]
  expect_silent(cl <- ilm_cluster(d, k = 2L, seed = 1))
  expect_identical(cl$k, 2L)
})

test_that("the missing-value note states the row cost, not just the column %", {
  ## 12% missing in each of three columns is a third of the ROWS, and nobody
  ## predicts that from "12% missing" -- which is the whole point of saying it
  skip_if_not_installed("PCAmixdata")
  d <- mixed_df(300)
  set.seed(2)
  for (v in c("n1", "n2", "n3")) d[[v]][sample(300, 36)] <- NA
  expect_message(ilm_cluster(d, k = 2L, seed = 1), "of rows are complete")
  expect_message(ilm_cluster(d, k = 2L, seed = 1), "ilm_impute\\(\\)")
  ## and no note at all when nothing is missing
  expect_message(ilm_cluster(mixed_df(), k = 2L, seed = 1), "reduced to")
  m <- utils::capture.output(
    invisible(ilm_cluster(mixed_df(), k = 2L, seed = 1)), type = "message")
  expect_false(any(grepl("complete", m)))
})

test_that("a k chosen at the edge of the search is reported as such", {
  ## six separated clusters, allowed to look for three: the curve had not
  ## turned, so the answer is where the search stopped
  set.seed(9); g <- rep(1:6, each = 50)
  m <- cbind(c(0, 6, 12, 18, 24, 30)[g] + rnorm(300, 0, .4), rnorm(300, 0, .4))
  expect_warning(ilm_cluster(m, k_max = 3L, B = 15, seed = 1),
                 "largest value searched")
  ## given room, it finds them and says nothing
  expect_silent(cl <- ilm_cluster(m, k_max = 10L, B = 15, seed = 1))
  expect_identical(cl$k, 6L)
})

test_that("ilm_anomaly(method = 'iforest') uses the categorical columns", {
  skip_if_not_installed("isotree")
  ## the anomaly is a category PAIRING that never otherwise occurs; no numeric
  ## value is unusual, so a numeric-only method is at chance here
  set.seed(5); n <- 400
  d <- data.frame(n1 = rnorm(n), n2 = rnorm(n), n3 = rnorm(n), n4 = rnorm(n))
  d$f2 <- factor(sample(c("high", "low"), n, TRUE))
  d$f1 <- factor(sample(c("p", "q", "r"), n, TRUE))
  i <- sample(n, 16)
  d$f1[i] <- "p"; d$f2[i] <- "high"
  d$f1[-i][d$f2[-i] == "high"] <- sample(c("q", "r"),
                                         sum(d$f2[-i] == "high"), TRUE)
  a <- ilm_anomaly(d, method = "iforest", alpha = 0.06, seed = 1)
  expect_s3_class(a, "ilm_anomaly")
  expect_identical(attr(a, "method"), "iforest")
  expect_true(all(is.na(a$p)))                 # no null, so no p-value
  auc <- local({
    r <- rank(a$score[order(a$row)]); lab <- seq_len(n) %in% i
    (sum(r[lab]) - sum(lab) * (sum(lab) + 1) / 2) / (sum(lab) * sum(!lab))
  })
  expect_gt(auc, 0.75)                         # measured about 0.89
  ## the reconstruction is at chance on the same data, which is why this exists
  b <- suppressMessages(suppressWarnings(ilm_anomaly(d, B = 20, seed = 1)))
  auc_b <- local({
    r <- rank(b$score[order(b$row)]); lab <- seq_len(n) %in% i
    (sum(r[lab]) - sum(lab) * (sum(lab) + 1) / 2) / (sum(lab) * sum(!lab))
  })
  expect_lt(auc_b, 0.70)
})

test_that("the iforest result is sorted and names a driver, as the other does", {
  skip_if_not_installed("isotree")
  set.seed(6); n <- 200
  d <- data.frame(a = rnorm(n), b = rnorm(n), c = rnorm(n),
                  f = factor(sample(c("x", "y"), n, TRUE)))
  d$a[1:6] <- 9                                 # driver should be `a`
  a <- ilm_anomaly(d, method = "iforest", seed = 1)
  expect_false(is.unsorted(rev(a$score)))       # most anomalous first
  expect_true(all(a$driver %in% names(d)))
  expect_identical(unique(a$driver[a$row %in% 1:6]), "a")
})

test_that("iforest refuses clearly when isotree is absent", {
  skip_if(requireNamespace("isotree", quietly = TRUE),
          "isotree is installed, so the guard cannot fire")
  expect_error(ilm_anomaly(mixed_df(), method = "iforest"), "isotree")
})

test_that("ilm_var_contrib() finds a variable the clustering ignored", {
  skip_if_not_installed("PCAmixdata")
  set.seed(21); n <- 300; g <- sample(1:3, n, TRUE)
  d <- data.frame(n1 = rnorm(n, 2.2 * (g - 2)), n2 = rnorm(n, 2.2 * (g == 3)),
                  n3 = rnorm(n))
  p <- suppressMessages(suppressWarnings(ilm_profile(d, k = 3)))
  v <- ilm_var_contrib(p, d, B = 99, seed = 1)
  expect_s3_class(v, "ilm_var_contrib")
  expect_identical(nrow(v), 3L)
  expect_false(is.unsorted(rev(v$separation)))  # strongest first
  ## the informative columns outrank the noise one
  expect_true(match("n1", v$variable) < match("n3", v$variable))
})

test_that("a clustering that is one variable's levels is called out", {
  ## THE dangerous case: from inside the clustering a variable that defines it
  ## looks like the best variable, and every informative column then scores
  ## low against it -- so "drop the low ones" would be exactly backwards
  skip_if_not_installed("PCAmixdata")
  ## made deterministic: f1 is the ONLY structure in the data, so FAMD has
  ## nothing else to find and the three clusters can only be its levels
  set.seed(11); n <- 300
  d <- data.frame(n1 = rnorm(n), n2 = rnorm(n), n3 = rnorm(n),
                  f1 = factor(rep(c("a", "b", "c"), length.out = n)))
  p <- suppressMessages(suppressWarnings(ilm_profile(d, k = 3)))
  v <- ilm_var_contrib(p, d, B = 99, seed = 1)
  expect_identical(attr(v, "dominated_by"), "f1")
  expect_identical(v$verdict[1], "defines the clusters")
  out <- utils::capture.output(print(v))
  expect_true(any(grepl("re-labelling", out)))
  expect_true(any(grepl("do NOT read the bottom", out)))
})

test_that("mismatched row counts are refused rather than silently misaligned", {
  skip_if_not_installed("PCAmixdata")
  d <- mixed_df(200)
  p <- suppressMessages(suppressWarnings(ilm_profile(d, k = 2)))
  expect_error(ilm_var_contrib(p, d[1:100, ], B = 9), "labels and `data` has")
})

test_that("ilm_profile() runs the variable check itself", {
  ## The point of moving it into the path: a diagnostic nobody calls is a
  ## complaint. On the case that scored 0.006 the front door now says so by
  ## itself, without the user knowing ilm_var_contrib() exists.
  skip_if_not_installed("PCAmixdata")
  set.seed(11); n <- 300
  d <- data.frame(n1 = rnorm(n), n2 = rnorm(n), n3 = rnorm(n),
                  f1 = factor(rep(c("a", "b", "c"), length.out = n)))
  expect_warning(p <- ilm_profile(d, k = 3), "re-labelling of `f1`")
  expect_s3_class(p$var_contrib, "ilm_var_contrib")
  out <- utils::capture.output(print(p))
  expect_true(any(grepl("re-labelling of `f1`", out)))
  expect_true(any(grepl("ilm_var_contrib()", out, fixed = TRUE)))
})

test_that("the variable check can be switched off", {
  skip_if_not_installed("PCAmixdata")
  set.seed(11); n <- 200
  d <- data.frame(n1 = rnorm(n), n2 = rnorm(n), n3 = rnorm(n),
                  f1 = factor(rep(c("a", "b", "c"), length.out = n)))
  p <- suppressWarnings(suppressMessages(
    ilm_profile(d, k = 3, var_contrib = FALSE)))
  expect_null(p$var_contrib)
  ## and with it off there is no warning to suppress
  expect_silent(suppressMessages(ilm_profile(d, k = 3, var_contrib = FALSE)))
})

test_that("a clean clustering gets no complaint", {
  skip_if_not_installed("PCAmixdata")
  set.seed(21); n <- 300; g <- sample(1:3, n, TRUE)
  d <- data.frame(n1 = rnorm(n, 2.4 * (g - 2)), n2 = rnorm(n, 2.4 * (g == 3)))
  p <- suppressMessages(ilm_profile(d, k = 3))
  expect_null(attr(p$var_contrib, "dominated_by"))
  out <- utils::capture.output(print(p))
  expect_false(any(grepl("re-labelling", out)))
})
