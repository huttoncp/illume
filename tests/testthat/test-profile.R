## The profiling set leans on PCAmixdata and cluster, both in Suggests.
skip_profile <- function() {
  skip_if_not_installed("PCAmixdata")
  skip_if_not_installed("cluster")
}

mixed_df <- function() {
  d <- mtcars
  d$cyl <- factor(d$cyl); d$am <- factor(d$am)
  d
}

test_that("the method is chosen from the column types", {
  skip_profile()
  expect_equal(ilm_reduce(mtcars)$method, "pca")
  expect_equal(ilm_reduce(mixed_df())$method, "famd")
  cats <- data.frame(a = factor(rep(letters[1:3], 20)),
                     b = factor(rep(LETTERS[1:2], 30)),
                     c = factor(sample(c("x", "y"), 60, TRUE)))
  expect_equal(ilm_reduce(cats)$method, "mca")
})

test_that("a reduction returns coordinates, loadings and eigenvalues", {
  skip_profile()
  r <- ilm_reduce(mtcars, ndim = 4)
  expect_s3_class(r, "ilm_reduce")
  expect_equal(r$n, nrow(mtcars))
  expect_equal(r$ndim, 4L)
  ## one coordinate column per retained dimension, plus the row id
  expect_equal(ncol(r$ind_coord), 5L)
  expect_equal(r$ind_coord$row_id, seq_len(nrow(mtcars)))
  expect_true(all(paste0("dim", 1:4) %in% names(r$ind_coord)))
  ## squared loadings live on 0..1 and every variable appears on every dim
  expect_true(all(r$var_contrib$sqload >= 0 & r$var_contrib$sqload <= 1.001))
  expect_setequal(unique(r$var_contrib$variable), names(mtcars))
  ## variance explained is a non-decreasing cumulative percentage
  expect_true(all(diff(r$eig$cum_pct_var) >= -1e-8))
  expect_lte(max(r$eig$cum_pct_var), 100.001)
  expect_output(print(r), "method = pca")
})

test_that("a reduction refuses what it cannot use", {
  skip_profile()
  expect_error(ilm_reduce("nope"), "must be a data frame")
  expect_error(ilm_reduce(mtcars, cols = "nosuch"), "not found")
  ## dates are usable by none of the three methods, so they are dropped aloud
  d <- mtcars; d$when <- as.Date("2024-01-01") + seq_len(nrow(d))
  expect_message(ilm_reduce(d), "neither numeric nor categorical")
  expect_error(ilm_reduce(data.frame(when = as.Date("2024-01-01") + 1:5)),
               "no numeric or categorical columns")
})

test_that("clustering reports size, stability and per-row ambiguity", {
  skip_profile()
  cl <- ilm_cluster(ilm_reduce(mtcars), k = 3, B = 20, seed = 1)
  expect_s3_class(cl, "ilm_cluster")
  expect_equal(cl$k, 3)
  expect_null(cl$gap)                       # k was given, so no search
  expect_equal(nrow(cl$clusters), 3L)
  expect_equal(sum(cl$clusters$size), nrow(mtcars))
  expect_equal(sum(cl$clusters$pct), 100, tolerance = 0.2)
  expect_true(all(cl$clusters$jaccard >= 0 & cl$clusters$jaccard <= 1))
  expect_setequal(as.character(cl$clusters$stability),
                  intersect(c("unstable", "moderate", "stable"),
                            as.character(cl$clusters$stability)))
  expect_equal(nrow(cl$ind_cluster), nrow(mtcars))
  expect_true(all(cl$ind_cluster$silhouette >= -1 &
                  cl$ind_cluster$silhouette <= 1))
  ## anomalous is the union of the two flags, which are different things
  expect_equal(cl$ind_cluster$is_anomalous,
               cl$ind_cluster$is_small_cluster | cl$ind_cluster$is_ambiguous)
})

test_that("well-separated clusters come back stable, and k is found", {
  skip_profile()
  ## three blobs far apart: any honest method should find them and call them
  ## stable. This is the case that caught the Jaccard universe bug, where a
  ## perfect clustering scored ~0.63 because the original cluster was compared
  ## against rows a bootstrap cannot draw.
  set.seed(7)
  co <- rbind(cbind(rnorm(40, -8), rnorm(40, -8)),
              cbind(rnorm(40, 8), rnorm(40, -8)),
              cbind(rnorm(40, 0), rnorm(40, 9)))
  cl <- ilm_cluster(co, k_max = 6, B = 30, seed = 1)
  expect_equal(cl$k, 3)
  expect_false(is.null(cl$gap))
  expect_true(all(cl$clusters$jaccard > 0.9))
  expect_true(all(as.character(cl$clusters$stability) == "stable"))
  expect_true(all(cl$clusters$mean_silhouette > 0.7))
  expect_equal(sum(cl$ind_cluster$is_ambiguous), 0L)
})

test_that("clustering refuses bad input and handles k = 1", {
  skip_profile()
  expect_error(ilm_cluster("nope"), "must be an ilm_reduce")
  ## A categorical column no longer refuses: it is reduced to coordinates on
  ## the way through, because arriving with raw mixed data is the ordinary way
  ## to reach this function rather than a mistake. A MATRIX still refuses --
  ## a matrix has one type, so there is nothing to reduce.
  expect_message(ilm_cluster(data.frame(a = letters[1:5])), "reduced to")
  expect_error(ilm_cluster(matrix(letters[1:10], 5, 2)), "must be numeric")
  ## the silhouette is undefined for a single cluster and is not invented
  one <- ilm_cluster(ilm_reduce(mtcars), k = 1, B = 5, seed = 1)
  expect_equal(one$k, 1)
  expect_true(all(is.na(one$ind_cluster$silhouette)))
  expect_equal(sum(one$ind_cluster$is_ambiguous), 0L)
})

test_that("a profile names what distinguishes each cluster", {
  skip_profile()
  p <- ilm_profile(mixed_df(), k = 3, B = 20, seed = 1)
  expect_s3_class(p, "ilm_profile")
  expect_length(p$summary, 3L)
  expect_true(all(nzchar(p$summary)))
  expect_true(all(grepl("^Cluster [0-9]", p$summary)))
  ## every characterising dimension was actually retained, and cleared the bar
  expect_true(all(p$characterization$dim <= p$reduce$ndim))
  expect_true(all(abs(p$characterization$vtest) >= 1.96))
  ## a numeric variable gets a direction; the threshold controls how many rows
  expect_match(paste(p$characterization$top_variables, collapse = " "),
               "high |low ")
  strict <- ilm_profile(mixed_df(), k = 3, B = 20, seed = 1,
                        vtest_threshold = 50)
  expect_lte(nrow(strict$characterization), nrow(p$characterization))
})

test_that("missingness profiling recovers which columns go missing together", {
  skip_profile()
  ## a and b are lost as a block; c goes on its own; e never goes at all
  set.seed(2); n <- 300
  d <- data.frame(a = rnorm(n), b = rnorm(n), c = rnorm(n), e = rnorm(n))
  blk <- sample(n, 90); d$a[blk] <- NA; d$b[blk] <- NA
  d$c[sample(n, 60)] <- NA
  expect_message(r <- ilm_reduce_na(d), "never varies")
  expect_s3_class(r, "ilm_reduce_na")
  expect_s3_class(r, "ilm_reduce")
  expect_equal(r$method, "mca")             # indicators are all categorical
  expect_false("e" %in% r$var_contrib$variable)

  ## the block shares a dimension; the independent column takes another
  d1 <- r$var_contrib[r$var_contrib$dim == 1, ]
  expect_gt(d1$sqload[d1$variable == "a"], 0.8)
  expect_gt(d1$sqload[d1$variable == "b"], 0.8)
  expect_lt(d1$sqload[d1$variable == "c"], 0.2)
  d2 <- r$var_contrib[r$var_contrib$dim == 2, ]
  expect_gt(d2$sqload[d2$variable == "c"], 0.8)
})

test_that("the missingness pipeline runs end to end and is type-checked", {
  skip_profile()
  p <- suppressMessages(ilm_profile_na(airquality, k_max = 4, B = 20, seed = 1))
  expect_s3_class(p, "ilm_profile_na")
  expect_s3_class(p, "ilm_profile")
  expect_gt(length(p$summary), 0L)
  expect_match(paste(p$summary, collapse = " "), "missing: ")
  ## no cluster names the same variable set twice in one sentence
  expect_false(any(grepl("(missing: Ozone, missing: Solar.R).*\\1", p$summary)))

  expect_error(ilm_cluster_na(ilm_reduce(mtcars)), "must be an ilm_reduce_na")
  expect_error(ilm_reduce_na(mtcars), "at least 2 columns")
  expect_error(ilm_reduce_na("nope"), "must be a data frame")
})

test_that("the profiling plots draw and check their input", {
  skip_profile()
  r <- ilm_reduce(mtcars)
  cl <- ilm_cluster(r, k_max = 5, B = 20, seed = 1)
  p <- ilm_profile(mtcars, k = 3, B = 20, seed = 1)
  pf <- file.path(tempdir(), "prof.png")
  grDevices::png(pf, width = 600, height = 400); on.exit(unlink(pf))
  expect_silent(ilm_plot_reduce(r))
  expect_silent(ilm_plot_reduce_scree(r))
  expect_silent(ilm_plot_reduce_contrib(r, dim = 1))
  expect_silent(ilm_plot_cluster(cl))
  expect_silent(ilm_plot_cluster_gap(cl))
  expect_silent(ilm_plot_profile(p))
  grDevices::dev.off()
  expect_true(file.exists(pf))

  expect_error(ilm_plot_reduce(1), "must be an ilm_reduce")
  expect_error(ilm_plot_reduce(r, dims = 1), "exactly 2")
  expect_error(ilm_plot_reduce(r, dims = c(1, 99)), "not retained")
  ## a bare `dim` here would resolve to var_contrib's own column and match
  ## everything, so an out-of-range dimension must still error
  expect_error(ilm_plot_reduce_contrib(r, dim = 99), "not retained")
  expect_error(ilm_plot_cluster_gap(ilm_cluster(r, k = 2, B = 5, seed = 1)),
               "no gap statistic")
  expect_error(ilm_plot_profile(1), "must be an ilm_profile")
  expect_error(ilm_plot_reduce_na(r), "must be an ilm_reduce_na")
  expect_error(ilm_plot_cluster_na(cl), "must be an ilm_cluster_na")
})
