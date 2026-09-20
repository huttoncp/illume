make_rd <- function(n = 2000L, jump = 0.8, slope = 0.5, sd = 0.5, seed = 1L) {
  set.seed(seed)
  r <- runif(n, -1, 1)
  data.frame(r = r, y = slope * r + jump * (r >= 0) + rnorm(n, 0, sd))
}

test_that("the jump at the cutoff recovers the truth", {
  f <- ilm_rdd(make_rd(), "y", "r", cutoff = 0, verbose = FALSE)
  expect_s3_class(f, "ilm_rdd")
  expect_lt(f$jump$lower, 0.8); expect_gt(f$jump$upper, 0.8)
  expect_equal(f$family, "gaussian")
  expect_equal(f$poly, 1L)
  expect_gt(f$n_below, 20L); expect_gt(f$n_above, 20L)
  ## only points within the bandwidth are used
  expect_lte(f$n_below + f$n_above, nrow(f$data))
})

test_that("no jump is reported as no jump", {
  f <- ilm_rdd(make_rd(jump = 0, seed = 2L), "y", "r", cutoff = 0,
               verbose = FALSE)
  expect_lt(f$jump$lower, 0); expect_gt(f$jump$upper, 0)
  expect_gt(f$jump$p_value, 0.05)
})

test_that("the bandwidth sweep is reported and includes the chosen one", {
  f <- ilm_rdd(make_rd(), "y", "r", cutoff = 0, verbose = FALSE)
  bw <- f$bandwidth
  expect_false(is.null(bw))
  expect_true(1 %in% bw$multiplier)
  expect_equal(bw$estimate[bw$multiplier == 1], f$jump$estimate)
  ## a wider window uses more points and gives a tighter interval
  expect_true(all(diff(bw$n) > 0))
  expect_lt(bw$se[nrow(bw)], bw$se[1])
  ## a supplied bandwidth is used as given
  expect_equal(ilm_rdd(make_rd(), "y", "r", 0, h = 0.3, verbose = FALSE)$h, 0.3)
})

test_that("the density check is about continuity, not symmetry", {
  ## A SLOPED but perfectly smooth density must not be flagged: this is the
  ## case a count-split test gets wrong, firing on three quarters of clean
  ## sloped data where the local linear density fires on about one in twenty.
  set.seed(11)
  sl <- data.frame(r = rnorm(3000, 0.6))
  sl$y <- 0.5 * sl$r + 0.8 * (sl$r >= 0) + rnorm(3000, 0, 0.5)
  expect_equal(ilm_rdd(sl, "y", "r", cutoff = 0, verbose = FALSE)$density$status,
               "OK")

  ## and a manipulated one must be
  set.seed(3); n <- 2000
  r <- runif(n, -1, 1)
  j <- which(r > -0.1 & r < 0)
  mv <- sample(j, floor(length(j) * 0.7)); r[mv] <- abs(r[mv])
  d <- data.frame(r = r, y = 0.5 * r + 0.8 * (r >= 0) + rnorm(n, 0, 0.5))
  f <- ilm_rdd(d, "y", "r", cutoff = 0, verbose = FALSE)
  expect_equal(f$density$status, "FAIL")
  expect_lt(f$density$p_value, 0.01)
  expect_output(print(f), "density at the cutoff: FAIL")
})

test_that("covariate balance catches a covariate that jumps", {
  set.seed(5); n <- 2000
  r <- runif(n, -1, 1)
  d <- data.frame(r = r, good = rnorm(n), bad = rnorm(n) + 0.9 * (r >= 0))
  d$y <- 0.5 * r + 0.8 * (r >= 0) + rnorm(n, 0, 0.5)
  f <- ilm_rdd(d, "y", "r", cutoff = 0, covariates = c("good", "bad"),
               verbose = FALSE)
  bal <- f$balance
  expect_equal(nrow(bal), 2L)
  expect_equal(bal$status[bal$covariate == "good"], "OK")
  expect_equal(bal$status[bal$covariate == "bad"], "FAIL")
})

test_that("placebo cutoffs are null where nothing happens", {
  f <- ilm_rdd(make_rd(), "y", "r", cutoff = 0, verbose = FALSE)
  expect_false(is.null(f$placebo))
  expect_gt(nrow(f$placebo), 1L)
  ## the placebos sit either side of the real cutoff, never on it
  expect_true(all(f$placebo$cutoff != 0))
  expect_lte(sum(f$placebo$status == "FAIL"), 1L)
})

test_that("fuzzy assignment is detected and not reported as sharp", {
  set.seed(4); n <- 2000
  r <- runif(n, -1, 1)
  tr <- rbinom(n, 1, ifelse(r >= 0, 0.75, 0.15))
  d <- data.frame(r = r, tr = tr, y = 0.5 * r + 0.8 * tr + rnorm(n, 0, 0.5))
  f <- ilm_rdd(d, "y", "r", cutoff = 0, treatment = "tr", verbose = FALSE)
  expect_false(f$sharp)
  expect_lt(f$compliance[["below"]], 0.3)
  expect_lt(f$compliance[["above"]], 0.9)
  expect_output(print(f), "FUZZY")

  ## a genuinely sharp design is labelled sharp
  d2 <- make_rd(); d2$tr <- as.integer(d2$r >= 0)
  expect_true(ilm_rdd(d2, "y", "r", cutoff = 0, treatment = "tr",
                      verbose = FALSE)$sharp)
})

test_that("rdd refuses what it cannot do and warns about high polynomials", {
  d <- make_rd()
  expect_error(ilm_rdd(d, "y", "nope", verbose = FALSE), "not found")
  expect_error(ilm_rdd("nope", "y", "r", verbose = FALSE), "must be a data frame")
  expect_error(ilm_rdd(d, "y", "r", cutoff = 99, verbose = FALSE),
               "too few observations")
  expect_error(ilm_rdd(d, "y", "r", poly = 0, verbose = FALSE),
               "at least 1")
  expect_warning(ilm_rdd(d, "y", "r", cutoff = 0, poly = 3, verbose = FALSE),
                 "manufacture a discontinuity")
  df <- d; df$f <- factor("a")
  expect_error(ilm_rdd(df, "y", "f", verbose = FALSE), "must be numeric")
  expect_error(ilm_rdd(d, "y", "r", cutoff = 0, kernel = "nope",
                       verbose = FALSE), "unknown `kernel`")
})

test_that("the rd plot draws", {
  f <- ilm_rdd(make_rd(), "y", "r", cutoff = 0, verbose = FALSE)
  pf <- file.path(tempdir(), "rdd.png")
  grDevices::png(pf, width = 600, height = 400); on.exit(unlink(pf))
  expect_silent(ilm_plot_rdd(f))
  grDevices::dev.off()
  expect_true(file.exists(pf))
  expect_error(ilm_plot_rdd(1), "must be an ilm_rdd")
})
