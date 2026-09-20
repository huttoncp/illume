# The residual variogram: correlation as a function of separation.
#
# This exists because ilm_check_ar() matches pairs at EXACT lags, which finds
# almost nothing once times are irregular. Measured on the design below, exact
# matching found 53 pairs at lag 1, five at lag 2 and none beyond; binning by
# distance uses all 900. Fitting ilm_car1() and then having no way to check it
# would be exactly the gap this package is meant not to have.

vg_gen <- function(seed, kind = "ar", rho = 0.85, ng = 60, nt = 6, tmax = 30) {
  set.seed(seed)
  tl <- lapply(seq_len(ng), function(i) sort(sample.int(tmax, nt)))
  d <- data.frame(id = factor(rep(seq_len(ng), each = nt)),
                  t = unlist(tl), x = stats::rnorm(ng * nt))
  u <- unlist(lapply(tl, function(tt) {
    if (kind == "ar") {
      z <- numeric(length(tt)); z[1] <- stats::rnorm(1)
      for (k in seq_along(tt)[-1]) {
        p <- rho^(tt[k] - tt[k - 1])
        z[k] <- p * z[k - 1] + stats::rnorm(1, 0, sqrt(1 - p^2))
      }
      z
    } else if (kind == "cycle") 1.4 * sin(2 * pi * tt / 10)
    else rep(stats::rnorm(1), length(tt))
  }))
  b <- stats::rnorm(ng, 0, 0.3)[as.integer(d$id)]
  d$y <- 0.5 + 0.8 * d$x + b + 1.2 * u + stats::rnorm(ng * nt, 0, 0.6)
  d
}

## ---- the pair index --------------------------------------------------------

test_that("every within-group pair is found and none across groups", {
  g <- c(1, 1, 1, 2, 2)
  t <- c(0, 4, 9, 2, 5)
  pr <- illume:::ilm_pair_index(g, t)
  # 3 choose 2 plus 2 choose 2 = 4
  expect_equal(pr$n_total, 4L)
  expect_setequal(pr$d, c(4, 9, 5, 3))
  expect_true(all(g[pr$i] == g[pr$j]))
  expect_true(all(pr$d > 0))
})

test_that("the pair count is capped rather than left to explode", {
  g <- rep(1L, 200); t <- seq_len(200)
  pr <- illume:::ilm_pair_index(g, t, max_pairs = 500L)
  expect_equal(pr$n_total, 19900L)      # 200 choose 2
  expect_length(pr$i, 500L)
  expect_length(pr$d, 500L)
  expect_error(illume:::ilm_pair_index(1:5, 1:5), "no pairs")
})

test_that("distance binning finds pairs where exact lag matching does not", {
  d <- vg_gen(3)
  ex <- illume:::ilm_resid_acf(stats::rnorm(nrow(d)), as.integer(d$id), d$t,
                               maxlag = 6L)
  np <- attr(ex, "n_pairs")
  expect_lt(sum(np[-1]), 30L)           # essentially nothing beyond lag 1
  pr <- illume:::ilm_pair_index(as.integer(d$id), d$t)
  expect_equal(pr$n_total, 900L)        # 60 groups x choose(6, 2)
})

## ---- the verdicts ----------------------------------------------------------

test_that("an omitted correlation is detected and CAR(1) is named", {
  d <- vg_gen(3, "ar")
  f <- ilm_model(y ~ x + (1 | id), data = d, family = "gaussian",
                 verbose = FALSE)
  v <- suppressWarnings(ilm_variogram(f, d$t, d$id, breaks = 6L, B = 60L,
                                      plot = FALSE, verbose = FALSE))
  expect_true(any(v$table$status %in% c("WARN", "FAIL")))
  expect_gt(v$table$estimate[1], 0.1)   # strongest at the shortest separation
  expect_match(illume:::ilm_variogram_advice(v), "ilm_car1", fixed = TRUE)
})

test_that("fitting the structure clears the diagnosis that asked for it", {
  d <- vg_gen(3, "ar")
  f <- ilm_model(y ~ x + (1 | id), data = d, family = "gaussian",
                 ar = suppressWarnings(ilm_car1(d$t, d$id, verbose = FALSE)),
                 verbose = FALSE)
  v <- suppressWarnings(ilm_variogram(f, d$t, d$id, breaks = 6L, B = 60L,
                                      plot = FALSE, verbose = FALSE))
  # this only passes because fitted() includes the AR latent; without that the
  # residual still carries the whole process and every near bin fires
  expect_lt(sum(v$table$status != "OK"), 2L)
})

test_that("the three shapes get three different remedies", {
  adv <- illume:::ilm_variogram_advice
  mk <- function(est, null = rep(-0.05, length(est)))
    list(type = "correlation",
         table = data.frame(estimate = est, null_mean = null))
  # decaying from the shortest separation: autoregression
  expect_match(adv(mk(c(0.40, 0.18, 0.02, -0.06, -0.10))), "ilm_car1",
               fixed = TRUE)
  # the same level everywhere: a group effect, not a decaying one
  expect_match(adv(mk(rep(0.75, 5))), "group effect")
  # alternating: a cycle, and an AR term cannot represent one
  expect_match(adv(mk(c(0.20, -0.55, 0.12, 0.42, -0.40))), "cycle")
  expect_match(adv(mk(c(0.20, -0.55, 0.12, 0.42, -0.40))), "ilm_fourier",
               fixed = TRUE)
  # too little to go on
  expect_match(adv(mk(c(NA_real_, NA_real_, 0.1))), "add a correlation over time")
})

test_that("a flat group effect is not mistaken for autocorrelation", {
  d <- vg_gen(5, "flat")
  f <- ilm_model(y ~ x, data = d, family = "gaussian", verbose = FALSE)
  v <- suppressWarnings(ilm_variogram(f, d$t, d$id, breaks = 5L, B = 40L,
                                      plot = FALSE, verbose = FALSE))
  expect_true(all(v$table$status %in% c("WARN", "FAIL")))
  expect_match(illume:::ilm_variogram_advice(v), "group effect")
})

test_that("a cycle is sent to the seasonal terms, not to CAR(1)", {
  d <- vg_gen(3, "cycle")
  f <- ilm_model(y ~ x + (1 | id), data = d, family = "gaussian",
                 verbose = FALSE)
  v <- suppressWarnings(ilm_variogram(f, d$t, d$id, breaks = 6L, B = 40L,
                                      plot = FALSE, verbose = FALSE))
  a <- illume:::ilm_variogram_advice(v)
  expect_match(a, "cycle")
  expect_false(grepl("ilm_car1", a, fixed = TRUE))
})

## ---- presentation ----------------------------------------------------------

test_that("semivariance is the same information turned over", {
  d <- vg_gen(3, "ar")
  f <- ilm_model(y ~ x + (1 | id), data = d, family = "gaussian",
                 verbose = FALSE)
  a <- suppressWarnings(ilm_variogram(f, d$t, d$id, breaks = 4L, B = 30L,
                                      plot = FALSE, verbose = FALSE))
  b <- suppressWarnings(ilm_variogram(f, d$t, d$id, breaks = 4L, B = 30L,
                                      type = "semivariance", plot = FALSE,
                                      verbose = FALSE))
  expect_equal(b$table$estimate, round(1 - a$table$estimate, 4))
  # the interval reverses with the sign, so lo and hi swap
  expect_equal(b$table$lo, round(1 - a$table$hi, 4))
  expect_equal(b$table$hi, round(1 - a$table$lo, 4))
  # the test is the same test either way
  expect_equal(b$table$p, a$table$p)
  expect_equal(b$table$status, a$table$status)
})

test_that("the plot renders and restores the device", {
  d <- vg_gen(3, "ar")
  f <- ilm_model(y ~ x + (1 | id), data = d, family = "gaussian",
                 verbose = FALSE)
  v <- suppressWarnings(ilm_variogram(f, d$t, d$id, breaks = 5L, B = 30L,
                                      plot = FALSE, verbose = FALSE))
  ff <- tempfile(fileext = ".png")
  grDevices::png(ff, width = 700, height = 450)
  on.exit({grDevices::dev.off(); unlink(ff)}, add = TRUE)
  before <- graphics::par(c("mar", "mfrow"))
  ilm_plot_variogram(v)
  expect_identical(graphics::par(c("mar", "mfrow")), before)
  ilm_plot_variogram(v, colour = "navy", fill = "#EEF", alpha = 0.6, size = 1.3)
  grDevices::dev.off()
  on.exit(unlink(ff), add = FALSE)
  expect_gt(file.size(ff), 4000)
})

test_that("the printed report names the remedy", {
  d <- vg_gen(3, "ar")
  f <- ilm_model(y ~ x + (1 | id), data = d, family = "gaussian",
                 verbose = FALSE)
  out <- capture.output(suppressWarnings(
    ilm_variogram(f, d$t, d$id, breaks = 5L, B = 40L, plot = FALSE,
                  verbose = TRUE)))
  expect_true(any(grepl("outside the envelope", out)))
  expect_true(any(grepl("ilm_car1", out, fixed = TRUE)))
  expect_true(any(grepl("bins are not independent", out)))
})

## ---- informative failure ---------------------------------------------------

test_that("misspecification is named", {
  d <- vg_gen(3, "ar", ng = 20, nt = 4, tmax = 15)
  f <- ilm_model(y ~ x + (1 | id), data = d, family = "gaussian",
                 verbose = FALSE)
  expect_error(ilm_variogram(d, d$t, d$id), "must be a fitted ilm_model")
  expect_error(ilm_variogram(f), "both required")
  expect_error(ilm_variogram(f, d$t[-1], d$id), "but the model has")
  expect_error(ilm_variogram(f, factor(d$t), d$id), "is a factor")
  expect_error(ilm_variogram(f, paste0("w", d$t), d$id), "must be numeric")
  expect_error(ilm_variogram(f, d$t, d$id, type = "nope"), "Options are")
  expect_error(ilm_variogram(f, d$t, d$id, B = 1), "at least 2")
  expect_error(ilm_variogram(f, d$t, d$id, breaks = 1), "at least 2 bins")
  expect_error(ilm_plot_variogram(data.frame(x = 1)), "returned by ilm_variogram")
  expect_warning(ilm_variogram(f, d$t, d$id, breaks = 3L, B = 20L,
                               plot = FALSE, verbose = FALSE),
                 "FAIL verdict is unreachable")
})

test_that("too few successful refits gives no verdict rather than a wrong one", {
  d <- vg_gen(3, "ar", ng = 20, nt = 4, tmax = 15)
  f <- ilm_model(y ~ x + (1 | id), data = d, family = "gaussian",
                 verbose = FALSE)
  v <- suppressWarnings(ilm_variogram(f, d$t, d$id, breaks = 3L, B = 5L,
                                      plot = FALSE, verbose = FALSE))
  expect_true(all(v$table$status == "INCONCLUSIVE"))
})
