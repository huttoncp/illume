## ilm_plot_anomaly(): four views of a scan, and the verdict the "scores" view
## reads off its reference band.

## A k-dimensional shared structure with rows pushed OFF it (orthogonal to the
## structure, so no single column need be extreme), or with t-distributed
## noise, whose heavy tail is a continuum rather than a handful of outliers.
## Realistic column names: tinyplot's legend deparsing breaks on long labels
## and hides behind short ones.
pa_data <- function(n = 400L, p = 8L, k = 2L, seed = 1L, contam = 0,
                    sd_e = 0.4, df = Inf) {
  set.seed(seed)
  F <- matrix(rnorm(n * k), n, k)
  L <- matrix(rnorm(k * p), k, p)
  E <- if (is.finite(df))
    matrix(rt(n * p, df) * sd_e / sqrt(df / (df - 2)), n, p)
  else matrix(rnorm(n * p, 0, sd_e), n, p)
  X <- F %*% L + E
  bad <- if (contam > 0) sample(n, max(1L, round(contam * n))) else integer(0)
  if (length(bad)) {
    P <- t(L) %*% solve(L %*% t(L)) %*% L
    for (i in bad) {
      v <- rnorm(p); v <- as.numeric(v - P %*% v)
      X[i, ] <- X[i, ] + 3.2 * v / sqrt(sum(v^2)) * sqrt(p) * sd_e
    }
  }
  d <- as.data.frame(X)
  names(d) <- c("north", "central", "south", "east", "west", "coastal",
                "interior", "highland", "lowland", "island", "valley",
                "plateau")[seq_len(p)]
  d
}

## the help-page example: three rows ordinary in every column that break the
## pattern the other rows make
pa_example <- function() {
  set.seed(1)
  n <- 300
  f <- rnorm(n)
  d <- data.frame(north = f + rnorm(n, 0, .3), central = 2 * f + rnorm(n, 0, .3),
                  south = -f + rnorm(n, 0, .3))
  d[1:3, ] <- rbind(c(1.2, -2.4, 1.2), c(-1, 2, 1), c(0.8, 1.6, 0.8))
  d
}

pa_sorted <- function(a) {
  d <- as.data.frame(a); class(d) <- "data.frame"
  d <- d[order(-d$score), , drop = FALSE]
  row.names(d) <- NULL
  d
}

test_that("every view draws, and says nothing when there is nothing to say", {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off(), add = TRUE)
  d <- pa_example()
  a <- ilm_anomaly(d, progress = FALSE)
  expect_equal(sum(a$flag), 3L)
  expect_silent(ilm_plot_anomaly(a))
  expect_null(ilm_plot_anomaly(a))
  expect_silent(ilm_plot_anomaly(a, "drivers"))
  expect_silent(ilm_plot_anomaly(a, "row"))
  expect_silent(ilm_plot_anomaly(a, "row", row = 2))
  expect_silent(ilm_plot_anomaly(a, "scores", top_n = 10))
  expect_silent(ilm_plot_anomaly(a, pch = "open circle"))
  expect_error(ilm_plot_anomaly(a, pch = "blob"), "plotting character")
  skip_if_not_installed("PCAmixdata")
  expect_silent(ilm_plot_anomaly(a, "map"))
})

test_that("the band restarts at the line: j ranks past it read the null at j * n / m", {
  ## a null curve that falls linearly from 100 at rank 1 to 1 at rank 100
  n <- 100L
  mid <- seq(100, 1, length.out = n)
  nq <- cbind(lo = mid - 1, mid = mid, hi = mid + 1)
  ## left of the line, the null's own top ranks
  expect_equal(ilm_anom_band(nq, 10L, c(1, 5, 10), shifted = FALSE)[, "mid"],
               mid[c(1, 5, 10)])
  ## right of it, the unflagged rows on their own: the first of the 90 reads
  ## rank 100 / 90 of the null, the 45th reads rank 50
  b <- ilm_anom_band(nq, 10L, c(11, 55), shifted = TRUE)
  expect_equal(b[, "mid"], stats::approx(seq_len(n), mid,
                                         xout = c(1, 45) * n / 90)$y)
  ## a rank before the first reads as the first
  expect_equal(unname(ilm_anom_band(nq, 10L, 10.5, shifted = TRUE)[, "mid"]),
               mid[1])

  ## Ten outliers followed by a remainder that IS the null for 90 rows sit
  ## apart. The same remainder read against the band without the restart
  ## would look like a continuum -- which is why it restarts.
  rest <- stats::approx(seq_len(n), mid, xout = seq_len(90) * n / 90)$y
  s <- c(rep(500, 10), rest)
  sh <- ilm_anom_shape(s, nq, 10L)
  expect_equal(sh$shape, "separate")
  expect_equal(sh$J, 20L)
  expect_gt(mean(s[10 + seq_len(20)] > nq[10 + seq_len(20), "hi"]), 0.75)
  ## a remainder well above the band all the way along is a continuum
  expect_equal(ilm_anom_shape(c(rep(500, 10), rest * 1.5), nq, 10L)$shape,
               "continuum")
  ## and one above it for about half the window is not clean either way
  half <- rest; half[seq_len(10)] <- half[seq_len(10)] * 1.5
  expect_equal(ilm_anom_shape(c(rep(500, 10), half), nq, 10L)$shape, "unclear")
  ## nothing flagged, or no reference: nothing is read off the band
  expect_true(is.na(ilm_anom_shape(mid, nq, 0L)$shape))
  expect_true(is.na(ilm_anom_shape(s, NULL, 10L)$shape))
})

test_that("planted anomalies stand clear and heavy tails do not", {
  ## Rates across scans rather than one scan each, as the calibration in the
  ## help page was measured: planted anomalies were called a continuum in 0
  ## of 200 scans there, and t(5) noise at 1000 rows in 58 of 60.
  planted <- vapply(1:3, function(s) {
    a <- suppressMessages(ilm_anomaly(pa_data(seed = 20L + s, contam = 0.05),
                                      seed = s, progress = FALSE))
    ilm_anom_shape(a$score, attr(a, "null_curve"), sum(a$flag))$shape
  }, "")
  expect_false(any(planted == "continuum"))
  expect_gte(sum(planted == "separate"), 2L)
  heavy <- vapply(1:3, function(s) {
    a <- suppressMessages(ilm_anomaly(pa_data(n = 1000L, p = 12L, k = 3L,
                                              seed = 40L + s, df = 5),
                                      seed = s, progress = FALSE))
    ilm_anom_shape(a$score, attr(a, "null_curve"), sum(a$flag))$shape
  }, "")
  expect_gte(sum(heavy == "continuum"), 2L)
})

test_that("no break at the line names the remedy; a clean scan reads nothing", {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off(), add = TRUE)
  a <- suppressMessages(ilm_anomaly(pa_data(n = 1000L, p = 12L, k = 3L,
                                            seed = 41L, df = 5),
                                    seed = 1L, progress = FALSE))
  expect_equal(ilm_anom_shape(a$score, attr(a, "null_curve"),
                              sum(a$flag))$shape, "continuum")
  expect_message(ilm_plot_anomaly(a), "fixed share")
  expect_message(ilm_plot_anomaly(a), "ilm_profile")
  ## the words on the plot come from the same place
  tx <- ilm_anom_scores_text(FALSE, sum(a$flag),
                             ilm_anom_shape(a$score, attr(a, "null_curve"),
                                            sum(a$flag)), 0.05)
  expect_match(tx$sub, "No break", fixed = TRUE)
  expect_lte(nchar(tx$sub), 70L)

  ## nothing flagged: no line, and no verdict read off the band
  set.seed(3)
  n <- 300; f <- rnorm(n)
  clean <- data.frame(north = f + rnorm(n, 0, .3),
                      central = 2 * f + rnorm(n, 0, .3),
                      south = -f + rnorm(n, 0, .3))
  b <- ilm_anomaly(clean, progress = FALSE)
  expect_equal(sum(b$flag), 0L)
  expect_match(ilm_anom_scores_text(FALSE, 0L, NULL, 0.05)$sub,
               "Nothing flagged", fixed = TRUE)
  expect_silent(ilm_plot_anomaly(b))
  expect_error(ilm_plot_anomaly(b, "drivers"), "no rows were flagged")
})

test_that("one column carrying the flags is called a problem in that column", {
  ## Straight from the counts. A column that drives 8 of 10 flags but also
  ## drives half of the unflagged rows is not singled out -- in clean data a
  ## column the structure explains poorly drives more rows for that reason
  ## alone. The same 8 of 10 from a column that drives 5% of the others is.
  mk <- function(base) data.frame(
    flag = rep(c(TRUE, FALSE), c(10, 390)),
    driver = c(rep("central", 8), rep("south", 2),
               rep(c("central", "north"), round(390 * c(base, 1 - base)))),
    stringsAsFactors = FALSE)
  cols <- c("north", "central", "south")
  expect_false(ilm_anom_drivers_table(mk(0.52), cols)$one_column)
  v <- ilm_anom_drivers_table(mk(0.05), cols)
  expect_true(v$one_column)
  expect_equal(v$table$column[1], "central")
  expect_equal(v$table$flagged[1], 8)
  expect_match(v$message, "ilm_recode_errors()", fixed = TRUE)
  ## too few flags to read anything into
  few <- data.frame(flag = rep(c(TRUE, FALSE), c(3, 97)),
                    driver = c(rep("central", 3), rep("north", 97)))
  expect_false(ilm_anom_drivers_table(few, cols)$one_column)

  ## and end to end: a misplaced decimal in 15 rows of one column
  grDevices::pdf(NULL); on.exit(grDevices::dev.off(), add = TRUE)
  d <- pa_data(seed = 8L)
  set.seed(2); i <- sample(nrow(d), 15)
  d$central[i] <- d$central[i] * 10
  a <- suppressMessages(ilm_anomaly(d, seed = 1L, progress = FALSE))
  expect_message(ilm_plot_anomaly(a, "drivers"), "one column, `central`")
})

test_that("the row view sets z against its column and the residual against its own size", {
  d <- pa_example()
  a <- ilm_anomaly(d, progress = FALSE)
  v <- ilm_anom_row_values(a, pa_sorted(a), 1L, d)
  ## z is the value against its own column, exactly as the scan standardised it
  z <- scale(as.matrix(d))[1, ]
  expect_equal(v$values$z[match(names(z), v$values$column)], unname(z))
  ## the residual, unscaled, is what the score is made of
  R <- attr(a, "residual")
  expect_equal(sum(R[1, ]^2), a$score[a$row == 1])
  expect_equal(v$values$residual[match(colnames(R), v$values$column)],
               unname(R[1, ] / apply(R, 2, stats::mad)))
  ## ordinary in every column, far off the structure: the case this exists for
  expect_lt(max(abs(v$values$z)), 3)
  expect_gte(max(abs(v$values$residual)), 3)
  expect_match(v$sub, "Odd only as a combination", fixed = TRUE)
  ## a value extreme in its own column is said to be one
  d2 <- d; d2$south[10] <- 8
  a2 <- ilm_anomaly(d2, progress = FALSE)
  v2 <- ilm_anom_row_values(a2, pa_sorted(a2), 10L, d2)
  expect_match(v2$sub, "`south` is extreme on its own", fixed = TRUE)
  ## the default row is the highest-scoring one
  expect_equal(ilm_anom_row_values(a, pa_sorted(a), NULL, d)$row,
               pa_sorted(a)$row[1])
})

test_that("rows that cannot be shown are refused with the reason", {
  d <- pa_example()
  d$north[5] <- NA
  a <- suppressMessages(ilm_anomaly(d, progress = FALSE))
  s <- pa_sorted(a)
  expect_error(ilm_anom_row_values(a, s, 5L, d), "was not scored")
  expect_error(ilm_anom_row_values(a, s, 9999L, d), "there is no row 9999")
  expect_error(ilm_anom_row_values(a, s, 1.5, d), "single row number")
  expect_error(ilm_anom_row_values(a, s, c(1, 2), d), "single row number")
  ## a frame that does not match the scan is caught rather than misaligned
  expect_error(ilm_anom_row_values(a, s, 1L, d[-(1:10), ]), "does not line up")
})

test_that("views that need the data say how to supply it", {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off(), add = TRUE)
  d <- pa_example()
  a <- ilm_anomaly(d, keep_data = FALSE, progress = FALSE)
  expect_error(ilm_plot_anomaly(a, "row"), "keep_data = TRUE")
  expect_silent(ilm_plot_anomaly(a, "row", data = d))
  expect_error(ilm_plot_anomaly(a, "row", data = d[1:2]), "no column south")
  expect_error(ilm_plot_anomaly(a, "row", data = d[1:50, ]), "refers to row")
  skip_if_not_installed("PCAmixdata")
  expect_error(ilm_plot_anomaly(a, "map"), "keep_data = TRUE")
  expect_silent(ilm_plot_anomaly(a, "map", data = d))
})

test_that("bad input is refused with the reason", {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off(), add = TRUE)
  d <- pa_example()
  a <- ilm_anomaly(d, progress = FALSE)
  expect_error(ilm_plot_anomaly(d), "must be an ilm_anomaly\\(\\) result; it is data.frame")
  expect_error(ilm_plot_anomaly(a, "nope"), "should be one of")
  expect_error(ilm_plot_anomaly(a, top_n = 0), "`top_n`")
  expect_error(ilm_plot_anomaly(a, top_n = c(5, 10)), "`top_n`")
  ## `[` keeps the attributes, so a slice of a scan still carries the whole
  ## scan's reference, and its ranks would be read against the wrong one
  expect_error(ilm_plot_anomaly(a[1:20, ]), "part of a scan")
})

test_that("inside a multi-panel layout the legend is dropped, not the figure", {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off(), add = TRUE)
  a <- ilm_anomaly(pa_example(), progress = FALSE)
  op <- graphics::par(mfrow = c(2, 2)); on.exit(graphics::par(op), add = TRUE)
  expect_silent(ilm_plot_anomaly(a))
  expect_silent(ilm_plot_anomaly(a, "drivers"))
  expect_message(ilm_plot_anomaly(a, "row"), "legend dropped")
  ## and the layout survives the plots that change the margins
  expect_equal(graphics::par("mfrow"), c(2L, 2L))
})

test_that("an isolation forest gets no band and no row view, and says why", {
  skip_if_not_installed("isotree")
  grDevices::pdf(NULL); on.exit(grDevices::dev.off(), add = TRUE)
  d <- pa_data(seed = 3L, contam = 0.03)
  d$zone <- factor(sample(c("urban", "rural", "remote"), nrow(d), TRUE))
  a <- suppressMessages(ilm_anomaly(d, method = "iforest", seed = 1L))
  expect_match(ilm_anom_scores_text(TRUE, sum(a$flag), NULL, 0.05)$sub,
               "not a test", fixed = TRUE)
  expect_silent(ilm_plot_anomaly(a))
  expect_silent(ilm_plot_anomaly(a, "drivers"))
  expect_error(ilm_plot_anomaly(a, "row"), "method = \"reconstruction\"",
               fixed = TRUE)
  skip_if_not_installed("PCAmixdata")
  expect_silent(ilm_plot_anomaly(a, "map"))
})

test_that("a tibble scan draws the views that go back to the data", {
  skip_if_not_installed("tibble")
  grDevices::pdf(NULL); on.exit(grDevices::dev.off(), add = TRUE)
  d <- tibble::as_tibble(pa_example())
  a <- ilm_anomaly(d, progress = FALSE)
  expect_silent(ilm_plot_anomaly(a, "row"))
  skip_if_not_installed("PCAmixdata")
  expect_silent(ilm_plot_anomaly(a, "map"))
})
