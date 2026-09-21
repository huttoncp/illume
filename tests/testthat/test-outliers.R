test_that("the iqr rule is the one a boxplot draws", {
  ## base R's boxplot.stats() applies the same 1.5 IQR fence, so the two must
  ## flag the same values. An independent implementation is the only thing that
  ## catches a fence applied at the wrong width.
  for (v in list(mtcars$hp, mtcars$wt, airquality$Ozone, rexp(300))) {
    mine <- ilm_outliers(v)
    flagged <- sort(unique(mine$value[!is.na(mine$is_outlier) & mine$is_outlier]))
    theirs <- sort(unique(grDevices::boxplot.stats(v)$out))
    expect_equal(flagged, theirs)
  }
})

test_that("a lower threshold flags at least as much, never less", {
  ## the fence is applied exactly once. An earlier form subtracted the
  ## threshold into the fence AND compared the score against it, doubling the
  ## width, so lowering the threshold flagged FEWER values.
  y <- c(rnorm(200), 8, -9, 12)
  n <- vapply(c(3, 2, 1.5, 1, 0.5), function(th)
    sum(ilm_outliers(y, threshold = th)$is_outlier, na.rm = TRUE), 1L)
  expect_true(all(diff(n) >= 0))
  expect_gt(n[length(n)], n[1])
})

test_that("each rule scores in its own units and defaults sensibly", {
  set.seed(1); y <- c(rnorm(200), 10)
  for (m in c("iqr", "mad", "zscore")) {
    r <- ilm_outliers(y, method = m)
    expect_equal(nrow(r), length(y))
    expect_equal(unique(r$method), m)
    expect_true(all(r$score >= 0, na.rm = TRUE))
    expect_true(r$is_outlier[length(y)])          # the planted value
  }
  expect_equal(unique(ilm_outliers(y, "iqr")$threshold), 1.5)
  expect_equal(unique(ilm_outliers(y, "mad")$threshold), 3.5)
  expect_equal(unique(ilm_outliers(y, "zscore")$threshold), 3)

  ## a single huge value inflates the sd and so hides itself under zscore,
  ## where the median and MAD are not dragged around by it
  z <- c(rnorm(30), 100)
  expect_gt(max(ilm_outliers(z, "mad")$score),
            max(ilm_outliers(z, "zscore")$score))
})

test_that("degenerate and missing input are handled rather than crashed on", {
  ## no spread means nothing is unusual, not a division by zero
  flat <- ilm_outliers(rep(3, 20))
  expect_true(all(flat$score == 0))
  expect_false(any(flat$is_outlier))
  ## one row per input value, with NA carried through
  y <- c(1, 2, NA, 3, 50)
  r <- ilm_outliers(y)
  expect_equal(nrow(r), 5L)
  expect_true(is.na(r$is_outlier[3]))
  expect_true(is.na(r$score[3]))
  expect_error(ilm_outliers(letters), "must be numeric")
  expect_error(ilm_outliers(c(1, NA)), "at least 2 non-missing")
})

test_that("outliers_all covers every numeric column and keeps the row", {
  r <- ilm_outliers_all(mtcars)
  expect_true(all(r$is_outlier))                 # flagged_only by default
  expect_true(all(r$variable %in% names(mtcars)))
  expect_true(all(r$row_id >= 1 & r$row_id <= nrow(mtcars)))
  ## row_id and value point at the same cell in the original data
  for (i in seq_len(min(nrow(r), 10L)))
    expect_equal(r$value[i], mtcars[[r$variable[i]]][r$row_id[i]])
  ## and each column's flags agree with the single-vector function
  for (v in unique(r$variable)) {
    got <- sort(r$value[r$variable == v])
    want <- sort(mtcars[[v]][which(ilm_outliers(mtcars[[v]])$is_outlier)])
    expect_equal(got, want)
  }
  ## flagged_only = FALSE returns everything scored
  allr <- ilm_outliers_all(mtcars, flagged_only = FALSE)
  expect_equal(nrow(allr), nrow(mtcars) * sum(vapply(mtcars, is.numeric, TRUE)))
})

test_that("grouping changes what counts as unusual", {
  ## a value ordinary for its own group can be extreme against the pooled
  ## distribution, so flagging without `by` on grouped data partly rediscovers
  ## the groups
  g <- ilm_outliers_all(mtcars, by = "cyl")
  expect_true("cyl" %in% names(g))
  expect_equal(names(g)[1], "cyl")
  expect_false(identical(nrow(g), nrow(ilm_outliers_all(mtcars))))
  ## every flag is judged within its own group
  for (i in seq_len(min(nrow(g), 6L))) {
    sub <- mtcars[[g$variable[i]]][mtcars$cyl == g$cyl[i]]
    expect_true(g$value[i] %in% sub)
  }
})

test_that("outliers_all refuses what it cannot do", {
  expect_error(ilm_outliers_all("nope"), "must be a data frame")
  expect_error(ilm_outliers_all(mtcars, cols = "nosuch"), "not found")
  expect_error(ilm_outliers_all(data.frame(a = letters[1:5])),
               "no numeric columns")
  ## nothing flagged is an empty table with the right shape, not an error
  e <- ilm_outliers_all(data.frame(a = rep(1, 10), b = rep(2, 10)))
  expect_equal(nrow(e), 0L)
  expect_true(all(c("row_id", "variable", "value", "score", "is_outlier")
                  %in% names(e)))
})

test_that("the EDA plots draw and check their columns", {
  pf <- file.path(tempdir(), "eda.png")
  grDevices::png(pf, width = 700, height = 500); on.exit({
    grDevices::dev.off(); unlink(pf) }, add = TRUE)
  expect_silent(ilm_plot_histogram(mtcars, "mpg"))
  expect_silent(ilm_plot_density(mtcars, "mpg", by = "cyl"))
  expect_silent(ilm_plot_box(mtcars, "mpg", x = "cyl"))
  expect_silent(ilm_plot_violin(mtcars, "mpg", x = "cyl"))
  expect_silent(ilm_plot_bar(mtcars, "cyl"))
  expect_silent(ilm_plot_line(data.frame(t = 1:20, v = 1:20), "v", "t"))
  expect_silent(ilm_plot_stat_error(mtcars, "mpg", "cyl"))
  expect_silent(ilm_plot_stat_error(mtcars, "mpg", "cyl", stat = "median"))
  expect_silent(ilm_plot_na_all(airquality))
  expect_silent(ilm_plot_na(airquality, "Ozone", by = "Month"))
  expect_silent(ilm_plot_var(mtcars, "mpg"))
  expect_silent(ilm_plot_var_all(mtcars, cols = c("mpg", "cyl")))
  expect_silent(ilm_plot_c(ilm_plot_bar(mtcars, "cyl"),
                           ilm_plot_bar(mtcars, "am"), nrow = 1))

  expect_error(ilm_plot_histogram("nope", "mpg"), "must be a data frame")
  expect_error(ilm_plot_histogram(mtcars, "nosuch"), "column not found")
  expect_error(ilm_plot_histogram(mtcars, c("a", "b")), "single column name")
  expect_error(ilm_plot_na(airquality, "Ozone"), "`by` is required")
  expect_error(ilm_plot_var_pairs(mtcars, cols = "mpg"), "at least 2 columns")
  expect_error(ilm_plot_c(nrow = 1), "at least one plot expression")
})

test_that("plot_var picks the geometry from the column types", {
  expect_message(ilm_plot_var(mtcars, "mpg", verbose = TRUE), "continuous")
  expect_message(ilm_plot_var(mtcars, "mpg", "wt", verbose = TRUE), "scatter")
  ## mtcars$cyl is stored as a NUMBER, so it is continuous until made a factor
  expect_message(ilm_plot_var(mtcars, "mpg", "cyl", verbose = TRUE), "scatter")
  d <- mtcars; d$cyl <- factor(d$cyl); d$am <- factor(d$am)
  expect_message(ilm_plot_var(d, "mpg", "cyl", verbose = TRUE), "boxplot")
  expect_message(ilm_plot_var(d, "cyl", verbose = TRUE), "categorical")
  expect_message(ilm_plot_var(d, "cyl", "am", verbose = TRUE), "grouped bar")
})

test_that("boot_diff keeps its draws so the distribution can be drawn", {
  d <- ilm_sim()
  b <- ilm_boot_diff(d, "score", "grp", R = 200, seed = 1)
  dr <- attr(b, "draws")
  expect_false(is.null(dr))
  expect_equal(ncol(dr), nrow(b))
  expect_lte(nrow(dr), 200L)
  ## p_superiority is the share of replicates above zero, and agrees with them
  expect_equal(b$p_superiority, unname(round(colMeans(dr > 0), 4)))
  ## the draws are labelled by the comparison they belong to
  expect_equal(colnames(dr), paste(b$from, b$to, sep = " -> "))
  ## and it lines up with the sign of the estimate
  expect_true(all((b$p_superiority > 0.5) == (b$observed > 0)))

  pf <- file.path(tempdir(), "bd.png")
  grDevices::png(pf, width = 600, height = 400); on.exit({
    grDevices::dev.off(); unlink(pf) }, add = TRUE)
  expect_silent(ilm_plot_boot_diff(b, row = 1))
  expect_silent(ilm_plot_boot_diff(b, row = nrow(b), type = "histogram"))
  expect_error(ilm_plot_boot_diff(b, row = nrow(b) + 1L), "must be one of")
  expect_error(ilm_plot_boot_diff(mtcars), "must be an ilm_boot_diff")
})
