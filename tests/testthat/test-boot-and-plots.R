# Bootstrap intervals, missingness, and the plot layer.
#
# The plot tests assert on RENDERED OUTPUT, not on the absence of an error.
# Two real bugs in this layer produced completely blank images while raising
# nothing at all, so "it ran" is not evidence that anything was drawn.

png_bytes <- function(expr) {
  f <- tempfile(fileext = ".png")
  grDevices::png(f, width = 600, height = 450)
  on.exit(grDevices::dev.off(), add = TRUE)
  suppressMessages(suppressWarnings(force(expr)))
  grDevices::dev.off()
  on.exit()
  file.size(f)
}
BLANK <- 2000     # an empty device lands near 500 bytes; a real plot far above

## ---- bootstrap -------------------------------------------------------------

test_that("bootstrap intervals bracket the observed statistic", {
  set.seed(1)
  x <- stats::rlnorm(300)
  for (ct in c("percentile", "bca", "normal", "basic")) {
    r <- ilm_boot_ci(x, stat = "mean", R = 300, ci_type = ct, seed = 1)
    expect_lt(r$lower, r$observed, label = ct)
    expect_gt(r$upper, r$observed, label = ct)
    expect_equal(r$ci_type, ct)
  }
})

test_that("BCa agrees with boot::boot.ci", {
  skip_if_not_installed("boot")
  set.seed(9)
  x <- stats::rlnorm(200, 0, 0.7)
  mine <- ilm_boot_ci(x, stat = "mean", R = 3000, ci_type = "bca", seed = 9)
  set.seed(9)
  bb <- boot::boot(x, function(d, i) mean(d[i]), R = 3000)
  ref <- boot::boot.ci(bb, type = "bca", conf = 0.95)$bca[4:5]
  # different RNG streams, so agreement is to within Monte Carlo error
  expect_equal(mine$lower, ref[1], tolerance = 0.05)
  expect_equal(mine$upper, ref[2], tolerance = 0.05)
})

test_that("bootstrap grouping and argument checks behave", {
  d <- ilm_sim(n_id = 20)
  r <- ilm_boot_ci(d, "score", by = "grp", R = 100, seed = 1)
  expect_true("grp" %in% names(r))
  expect_gt(nrow(r), 1L)
  expect_error(ilm_boot_ci(d, "score", ci_type = "student"), "Options are")
  expect_error(ilm_boot_ci(d, "score", conf = 2), "between 0 and 1")
  expect_error(ilm_boot_ci(d, "grp"), "must be numeric")
  expect_error(ilm_boot_ci(d, "nope"), "not found")
})

test_that("boot_diff reports the difference and whether it excludes zero", {
  d <- ilm_sim()
  d2 <- d[d$grp %in% c("alpha", "beta"), ]
  d2$grp <- factor(d2$grp)
  r <- ilm_boot_diff(d2, "score", "grp", R = 300, seed = 1)
  expect_equal(nrow(r), 1L)
  expect_lt(r$lower, r$observed)
  expect_gt(r$upper, r$observed)
  expect_type(r$excludes_zero, "logical")
  expect_error(ilm_boot_diff(d, "score", "grp"), "exactly 2 levels")
})

test_that("describe_na counts missing values and sorts by them", {
  d <- ilm_sim(n_id = 20)
  r <- ilm_describe_na(d, "lab_value")
  expect_equal(r$na, sum(is.na(d$lab_value)))
  expect_equal(r$obs, nrow(d))
  a <- ilm_describe_na_all(d)
  expect_equal(a$variable[1], "lab_value")   # most missing first
  expect_equal(nrow(a), ncol(d))
})

## ---- plots -----------------------------------------------------------------

test_that("the geom table and the argument checks agree", {
  sp <- ilm_geom_spec()
  expect_setequal(sp$geom, c("histogram", "density", "bar", "box", "violin",
                             "point", "bin2d", "line", "spine"))
  need <- sp$geom[sp$needs_y]
  d <- ilm_sim(n_id = 20)
  for (g in need)
    expect_error(ilm_plot(d, "score", geom = g), "needs both", info = g)
})

test_that("geom auto-selection follows the data types", {
  d <- ilm_sim(n_id = 20)
  expect_equal(ilm_pick_geom(d$score)$geom, "histogram")
  expect_equal(ilm_pick_geom(d$grp)$geom, "bar")
  expect_equal(ilm_pick_geom(d$date)$geom, "line")
  expect_equal(ilm_pick_geom(d$score, d$income)$geom, "point")
  expect_equal(ilm_pick_geom(d$grp, d$score)$geom, "box")
  expect_equal(ilm_pick_geom(d$grp, d$site)$geom, "spine")
})

test_that("a scatter becomes a binned density above n_max", {
  set.seed(1)
  x <- stats::rnorm(20000); y <- stats::rnorm(20000)
  p <- ilm_pick_geom(x, y, n_max = 5000)
  expect_equal(p$geom, "bin2d")
  expect_match(p$reason, "exceeds n_max")
  expect_equal(ilm_pick_geom(x, y, n_max = 50000)$geom, "point")
})

test_that("every geom actually draws something", {
  d <- ilm_sim(n_id = 25)
  big <- data.frame(a = stats::rnorm(20000), b = stats::rnorm(20000))
  expect_gt(png_bytes(ilm_plot(d, "score")), BLANK)
  expect_gt(png_bytes(ilm_plot(d, "grp")), BLANK)
  expect_gt(png_bytes(ilm_plot(d, "grp", "score")), BLANK)
  expect_gt(png_bytes(ilm_plot(d, "grp", "score", geom = "violin")), BLANK)
  expect_gt(png_bytes(ilm_plot(d, "score", "income")), BLANK)
  expect_gt(png_bytes(ilm_plot(d, "grp", "site")), BLANK)
  expect_gt(png_bytes(ilm_plot(d, "date", "score")), BLANK)
  expect_gt(png_bytes(ilm_plot(d, "score", geom = "density")), BLANK)
  expect_gt(png_bytes(ilm_plot(big, "a", "b")), BLANK)
})

test_that("a legend does not blank the plot", {
  # passing vectors through do.call() rendered nothing, silently; the formula
  # interface does not
  d <- ilm_sim(n_id = 25)
  expect_gt(png_bytes(ilm_plot(d, "score", "income", by = "grp")), BLANK)
  expect_gt(png_bytes(ilm_plot(d, "score", geom = "density", by = "grp")), BLANK)
})

test_that("a legend inside a multi-panel layout does not blank the figure", {
  # tinyplot reserves legend space by altering the device layout
  d <- ilm_sim(n_id = 25)
  sz <- png_bytes({
    op <- graphics::par(mfrow = c(1, 2))
    on.exit(graphics::par(op), add = TRUE)
    ilm_plot(d, "score")
    ilm_plot(d, "score", "income", by = "grp")
  })
  expect_gt(sz, BLANK)
})

test_that("aesthetics and themes draw", {
  d <- ilm_sim(n_id = 25)
  expect_gt(png_bytes(ilm_plot(d, "score", fill = "steelblue", alpha = 0.5)), BLANK)
  expect_gt(png_bytes(ilm_plot(d, "score", "income", size = 1.4, colour = "darkred")), BLANK)
  expect_gt(png_bytes(ilm_plot(d, "grp", "score", theme = "clean")), BLANK)
  expect_gt(png_bytes(ilm_plot(d, "score", "income", by = "grp", palette = "Dark 2")), BLANK)
})

test_that("plot arguments are validated with the options named", {
  d <- ilm_sim(n_id = 20)
  expect_error(ilm_plot(d, "score", geom = "scatter"), "Options are")
  expect_error(ilm_plot(d, "score", geom = "scatter"), "would use")
  expect_error(ilm_plot(d, "nope"), "not found")
  expect_error(ilm_plot(d, "grp", geom = "histogram"), "needs `x` to be numeric")
  expect_error(ilm_plot(d, "score", theme = "ggplot"), "unknown `theme`")
  expect_error(ilm_plot(d, "score", alpha = 3), "between 0")
  expect_error(ilm_plot(d, "score", size = -1), "positive")
  expect_error(ilm_plot(d, "score", colour = "red", color = "blue"), "use one")
  # png_bytes() suppresses warnings, so this one is checked on its own device
  f <- tempfile(fileext = ".png"); grDevices::png(f)
  expect_warning(suppressMessages(ilm_plot(d, "score", "income", geom = "histogram")),
                 "is ignored")
  grDevices::dev.off()
})

test_that("plot_all and plot_missing draw", {
  d <- ilm_sim(n_id = 25)
  expect_gt(png_bytes(ilm_plot_all(d, class = "numeric", max_panels = 4)), BLANK)
  expect_gt(png_bytes(ilm_plot_missing(d)), BLANK)
  expect_error(ilm_plot_all(d, class = "numerical"), "Options are")
})

## ---- model plots -----------------------------------------------------------

test_that("model plots draw and validate their arguments", {
  set.seed(1)
  d <- ilm_sim(n_id = 30)
  f <- ilm_model(score ~ income + grp + (1 | id), data = d,
                 family = "gaussian", verbose = FALSE)
  expect_gt(png_bytes(ilm_plot_model(f, "coef")), BLANK)
  expect_gt(png_bytes(ilm_plot_model(f, "random")), BLANK)
  expect_gt(png_bytes(ilm_plot_model(f, "effect", term = "income")), BLANK)
  expect_gt(png_bytes(plot(f)), BLANK)
  expect_error(ilm_plot_model(f, "coeff"), "Options are")
  expect_error(ilm_plot_model(f, "effect", term = "nope"), "Available")
  expect_error(ilm_plot_model(d, "coef"), "must be a fitted ilm_model")
})

test_that("the coefficient plot returns the interval table", {
  set.seed(1)
  d <- ilm_sim(n_id = 30)
  f <- ilm_model(score ~ income + (1 | id), data = d, family = "gaussian",
                 verbose = FALSE)
  ## drawn inside a device: a plot call with none open writes Rplots.pdf into
  ## the test directory
  ff <- tempfile(fileext = ".png"); grDevices::png(ff)
  tb <- suppressMessages(ilm_plot_model(f, "coef"))
  grDevices::dev.off()
  expect_true(all(c("term", "estimate", "lower", "upper", "excludes_zero") %in%
                    names(tb)))
  expect_true(all(tb$lower < tb$estimate & tb$estimate < tb$upper))
})

## ---- effect plots under interactions ---------------------------------------

test_that("an effect plot with an interaction draws one curve per partner", {
  # Before this, the plot held the interacting partner at its most common
  # level, drew one of several quite different slopes, and said nothing.
  set.seed(11)
  n <- 600
  d <- data.frame(id = factor(rep(1:60, each = 10)), x = stats::rnorm(n),
                  z = stats::rnorm(n),
                  g = factor(sample(c("a", "b", "c"), n, TRUE)))
  b <- stats::rnorm(60, 0, 0.5)[as.integer(d$id)]
  d$y <- 0.5 + 0.8 * d$x + 0.3 * d$z +
    c(a = 0, b = 0.9, c = -0.6)[as.character(d$g)] +
    ifelse(d$g == "b", 1, 0) * d$x + 0.7 * d$x * d$z + b +
    stats::rnorm(n, 0, 0.8)

  ff <- tempfile(fileext = ".png")
  grDevices::png(ff, width = 700, height = 450)
  on.exit({grDevices::dev.off(); unlink(ff)}, add = TRUE)

  # factor partner: one curve per level
  f1 <- ilm_model(y ~ g * x + z + (1 | id), d, family = "gaussian",
                  verbose = FALSE)
  r1 <- ilm_plot_model(f1, "effect", term = "x")
  expect_setequal(unique(r1$by), c("a", "b", "c"))
  # and the curves genuinely differ: b was generated with an extra slope of 1
  sl <- tapply(seq_len(nrow(r1)), r1$by, function(i)
    stats::coef(stats::lm(r1$fit[i] ~ as.numeric(r1$value[i])))[2])
  expect_gt(sl[["b"]] - sl[["a"]], 0.6)

  # numeric partner: quartiles, labelled with the value
  f2 <- ilm_model(y ~ x * z + (1 | id), d, family = "gaussian", verbose = FALSE)
  r2 <- ilm_plot_model(f2, "effect", term = "x")
  expect_length(unique(r2$by), 3L)
  expect_true(all(grepl("^z = ", unique(r2$by))))

  # no interaction: one curve, as before
  f3 <- ilm_model(y ~ x + z + (1 | id), d, family = "gaussian", verbose = FALSE)
  r3 <- ilm_plot_model(f3, "effect", term = "x")
  expect_identical(unique(r3$by), "")
})

test_that("the interaction partner is chosen sensibly", {
  mf <- data.frame(x = stats::rnorm(50), z = stats::rnorm(50),
                   g = factor(sample(letters[1:9], 50, TRUE)))
  pick <- illume:::ilm_effect_partner
  expect_null(pick("x", c("x", "z"), mf))
  expect_equal(pick("x", c("x", "z", "x:z"), mf)$var, "z")
  expect_length(pick("x", c("x", "z", "x:z"), mf)$values, 3L)
  expect_equal(pick("z", c("x", "z", "x:z"), mf)$var, "x")
  # a factor with more levels than fit in a legend is capped, not dropped
  p <- pick("x", c("x", "g", "x:g"), mf, max_levels = 4L)
  expect_length(p$values, 4L)
  # a partner that is not in the model frame is not offered
  expect_null(pick("x", c("x", "x:nope"), mf))
})
