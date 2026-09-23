# Model plots: the part of the plot layer that needs a fitted model. The rest
# of the plot layer is illumex's.
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
