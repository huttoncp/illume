# Residual autocorrelation: the statistic, the envelope, and the picture.
#
# The envelope is the whole point. Within-group residuals from a mixed model
# have a NEGATIVELY biased autocorrelation at every lag, so the symmetric
# 2/sqrt(n) band that acf() draws flags correlation that is not there. One of
# the tests below asserts that bias exists, because if it ever stopped
# existing the expensive refit machinery would no longer be earning its keep.

panel_ar <- function(rho, seed, ng = 40, nt = 8) {
  set.seed(seed)
  n <- ng * nt
  d <- data.frame(id = factor(rep(seq_len(ng), each = nt)),
                  t = rep(seq_len(nt), ng), x = stats::rnorm(n))
  b <- stats::rnorm(ng, 0, 0.6)[as.integer(d$id)]
  e <- unlist(lapply(seq_len(ng), function(i) {
    z <- numeric(nt); z[1] <- stats::rnorm(1)
    for (k in seq_len(nt)[-1])
      z[k] <- rho * z[k - 1] + stats::rnorm(1, 0, sqrt(1 - rho^2))
    z
  }))
  d$y <- 0.5 + 0.8 * d$x + b + e
  d$cnt <- stats::rpois(n, exp(pmin(0.3 + 0.5 * d$x + b + e, 3)))
  d$cat <- factor(c("a", "b", "c")[
    pmax(1, pmin(3, round(2 + (d$y - mean(d$y)) / stats::sd(d$y))))])
  d
}

gfit <- function(d, resp = "y", fam = "gaussian")
  ilm_model(stats::as.formula(paste(resp, "~ x + (1 | id)")), data = d,
            family = fam, verbose = FALSE)

## ---- Durbin-Levinson -------------------------------------------------------

test_that("the partial autocorrelation matches stats::pacf", {
  set.seed(4)
  x <- as.numeric(stats::arima.sim(list(ar = c(0.6, -0.3)), 500))
  a <- as.numeric(stats::acf(x, lag.max = 10, plot = FALSE)$acf)[-1]
  expect_equal(illume:::ilm_pacf_dl(a),
               as.numeric(stats::pacf(x, lag.max = 10, plot = FALSE)$acf),
               tolerance = 1e-8)
})

test_that("an invalid autocorrelation sequence stops the recursion", {
  # |r1| >= 1 cannot come from a real sequence, and the recursion divides by
  # 1 - r1^2 at the very next step
  expect_true(is.na(illume:::ilm_pacf_dl(c(1, 0.5, 0.2))[2]))
  expect_equal(illume:::ilm_pacf_dl(c(1, 0.5, 0.2))[1], 1)
  # a gap propagates rather than being skipped over
  expect_true(all(is.na(illume:::ilm_pacf_dl(c(0.5, NA, 0.2))[2:3])))
  expect_true(all(is.na(illume:::ilm_pacf_dl(c(NA, 0.5)))))
})

## ---- detection -------------------------------------------------------------

test_that("independent errors are not flagged", {
  f <- gfit(panel_ar(0, 21))
  r <- ilm_check_ar(f, panel_ar(0, 21)$t, panel_ar(0, 21)$id,
                    maxlag = 4L, B = 40L, verbose = FALSE)
  expect_false(any(r$table$status == "FAIL"))
  expect_false(any(r$pacf$status == "FAIL"))
})

test_that("AR(1) errors are detected, and the partial panel points at lag 1", {
  d <- panel_ar(0.7, 22)
  r <- ilm_check_ar(gfit(d), d$t, d$id, maxlag = 4L, B = 40L, verbose = FALSE)
  expect_true(r$table$status[1] %in% c("WARN", "FAIL"))
  expect_gt(r$table$estimate[1], 0.3)
  # the partial autocorrelation at lag 1 is the largest, which is the AR(1)
  # signature and what tells the user the order
  expect_equal(which.max(abs(r$pacf$estimate)), 1L)
})

test_that("the envelope is not centred on zero, which is why it exists", {
  d <- panel_ar(0, 23)
  r <- ilm_check_ar(gfit(d), d$t, d$id, maxlag = 4L, B = 40L, verbose = FALSE)
  # residuals within a group sum to roughly zero, so pairs of them correlate
  # negatively even under a correct model; a band symmetric about zero would
  # call every one of these lags significant
  expect_lt(max(r$table$null_mean), -0.01)
  expect_true(all(r$table$hi < 0.25))
})

test_that("the band drawn is the band the verdict is made on", {
  d <- panel_ar(0.5, 24)
  r <- ilm_check_ar(gfit(d), d$t, d$id, maxlag = 4L, B = 60L, verbose = FALSE)
  tab <- r$table
  outside <- tab$estimate < tab$lo | tab$estimate > tab$hi
  # away from the 0.05 boundary the two must agree exactly, or the picture is
  # telling the user something different from the table
  clear <- is.finite(tab$p) & (tab$p < 0.02 | tab$p > 0.15)
  expect_equal(outside[clear], tab$p[clear] <= 0.05)
})

## ---- families --------------------------------------------------------------

test_that("the check works for every family, not just the multinomial", {
  d <- panel_ar(0.6, 25)
  for (spec in list(c("y", "gaussian"), c("cnt", "poisson"),
                    c("cat", "multinomial"))) {
    r <- ilm_check_ar(gfit(d, spec[1], spec[2]), d$t, d$id, maxlag = 3L,
                      B = 25L, verbose = FALSE)

    expect_equal(nrow(r$table), 3L, label = spec[2])
    expect_true(all(is.finite(r$table$estimate)), info = spec[2])
    # one residual series for a univariate family, one per category otherwise
    expect_equal(ncol(r$acf_by_component),
                 if (spec[2] == "multinomial") 3L else 1L, label = spec[2])
  }
})

## ---- the plot --------------------------------------------------------------

draws <- function(expr) {
  ff <- tempfile(fileext = ".png")
  grDevices::png(ff, width = 900, height = 450)
  force(expr)
  grDevices::dev.off()
  file.size(ff)
}

test_that("both panels render, and the envelope can be reused", {
  d <- panel_ar(0.6, 26)
  r <- ilm_check_ar(gfit(d), d$t, d$id, maxlag = 4L, B = 30L, verbose = FALSE)
  # a blank device is about 500 bytes
  expect_gt(draws(ilm_plot_acf(r)), 4000)
  expect_gt(draws(ilm_plot_acf(r, which = "acf")), 4000)
  expect_gt(draws(ilm_plot_acf(r, which = "pacf")), 4000)
  # passing the finished check back in must not refit anything
  s1 <- ilm_plot_acf(r, which = "acf")
  expect_identical(s1$acf, r$table)
})

test_that("appearance arguments are honoured and par is restored", {
  d <- panel_ar(0.6, 27)
  r <- ilm_check_ar(gfit(d), d$t, d$id, maxlag = 3L, B = 25L, verbose = FALSE)
  ff <- tempfile(fileext = ".png")
  grDevices::png(ff, width = 900, height = 450)
  on.exit({grDevices::dev.off(); unlink(ff)}, add = TRUE)
  before <- graphics::par(c("mfrow", "mar"))
  ilm_plot_acf(r, colour = "navy", fill = "#FDD", alpha = 0.6, size = 1.4)
  # two panels means touching mfrow and mar; a diagnostic that leaves the
  # device reconfigured breaks whatever the user draws next
  expect_identical(graphics::par(c("mfrow", "mar")), before)
})

test_that("the model plot interface reaches it", {
  d <- panel_ar(0.6, 28)
  f <- gfit(d)
  # B is small here to keep the test quick, which the warning is right to
  # object to
  expect_gt(draws(suppressWarnings(
    ilm_plot_model(f, "acf", time = d$t, group = d$id, maxlag = 3L, B = 25L))),
    4000)
  expect_gt(draws(suppressWarnings(
    plot(f, "acf", time = d$t, group = d$id, maxlag = 3L, B = 25L))), 4000)
})

## ---- informative failure ---------------------------------------------------

test_that("misspecification is named rather than guessed at", {
  d <- panel_ar(0.3, 29)
  f <- gfit(d)
  expect_error(ilm_check_ar(d, d$t, d$id), "must be a fitted ilm_model")
  expect_error(ilm_check_ar(f), "both required")
  expect_error(ilm_check_ar(f, d$t[-1], d$id), "but the model has")
  expect_error(ilm_check_ar(f, d$t, d$id[-1]), "but the model has")
  expect_error(ilm_check_ar(f, d$t, d$id, maxlag = 50L),
               "longest series spans only")
  expect_error(ilm_check_ar(f, d$t, d$id, maxlag = 0L),
               "positive whole number")
  expect_error(ilm_check_ar(f, paste0("w", d$t), d$id, maxlag = 3L),
               "must be numeric")
  expect_error(ilm_check_ar(f, factor(d$t), d$id, maxlag = 3L),
               "is a factor")
  expect_error(ilm_check_ar(f, d$t, d$id, maxlag = 3L, B = 1L),
               "at least 2")
  expect_error(ilm_plot_acf(f, d$t, d$id, which = "lines"), "Options are")
  expect_error(ilm_plot_acf("nope"), "or the result of ilm_check_ar")
})

test_that("conditions the pairing cannot handle are warned about", {
  d <- panel_ar(0.3, 30)
  f <- gfit(d)
  expect_warning(ilm_check_ar(f, d$t + 0.5, d$id, maxlag = 2L, B = 5L,
                              verbose = FALSE),
                 "not whole-numbered")
  dd <- d; dd$t[2] <- dd$t[1]
  expect_warning(ilm_check_ar(f, dd$t, dd$id, maxlag = 2L, B = 5L,
                              verbose = FALSE),
                 "combinations are repeated")
  expect_warning(ilm_plot_acf(f, d$t, d$id, maxlag = 2L, B = 20L),
                 "FAIL verdict is unreachable")
})

test_that("too few successful refits gives no verdict rather than a wrong one", {
  d <- panel_ar(0.6, 31)
  r <- ilm_check_ar(gfit(d), d$t, d$id, maxlag = 3L, B = 5L, verbose = FALSE)
  expect_true(all(r$table$status == "INCONCLUSIVE"))
  expect_gt(draws(ilm_plot_acf(r)), 4000)
})

test_that("the printed report names the remedy", {
  d <- panel_ar(0.7, 32)
  out <- capture.output(
    ilm_check_ar(gfit(d), d$t, d$id, maxlag = 3L, B = 40L, verbose = TRUE))
  expect_true(any(grepl("outside the envelope", out)))
  expect_true(any(grepl("AR(1) term over time within group", out, fixed = TRUE)))
  expect_true(any(grepl("lags are not independent", out)))
})
