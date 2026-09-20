# Cycles: telling one apart from autocorrelation, and fitting it.
#
# The reason this exists: before it, ilm_check_ar() answered a pure 12-month
# cycle with "add an AR(1) term over time within group". An AR term cannot
# represent a cycle at all, so that is not a weak suggestion, it is a wrong
# one, and the user would have followed it.

## ---- the discriminator, on hand-made sequences -----------------------------

test_that("a returning peak is read as a cycle and a decay is not", {
  hint <- illume:::ilm_period_hint
  # sinusoidal: down through zero, back up to a peak at lag 12
  seas <- cos(2 * pi * seq_len(14) / 12)
  expect_equal(hint(seas, which(abs(seas) > 0.2), 14L), 12L)
  # autoregressive with the group mean absorbed: decays through zero and stays
  ar <- c(0.54, 0.14, -0.15, -0.34, -0.44)
  expect_null(hint(ar, 1:5, 5L))
  # plain geometric decay, never negative
  expect_null(hint(0.7^(1:8), 1:8, 8L))
  # a single flagged lag is not a pattern
  expect_null(hint(seas, 12L, 14L))
  # alternation is not seasonality
  expect_null(hint(c(-0.4, 0.4, -0.4, 0.4), 1:4, 4L))
})

test_that("a cycle longer than maxlag is reported as unfinished", {
  trunc <- illume:::ilm_period_truncated
  # still climbing back at the last lag: there is more to see
  expect_true(trunc(cos(2 * pi * seq_len(8) / 12), 8L))
  # already decayed and flat
  expect_false(trunc(c(0.54, 0.14, -0.15, -0.34, -0.44), 5L))
  expect_false(trunc(c(NA_real_, NA_real_), 2L))
})

## ---- end to end ------------------------------------------------------------

panel_seas <- function(seed, amp = 1.1, ng = 30, nt = 24) {
  set.seed(seed)
  n <- ng * nt
  d <- data.frame(id = factor(rep(seq_len(ng), each = nt)),
                  t = rep(seq_len(nt), ng), x = stats::rnorm(n))
  b <- stats::rnorm(ng, 0, 0.6)[as.integer(d$id)]
  d$y <- 0.5 + 0.8 * d$x + b + amp * sin(2 * pi * d$t / 12) +
    0.5 * cos(2 * pi * d$t / 12) + stats::rnorm(n, 0, 0.8)
  d
}

test_that("a seasonal residual is named as a cycle, not as autocorrelation", {
  d <- panel_seas(41)
  f <- ilm_model(y ~ x + (1 | id), data = d, family = "gaussian",
                 verbose = FALSE)
  r <- ilm_check_ar(f, d$t, d$id, maxlag = 14L, B = 25L, verbose = FALSE)
  ad <- illume:::ilm_ar_advice(r$envelope, "acf")
  expect_match(ad$why, "cycle of period 12")
  expect_match(ad$fix, "ilm_fourier", fixed = TRUE)
  expect_false(grepl("add an AR(1) term", ad$fix, fixed = TRUE))
  expect_match(ad$fix_short, "not AR")
})

test_that("an autoregressive residual still gets the autoregressive remedy", {
  set.seed(42)
  ng <- 40; nt <- 8
  d <- data.frame(id = factor(rep(seq_len(ng), each = nt)),
                  t = rep(seq_len(nt), ng), x = stats::rnorm(ng * nt))
  b <- stats::rnorm(ng, 0, 0.6)[as.integer(d$id)]
  e <- unlist(lapply(seq_len(ng), function(i) {
    z <- numeric(nt); z[1] <- stats::rnorm(1)
    for (k in seq_len(nt)[-1])
      z[k] <- 0.7 * z[k - 1] + stats::rnorm(1, 0, sqrt(1 - 0.49))
    z
  }))
  d$y <- 0.5 + 0.8 * d$x + b + e
  f <- ilm_model(y ~ x + (1 | id), data = d, family = "gaussian",
                 verbose = FALSE)
  ad <- illume:::ilm_ar_advice(
    ilm_check_ar(f, d$t, d$id, maxlag = 4L, B = 25L, verbose = FALSE)$envelope,
    "acf")
  expect_match(ad$fix, "AR(1) term over time within group", fixed = TRUE)
  expect_false(grepl("cycle of period", ad$why))
})

## ---- the remedy ------------------------------------------------------------

test_that("Fourier columns are what they claim to be", {
  m <- ilm_fourier(1:12, period = 12, K = 2)
  expect_equal(colnames(m), c("sin1", "cos1", "sin2", "cos2"))
  expect_equal(attr(m, "period"), 12)
  # orthogonal over exactly one cycle, which is why they behave in a model
  expect_equal(crossprod(m), diag(c(6, 6, 6, 6)), tolerance = 1e-10,
               ignore_attr = TRUE)
  # and they repeat with the period
  expect_equal(ilm_fourier(1:6, 12), ilm_fourier(13:18, 12),
               tolerance = 1e-10, ignore_attr = TRUE)
  # at 2K == period the sine is zero at every whole time, so it is dropped
  # rather than carried as an empty column
  expect_equal(colnames(ilm_fourier(1:4, period = 4, K = 2)),
               c("sin1", "cos1", "cos2"))
  expect_true(all(is.na(ilm_fourier(c(1, NA, 3), 12))[2, ]))
})

test_that("a model with Fourier terms recovers the cycle and predicts it", {
  d <- panel_seas(43)
  f <- ilm_model(y ~ x + ilm_fourier(t, 12) + (1 | id), data = d,
                 family = "gaussian", verbose = FALSE)
  b <- coef(f)
  expect_equal(unname(b[["x"]]), 0.8, tolerance = 0.12)
  # measured over 20 replicates these average 1.098 and 0.494 with standard
  # deviations of 0.034 and 0.046, so a single fit needs room to be a couple
  # of deviations out without failing
  expect_equal(unname(b[["ilm_fourier(t, 12)sin1"]]), 1.1, tolerance = 0.15)
  expect_equal(unname(b[["ilm_fourier(t, 12)cos1"]]), 0.5, tolerance = 0.35)
  # one term, two degrees of freedom, not two separate terms
  expect_true("ilm_fourier(t, 12)" %in% f$term_labels)
  expect_equal(ilm_anova(f, test = "Wald")$Df[2], 2L)
  # unlike poly(), it carries no fitted state, so times beyond the data are
  # predicted by the cycle rather than extrapolated wrongly
  nd <- data.frame(x = 0, t = c(3, 15, 27),
                   id = factor(1, levels = levels(d$id)))
  p <- as.numeric(predict(f, newdata = nd))
  expect_equal(p[1], p[2], tolerance = 1e-8)
  expect_equal(p[1], p[3], tolerance = 1e-8)
})

test_that("adding the cycle clears the diagnosis that asked for it", {
  d <- panel_seas(44)
  f <- ilm_model(y ~ x + ilm_fourier(t, 12) + (1 | id), data = d,
                 family = "gaussian", verbose = FALSE)
  r <- ilm_check_ar(f, d$t, d$id, maxlag = 14L, B = 30L, verbose = FALSE)
  # a couple of the fourteen may still trip at alpha = 0.05 by chance; a
  # cycle that was not removed trips almost all of them
  expect_lt(sum(r$table$status != "OK"), 4L)
  expect_null(illume:::ilm_period_hint(
    r$table$estimate, which(r$table$status != "OK"), 14L))
})

test_that("misspecified Fourier arguments are named", {
  expect_error(ilm_fourier(1:10), "`period` is required")
  expect_error(ilm_fourier(factor(1:10), 12), "is a factor")
  expect_error(ilm_fourier(letters[1:5], 12), "must be numeric")
  expect_error(ilm_fourier(1:10, period = 1), "greater than 1")
  expect_error(ilm_fourier(1:10, 12, K = 0), "at least 1")
  expect_error(ilm_fourier(1:10, 12, K = 9), "supports at most 6 harmonics")
  expect_error(ilm_fourier(1:10, 3, K = 2), "supports at most 1 harmonic:")
})

## ---- the cyclic spline -----------------------------------------------------

test_that("the cyclic basis closes on itself", {
  b <- ilm_cyclic(0:11, period = 12, df = 5)
  expect_equal(dim(b), c(12L, 5L))
  # the same point one cycle on
  expect_equal(ilm_cyclic(0:5, 12, 5), ilm_cyclic(12:17, 12, 5),
               tolerance = 1e-10, ignore_attr = TRUE)
  expect_equal(ilm_cyclic(-3:-1, 12, 5), ilm_cyclic(9:11, 12, 5),
               tolerance = 1e-10, ignore_attr = TRUE)
  # one column is dropped because a cyclic basis sums to one everywhere; what
  # is left must still be independent of the intercept
  expect_equal(qr(cbind(1, ilm_cyclic(1:12, 12, 5)))$rank, 6L)
  expect_true(all(is.na(ilm_cyclic(c(1, NA, 3), 12, 4))[2, ]))
})

test_that("the curve is smooth across the wrap, not merely continuous", {
  # a discontinuity in the second derivative is what an ordinary bs() on a
  # phase variable would leave behind, and it would show as a spike here
  set.seed(2)
  g <- seq(-3, 3, by = 0.01)
  f <- as.numeric(ilm_cyclic(g, 12, 5) %*% stats::rnorm(5))
  d2 <- diff(f, differences = 2)
  wrap <- which.min(abs(g[-(1:2)]))
  expect_lt(abs(d2[wrap]), max(abs(d2[-(wrap + (-3:3))])))
})

test_that("a cyclic spline fits a seasonal model", {
  d <- panel_seas(45)
  f <- ilm_model(y ~ x + ilm_cyclic(t, 12, 5) + (1 | id), data = d,
                 family = "gaussian", verbose = FALSE)
  expect_equal(unname(coef(f)[["x"]]), 0.8, tolerance = 0.12)
  expect_equal(ilm_anova(f, test = "Wald")$Df[2], 5L)
  # periodic prediction, same as the Fourier basis gives
  nd <- data.frame(x = 0, t = c(4, 16, 28),
                   id = factor(1, levels = levels(d$id)))
  p <- as.numeric(predict(f, newdata = nd))
  expect_equal(p[1], p[2], tolerance = 1e-8)
  expect_equal(p[1], p[3], tolerance = 1e-8)
})

test_that("misspecified spline arguments are named", {
  expect_error(ilm_cyclic(1:10), "`period` is required")
  expect_error(ilm_cyclic(factor(1:10), 12), "is a factor")
  expect_error(ilm_cyclic(letters[1:5], 12), "must be numeric")
  expect_error(ilm_cyclic(1:10, period = 0.5), "greater than 1")
  expect_error(ilm_cyclic(1:10, 12, df = 2), "at least 3")
  expect_error(ilm_cyclic(1:10, 12, df = 12), "at most 11 columns")
})

## ---- the formula environment ------------------------------------------------

test_that("a term can use a variable local to the calling function", {
  # lme4::nobars(), mgcv::interpret.gam() and reformulate() each hand back a
  # formula carrying an environment of their own. Without carrying the
  # caller's through, any term with a local argument -- ilm_fourier(t, 12, K),
  # ns(x, df = d), a locally defined helper -- fails to find it.
  d <- panel_seas(46)
  fit_with <- function(K) {
    ilm_model(y ~ x + ilm_fourier(t, 12, K) + (1 | id), data = d,
              family = "gaussian", verbose = FALSE)
  }
  f <- fit_with(2L)
  expect_equal(ncol(f$X), 6L)          # intercept + x + 4 Fourier columns
  expect_equal(unname(coef(f)[["x"]]), 0.8, tolerance = 0.12)

  local_period <- 12
  g <- ilm_model(y ~ x + ilm_cyclic(t, local_period, 4) + (1 | id), data = d,
                 family = "gaussian", verbose = FALSE)
  expect_equal(ilm_anova(g, test = "Wald")$Df[2], 4L)

  # and it still predicts, which needs the environment on the stored terms too
  nd <- data.frame(x = 0, t = 5, id = factor(1, levels = levels(d$id)))
  expect_true(is.finite(as.numeric(predict(f, newdata = nd))))
  expect_true(is.finite(as.numeric(predict(g, newdata = nd))))
})

test_that("terms that store fitted state predict on the basis they were fitted on", {
  # ns(), bs(), poly() and scale() all derive something from the data -- knots,
  # orthogonalising coefficients, a centre and a scale -- and predict() must
  # reuse it. Recomputing from the rows predict() happens to be handed makes
  # ns() error and poly() and scale() quietly answer the wrong number.
  d <- panel_seas(47)
  nd <- data.frame(x = c(-1, 0, 1), id = factor(1, levels = levels(d$id)))
  for (trm in c("splines::ns(x, df = 3)", "poly(x, 2)", "scale(x)")) {
    f <- ilm_model(stats::as.formula(paste("y ~", trm, "+ (1 | id)")), data = d,
                   family = "gaussian", verbose = FALSE)
    expect_false(is.null(attr(f$terms, "predvars")), info = trm)
    p1 <- as.numeric(predict(f, newdata = nd))
    p2 <- as.numeric(predict(f, newdata = nd[2, , drop = FALSE]))
    expect_equal(p1[2], p2, tolerance = 1e-8, label = trm)
  }
  f <- ilm_model(y ~ splines::ns(x, df = 3) + (1 | id), data = d,
                 family = "gaussian", verbose = FALSE)
  expect_equal(ilm_anova(f, test = "Wald")$Df[1], 3L)
})
