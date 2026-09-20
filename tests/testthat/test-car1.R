# CAR(1): a first-order autoregression over time that need not be evenly spaced.
#
# The load-bearing test here is that CAR(1) and AR(1) agree EXACTLY when every
# gap is one step. They are written as separate branches of the likelihood, and
# if the continuous one is wrong in its innovation variance or its log
# determinant, that is where it shows.

panel_even <- function(seed, ng = 40, nt = 8, rho = 0.7) {
  set.seed(seed)
  d <- data.frame(id = factor(rep(seq_len(ng), each = nt)),
                  t = rep(seq_len(nt), ng), x = stats::rnorm(ng * nt))
  u <- unlist(lapply(seq_len(ng), function(i) {
    z <- numeric(nt); z[1] <- stats::rnorm(1)
    for (k in 2:nt)
      z[k] <- rho * z[k - 1] + stats::rnorm(1, 0, sqrt(1 - rho^2))
    z
  }))
  d$y <- 0.5 + 0.8 * d$x + 1.2 * u + stats::rnorm(ng * nt, 0, 0.6)
  d
}

panel_irreg <- function(seed, rho = 0.85, ng = 60, nt = 6, tmax = 30) {
  set.seed(seed)
  tl <- lapply(seq_len(ng), function(i) sort(sample.int(tmax, nt)))
  d <- data.frame(id = factor(rep(seq_len(ng), each = nt)),
                  t = unlist(tl), x = stats::rnorm(ng * nt))
  u <- unlist(lapply(tl, function(tt) {
    z <- numeric(length(tt)); z[1] <- stats::rnorm(1)
    for (k in seq_along(tt)[-1]) {
      ph <- rho^(tt[k] - tt[k - 1])
      z[k] <- ph * z[k - 1] + stats::rnorm(1, 0, sqrt(1 - ph^2))
    }
    z
  }))
  b <- stats::rnorm(ng, 0, 0.5)[as.integer(d$id)]
  d$y <- 0.5 + 0.8 * d$x + b + 1.2 * u + stats::rnorm(ng * nt, 0, 0.6)
  d
}

## ---- the index arithmetic --------------------------------------------------

test_that("the continuous index maps observations to shared latent values", {
  d <- data.frame(id = rep(1:3, each = 4),
                  t = c(0, 2, 7, 15, 0, 3, 4, 20, 1, 6, 9, 11))
  s <- ilm_car1(d$t, d$id, verbose = FALSE)
  expect_equal(s$idx, 1:12)             # all distinct, in group-then-time order
  expect_equal(s$first, c(1L, 5L, 9L))  # one per group
  expect_equal(s$prev, c(1:3, 5:7, 9:11))
  expect_equal(s$gap, c(2, 5, 8, 3, 1, 16, 5, 3, 2))
  expect_equal(s$n_cell, 12L)
  # transitions never cross a group boundary
  expect_true(all(s$gap > 0))
})

test_that("observations at the same time in the same group share one latent", {
  d <- data.frame(id = c(1, 1, 1, 2, 2), t = c(5, 5, 9, 5, 9))
  s <- ilm_car1(d$t, d$id, verbose = FALSE)
  expect_equal(s$idx, c(1L, 1L, 2L, 3L, 4L))
  expect_equal(s$n_cell, 4L)
  expect_equal(s$gap, c(4, 4))
  expect_equal(s$obs_per_latent, 5 / 4)
})

test_that("the evenly spaced index spans missed steps", {
  d <- data.frame(id = rep(1:2, each = 3), t = c(1, 2, 4, 1, 3, 4))
  s <- ilm_ar1(d$t, d$id, verbose = FALSE)
  # a grid of 4 slots per group, so the chain continues across the gap
  expect_equal(s$Tt, 4L)
  expect_equal(s$idx, c(1L, 2L, 4L, 5L, 7L, 8L))
  expect_equal(s$n_cell, 8L)
  # a common step other than 1 is rescaled, not rejected
  s2 <- ilm_ar1(c(0, 5, 10, 0, 5, 10), rep(1:2, each = 3), verbose = FALSE)
  expect_equal(s2$step, 5)
  expect_equal(s2$Tt, 3L)
})

## ---- the load-bearing agreement --------------------------------------------

test_that("CAR(1) reduces exactly to AR(1) when every gap is one", {
  d <- panel_even(7)
  fa <- ilm_model(y ~ x + (1 | id), data = d, family = "gaussian",
                  ar = suppressWarnings(ilm_ar1(d$t, d$id)), verbose = FALSE)
  fc <- ilm_model(y ~ x + (1 | id), data = d, family = "gaussian",
                  ar = suppressWarnings(ilm_car1(d$t, d$id)), verbose = FALSE)
  # the two likelihoods are written separately; on equal spacing they describe
  # the same process and must agree to optimiser tolerance, not merely closely
  expect_equal(as.numeric(logLik(fa)), as.numeric(logLik(fc)), tolerance = 1e-6)
  expect_equal(fa$rho, fc$rho, tolerance = 1e-4)
  expect_equal(unname(coef(fa)), unname(coef(fc)), tolerance = 1e-5)
})

test_that("the bare list the fitter took before the constructors still works", {
  d <- panel_even(7)
  ng <- nlevels(d$id); nt <- max(d$t)
  f1 <- ilm_model(y ~ x + (1 | id), data = d, family = "gaussian",
                  ar = list(idx = (as.integer(d$id) - 1L) * nt + d$t,
                            n_group = ng, Tt = nt), verbose = FALSE)
  f2 <- ilm_model(y ~ x + (1 | id), data = d, family = "gaussian",
                  ar = suppressWarnings(ilm_ar1(d$t, d$id)), verbose = FALSE)
  expect_equal(f1$rho, f2$rho, tolerance = 1e-6)
  expect_equal(as.numeric(logLik(f1)), as.numeric(logLik(f2)), tolerance = 1e-8)
})

## ---- recovery --------------------------------------------------------------

test_that("CAR(1) recovers the correlation from irregular spacing", {
  d <- panel_irreg(3, rho = 0.85)
  f <- ilm_model(y ~ x + (1 | id), data = d, family = "gaussian",
                 ar = suppressWarnings(ilm_car1(d$t, d$id)), verbose = FALSE)
  # measured over 12 replicates of this design: mean 0.845, sd 0.035
  expect_equal(f$rho, 0.85, tolerance = 0.15)
  expect_equal(unname(coef(f)[["x"]]), 0.8, tolerance = 0.25)
  # the range is the same thing on a different scale
  expect_equal(f$ar_range, -1 / log(f$rho), tolerance = 1e-8)
})

test_that("the simulator walks the same chain the likelihood scores", {
  # every envelope diagnostic is built on ilm_simulate(); if it walked a
  # different process the envelopes would be calibrated against the wrong null
  d <- panel_irreg(4, rho = 0.9)
  f <- ilm_model(y ~ x + (1 | id), data = d, family = "gaussian",
                 ar = suppressWarnings(ilm_car1(d$t, d$id)), verbose = FALSE)
  ys <- ilm_simulate(f, 40, seed = 1)
  expect_equal(dim(ys), c(nrow(d), 40L))
  expect_true(all(is.finite(ys)))

  # correlation between consecutive observations of a unit should fall away as
  # the gap widens, which is the whole content of the structure
  s <- suppressWarnings(ilm_car1(d$t, d$id, verbose = FALSE))
  near <- s$rest[s$gap <= 2]; near_p <- s$prev[s$gap <= 2]
  far  <- s$rest[s$gap >= 8]; far_p  <- s$prev[s$gap >= 8]
  skip_if(length(near) < 20 || length(far) < 20)
  cell1 <- match(seq_len(s$n_cell), s$idx)      # one observation per cell
  cr <- function(a, b) {
    i <- cell1[a]; j <- cell1[b]
    ok <- is.finite(i) & is.finite(j)
    mean(vapply(seq_len(ncol(ys)), function(k)
      stats::cor(ys[i[ok], k], ys[j[ok], k]), 0))
  }
  expect_gt(cr(near, near_p), cr(far, far_p))
})

## ---- informative failure ---------------------------------------------------

test_that("irregular times through ilm_ar1 are named and redirected", {
  d <- data.frame(id = rep(1:3, each = 3), t = rep(c(0, 1, 2.5), 3))
  expect_error(ilm_ar1(d$t, d$id), "do not sit on a common grid")
  expect_error(ilm_ar1(d$t, d$id), "ilm_car1", fixed = TRUE)
  # a grid that is technically common but mostly empty is a warning, since the
  # empty slots are latent values with no data behind them
  d2 <- data.frame(id = rep(1:3, each = 3), t = rep(c(0, 1, 40), 3))
  w <- character(0)
  withCallingHandlers(ilm_ar1(d2$t, d2$id),
    warning = function(cnd) { w <<- c(w, conditionMessage(cnd))
                              invokeRestart("muffleWarning") })
  # the empty-grid warning is the actionable one; the thin-budget warning is
  # its consequence and both are worth saying
  expect_true(any(grepl("steps are ever observed", w)))
  expect_true(any(grepl("per latent", w)))
})

test_that("misspecified correlation structures are named", {
  g <- rep(1:3, each = 3); t <- rep(1:3, 3)
  for (f in list(ilm_ar1, ilm_car1)) {
    expect_error(f(t), "both required")
    expect_error(f(t, g[-1]), "same length")
    expect_error(f(factor(t), g), "is a factor")
    expect_error(f(letters[t], g), "must be numeric")
    expect_error(f(c(NA, t[-1]), g), "cannot contain missing values")
    expect_error(f(rep(1, 9), g), "only one distinct value")
  }
  # one time per unit is a random intercept, not a correlation over time
  expect_error(ilm_car1(1:3, 1:3), "random intercept")
  # a spec built from different rows than the model is fitted to
  d <- panel_even(1, ng = 5, nt = 4)
  s <- suppressWarnings(ilm_car1(d$t, d$id))
  s$idx <- s$idx[-1]
  expect_error(ilm_model(y ~ x + (1 | id), data = d, family = "gaussian",
                         ar = s, verbose = FALSE),
               "covers 19 observations but the model matrix has 20 rows")
  expect_error(ilm_model(y ~ x + (1 | id), data = d, family = "gaussian",
                         ar = list(nope = 1), verbose = FALSE),
               "must come from ilm_ar1")
})

test_that("a thin latent budget is warned about by the constructor", {
  d <- panel_irreg(1, ng = 20, nt = 5, tmax = 50)
  expect_warning(ilm_car1(d$t, d$id), "per latent")
  expect_warning(ilm_car1(d$t, d$id), "coarser grid")
  # and not warned about when observations genuinely share latent values
  d2 <- data.frame(id = rep(1:4, each = 6), t = rep(rep(1:3, each = 2), 4))
  expect_silent(ilm_car1(d2$t, d2$id))
})

test_that("printing a structure says what it is and how thin it is", {
  d <- panel_irreg(1, ng = 10, nt = 4, tmax = 20)
  out <- capture.output(print(suppressWarnings(ilm_car1(d$t, d$id))))
  expect_true(any(grepl("CAR(1), continuous time", out, fixed = TRUE)))
  expect_true(any(grepl("gaps: min", out)))
  expect_true(any(grepl("latent budget is thin", out)))
  out2 <- capture.output(print(suppressWarnings(ilm_ar1(d$t, d$id))))
  expect_true(any(grepl("AR(1), evenly spaced", out2, fixed = TRUE)))
  expect_true(any(grepl("step:", out2)))
})
