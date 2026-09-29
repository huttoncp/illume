## Draws of an AR latent held with its SD at zero.
##
## The fit holds such a block, and its standard errors say so; the draws
## held only the directions the Hessian called flat. With the SD at zero the
## correlation is unidentified too, and when only one of the two was among
## those directions the other was drawn from a variance in the thousands: in
## this case log SDs from -279 to 165, so a held latent's SD reached e^165 in
## the draws and every forecast built on them. The whole block is held now.

held_ar_fit <- function() {
  ## a weak correlation under small noise: the residual takes up the latent
  set.seed(1095029)
  d <- data.frame(g = factor("s1"), t = 1:48)
  d$y <- 5 + as.numeric(stats::arima.sim(list(ar = 0.1), 48, sd = 0.8 * sqrt(1 - 0.01))) +
    stats::rnorm(48, 0, 0.1)
  suppressMessages(suppressWarnings(
    ilm_model(reml = FALSE, y ~ 1, data = d, family = "gaussian", ar = ilm_ar1(~ t | g),
              verbose = FALSE)))
}

test_that("an AR latent held with its SD at zero is held whole in the draws", {
  f <- held_ar_fit()
  expect_true("ar" %in% f$hessian_held)
  expect_lt(sqrt(as.matrix(f$Sigma[["ar"]])[1, 1]), 1e-3)
  dr <- ilm_draws(f, nsim = 200, seed = 1)
  ar <- dr$draws[dr$map$block %in% c("lchol_ar", "rho_raw"), , drop = FALSE]
  expect_true(all(is.finite(dr$draws)))
  ## at the mode, every draw
  expect_equal(unname(ar), unname(matrix(dr$mode[dr$map$block %in% c("lchol_ar", "rho_raw")],
                                         nrow(ar), ncol(ar))), tolerance = 1e-12)
  expect_gte(dr$held$n, 2L)
  ## the rest is still drawn: the intercept and the residual SD vary
  expect_gt(stats::sd(dr$draws[dr$map$block == "beta", ]), 0)
  expect_gt(stats::sd(dr$draws[dr$map$block == "logdisp", ]), 0)
})

test_that("a fit made before the held directions were stored draws as it did", {
  f <- held_ar_fit()
  f$hessian_dirs <- NULL
  dr <- ilm_draws(f, nsim = 200, seed = 1)
  ar <- dr$draws[dr$map$block %in% c("lchol_ar", "rho_raw"), , drop = FALSE]
  expect_true(all(apply(ar, 1L, function(v) diff(range(v))) == 0))
  expect_true(all(is.finite(dr$draws)))
  ## and the two routes agree where the AR block is concerned
  dn <- ilm_draws(held_ar_fit(), nsim = 200, seed = 1)
  expect_equal(unname(ar), unname(dn$draws[dn$map$block %in% c("lchol_ar", "rho_raw"), , drop = FALSE]))
})

test_that("an AR block held for its correlation, with an identified SD, keeps its SD", {
  ## a healthy AR fit, held as one is when only its correlation is at an
  ## edge: the one flat direction along the correlation, the SD identified
  set.seed(7)
  d <- data.frame(g = factor(rep(1:10, each = 30)), t = rep(1:30, 10))
  d$y <- unlist(lapply(1:10, function(i) as.numeric(stats::arima.sim(list(ar = 0.6), 30)))) +
    stats::rnorm(300, 0, 0.5)
  f <- suppressMessages(suppressWarnings(
    ilm_model(reml = FALSE, y ~ 1, data = d, family = "gaussian", ar = ilm_ar1(~ t | g),
              verbose = FALSE)))
  expect_gt(sqrt(as.matrix(f$Sigma[["ar"]])[1, 1]), 0.1)
  pn <- names(f$opt$par)
  f$hessian_held <- "ar"
  f$hessian_dirs <- matrix(as.numeric(pn == "rho_raw"), ncol = 1L)
  expect_false(illume:::ilm_ar_sd_at_zero(f))
  dr <- ilm_draws(f, nsim = 200, seed = 1)
  ## the correlation held along its direction, the SD still drawn
  expect_equal(stats::sd(dr$draws[dr$map$block == "rho_raw", ]), 0, tolerance = 1e-10)
  expect_gt(stats::sd(dr$draws[dr$map$block == "lchol_ar", ]), 0.01)
  expect_identical(dr$held$n, 1L)
})

## A panel with no noise, whose fit stopped short of a stationary point
## (gradient 1.2e-2 on the slope) with sigma at 1e-4 of sd(y): the
## recomputed covariance is refused off a stationary point, so the flagged
## dispersion was not held, and its draws, from an SE of 14 on the log scale,
## reached e^42 -- while the checks said it was held.
stalled_fit <- function() {
  set.seed(1293099)
  G <- 20L
  d <- expand.grid(t = 1:12, g = factor(seq_len(G)))
  lat <- unlist(lapply(seq_len(G), function(i)
    as.numeric(stats::arima.sim(list(ar = 0.7), 12, sd = 0.8 * sqrt(1 - 0.49)))))
  d$x <- stats::rnorm(nrow(d))
  d$y <- 1 + 0.3 * d$x + stats::rnorm(G, 0, 0.5)[d$g] + lat + stats::rnorm(nrow(d), 0, 0)
  suppressMessages(suppressWarnings(
    ilm_model(reml = FALSE, y ~ x + (1 | g), data = d, family = "gaussian",
              ar = ilm_ar1(~ t | g), verbose = FALSE)))
}

test_that("a dispersion at its limit is held where the fit did not converge", {
  f <- stalled_fit()
  ck <- f$checks
  expect_identical(ck$status[ck$check == "gradient"], "FAIL")
  expect_true("dispersion" %in% f$hessian_held)
  ## what the checks say is true
  expect_match(ck$detail[ck$check == "dispersion_limit"], "held at its estimate",
               fixed = TRUE)
  ## the draws hold it, and say once that the fit did not converge
  expect_warning(dr <- ilm_draws(f, nsim = 200, seed = 1),
                 "did not reach a stationary point (optimizer and gradient checks FAIL)",
                 fixed = TRUE)
  ld <- dr$draws[dr$map$block == "logdisp", ]
  expect_equal(diff(range(ld)), 0)
  expect_true(all(is.finite(dr$draws)))
})

test_that("a converged fit draws without the warning", {
  expect_silent(ilm_draws(held_ar_fit(), nsim = 20, seed = 1))
})

test_that("a smooth's band beside a held AR is the band without it", {
  ## the joint draws behind a typical or population band now hold what the
  ## fit holds; the AR block at its zero SD drew from -114 to 135 there
  set.seed(1095029)
  d <- data.frame(g = factor("s1"), t = 1:48)
  d$x <- stats::runif(48, -2, 2)
  d$y <- 5 + sin(d$x) + as.numeric(stats::arima.sim(list(ar = 0.1), 48, sd = 0.8 * sqrt(1 - 0.01))) +
    stats::rnorm(48, 0, 0.1)
  fa <- suppressMessages(suppressWarnings(ilm_model(reml = FALSE, y ~ s(x), data = d, family = "gaussian",
          ar = ilm_ar1(~ t | g), verbose = FALSE)))
  f0 <- suppressMessages(suppressWarnings(ilm_model(reml = FALSE, y ~ s(x), data = d, family = "gaussian",
          verbose = FALSE)))
  expect_true("ar" %in% fa$hessian_held)
  jd <- illume:::ilm_joint_draws(fa, 500, 2)
  k <- jd$which %in% c("lchol_ar", "rho_raw")
  expect_equal(max(apply(jd$draws[k, , drop = FALSE], 1L, function(v) diff(range(v)))), 0)
  nd <- data.frame(x = seq(-2, 2, length.out = 9), t = 49, g = "s1")
  w <- function(f) { p <- suppressWarnings(predict(f, newdata = nd, groups = "typical",
                       interval = "confidence", nsim = 2000, seed = 1))
                     as.numeric(p$upper - p$lower) }
  r <- w(fa) / w(f0)
  expect_true(all(r > 0.9 & r < 1.1))
})
