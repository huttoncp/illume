## The gaussian hold lines as ruled, and the checks agreeing with the holds.
##
## sigma is held below 0.2 of sd(y), a gaussian AR(1)'s SD and a random
## effect's SD below 0.1 of it, each only where the likelihood is flat below
## (studies/scripts/sigma_limit.R, and re_sd_limit.R's fits rebuilt from their
## seeds). Each case below is a fit those studies recorded, rebuilt by the
## study's own generator from its seed.

ar_series <- function(Tn, rho, sd_marg)
  as.numeric(stats::arima.sim(list(ar = rho), Tn, sd = sd_marg * sqrt(1 - rho^2)))
q <- function(e) suppressMessages(suppressWarnings(e))
ck_of <- function(f, nm) f$checks[f$checks$check == nm, , drop = FALSE]

test_that("sigma short of zero is held below 0.2 of sd(y)", {
  ## 20 AR(1) series of 12 with little noise: sigma stops at 0.046 of sd(y),
  ## far above the old 1e-3 line, on a flat likelihood
  set.seed(1443467)
  d <- expand.grid(t = seq_len(12), g = factor(seq_len(20)))
  d$y <- unlist(lapply(seq_len(20), function(i) ar_series(12, 0.5, 1))) +
    stats::rnorm(nrow(d), 0, 0.1)
  f <- q(ilm_model(y ~ 1, data = d, family = "gaussian", ar = ilm_ar1(~ t | g),
                   verbose = FALSE))
  expect_lt(exp(f$opt$par[["logdisp"]]) / sd(d$y), 0.2)
  expect_true("dispersion" %in% f$hessian_held)
  expect_match(ck_of(f, "dispersion_limit")$detail, "held at its estimate", fixed = TRUE)
})

test_that("a gaussian AR(1) SD short of zero is held below 0.1 of sd(y)", {
  ## one series of 12 whose noise took up the latent: the AR SD stops at
  ## 0.014 of sd(y) on a flat likelihood, and sigma is well determined
  set.seed(1269247)
  d <- data.frame(g = factor("s1"), t = seq_len(12))
  d$y <- 5 + ar_series(12, 0.7, 0.8) + stats::rnorm(12, 0, 0.5)
  f <- q(ilm_model(y ~ 1, data = d, family = "gaussian", ar = ilm_ar1(~ t | g),
                   verbose = FALSE))
  expect_true("ar" %in% f$hessian_held)
  expect_false("dispersion" %in% f$hessian_held)
  ## and variance_boundary says so, rather than reading the stopping point
  vb <- ck_of(f, "variance_boundary")
  expect_identical(vb$status, "WARN")
  expect_match(vb$detail, "held at zero: ar (the optimiser stopped at", fixed = TRUE)
})

test_that("a random-intercept SD's line is relative to sd(y) for gaussian", {
  ## re_sd_limit's fit: 5 groups of 10, an intercept SD stopped at 0.022
  set.seed(328503)
  d <- data.frame(g = factor(rep(seq_len(5), each = 10L)))
  lat <- stats::rnorm(5, 0, 0.1)[d$g]
  d$x <- stats::rnorm(nrow(d))
  d$y <- 0.5 + 0.3 * d$x + lat + stats::rnorm(nrow(d), 0, 0.5)
  f1 <- q(ilm_model(y ~ x + (1 | g), data = d, family = "gaussian", verbose = FALSE))
  expect_true("g" %in% f1$hessian_held)
  ## the same data in units 1,000 times smaller: the likelihood's flatness
  ## is the same, and so is the hold; an absolute line would have missed it
  d2 <- d; d2$y <- 1000 * d$y
  f2 <- q(ilm_model(y ~ x + (1 | g), data = d2, family = "gaussian", verbose = FALSE))
  expect_true("g" %in% f2$hessian_held)
  for (f in list(f1, f2)) {
    vb <- ck_of(f, "variance_boundary")
    expect_identical(vb$status, "WARN")
    expect_match(vb$detail, "held at zero: g", fixed = TRUE)
  }
})

test_that("a held random intercept is at zero in variance_boundary (Poisson, beside an AR)", {
  ## re_sd_limit's fit: the intercept SD held at 0.0117, which the check read
  ## as OK while the hessian check said held
  set.seed(189252)
  G <- 20L; Tn <- 24L
  d <- expand.grid(t = seq_len(Tn), g = factor(seq_len(G)))
  ar <- unlist(lapply(seq_len(G), function(i) as.numeric(
    stats::arima.sim(list(ar = 0.8), Tn, sd = 0.5 * sqrt(1 - 0.64)))))
  lat <- ar + stats::rnorm(G, 0, 0.3)[d$g]
  d$x <- stats::rnorm(nrow(d))
  d$y <- stats::rpois(nrow(d), exp(log(3) + 0.3 * d$x + lat))
  f <- q(ilm_model(y ~ x + (1 | g), data = d, family = "poisson",
                   ar = ilm_ar1(~ t | g), verbose = FALSE))
  expect_true("g" %in% f$hessian_held)
  vb <- ck_of(f, "variance_boundary")
  expect_identical(vb$status, "WARN")
  expect_match(vb$detail, "held at zero: g", fixed = TRUE)
})

test_that("with sigma and an AR SD both flat, only the flatter is held", {
  ## the policy on its own: noise and latent interchangeable, each flat alone
  ## while their sum is determined, so holding both would set it to zero
  par <- c(beta = 1, logdisp = -3, lchol_ar = -3, rho_raw = 0.1)
  push <- c(lchol_ar = 4e-4)
  obj <- list(fn = function(p) {
    d <- p - par
    if (d[["lchol_ar"]] != 0) push[["lchol_ar"]] else 0
  })
  mk <- function(disp_push) list(
    flagged = "dispersion", force = character(0), maybe = "ar",
    blocks = list(dispersion = 2L, ar = 3:4), keep = 1L,
    push = list(dispersion = disp_push))
  ## sigma's push moves the objective less: sigma is held, the AR SD is not
  cb <- illume:::ilm_re_flat(obj, par, mk(1e-4))
  expect_true("dispersion" %in% cb$flagged)
  expect_false("ar" %in% cb$flagged)
  ## the AR SD's push moves it less: the AR SD is held, sigma is not
  cb <- illume:::ilm_re_flat(obj, par, mk(9e-4))
  expect_true("ar" %in% cb$flagged)
  expect_false("dispersion" %in% cb$flagged)
  expect_true(2L %in% cb$keep)
})
