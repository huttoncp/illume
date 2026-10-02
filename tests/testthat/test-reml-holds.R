## The hold rules on the default path. A gaussian mixed model is fitted by
## REML by default (Craig's item 249), and the hold rules' own tests pin
## maximum likelihood, where their constructed boundaries were found. These
## pin the same rules under REML, on boundaries REML reaches: REML estimates
## variances larger, and several of the ML cases do not reach a boundary
## under it at all (the hold-rule study, studies/findings/reml_hold.md).

q <- function(e) suppressMessages(suppressWarnings(e))
ck_of <- function(f, nm) f$checks[f$checks$check == nm, , drop = FALSE]

## a random intercept with no variation to find: y and x double-centred on
## two crossed factors, so every group's mean is the grand mean
zero_g_data <- function() {
  set.seed(1)
  n_g <- 20; per <- 10
  d <- data.frame(g = factor(rep(seq_len(n_g), each = per)), h = factor(rep(seq_len(per), n_g)),
                  x = stats::rnorm(n_g * per), e = stats::rnorm(n_g * per))
  dc <- function(v) v - stats::ave(v, d$g) - stats::ave(v, d$h) + mean(v)
  d$x <- dc(d$x); d$y <- 1 + 0.5 * d$x + dc(d$e)
  d
}

test_that("under REML a random intercept at zero is held, with the reduced model's SEs", {
  d <- zero_g_data()
  f <- q(ilm_model(y ~ x + (1 | g), data = d, family = "gaussian", verbose = FALSE))
  expect_true(isTRUE(f$reml))
  expect_true("g" %in% f$hessian_held)
  expect_identical(f$hessian_how, "boundary")
  ## variance_boundary warns: the SD reached zero, well below the hold line
  ## (an SD held above the old 1e-3 line would be named "held at zero")
  vb <- ck_of(f, "variance_boundary")
  expect_identical(vb$status, "WARN")
  expect_lt(sqrt(ilm_varcorr(f)$re$g[1, 1]), 1e-3)
  ## the gradient along the held direction is reported, not judged
  expect_match(ck_of(f, "gradient")$detail, "along the held g", fixed = TRUE)
  expect_identical(ck_of(f, "gradient")$status, "OK")
  ## rule C: the fixed effects' SEs are the reduced model's, which here is the
  ## linear model (REML's residual variance is RSS / (n - p), as lm's is)
  f0 <- q(ilm_model(y ~ x, data = d, family = "gaussian", verbose = FALSE))
  expect_equal(sqrt(diag(vcov(f)))[["x"]], sqrt(diag(vcov(f0)))[["x"]], tolerance = 1e-6)
  ## the line is relative to sd(y): the same data 1,000 times larger hold too
  d2 <- d; d2$y <- 1000 * d$y
  f2 <- q(ilm_model(y ~ x + (1 | g), data = d2, family = "gaussian", verbose = FALSE))
  expect_true(isTRUE(f2$reml) && "g" %in% f2$hessian_held)
  ## and the draws hold it: the held variance is at its estimate in every draw
  dr <- ilm_draws(f, nsim = 100, seed = 1)
  th <- dr$draws[rownames(dr$draws) == "theta", , drop = FALSE]
  expect_true(all(th == th[1]))
  expect_true(all(is.finite(dr$draws)))
})

test_that("under REML boundary = 'avoid' keeps that intercept off zero", {
  d <- zero_g_data()
  f <- q(ilm_model(y ~ x + (1 | g), data = d, family = "gaussian", boundary = "avoid",
                   verbose = FALSE))
  expect_true(isTRUE(f$reml))
  expect_false("g" %in% f$hessian_held)
  expect_gt(sqrt(ilm_varcorr(f)$re$g[1, 1]), 1e-3)
  f0 <- q(ilm_model(y ~ x + (1 | g), data = d, family = "gaussian", verbose = FALSE))
  ## the penalty moves the variance, and the slope hardly at all
  expect_equal(unname(coef(f)[["x"]]), unname(coef(f0)[["x"]]), tolerance = 1e-3)
})

test_that("under REML a residual SD at its limit is held, and its gradient not judged", {
  ## one AR(1) series of 12 whose latent took up the noise under REML: sigma
  ## runs to its limit, where under maximum likelihood the AR SD did
  ar_series <- function(Tn, rho, sd_marg)
    as.numeric(stats::arima.sim(list(ar = rho), Tn, sd = sd_marg * sqrt(1 - rho^2)))
  set.seed(1269247)
  d <- data.frame(g = factor("s1"), t = seq_len(12))
  d$y <- 5 + ar_series(12, 0.7, 0.8) + stats::rnorm(12, 0, 0.5)
  f <- q(ilm_model(y ~ 1, data = d, family = "gaussian", ar = ilm_ar1(~ t | g), verbose = FALSE))
  expect_true(isTRUE(f$reml))
  expect_true("dispersion" %in% f$hessian_held)
  expect_lt(exp(f$opt$par[["logdisp"]]) / sd(d$y), 0.2)
  dl <- ck_of(f, "dispersion_limit")
  expect_identical(dl$status, "BOUNDARY")
  expect_match(dl$detail, "held at its estimate", fixed = TRUE)
  expect_match(ck_of(f, "gradient")$detail, "along the held dispersion", fixed = TRUE)
  ## the draws hold it
  dr <- ilm_draws(f, nsim = 100, seed = 1)
  ld <- dr$draws[dr$map$block == "logdisp", ]
  expect_equal(diff(range(ld)), 0)
  expect_true(all(is.finite(dr$draws)))
})

test_that("under REML an AR block held for its correlation keeps its SD in the draws", {
  ## white noise in four series of 24: the AR correlation has nothing to find
  ## and sits at an edge, while its SD is identified
  set.seed(3)
  d <- expand.grid(t = seq_len(24), g = factor(seq_len(4)))
  d$x <- stats::rnorm(nrow(d)); d$y <- 1 + 0.5 * d$x + stats::rnorm(nrow(d))
  f <- q(ilm_model(y ~ x, data = d, family = "gaussian", ar = ilm_ar1(~ t | g), verbose = FALSE))
  expect_true(isTRUE(f$reml))
  expect_true("ar" %in% f$hessian_held)
  expect_gt(sqrt(as.matrix(f$Sigma[["ar"]])[1, 1]), 0.1)
  expect_false(illume:::ilm_ar_sd_at_zero(f))
  dr <- ilm_draws(f, nsim = 200, seed = 1)
  expect_true(all(is.finite(dr$draws)))
  expect_gt(stats::sd(dr$draws[dr$map$block == "lchol_ar", ]), 0.01)
})
