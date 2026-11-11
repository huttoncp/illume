## A factor-smooth interaction, s(x, id, bs = "fs"), checked against mgcv's
## own fit of the same formula: mgcv::gam(), with method = "REML" or "ML",
## and its results read through stats::fitted(), the fit's `sp` and
## `gcv.ubre` elements, and `reml.scale` (the REML or ML scale; under ML,
## `sig2` is a different estimate).
##
## illume fits the smooth as variance components, one per penalty, so each
## component should be the scale divided by the matching smoothing parameter,
## and illume's -logLik should be gam()'s REML or ML score at the optimum.

fs_data <- function(seed, n_id = 30, n_per = 10, pois = FALSE, unbal = FALSE) {
  set.seed(seed)
  sizes <- if (unbal) sample(3:18, n_id, replace = TRUE) else rep(n_per, n_id)
  id <- factor(rep(seq_len(n_id), times = sizes))
  x <- stats::runif(length(id))
  a <- stats::rnorm(n_id, 0, 0.6)[id]; b <- stats::rnorm(n_id, 0, 0.8)[id]
  w <- stats::rnorm(n_id, 0, 0.3)[id]
  eta <- sin(2 * pi * x) + a + b * (x - 0.5) + w * cos(3 * pi * x)
  y <- if (pois) stats::rpois(length(x), exp(0.5 + 0.5 * eta)) else
    eta + stats::rnorm(length(x), 0, 0.4)
  data.frame(y = y, x = x, id = id)
}

## illume's fit and gam()'s, side by side
fs_compare <- function(form, d, family = "gaussian", reml = TRUE) {
  gfam <- if (family == "gaussian") stats::gaussian() else stats::poisson()
  g <- suppressWarnings(mgcv::gam(form, data = d, family = gfam,
                                  method = if (reml) "REML" else "ML"))
  f <- suppressMessages(suppressWarnings(
    ilm_model(form, data = d, family = family, reml = reml, verbose = FALSE)))
  vc <- as.data.frame(nlme::VarCorr(f))
  s2 <- if (family == "gaussian") vc$vcov[vc$grp == "Residual"] else 1
  list(fitted = max(abs(as.numeric(ilm_fitted(f)) - as.numeric(stats::fitted(g)))),
       score = c(illume = -as.numeric(logLik(f)), gam = as.numeric(g$gcv.ubre)),
       log_sp = max(abs(log(s2 / vc$vcov[vc$grp != "Residual"]) - log(g$sp))),
       sigma = c(illume = sqrt(s2), gam = sqrt(g$reml.scale)))
}

test_that("a gaussian fs smooth is the model mgcv::gam() fits, by REML and by ML", {
  skip_on_cran()
  d <- fs_data(1)
  for (reml in c(TRUE, FALSE)) {
    r <- fs_compare(y ~ s(x, id, bs = "fs", k = 5), d, reml = reml)
    expect_lt(r$fitted, 1e-5)
    expect_equal(r$score[["illume"]], r$score[["gam"]], tolerance = 1e-7)
    expect_lt(r$log_sp, 1e-3)
    expect_equal(r$sigma[["illume"]], r$sigma[["gam"]], tolerance = 1e-5)
  }
})

test_that("it agrees with mgcv::gam() beside a common curve, unbalanced, and on a cr basis", {
  skip_on_cran()
  cases <- list(
    list(d = fs_data(5, unbal = TRUE),
         form = y ~ s(x) + s(x, id, bs = "fs", k = 5, m = 1)),
    list(d = fs_data(4, 12, 6), form = y ~ s(x, id, bs = "fs", k = 4)),
    list(d = fs_data(6),
         form = y ~ s(x, id, bs = "fs", k = 6, xt = list(bs = "cr"))))
  for (z in cases) {
    r <- fs_compare(z$form, z$d)
    expect_lt(r$fitted, 1e-5)
    expect_equal(r$score[["illume"]], r$score[["gam"]], tolerance = 1e-7)
    expect_lt(r$log_sp, 1e-3)
  }
})

test_that("a Poisson fs smooth by REML is the model mgcv::gam() fits", {
  skip_on_cran()
  r <- fs_compare(y ~ s(x, id, bs = "fs", k = 5), fs_data(2, pois = TRUE),
                  family = "poisson", reml = TRUE)
  expect_lt(r$fitted, 1e-3)
  expect_equal(r$score[["illume"]], r$score[["gam"]], tolerance = 1e-7)
  expect_lt(r$log_sp, 1e-3)
})
