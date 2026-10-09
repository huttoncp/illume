## A fit whose fixed-effect standard errors cannot be computed fails its
## Hessian check and reports them as NaN. A term held at its boundary used to
## let it pass: the check judged only the variances that were finite, and
## vcov() set every non-finite entry to 0, so the fit graded ok with standard
## errors of 0. Holding a term means zero uncertainty for the held covariance
## directions only (rule C, item 4), never for a fixed effect.

q <- function(e) suppressMessages(suppressWarnings(e))
ck_of <- function(f, nm) f$checks[f$checks$check == nm, , drop = FALSE]

## the earlier rescaling study's random-slope arm (rescale.R, R6, replicate
## 1): a slope on income in dollars at its boundary, where the fixed
## effects' covariance comes back NaN once the term is held
r6_data <- function(rep = 1) {
  set.seed(1e4 * 11 + rep); n <- 600; G <- 25
  d <- data.frame(g = factor(sample(G, n, TRUE)), income = round(stats::rnorm(n, 33000, 9000)),
                  tiny = stats::rnorm(n, 0, 1e-4), x = stats::runif(n))
  u <- stats::rnorm(G, 0, 0.4)[d$g]
  lin <- 0.00003 * (d$income - 33000) + 2000 * d$tiny
  lin <- lin + stats::rnorm(G, 0, 0.00001)[d$g] * (d$income - 33000)
  d$y <- stats::rpois(n, exp(0.4 + lin + u))
  d
}

## a random intercept with no variation to find (test-reml-holds.R): a held
## term whose fixed effects are healthy
zero_g_data <- function() {
  set.seed(1)
  n_g <- 20; per <- 10
  d <- data.frame(g = factor(rep(seq_len(n_g), each = per)), h = factor(rep(seq_len(per), n_g)),
                  x = stats::rnorm(n_g * per), e = stats::rnorm(n_g * per))
  dc <- function(v) v - stats::ave(v, d$g) - stats::ave(v, d$h) + mean(v)
  d$x <- dc(d$x); d$y <- 1 + 0.5 * d$x + dc(d$e)
  d
}

test_that("a held term with fixed-effect variances that cannot be formed fails, and says NaN", {
  d <- r6_data(1)
  f <- q(ilm_model(y ~ income + (1 + income | g), data = d, family = "poisson", verbose = FALSE))
  ## the case: the term is held, and the fixed effects' variances are not finite
  expect_true(length(f$hessian_held) > 0L)
  expect_false(all(is.finite(diag(f$sdr$cov.fixed)[1:2])))
  expect_identical(ck_of(f, "hessian")$status, "FAIL")
  expect_match(ck_of(f, "hessian")$detail, "non-finite or non-positive variances: TRUE",
               fixed = TRUE)
  expect_false(f$ok)
  expect_false(illume:::ilm_fixed_usable(f))
  se <- suppressWarnings(sqrt(diag(vcov(f))))
  expect_true(all(is.nan(se) | is.na(se)))
  expect_false(any(se == 0, na.rm = TRUE))
})

test_that("under REML a fixed-effect covariance that cannot be formed fails too", {
  d <- zero_g_data()
  ## the fixed effects' covariance under REML is the Schur complement of the
  ## joint precision; made to fail, it must fail the fit
  local_mocked_bindings(ilm_schur_vb = function(Hbb, ...) matrix(NaN, nrow(Hbb), nrow(Hbb)))
  f <- q(ilm_model(y ~ x + (1 | g), data = d, family = "gaussian", verbose = FALSE))
  expect_true(isTRUE(f$reml))
  expect_true("g" %in% f$hessian_held)
  expect_identical(ck_of(f, "hessian")$status, "FAIL")
  expect_false(f$ok)
  expect_false(illume:::ilm_fixed_usable(f))
  se <- suppressWarnings(sqrt(diag(vcov(f))))
  expect_false(any(se == 0, na.rm = TRUE))
})

test_that("a held term with healthy fixed effects is still a usable boundary fit", {
  d <- zero_g_data()
  for (reml in c(TRUE, FALSE)) {
    f <- q(ilm_model(y ~ x + (1 | g), data = d, family = "gaussian", reml = reml,
                     verbose = FALSE))
    expect_true("g" %in% f$hessian_held)
    expect_identical(ck_of(f, "hessian")$status, "BOUNDARY")
    expect_true(f$ok)
    expect_true(illume:::ilm_fixed_usable(f))
    se <- sqrt(diag(vcov(f)))
    expect_true(all(is.finite(se) & se > 0))
  }
})
