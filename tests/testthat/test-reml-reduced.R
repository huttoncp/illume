## Under REML the fixed effects are integrated out, so a count or yes/no
## model's only outer parameters are its variances. With every one of them at
## its boundary the fit is the plain model: all are held, the check line says
## so (WARN, never a clean pass), and group-level predictions say they are
## conditional on zero variance. Before this, such fits fell through unheld --
## one with no intervals, one graded usable whose region intervals covered a
## third of the regions, and a binomial (1 | area) at zero graded clean.

q <- function(e) suppressMessages(suppressWarnings(e))
hess_line <- function(f) f$checks[f$checks$check == "hessian", ]

test_that("a binomial REML fit whose only variance is at zero is the plain model", {
  set.seed(10); A <- 30; n <- 600
  d <- data.frame(area = factor(sample(A, n, TRUE)), x = rnorm(n))
  d$y <- rbinom(n, 1, plogis(-0.3 + 0.5 * d$x))      # no area variance at all
  f <- q(ilm_model(y ~ x + (1 | area), data = d, family = "binomial", reml = TRUE, verbose = FALSE))
  expect_lt(exp(f$opt$par[names(f$opt$par) == "theta"]), 1e-3)
  expect_identical(f$hessian_how, "reduced")
  expect_identical(f$hessian_held, "area")
  h <- hess_line(f)
  expect_identical(h$status, "WARN")
  expect_match(h$detail, "reduces to the plain model", fixed = TRUE)
  ## the plain model's coefficients and standard errors
  g <- stats::glm(y ~ x, data = d, family = stats::binomial())
  expect_equal(unname(coef(f)), unname(coef(g)), tolerance = 1e-6)
  expect_equal(unname(sqrt(diag(vcov(f)))), unname(sqrt(diag(vcov(g)))), tolerance = 1e-5)
  ## maximum likelihood on the same data holds the term the usual way
  f0 <- q(ilm_model(y ~ x + (1 | area), data = d, family = "binomial", verbose = FALSE))
  expect_identical(f0$hessian_how, "boundary")
})

test_that("a gaussian REML fit keeps its residual SD, and holds the usual way", {
  set.seed(10); A <- 30; n <- 600
  d <- data.frame(area = factor(sample(A, n, TRUE)), x = rnorm(n))
  d$y <- 1 + 0.5 * d$x + rnorm(n)
  f <- q(ilm_model(y ~ x + (1 | area), data = d, family = "gaussian", reml = TRUE, verbose = FALSE))
  expect_identical(f$hessian_how, "boundary")
  expect_identical(hess_line(f)$status, "BOUNDARY")
})

## Cases chosen clearly away from the edge, so the same branch runs on every
## platform: REML puts every variance at about e^-10 in these, far below the
## line (measured under RTMB 1.9 and 2.0 alike). Two replicates of the BYM
## study used here before held both terms on one machine and neither on
## others, because whether TMB calls the Hessian positive definite at a zero
## variance is float noise (its smallest eigenvalue is +/-1e-11 there); they
## are kept as fragile-band examples in the notes of item 269.
zero_area <- function(seed) {
  set.seed(seed)
  d <- data.frame(area = factor(rep(1:60, each = 10)), x = stats::rnorm(600))
  d$y <- stats::rbinom(600, 1, stats::plogis(-0.3 + 0.4 * d$x))   # no area variance at all
  d
}
zero_crossed <- function(seed) {
  set.seed(seed)
  d <- expand.grid(a = factor(1:30), b = factor(1:20))
  d$x <- stats::rnorm(nrow(d))
  d$y <- stats::rbinom(nrow(d), 1, stats::plogis(-0.3 + 0.4 * d$x))
  d
}

test_that("a REML fit whose every variance is at zero is held, and says so", {
  cases <- list(list(d = zero_area(4), fo = y ~ x + (1 | area), held = "area"),
                list(d = zero_area(23), fo = y ~ x + (1 | area), held = "area"),
                list(d = zero_crossed(20), fo = y ~ x + (1 | a) + (1 | b), held = c("a", "b")),
                list(d = zero_crossed(25), fo = y ~ x + (1 | a) + (1 | b), held = c("a", "b")))
  for (k in cases) {
    f <- q(ilm_model(k$fo, data = k$d, family = "binomial", reml = TRUE, verbose = FALSE))
    expect_identical(f$hessian_how, "reduced")
    expect_setequal(f$hessian_held, k$held)
    expect_identical(hess_line(f)$status, "WARN")
    expect_true(all(is.finite(sqrt(diag(vcov(f))))))
    ## group-level intervals exist, and are not printed silently
    expect_warning(p <- stats::predict(f, newdata = k$d[1:20, ], groups = "fitted",
                                       interval = "confidence", nsim = 100, seed = 1),
                   "conditional on those variances being zero")
    expect_true(is.list(p) && all(is.finite(p$lower)))
    ## maximum likelihood holds the terms the usual way
    f0 <- q(ilm_model(k$fo, data = k$d, family = "binomial", verbose = FALSE))
    expect_identical(f0$hessian_how, "boundary")
  }
})

test_that("both ways into the hold -- TMB's Hessian passed or failed -- reduce the fit", {
  ## Before the reduction, such a fit fell through unheld by two routes: TMB
  ## called its Hessian positive definite, or it did not. Which one a fit takes
  ## at a zero variance is float noise, so no data choose it on every platform;
  ## the hold's own arguments from a real fit are taken, and the function is
  ## called with TMB's verdict set each way.
  got <- NULL
  real <- ilm_hess_recover
  local_mocked_bindings(ilm_hess_recover = function(obj, opt, sdr, cb, joint) {
    if (is.null(got)) got <<- list(obj = obj, opt = opt, sdr = sdr, cb = cb, joint = joint)
    real(obj, opt, sdr, cb, joint)
  })
  f <- q(ilm_model(y ~ x + (1 | area), data = zero_area(4), family = "binomial",
                   reml = TRUE, verbose = FALSE))
  expect_identical(f$hessian_how, "reduced")
  expect_false(is.null(got))
  for (pd in c(TRUE, FALSE)) {
    s <- got$sdr; s$pdHess <- pd
    r <- real(got$obj, got$opt, s, got$cb, got$joint)
    expect_identical(r$how, "reduced", label = paste("pdHess", pd))
    expect_identical(r$held, "area")
    expect_true(all(is.na(r$sdr$cov.fixed)))
  }
})

test_that("no fit with every variance at its boundary grades as a clean pass", {
  ## Poisson counts with no area variance, 8 data sets, REML: every fit whose
  ## area SD is at zero is reduced and carries a WARN; at least one is
  at_zero <- 0L
  for (sd in 1:8) {
    set.seed(sd); A <- 25; n <- 500
    d <- data.frame(area = factor(sample(A, n, TRUE)), x = rnorm(n))
    d$y <- rpois(n, exp(0.5 + 0.3 * d$x))
    f <- q(ilm_model(y ~ x + (1 | area), data = d, family = "poisson", reml = TRUE, verbose = FALSE))
    if (exp(f$opt$par[names(f$opt$par) == "theta"]) < 1e-3) {
      at_zero <- at_zero + 1L
      expect_identical(f$hessian_how, "reduced")
      expect_false(hess_line(f)$status %in% c("OK", "BOUNDARY"))
    }
  }
  expect_gt(at_zero, 0L)
})
