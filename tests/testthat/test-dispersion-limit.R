# A dispersion at its unbounded limit is held, as a covariance at its
# boundary is.

## Counts around an AR(1) latent at one observation per cell, negative
## binomial with size 5: the latent can take up the overdispersion, and for
## seeds 5 and 6 it takes up all of it, so k runs to infinity. Seed 2 keeps a
## finite k. Reported by another agent, whose draws of log k spanned -743 to
## +788 there.
dl_counts <- function(seed) {
  set.seed(seed); G <- 20; Tn <- 16
  d <- expand.grid(t = seq_len(Tn), g = factor(seq_len(G)))
  d$x <- stats::rnorm(nrow(d))
  lat <- unlist(lapply(seq_len(G), function(i) as.numeric(
    stats::arima.sim(list(ar = 0.7), Tn, sd = 0.8 * sqrt(1 - 0.49)))))
  d$y <- stats::rnbinom(nrow(d), size = 5, mu = exp(0.5 + 0.3 * d$x + lat))
  d
}
dl_fit <- function(d, family = "nbinom")
  suppressWarnings(ilm_model(y ~ x + (1 | g), data = d, family = family,
                             ar = ilm_ar1(~ t | g), verbose = FALSE))

test_that("a negative binomial whose k ran to infinity holds k", {
  d <- dl_counts(5)
  expect_message(f <- dl_fit(d), "k has run to its limit")
  expect_true("dispersion" %in% f$hessian_held)
  expect_identical(f$hessian_how, "boundary")
  ck <- f$checks[f$checks$check == "dispersion_limit", ]
  expect_identical(ck$status, "BOUNDARY")
  expect_match(ck$cause, "no overdispersion beyond what the model's other terms")
  ## the fit is the Poisson model: the same standard errors
  fp <- suppressMessages(dl_fit(d, "poisson"))
  expect_equal(sqrt(diag(vcov(f))), sqrt(diag(vcov(fp))), tolerance = 1e-4)
  ## k has no standard error, and the draws leave it where it is
  pn <- names(f$opt$par)
  expect_true(is.na(diag(f$sdr$cov.fixed)[pn == "logdisp"]))
  dr <- ilm_draws(f, nsim = 50, seed = 1, natural = FALSE)
  ld <- dr$draws[rownames(dr$draws) == "logdisp", ]
  expect_true(all(ld == ld[1]))
  expect_true(all(is.finite(dr$draws)))
  ## said where a reader looks
  expect_warning(vcov(f, full = TRUE), "k has run to its limit")
  r <- ilm_remedies(f)
  i <- which(r$check == "dispersion_limit")
  expect_length(i, 1L)
  expect_match(r$change[i], "poisson")
})

test_that("summary() says what a dispersion at its limit means", {
  f <- suppressMessages(dl_fit(dl_counts(5)))
  out <- paste(capture.output(summary(f)), collapse = " ")
  expect_match(out, "dispersion_limit", fixed = TRUE)
  expect_match(out, "family = \"poisson\" is the simpler", fixed = TRUE)
})

test_that("a negative binomial with a finite k is left as it was", {
  f <- suppressMessages(dl_fit(dl_counts(2)))
  expect_false("dispersion" %in% f$hessian_held)
  expect_false("dispersion_limit" %in% f$checks$check)
  pn <- names(f$opt$par)
  expect_true(is.finite(diag(f$sdr$cov.fixed)[pn == "logdisp"]))
})

## The study's beta generator: an AR(1) latent at one observation per cell,
## beta noise of precision phi around it (studies/scripts/dispersion_limit.R)
dl_beta <- function(seed, phi) {
  set.seed(seed); G <- 20; Tn <- 16
  d <- expand.grid(t = seq_len(Tn), g = factor(seq_len(G)))
  lat <- unlist(lapply(seq_len(G), function(i) as.numeric(
    stats::arima.sim(list(ar = 0.7), Tn, sd = 0.8 * sqrt(1 - 0.49)))))
  d$x <- stats::rnorm(nrow(d))
  mu <- stats::plogis(0.5 + 0.3 * d$x + lat - 0.5)
  d$y <- pmin(pmax(stats::rbeta(nrow(d), mu * phi, (1 - mu) * phi), 1e-6),
              1 - 1e-6)
  d
}

test_that("a beta's phi is held at its limit, and not at a curved optimum", {
  ## flat beyond phi = 6e7: held
  f <- suppressMessages(dl_fit(dl_beta(71272, 200), "beta"))
  expect_true("dispersion" %in% f$hessian_held)
  expect_match(f$checks$cause[f$checks$check == "dispersion_limit"],
               "other terms carry all the variation")
  ## phi near 1e12, but at a curved optimum, the latent interpolating nearly
  ## noiseless data: past the line, not flat, so not held
  f2 <- suppressMessages(dl_fit(dl_beta(63362, 1e4), "beta"))
  expect_lt(exp(-f2$opt$par[names(f2$opt$par) == "logdisp"] / 2), 1e-2)
  expect_false("dispersion" %in% f2$hessian_held)
})

test_that("the hold reads the same whatever the covariate's units, at both edges", {
  ## The rule reads logdisp and the variance parameters, which rescaling never
  ## touches, on the tape the fit ran on. x in large units is rescaled inside
  ## the fit; x / sd(x) by hand is inside the band and left alone. At the
  ## limit (k run off, the push of 3 moving the objective by about 1e-7) both
  ## are held; with a finite k neither is flagged, and both land where the
  ## fit on x as given does. (A beta near the flatness tolerance is not used
  ## here: there the push sits at the objective's own noise, and columns
  ## equal to 4e-16 already stop the optimiser in different places.)
  push <- function(f) {
    p <- f$opt$par; i <- names(p) == "logdisp"; p2 <- p; p2[i] <- p2[i] + 3
    v <- f$obj$fn(p2) - f$obj$fn(p); invisible(f$obj$fn(p)); v
  }
  for (case in list(list(5, TRUE), list(2, FALSE))) {
    d <- dl_counts(case[[1]])
    dh <- d; dh$x <- d$x / stats::sd(d$x)
    db <- d; db$x <- d$x * 1e4
    f0 <- suppressMessages(dl_fit(d))
    fh <- suppressMessages(dl_fit(dh))
    fb <- suppressMessages(dl_fit(db))
    expect_null(fh$rescale)
    expect_true(fb$rescale$x$any)
    for (f in list(f0, fh, fb)) {
      expect_identical("dispersion" %in% f$hessian_held, case[[2]])
      expect_identical("dispersion_limit" %in% f$checks$check, case[[2]])
      expect_equal(f$opt$objective, f0$opt$objective, tolerance = 1e-8)
    }
    if (case[[2]]) expect_true(all(abs(c(push(fh), push(fb))) < 1e-5))
    else expect_equal(push(fb), push(fh), tolerance = 1e-6)
  }
})

## The gaussian study's AR(1) design at one observation per cell: the latent
## process can take up all the noise (studies/scripts/dispersion_limit_gaussian.R)
dl_gauss <- function(seed, noise = 0.5) {
  set.seed(seed); G <- 20; Tn <- 16
  d <- expand.grid(t = seq_len(Tn), g = factor(seq_len(G)))
  lat <- unlist(lapply(seq_len(G), function(i) as.numeric(
    stats::arima.sim(list(ar = 0.7), Tn, sd = 0.8 * sqrt(1 - 0.49)))))
  d$x <- stats::rnorm(nrow(d))
  d$y <- 0.5 + 0.3 * d$x + lat + stats::rnorm(nrow(d), 0, noise)
  d
}

test_that("a gaussian residual SD at zero is held, and a healthy one is not", {
  expect_message(f <- dl_fit(dl_gauss(6027), "gaussian"),
                 "residual SD has run to zero")
  expect_true("dispersion" %in% f$hessian_held)
  ck <- f$checks[f$checks$check == "dispersion_limit", ]
  expect_identical(ck$status, "BOUNDARY")
  expect_match(ck$suggestion, "^coarsen")
  dr <- ilm_draws(f, nsim = 50, seed = 1, natural = FALSE)
  ld <- dr$draws[rownames(dr$draws) == "logdisp", ]
  expect_true(all(ld == ld[1]))
  f2 <- suppressMessages(dl_fit(dl_gauss(6008), "gaussian"))
  expect_false("dispersion" %in% f2$hessian_held)
})
