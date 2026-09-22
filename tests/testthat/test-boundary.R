## A random-effect covariance at the edge of its range -- a variance of zero,
## a correlation of +/-1 -- leaves the likelihood flat in that direction. The
## fixed effects are still identified, and ilm_model() reports them, with the
## verdict saying which parts of the fit can be trusted.
##
## Measured on the messy-data regime with a true random-effect SD of 0.05
## (studies/scripts/messy_compare.R, its own seeds, 400 fits): usable fits
## rose from 63.2% to 97.5%, the newly usable ones covered at 0.944, and their
## standard errors were 1.011 times those of the same model without the random
## term -- which is what a variance of zero says they should be.

bnd_data <- function(seed) {
  set.seed(seed); n <- 300
  d <- data.frame(x = rnorm(n),
                  g = factor(sample(c("north", "central", "south"), n, TRUE)),
                  site = factor(sample(30, n, TRUE)))
  ## no site effect at all
  eta <- cbind(0.4 * d$x, -0.3 * d$x + 0.5 * (d$g == "south"))
  P <- exp(cbind(eta, 0)); P <- P / rowSums(P)
  d$y <- factor(apply(P, 1, function(p) sample(c("lo", "mid", "hi"), 1, prob = p)),
                levels = c("lo", "mid", "hi"))
  d
}

test_that("a variance at zero is held, and the fixed effects are those of the model without it", {
  d <- bnd_data(4)
  f <- ilm_model(y ~ x + g + (1 | site), data = d, family = "multinomial",
                 verbose = FALSE)
  f0 <- ilm_model(y ~ x + g, data = d, family = "multinomial", verbose = FALSE)
  expect_equal(f$hessian_how, "boundary")
  expect_equal(f$hessian_held, "site")
  expect_false(isTRUE(f$sdr$pdHess))            # the record stays honest
  expect_true(f$ok)
  expect_equal(f$checks$status[f$checks$check == "hessian"], "BOUNDARY")
  ## the fixed effects: usable, silent, and those of the model without the term
  expect_silent(V <- vcov(f))
  expect_true(all(is.finite(V)))
  r <- sqrt(diag(V)) / sqrt(diag(vcov(f0)))
  expect_true(all(r > 0.97 & r < 1.05))
  expect_lt(max(abs(coef(f) - coef(f0))), 0.05)
  expect_silent(suppressMessages(ilm_anova(f)))
  ## the covariance parameters are held: asking for them says so
  expect_warning(Vf <- vcov(f, full = TRUE), "held at its estimate")
  expect_true(all(is.finite(Vf)))
  ## and the verdict says, where the numbers are read, what to trust
  out <- paste(utils::capture.output(summary(f)), collapse = "\n")
  expect_match(out, "What can be trusted", fixed = TRUE)
  expect_match(out, "held at its", fixed = TRUE)
  expect_match(paste(utils::capture.output(print(f)), collapse = "\n"),
               "[BOUNDARY: site at a boundary; fixed effects usable]",
               fixed = TRUE)
})

test_that("a covariance at a correlation of -1 is a boundary, not a failure", {
  ## the Hessian is positive definite here, so nothing needs holding; the
  ## rank-deficient covariance used to be graded FAIL all the same, with
  ## "the standard errors above are not usable" beneath standard errors that
  ## were within 4% of the model without the term
  d <- bnd_data(1)
  f <- ilm_model(y ~ x + g + (1 | site), data = d, family = "multinomial",
                 verbose = FALSE)
  expect_equal(f$hessian_how, "tmb")
  expect_true(f$ok)
  expect_equal(f$checks$status[f$checks$check == "sigma_rank[site]"], "BOUNDARY")
  expect_true("site" %in% f$boundary_terms)
  out <- paste(utils::capture.output(summary(f)), collapse = "\n")
  expect_match(out, "What can be trusted", fixed = TRUE)
  expect_match(out, "positive definite", fixed = TRUE)
})

test_that("fixed effects that are not identified are not rescued", {
  ## aliased columns: the singular direction is in the fixed effects
  ## themselves, and no holding of a covariance can make them estimable
  set.seed(2); n <- 200
  d <- data.frame(x1 = rnorm(n), site = factor(sample(10, n, TRUE)))
  d$x2 <- 2 * d$x1
  d$y <- rpois(n, exp(0.3 + 0.2 * d$x1))
  f <- suppressWarnings(ilm_model(y ~ x1 + x2 + (1 | site), data = d,
                                  family = "poisson", verbose = FALSE))
  expect_false(f$hessian_how %in% c("recomputed", "boundary"))
  expect_false(f$ok)
  expect_warning(vcov(f), "not positive definite")
})

test_that("the recomputed Hessian agrees with the exact one where both exist", {
  ## a fixed-effects fit has an exact Hessian through automatic
  ## differentiation, so the finite-difference one can be checked against it
  set.seed(5); n <- 200
  d <- data.frame(x = rnorm(n)); d$y <- rpois(n, exp(0.2 + 0.4 * d$x))
  f <- ilm_model(y ~ x, data = d, family = "poisson", verbose = FALSE)
  H <- ilm_hessian(function(p) as.numeric(f$obj$gr(p)), f$opt$par)
  expect_equal(H, unname(f$obj$he(f$opt$par)), tolerance = 1e-6)
})
