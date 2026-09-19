# Several functions refit the model internally: the LRT, the parametric
# bootstrap, the simulation-calibrated diagnostics and the null-model
# log-likelihood behind pseudo-R-squared.  lum_fit() defaults to the
# multinomial family, so an internal refit that forgets to pass the fitted
# family is silently a DIFFERENT model -- and because the refit is wrapped in
# try(), the failure surfaces as an empty table rather than an error.
#
# These tests pin the invariant: an internal refit must inherit the family and
# the weights of the model it came from.

sim_fam <- function(fam, seed = 1, n = 400, ncl = 25) {
  set.seed(seed)
  dd <- data.frame(g = factor(rep(seq_len(ncl), each = n / ncl)))
  dd$x1  <- stats::rnorm(n)
  dd$grp <- factor(sample(c("a", "b", "c"), n, TRUE))
  b   <- stats::rnorm(ncl, 0, 0.6)[as.integer(dd$g)]
  eta <- 0.3 + 0.4 * dd$x1 + b
  dd$y <- switch(fam,
    gaussian = eta + stats::rnorm(n),
    poisson  = stats::rpois(n, exp(pmin(eta, 5))),
    binomial = stats::rbinom(n, 1, 1 / (1 + exp(-eta))))
  dd
}

test_that("the LRT returns finite statistics for every univariate family", {
  for (fam in c("gaussian", "poisson", "binomial")) {
    dd <- sim_fam(fam)
    f <- lum_model(y ~ x1 + grp + (1 | g), data = dd, family = fam,
                   verbose = FALSE)
    a <- suppressWarnings(lum_anova(f, type = 3, test = "LRT"))
    pc <- intersect(c("Pr(>Chisq)", "Pr(>F)"), names(a))[1]
    expect_true(all(is.finite(a[["Chisq"]])),
                info = paste("LRT statistic not finite for", fam))
    expect_true(all(is.finite(a[[pc]])),
                info = paste("LRT p-value not finite for", fam))
  }
})

test_that("the LRT statistic matches an explicit pair of refits", {
  # the sharpest check: if the reduced model were fitted under the wrong
  # family, its log-likelihood would not be comparable and this would not match
  dd <- sim_fam("poisson", seed = 3)
  full <- lum_model(y ~ x1 + grp + (1 | g), data = dd, family = "poisson",
                    verbose = FALSE)
  red  <- lum_model(y ~ grp + (1 | g), data = dd, family = "poisson",
                    verbose = FALSE)
  manual <- 2 * (as.numeric(logLik(full)) - as.numeric(logLik(red)))
  a <- suppressWarnings(lum_anova(full, type = 3, test = "LRT"))
  expect_equal(a[["x1", "Chisq"]], manual, tolerance = 1e-3)
})

test_that("an internal refit inherits weights as well as family", {
  dd <- sim_fam("gaussian", seed = 4)
  dd$w <- rep(c(1, 3), length.out = nrow(dd))
  full <- lum_model(y ~ x1 + grp + (1 | g), data = dd, family = "gaussian",
                    weights = w, verbose = FALSE)
  red  <- lum_model(y ~ grp + (1 | g), data = dd, family = "gaussian",
                    weights = w, verbose = FALSE)
  manual <- 2 * (as.numeric(logLik(full)) - as.numeric(logLik(red)))
  a <- suppressWarnings(lum_anova(full, type = 3, test = "LRT"))
  expect_equal(a[["x1", "Chisq"]], manual, tolerance = 1e-3)
})

test_that("simulation-calibrated residual tests run for a univariate family", {
  dd <- sim_fam("binomial", seed = 5)
  f <- lum_model(y ~ x1 + grp + (1 | g), data = dd, family = "binomial",
                 verbose = FALSE)
  r <- suppressWarnings(lum_rqr_test(f, B = 5L, seed = 1L, verbose = FALSE))
  expect_true(is.finite(r$obs))
  # the simulated null must not be entirely NA: that is what a wrong-family
  # refit would produce, since every refit would fail and be discarded
  expect_true(sum(is.finite(r$null)) > 0L)
})

test_that("the null-model log-likelihood is finite for a univariate family", {
  dd <- sim_fam("poisson", seed = 6)
  f <- lum_model(y ~ x1 + grp + (1 | g), data = dd, family = "poisson",
                 verbose = FALSE)
  expect_true(is.finite(lum_null_ll(f)))
})

test_that("the parametric bootstrap runs for a univariate family", {
  dd <- sim_fam("gaussian", seed = 7)
  f <- lum_model(y ~ x1 + grp + (1 | g), data = dd, family = "gaussian",
                 verbose = FALSE)
  pb <- suppressWarnings(lum_pb_lrt(f, term = "x1", B = 10L, seed = 1L,
                                    verbose = FALSE))
  expect_true(is.finite(pb$p_boot))
  expect_equal(pb$n_ok, 10L)
  expect_true(all(is.finite(pb$null)))
})
