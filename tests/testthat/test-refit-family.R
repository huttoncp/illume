# Several functions refit the model internally: the LRT, the parametric
# bootstrap, the simulation-calibrated diagnostics and the null-model
# log-likelihood behind pseudo-R-squared.  ilm_fit() defaults to the
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
    f <- ilm_model(y ~ x1 + grp + (1 | g), data = dd, family = fam,
                   verbose = FALSE)
    a <- suppressWarnings(ilm_anova(f, type = 3, test = "LRT"))
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
  full <- ilm_model(y ~ x1 + grp + (1 | g), data = dd, family = "poisson",
                    verbose = FALSE)
  red  <- ilm_model(y ~ grp + (1 | g), data = dd, family = "poisson",
                    verbose = FALSE)
  manual <- 2 * (as.numeric(logLik(full)) - as.numeric(logLik(red)))
  a <- suppressWarnings(ilm_anova(full, type = 3, test = "LRT"))
  expect_equal(a[["x1", "Chisq"]], manual, tolerance = 1e-3)
})

test_that("an internal refit inherits weights as well as family", {
  dd <- sim_fam("gaussian", seed = 4)
  dd$w <- rep(c(1, 3), length.out = nrow(dd))
  full <- ilm_model(y ~ x1 + grp + (1 | g), data = dd, family = "gaussian",
                    weights = w, verbose = FALSE)
  red  <- ilm_model(y ~ grp + (1 | g), data = dd, family = "gaussian",
                    weights = w, verbose = FALSE)
  manual <- 2 * (as.numeric(logLik(full)) - as.numeric(logLik(red)))
  a <- suppressWarnings(ilm_anova(full, type = 3, test = "LRT"))
  expect_equal(a[["x1", "Chisq"]], manual, tolerance = 1e-3)
})

test_that("simulation-calibrated residual tests run for a univariate family", {
  dd <- sim_fam("binomial", seed = 5)
  f <- ilm_model(y ~ x1 + grp + (1 | g), data = dd, family = "binomial",
                 verbose = FALSE)
  r <- suppressWarnings(ilm_rqr_test(f, B = 5L, seed = 1L, verbose = FALSE))
  expect_true(is.finite(r$obs))
  # the simulated null must not be entirely NA: that is what a wrong-family
  # refit would produce, since every refit would fail and be discarded
  expect_true(sum(is.finite(r$null)) > 0L)
})

test_that("the null-model log-likelihood is finite for a univariate family", {
  dd <- sim_fam("poisson", seed = 6)
  f <- ilm_model(y ~ x1 + grp + (1 | g), data = dd, family = "poisson",
                 verbose = FALSE)
  expect_true(is.finite(ilm_null_ll(f)))
})

test_that("the parametric bootstrap runs for a univariate family", {
  dd <- sim_fam("gaussian", seed = 7)
  f <- ilm_model(y ~ x1 + grp + (1 | g), data = dd, family = "gaussian",
                 verbose = FALSE)
  pb <- suppressWarnings(ilm_pb_lrt(f, term = "x1", B = 10L, seed = 1L,
                                    verbose = FALSE))
  expect_true(is.finite(pb$p_boot))
  expect_equal(pb$n_ok, 10L)
  expect_true(all(is.finite(pb$null)))
})

## ---- a reduced fit has to be the SAME model --------------------------------
##
## Every internal refit rebuilds the model with a different response or a
## reduced design. ilm_refit_drops() hands its workers a STUB rather than the
## whole object, and a field missing from that stub is a field the refit
## silently does without. It carried the family and the weights and nothing
## else, so for a censored model the reduced likelihood was computed as if
## nothing were censored -- two likelihoods on different scales, differenced.
## The expanded coverage study caught it: the null rejection rate was 1.000.

test_that("a likelihood-ratio test holds its size for every structure", {
  # measured at 200 replicates: 0.045, 0.050, 0.035, 0.055 against 0.05
  set.seed(30101)
  n <- 500
  mk <- function(kind) {
    d <- data.frame(x1 = stats::rnorm(n),
                    grp = factor(sample(c("a", "b", "c"), n, TRUE)),
                    s = factor(sample(c("lo", "hi"), n, TRUE)))
    eta <- 0.3 + 0.5 * (d$grp == "b") - 0.2 * (d$grp == "c")
    switch(kind,
      tobit = {
        ys <- eta + stats::rnorm(n); up <- unname(stats::quantile(ys, 0.8))
        d$y <- pmin(ys, up)
        ilm_model(y ~ x1 + grp, data = d, family = "gaussian",
                  censor = ilm_censor(d$y, upper = up), verbose = FALSE)
      },
      aft = {
        tt <- exp(eta + 0.7 * log(stats::rexp(n)))
        q <- unname(stats::quantile(tt, 0.65))
        ct <- pmin(q, stats::rexp(n, 1 / (3 * q))); d$y <- pmin(tt, ct)
        ilm_model(y ~ x1 + grp, data = d, family = "weibull",
                  censor = ilm_surv(d$y, as.integer(tt <= ct)), verbose = FALSE)
      },
      disp = {
        d$y <- eta + stats::rnorm(n, 0, ifelse(d$s == "hi", 2, 0.5))
        ilm_model(y ~ x1 + grp, data = d, family = "gaussian",
                  dispformula = ~ s, verbose = FALSE)
      },
      rp = {
        e <- -log(stats::runif(n)); lp <- eta - 0.3
        t1 <- e / (0.15 * exp(lp))
        tt <- pmax(ifelse(t1 <= 3, t1,
                          3 + (e - 0.45 * exp(lp)) / (0.9 * exp(lp))), 1e-4)
        ct <- stats::rexp(n, 1 / (2 * stats::median(tt))); d$y <- pmin(tt, ct)
        ilm_model(y ~ x1 + grp, data = d, family = "rp", rp_df = 3,
                  censor = ilm_surv(d$y, as.integer(tt <= ct)), verbose = FALSE)
      })
  }
  for (kind in c("tobit", "aft", "disp", "rp")) {
    f <- mk(kind)
    a <- suppressWarnings(ilm_anova(f, test = "LRT"))
    w <- suppressWarnings(ilm_anova(f))
    # x1 has no effect, so the two tests should broadly agree and neither
    # should be extreme. A reduced fit that dropped the structure gave a
    # statistic in the hundreds.
    expect_true(is.finite(a[["x1", "Chisq"]]), info = kind)
    expect_lt(a[["x1", "Chisq"]], 10, label = paste(kind, "LRT chisq"))
    expect_equal(a[["x1", "Chisq"]], w[["x1", "Chisq"]], tolerance = 0.25,
                 label = kind)
  }
})

test_that("the reduced fit keeps the censoring it was fitted with", {
  # the mechanism directly: the same reduced design fitted with and without
  # the censoring gives two different likelihoods, and differencing the wrong
  # one against the full model is what produced the 100% rejection rate
  set.seed(30202)
  n <- 400
  d <- data.frame(x1 = stats::rnorm(n))
  ys <- 0.4 * d$x1 + stats::rnorm(n)
  up <- unname(stats::quantile(ys, 0.75)); d$y <- pmin(ys, up)
  cs <- ilm_censor(d$y, upper = up)
  f <- ilm_model(y ~ x1, data = d, family = "gaussian", censor = cs,
                 verbose = FALSE)
  keep <- is.na(f$assign) | !(f$assign %in% 1L)
  with_c <- illume:::ilm_refit_like(f, X = f$X[, keep, drop = FALSE],
                                    keep = keep)
  without <- illume:::ilm_fit(f$X[, keep, drop = FALSE], f$y, f$J,
                              illume:::ilm_re_list_of(f), f$re_struct, f$ar,
                              family = f$family, verbose = FALSE, restarts = 1L)
  expect_false(isTRUE(all.equal(as.numeric(logLik(with_c)),
                                as.numeric(logLik(without)))))
  # and the one the refit uses is the censored one
  expect_equal(with_c$n_censored, f$n_censored)
})
