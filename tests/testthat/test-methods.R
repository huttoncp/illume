test_that("coef and vcov names agree exactly", {
  # every downstream package keys off this contract
  fit <- fit_basic()
  b <- coef(fit); V <- vcov(fit)
  expect_equal(length(b), ncol(V))
  expect_identical(names(b), colnames(V))
  expect_identical(names(b), rownames(V))
})

test_that("full = TRUE adds the covariance parameters", {
  fit <- fit_basic()
  bf <- coef(fit, full = TRUE); Vf <- vcov(fit, full = TRUE)
  expect_gt(length(bf), length(coef(fit)))
  expect_identical(names(bf), colnames(Vf))
  expect_true(any(grepl("^subj:", names(bf))))
})

test_that("vcov warns when the Hessian is unusable", {
  fit <- fit_basic()
  bad <- fit; bad$sdr$pdHess <- FALSE
  expect_warning(vcov(bad), "not positive definite")
})

test_that("fixef returns a predictors-by-categories matrix", {
  fit <- fit_basic()
  B <- fixef(fit)
  expect_equal(dim(B), c(ncol(fit$X), fit$C))
  expect_equal(colnames(B), fit$ylevels[seq_len(fit$C)])
  expect_equal(as.vector(B), unname(coef(fit)))
  expect_equal(dim(ilm_se_fixef(fit)), dim(B))
})

test_that("logLik carries df and nobs, and AIC/BIC follow from them", {
  fit <- fit_basic()
  ll <- logLik(fit)
  expect_s3_class(ll, "logLik")
  expect_equal(attr(ll, "df"), length(fit$opt$par))  # fixed + covariance
  expect_equal(attr(ll, "nobs"), nobs(fit))
  expect_equal(AIC(fit), -2 * as.numeric(ll) + 2 * attr(ll, "df"))
  expect_gt(BIC(fit), AIC(fit))                      # log(n) > 2 here
})

test_that("BIC's two sample-size conventions differ", {
  fit <- fit_basic()
  expect_false(isTRUE(all.equal(BIC(fit), BIC(fit, n = "groups"))))
  expect_lt(BIC(fit, n = "groups"), BIC(fit))        # fewer groups than rows
})

test_that("ilm_coef_table reports a sensible Wald table", {
  fit <- fit_basic()
  ct <- ilm_coef_table(fit)
  expect_equal(nrow(ct), length(coef(fit)))
  expect_true(all(ct[["Std. Error"]] > 0))
  expect_true(all(ct[["Pr(>|z|)"]] >= 0 & ct[["Pr(>|z|)"]] <= 1))
  expect_equal(ct[["z value"]], ct$Estimate / ct[["Std. Error"]])
})

test_that("print and summary run and surface the checks", {
  fit <- fit_basic()
  expect_output(print(fit), "ilm_model fit")
  out <- capture.output(print(summary(fit)))
  expect_true(any(grepl("Model checks", out)))
  expect_true(any(grepl("SUM-TO-ZERO", out)))   # the contrast warning
})

## ---- the log-likelihood is on the same scale as everyone else's ------------

test_that("logLik agrees with lme4 and nlme, constants included", {
  # This is pinned to external implementations on purpose. The random-effect
  # priors are written by hand as 0.5 u'Sigma^-1 u + 0.5 log|Sigma|, which
  # omits the (1/2) log(2*pi) per latent scalar; TMB's Laplace step then
  # subtracts (q/2) log(2*pi) of its own, so the reported objective was the
  # true negative log-likelihood minus a constant that GROWS with the number
  # of latent values. Every extra latent bought about 0.92 log-likelihood
  # units for free and AIC preferred the bigger random structure.
  skip_if_not_installed("lme4")
  set.seed(4)
  ng <- 50; nt <- 8; n <- ng * nt
  d <- data.frame(id = factor(rep(seq_len(ng), each = nt)),
                  t = rep(seq_len(nt), ng), x = stats::rnorm(n))
  b <- stats::rnorm(ng, 0, 0.5)[as.integer(d$id)]
  d$y <- 0.5 + 0.8 * d$x + b + stats::rnorm(n, 0, 0.9)

  f <- ilm_model(y ~ x + (1 | id), data = d, family = "gaussian",
                 verbose = FALSE)
  l <- lme4::lmer(y ~ x + (1 | id), data = d, REML = FALSE)
  expect_equal(as.numeric(logLik(f)), as.numeric(logLik(l)), tolerance = 1e-4)
  expect_equal(unname(AIC(f)), unname(AIC(l)), tolerance = 1e-3)

  # the correction is exactly (q/2) log(2*pi) for q integrated scalars
  expect_equal(f$n_integrated, ng)
  expect_equal(as.numeric(logLik(f)),
               -f$opt$objective - (ng / 2) * log(2 * pi), tolerance = 1e-10)
})

test_that("the latent AR process matches the marginal form nlme fits", {
  # illume's AR is a latent process PLUS an independent residual, so the
  # covariance it implies marginally is corExp with a nugget -- not corAR1,
  # which has no nugget and is a different model.
  skip_if_not_installed("nlme")
  set.seed(4)
  ng <- 50; nt <- 8; n <- ng * nt
  d <- data.frame(id = factor(rep(seq_len(ng), each = nt)),
                  t = rep(seq_len(nt), ng), x = stats::rnorm(n))
  u <- unlist(lapply(seq_len(ng), function(i) {
    z <- numeric(nt); z[1] <- stats::rnorm(1)
    for (k in 2:nt) z[k] <- 0.7 * z[k - 1] + stats::rnorm(1, 0, sqrt(1 - 0.49))
    z
  }))
  b <- stats::rnorm(ng, 0, 0.5)[as.integer(d$id)]
  d$y <- 0.5 + 0.8 * d$x + b + 1.2 * u + stats::rnorm(n, 0, 0.6)

  f <- ilm_model(y ~ x + (1 | id), data = d, family = "gaussian",
                 ar = suppressWarnings(ilm_ar1(d$t, d$id)), verbose = FALSE)
  g <- nlme::lme(y ~ x, random = ~ 1 | id, data = d, method = "ML",
                 correlation = nlme::corExp(form = ~ t | id, nugget = TRUE),
                 control = nlme::lmeControl(opt = "optim", msMaxIter = 400))
  expect_equal(as.numeric(logLik(f)), as.numeric(logLik(g)), tolerance = 1e-3)
  expect_equal(unname(coef(f)[["x"]]), unname(nlme::fixef(g)[["x"]]),
               tolerance = 1e-3)
  # same number of estimated parameters, so AIC must agree too
  expect_equal(attr(logLik(f), "df"), attr(logLik(g), "df"))
})

test_that("AIC no longer rewards a model for carrying more latent values", {
  # the symptom: adding a correlation structure to data with none used to
  # improve AIC by hundreds, because the extra latents came with a bonus
  skip_if_not_installed("nlme")
  set.seed(9)
  ng <- 40; nt <- 6; n <- ng * nt
  d <- data.frame(id = factor(rep(seq_len(ng), each = nt)),
                  t = rep(seq_len(nt), ng), x = stats::rnorm(n))
  b <- stats::rnorm(ng, 0, 0.5)[as.integer(d$id)]
  d$y <- 0.5 + 0.8 * d$x + b + stats::rnorm(n, 0, 0.9)   # no autocorrelation

  f0 <- ilm_model(y ~ x + (1 | id), data = d, family = "gaussian",
                  verbose = FALSE)
  f1 <- ilm_model(y ~ x + (1 | id), data = d, family = "gaussian",
                  ar = suppressWarnings(ilm_ar1(d$t, d$id)), verbose = FALSE)
  expect_gt(f1$n_integrated, f0$n_integrated)
  # a structure the data do not support must not look hugely better
  expect_lt(as.numeric(logLik(f1)) - as.numeric(logLik(f0)), 8)
})
