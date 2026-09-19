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
  expect_equal(dim(lum_se_fixef(fit)), dim(B))
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

test_that("lum_coef_table reports a sensible Wald table", {
  fit <- fit_basic()
  ct <- lum_coef_table(fit)
  expect_equal(nrow(ct), length(coef(fit)))
  expect_true(all(ct[["Std. Error"]] > 0))
  expect_true(all(ct[["Pr(>|z|)"]] >= 0 & ct[["Pr(>|z|)"]] <= 1))
  expect_equal(ct[["z value"]], ct$Estimate / ct[["Std. Error"]])
})

test_that("print and summary run and surface the checks", {
  fit <- fit_basic()
  expect_output(print(fit), "lum_model fit")
  out <- capture.output(print(summary(fit)))
  expect_true(any(grepl("Model checks", out)))
  expect_true(any(grepl("SUM-TO-ZERO", out)))   # the contrast warning
})
