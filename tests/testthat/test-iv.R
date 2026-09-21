## closed-form 2SLS, written out separately so the test has its own reference
ref_2sls <- function(y, X, Z) {
  Pz <- Z %*% solve(crossprod(Z), t(Z))
  A <- t(X) %*% Pz %*% X
  b <- as.numeric(solve(A, t(X) %*% Pz %*% y))
  u <- y - X %*% b
  s2 <- sum(u^2) / (length(y) - ncol(X))
  list(beta = b, se = sqrt(diag(s2 * solve(A))))
}

iv_data <- function(n = 800L, seed = 1L, strength = 0.9) {
  set.seed(seed)
  z1 <- rnorm(n); z2 <- rnorm(n); w <- rnorm(n); u <- rnorm(n)
  x <- strength * z1 + 0.6 * z2 + 0.4 * w + u + rnorm(n)
  y <- 1 + 0.5 * x + 0.3 * w + u + rnorm(n)
  data.frame(y = y, x = x, w = w, z1 = z1, z2 = z2)
}

test_that("2SLS matches the closed form exactly", {
  d <- iv_data()
  f <- ilm_iv(y ~ x + w | z1 + z2 + w, data = d)
  X <- stats::model.matrix(~ x + w, d)
  Z <- stats::model.matrix(~ z1 + z2 + w, d)
  r <- ref_2sls(d$y, X, Z)
  expect_equal(unname(coef(f)), unname(r$beta), tolerance = 1e-10)
  expect_equal(unname(f$se), unname(r$se), tolerance = 1e-10)
  expect_equal(unname(sqrt(diag(vcov(f)))), unname(r$se), tolerance = 1e-10)
  ## and it removes the bias least squares has here
  b_ols <- unname(coef(stats::lm(y ~ x + w, data = d))[2])
  expect_gt(abs(b_ols - 0.5), abs(unname(coef(f))[2] - 0.5))
  ## Consistency is tested at a size where one sample settles it, rather than
  ## by putting a band around a single draw at n = 800: over 200 samples the
  ## bias there is +0.003 against a standard deviation of 0.048, so an
  ## individual estimate two standard errors from the truth is an ordinary
  ## event and not evidence of anything.
  big <- ilm_iv(y ~ x + w | z1 + z2 + w, data = iv_data(30000L, seed = 3L))
  expect_equal(unname(coef(big))[2], 0.5, tolerance = 0.02)
  expect_equal(unname(coef(big))[3], 0.3, tolerance = 0.06)
})

test_that("the standard errors are not the ones lm() twice would give", {
  d <- iv_data()
  f <- ilm_iv(y ~ x + w | z1 + z2 + w, data = d)
  s1 <- stats::lm(x ~ z1 + z2 + w, data = d)
  d2 <- transform(d, xhat = stats::fitted(s1))
  s2 <- stats::lm(y ~ xhat + w, data = d2)
  ## the point estimates ARE the same -- that is why the mistake survives
  expect_equal(unname(coef(f)), unname(coef(s2)), tolerance = 1e-8)
  ## the standard errors are not, because the naive residual is u + v * beta
  ## and u and v are correlated by construction
  expect_false(isTRUE(all.equal(unname(f$se),
                                unname(sqrt(diag(vcov(s2)))))))
  ## the structural residual uses x as observed, not its fitted value
  expect_equal(f$residuals,
               d$y - as.numeric(stats::model.matrix(~ x + w, d) %*% coef(f)),
               tolerance = 1e-10)
})

test_that("the first-stage F and the two thresholds are reported", {
  d <- iv_data()
  f <- ilm_iv(y ~ x + w | z1 + z2 + w, data = d)
  fs <- f$first_stage$x
  expect_gt(fs$F, 100)
  expect_equal(fs$df1, 2L)              # two EXCLUDED instruments
  ## the F is on the excluded instruments only, so it matches an explicit
  ## comparison of the two first-stage models
  m1 <- stats::lm(x ~ z1 + z2 + w, data = d)
  m0 <- stats::lm(x ~ w, data = d)
  a <- stats::anova(m0, m1)
  expect_equal(fs$F, a$F[2], tolerance = 1e-8)

  ## a weak instrument is reported as weak
  set.seed(2); n <- 400
  z <- rnorm(n); u <- rnorm(n)
  x <- 0.06 * z + u + rnorm(n)
  dw <- data.frame(y = 1 + 0.5 * x + u + rnorm(n), x = x, z = z)
  fw <- ilm_iv(y ~ x | z, data = dw)
  expect_lt(fw$first_stage$x$F, 10)
  expect_output(print(fw), "NOT met")
  expect_output(print(fw), "ilm_iv_ar")
})

test_that("Durbin-Wu-Hausman spots an exogenous regressor", {
  ## endogenous: shared error
  d <- iv_data()
  expect_lt(ilm_iv(y ~ x + w | z1 + z2 + w, data = d)$dwh$p, 1e-6)
  ## exogenous: none
  set.seed(5); n <- 800
  z <- rnorm(n); x <- 0.9 * z + rnorm(n)
  d2 <- data.frame(y = 1 + 0.5 * x + rnorm(n), x = x, z = z)
  f2 <- ilm_iv(y ~ x | z, data = d2)
  expect_gt(f2$dwh$p, 0.05)
  expect_output(print(f2), "No evidence of endogeneity")
  ## and instrumenting an exogenous regressor costs precision for nothing
  expect_gt(unname(f2$se)[2],
            unname(sqrt(diag(vcov(stats::lm(y ~ x, data = d2)))))[2])
})

test_that("Sargan catches an instrument that acts on the outcome directly", {
  set.seed(6); n <- 1500
  za <- rnorm(n); zb <- rnorm(n); u <- rnorm(n)
  x <- 0.8 * za + 0.8 * zb + u + rnorm(n)
  ## zb belongs in the structural equation, so it is not a valid instrument
  y <- 1 + 0.5 * x + 0.7 * zb + u + rnorm(n)
  bad <- ilm_iv(y ~ x | za + zb, data = data.frame(y = y, x = x, za = za,
                                                   zb = zb))
  expect_equal(bad$sargan$df, 1L)
  expect_lt(bad$sargan$p, 1e-6)
  ## valid ones pass
  good <- ilm_iv(y ~ x + w | z1 + z2 + w, data = iv_data())
  expect_gt(good$sargan$p, 0.05)
  ## exactly identified: nothing to test, and it says so by returning df 0
  just <- ilm_iv(y ~ x + w | z1 + w, data = iv_data())
  expect_equal(just$sargan$df, 0L)
  expect_true(is.na(just$sargan$p))
})

test_that("the Anderson-Rubin set is unbounded when the instrument is weak", {
  set.seed(2); n <- 400
  z <- rnorm(n); u <- rnorm(n)
  x <- 0.06 * z + u + rnorm(n)
  dw <- data.frame(y = 1 + 0.5 * x + u + rnorm(n), x = x, z = z)
  fw <- ilm_iv(y ~ x | z, data = dw)
  ar <- ilm_iv_ar(fw)
  expect_s3_class(ar, "ilm_iv_ar")
  expect_false(ar$bounded)
  expect_true(is.infinite(ar$lower) || is.infinite(ar$upper))
  expect_output(print(ar), "UNBOUNDED")
  ## a Wald interval in the same situation is finite and tidy, and wrong
  wald <- unname(coef(fw))[2] + c(-1.96, 1.96) * unname(fw$se)[2]
  expect_true(all(is.finite(wald)))

  ## with a strong instrument the set is bounded and close to the Wald one
  d <- iv_data()
  fs <- ilm_iv(y ~ x + w | z1 + z2 + w, data = d)
  as_ <- ilm_iv_ar(fs)
  expect_true(as_$bounded)
  expect_false(as_$empty)
  b <- unname(coef(fs))[2]; s <- unname(fs$se)[2]
  expect_lt(abs(as_$lower - (b - 1.96 * s)), 0.1)
  expect_lt(abs(as_$upper - (b + 1.96 * s)), 0.1)
  ## the true value is inside it
  expect_lt(as_$lower, 0.5); expect_gt(as_$upper, 0.5)
})

test_that("robust and clustered standard errors are available", {
  set.seed(8); G <- 50L; ni <- 12L; n <- G * ni
  g <- factor(rep(seq_len(G), each = ni))
  z <- rnorm(n); u <- rnorm(n) + rep(rnorm(G, 0, 1.2), each = ni)
  x <- 0.9 * z + u + rnorm(n)
  d <- data.frame(y = 1 + 0.5 * x + u + rnorm(n), x = x, z = z, g = g)
  f0 <- ilm_iv(y ~ x | z, data = d)
  fr <- ilm_iv(y ~ x | z, data = d, robust = TRUE)
  fc <- ilm_iv(y ~ x | z, data = d, cluster = ~ g)
  expect_equal(unname(coef(f0)), unname(coef(fc)))   # only the variance moves
  expect_match(fc$vcov_type, "cluster-robust")
  expect_match(fc$vcov_type, "50 clusters")
  expect_match(fr$vcov_type, "heteroskedasticity")
  ## the cluster effect is in the intercept here, so that error grows most
  expect_gt(unname(fc$se)[1] / unname(f0$se)[1], 1.5)
})

test_that("a formula that does not say what it means is refused", {
  d <- iv_data(300L)
  expect_error(ilm_iv(y ~ x + w, data = d), "needs a `\\|`")
  ## w left out of the instrument list means w is endogenous too
  expect_error(ilm_iv(y ~ x + w | z1 + w, data = d), NA)   # this one is fine
  expect_error(ilm_iv(y ~ x + w | z1, data = d),
               "fewer instruments than endogenous regressors")
  ## nothing to instrument
  expect_error(ilm_iv(y ~ x | x, data = d), "nothing is being instrumented")
  ## two copies of the same instrument
  d2 <- transform(d, z3 = d$z1)
  expect_error(ilm_iv(y ~ x | z1 + z3, data = d2), "rank deficient")
  expect_error(ilm_iv_ar(mtcars), "must be an ilm_iv")
})

test_that("missing values drop the same rows from both sides", {
  d <- iv_data(400L)
  d$z1[c(2, 40)] <- NA
  d$w[7] <- NA
  f <- ilm_iv(y ~ x + w | z1 + z2 + w, data = d)
  expect_equal(f$n, 397L)
  expect_equal(nrow(f$X), nrow(f$Z))
  expect_equal(length(f$residuals), 397L)
})
