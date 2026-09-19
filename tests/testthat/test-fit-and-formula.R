test_that("a basic model fits and passes its own checks", {
  fit <- fit_basic()
  expect_s3_class(fit, "lum_model")
  expect_equal(fit$opt$convergence, 0L)
  expect_true(isTRUE(fit$sdr$pdHess))
  expect_true(fit$ok)
  expect_lt(max(abs(fit$obj$gr(fit$opt$par))), 1e-2)
})

test_that("parameter names line up with the parameter vector", {
  # this guard exists because a mismatch would silently mislabel every
  # coefficient, which is worse than having no names at all
  fit <- fit_basic()
  expect_equal(length(fit$pnames), length(fit$obj$par))
  expect_false(any(duplicated(fit$pnames)))
  expect_true(all(grepl(":", fit$pnames[seq_len(ncol(fit$X) * fit$C)])))
})

test_that("fixed effects are recovered on simulated data", {
  dd <- sim_mlmm(seed = 4, n_subj = 60, per = 25)
  fit <- lum_model(y ~ x1 + grp + (1 | subj), data = dd, family = "multinomial", verbose = FALSE)
  bt <- attr(dd, "beta_true")
  expect_lt(max(abs(fixef(fit) - bt)), 0.45)   # generous: small simulation
})

test_that("fewer than three categories is refused", {
  dd <- sim_mlmm(seed = 2)
  dd$y2 <- factor(ifelse(dd$y == "c1", "a", "b"))
  expect_error(lum_model(y2 ~ x1 + (1 | subj), data = dd, family = "multinomial", verbose = FALSE),
               "at least 3")
})

test_that("a model with no random terms is allowed", {
  # fixed-effects-only models were once refused; they are now the simple end of
  # the same engine.  See test-fixed-only.R for the lm() equivalence.
  dd <- sim_mlmm(seed = 2)
  f <- lum_model(y ~ x1, data = dd, family = "multinomial", verbose = FALSE)
  expect_s3_class(f, "lum_model")
  expect_equal(length(f$re), 0L)
})

test_that("formula interface parses intercepts, slopes and smooths", {
  dd <- sim_mlmm(seed = 5, n_subj = 30, per = 14, with_time = TRUE)
  f1 <- lum_model(y ~ x1 + (1 | subj), data = dd, family = "multinomial", verbose = FALSE)
  expect_equal(names(f1$re), "subj")
  expect_equal(unname(f1$dk), 1L)

  f2 <- lum_model(y ~ x1 + time + (1 + time | subj), data = dd, family = "multinomial", verbose = FALSE)
  expect_equal(unname(f2$dk), 2L)                 # intercept + slope
  expect_true(!is.null(f2$Sigma_d$subj))
  expect_equal(f2$Sigma_d$subj[1, 1], 1)          # scale fixed, see lum_mkLd()

  skip_if_not_installed("mgcv")
  f3 <- lum_model(y ~ x1 + s(xs, k = 6) + (1 | subj), data = dd, family = "multinomial", verbose = FALSE)
  expect_true(any(vapply(f3$re, function(e) e$kind == "basis", TRUE)))
  expect_true(length(f3$smooths) >= 1L)
  expect_true(!is.null(f3$jointPrecision))        # auto-enabled with smooths

  # a namespaced smooth call is a common mistake and must fail clearly
  expect_error(lum_model(y ~ x1 + mgcv::s(xs, k = 6) + (1 | subj), data = dd, family = "multinomial",
                    verbose = FALSE), "unqualified")
})

test_that("the formula fit stores what the ecosystem needs", {
  fit <- fit_basic()
  expect_false(is.null(fit$call))
  expect_false(is.null(fit$terms))
  expect_false(is.null(fit$xlev))
  expect_false(is.null(fit$model))
  expect_equal(length(fit$assign), ncol(fit$X))
  expect_s3_class(terms(fit), "terms")
  expect_s3_class(formula(fit), "formula")
})

test_that("category labels come from the outcome factor", {
  fit <- fit_basic()
  expect_equal(fit$ylevels, c("c1", "c2", "c3"))
  expect_true(all(grepl("^c1:|^c2:", names(coef(fit)))))
})
