## Restricted maximum likelihood.
##
## Integrating the fixed effects out under a flat prior IS the restricted
## likelihood, and for a linear-gaussian model the Laplace approximation to
## that integral is exact -- so these are equality checks against lme4, not
## approximation checks.

reml_data <- function(ns, nt, seed = 42) {
  set.seed(seed)
  d <- expand.grid(id = factor(seq_len(ns)), time = factor(seq_len(nt)))
  d$grp <- factor(rep(c("ctl", "trt"), each = ns / 2))[as.integer(d$id)]
  d$x <- rnorm(nrow(d))
  u <- rnorm(ns, 0, 1.2)
  d$y <- 2 + 0.6 * as.integer(d$time) + 0.8 * (d$grp == "trt") + 0.3 * d$x +
    u[as.integer(d$id)] + rnorm(nrow(d), 0, 1)
  d
}

test_that("REML reproduces lme4's restricted fit", {
  skip_if_not_installed("lme4")
  for (cfg in list(list(30, 4, y ~ grp * time + x + (1 | id)),
                   list(12, 3, y ~ grp + time + (1 | id)),
                   list(60, 5, y ~ grp + x + (1 | id)))) {
    d <- reml_data(cfg[[1]], cfg[[2]])
    f <- ilm_model(cfg[[3]], data = d, family = "gaussian", reml = TRUE,
                   verbose = FALSE)
    m <- lme4::lmer(cfg[[3]], data = d, REML = TRUE)
    expect_true(isTRUE(f$reml))
    expect_equal(unname(coef(f)), unname(lme4::fixef(m)), tolerance = 1e-6)
    expect_equal(unname(vcov(f)), unname(as.matrix(vcov(m))), tolerance = 1e-5)
    expect_equal(sqrt(unname(f$Sigma[[1]][1, 1])),
                 as.data.frame(lme4::VarCorr(m))$sdcor[1], tolerance = 1e-5)
    expect_equal(unname(f$dispersion), stats::sigma(m), tolerance = 1e-5)
  }
})

test_that("maximum likelihood stays the default", {
  skip_if_not_installed("lme4")
  d <- reml_data(30, 4)
  f <- ilm_model(y ~ grp + x + (1 | id), data = d, family = "gaussian",
                 verbose = FALSE)
  expect_false(isTRUE(f$reml))
  m <- lme4::lmer(y ~ grp + x + (1 | id), data = d, REML = FALSE)
  expect_equal(sqrt(unname(f$Sigma[[1]][1, 1])),
               as.data.frame(lme4::VarCorr(m))$sdcor[1], tolerance = 1e-5)
})

test_that("Satterthwaite on a REML fit matches lmerTest on a REML fit", {
  skip_if_not_installed("lmerTest")
  d <- reml_data(12, 3)
  f <- ilm_model(y ~ grp + time + (1 | id), data = d, family = "gaussian",
                 reml = TRUE, verbose = FALSE)
  m <- lmerTest::lmer(y ~ grp + time + (1 | id), data = d, REML = TRUE)
  lt <- coef(summary(m))
  for (k in seq_along(coef(f))) {
    l <- numeric(length(coef(f))); l[k] <- 1
    expect_equal(ilm_denom_df(f, l, method = "satterthwaite")$df,
                 unname(lt[k, "df"]), tolerance = 1e-3)
  }
})

test_that("REML on a fixed-effects gaussian model IS the n - p divisor", {
  ## Restricting the likelihood to contrasts orthogonal to X is what produces
  ## n - p, so REML here should land on lm() exactly -- and vcov() must not
  ## then apply its own n/(n - p) correction a second time.
  set.seed(7)
  d <- data.frame(x = rnorm(60), g = factor(rep(c("a", "b", "c"), 20)))
  d$y <- 0.5 * d$x + rnorm(60)
  f <- ilm_model(y ~ x + g, data = d, family = "gaussian", reml = TRUE,
                 verbose = FALSE)
  l <- lm(y ~ x + g, data = d)
  expect_equal(unname(sqrt(diag(vcov(f)))), unname(sqrt(diag(vcov(l)))),
               tolerance = 1e-7)
})

test_that("a likelihood-ratio test on a REML fit is refused, not fudged", {
  d <- reml_data(20, 4)
  f <- ilm_model(y ~ grp + x + (1 | id), data = d, family = "gaussian",
                 reml = TRUE, verbose = FALSE)
  expect_error(ilm_anova(f, test = "LRT"), "REML")
  expect_error(ilm_anova(f, test = "LRT"), "reml = FALSE")
  expect_error(ilm_pb_lrt(f, "grp", B = 2L), "REML")
  ## a Wald test compares nothing across structures and stays available
  expect_s3_class(ilm_anova(f), "data.frame")
})

test_that("REML works for non-gaussian families and matches glmmTMB", {
  ## glmmTMB does `if (REML) randomArg <- c(randomArg, "beta")` for every
  ## family, and so does this. For a linear model that integral is exact and
  ## the result is textbook REML; elsewhere the Laplace approximation to it is
  ## itself approximate, which `reml_exact` records.
  skip_if_not_installed("glmmTMB")
  set.seed(9)
  ncl <- 25; per <- 12; N <- ncl * per
  d <- data.frame(g = factor(rep(seq_len(ncl), each = per)), x = rnorm(N))
  d$grp <- factor(rep(c("a", "b"), each = N / 2))
  u <- rnorm(ncl, 0, 0.8)
  eta <- 0.3 + 0.5 * d$x + 0.4 * (d$grp == "b") + u[as.integer(d$g)]
  d$yp <- rpois(N, exp(pmin(eta, 4)))

  fi <- ilm_model(yp ~ x + grp + (1 | g), data = d, family = "poisson",
                  reml = TRUE, verbose = FALSE)
  fg <- glmmTMB::glmmTMB(yp ~ x + grp + (1 | g), data = d,
                         family = stats::poisson(), REML = TRUE)
  expect_equal(unname(coef(fi)), unname(glmmTMB::fixef(fg)$cond),
               tolerance = 1e-6)
  expect_equal(sqrt(unname(fi$Sigma[[1]][1, 1])),
               sqrt(unname(glmmTMB::VarCorr(fg)$cond[[1]][1, 1])),
               tolerance = 1e-6)
  ## the whole point: REML pulls the variance component up off ML's bias
  fm <- ilm_model(yp ~ x + grp + (1 | g), data = d, family = "poisson",
                  verbose = FALSE)
  expect_gt(fi$Sigma[[1]][1, 1], fm$Sigma[[1]][1, 1])
  ## and it does not claim to be exact where it is not
  expect_false(isTRUE(fi$reml_exact))
})

test_that("reml_exact separates the exact case from the approximate one", {
  d <- reml_data(20, 4)
  fg <- ilm_model(y ~ grp + (1 | id), data = d, family = "gaussian",
                  reml = TRUE, verbose = FALSE)
  expect_true(isTRUE(fg$reml_exact))
  d$yc <- pmax(round(exp(d$y / 3)), 0L)
  fp <- ilm_model(yc ~ grp + (1 | id), data = d, family = "poisson",
                  reml = TRUE, verbose = FALSE)
  expect_true(isTRUE(fp$reml))
  expect_false(isTRUE(fp$reml_exact))
  ## Satterthwaite is still a LINEAR mixed model result and still refuses
  expect_error(ilm_denom_df(fp, c(0, 1), method = "satterthwaite"),
               "LINEAR mixed models")
})

test_that("a sandwich is refused on a REML fit", {
  ## Under REML the coefficients are integrated out, so there is no
  ## per-observation score for a sandwich to sum. glmmTMB's estfun refuses the
  ## same combination.
  set.seed(4)
  d <- data.frame(x = rnorm(200), cl = factor(rep(1:20, 10)))
  d$y <- 0.4 * d$x + rnorm(200)
  f <- ilm_model(y ~ x, data = d, family = "gaussian", reml = TRUE,
                 verbose = FALSE)
  expect_error(ilm_robust(f, cluster = ~ cl), "REML")
  expect_error(ilm_robust(f, cluster = ~ cl), "reml = FALSE")
})

test_that("the variance components survive the layout rewrite", {
  ## The REML fit is reshaped to the ML parameter layout so that downstream
  ## code needs no branch. Reshaping before the covariance derivations ran
  ## recycled a short index over a long vector and silently corrupted every
  ## variance component, which looked like a plausible number.
  skip_if_not_installed("lme4")
  d <- reml_data(30, 4)
  f <- ilm_model(y ~ grp + x + (1 | id), data = d, family = "gaussian",
                 reml = TRUE, verbose = FALSE)
  expect_length(f$dispersion, 1L)
  expect_true(is.finite(f$dispersion))
  expect_lt(sqrt(f$Sigma[[1]][1, 1]), 10)     # not the corrupted 45.76
  expect_identical(names(coef(f)), f$pnames[seq_along(coef(f))])
})
