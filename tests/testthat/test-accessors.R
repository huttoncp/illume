# The parts of a fitted model: fixef() on nlme's generic, ilm_normal_expect(),
# ilm_ranef() and ilm_varcorr().

acc_data <- function(seed = 3) {
  set.seed(seed)
  d <- data.frame(id = factor(rep(sprintf("s%02d", 1:15), each = 8)),
                  t = rep(0:7, 15))
  d$x <- stats::rnorm(nrow(d))
  d$y <- 1 + 0.5 * d$x + stats::rnorm(15)[d$id] +
    stats::rnorm(15, 0, 0.3)[d$id] * d$t + stats::rnorm(nrow(d))
  d
}

## the fit's full parameter vector, each block centred BY NAME -- what the
## joint precision and a draw from it are laid out along
centre_by_name <- function(f) {
  env <- f$obj$env; rn <- names(env$par)
  pf <- f$opt$par; pr <- f$sdr$par.random
  mu <- rep(NA_real_, length(rn))
  for (nm in unique(rn))
    mu[rn == nm] <- if (nm %in% names(pr)) pr[names(pr) == nm]
                    else pf[names(pf) == nm]
  mu
}

test_that("fixef() reaches the method through nlme's and lme4's generic", {
  ## it used to be registered on a generic of illume's own, which neither
  ## package's fixef() dispatches to
  d <- acc_data()
  f <- ilm_model(y ~ x + (1 | id), data = d, family = "gaussian", verbose = FALSE)
  expect_equal(nlme::fixef(f), fixef.ilm_model(f))
  skip_if_not_installed("lme4")
  expect_equal(lme4::fixef(f), fixef.ilm_model(f))
})

test_that("ilm_normal_expect() is the quadrature predict() averages with", {
  ## against integrate(), through a logit, at spreads where a fixed rule drifts
  for (s in c(0.5, 2, 3.4)) for (e in c(-2, 0, 1.5)) {
    ref <- stats::integrate(function(z) stats::plogis(e + s * z) * stats::dnorm(z),
                            -Inf, Inf, rel.tol = 1e-13)$value
    expect_equal(ilm_normal_expect(stats::plogis, e, s), ref, tolerance = 1e-10)
  }
  ## vectorised and recycled, and a spread of zero is h itself
  v <- ilm_normal_expect(stats::plogis, c(-1, 0, 1), 1)
  expect_length(v, 3L)
  expect_equal(ilm_normal_expect(exp, 0.3, 0), exp(0.3))
  ## a log link's mean is exact: exp(eta + sd^2 / 2)
  expect_equal(ilm_normal_expect(exp, c(0.2, 1), c(0.5, 1.5)),
               exp(c(0.2, 1) + c(0.5, 1.5)^2 / 2), tolerance = 1e-12)
  ## `n` overrides the rule
  expect_false(isTRUE(all.equal(ilm_normal_expect(stats::plogis, 0.4, 3, n = 5),
                                ilm_normal_expect(stats::plogis, 0.4, 3))))
  expect_error(ilm_normal_expect(stats::plogis, 0, -1), "non-negative")
})

test_that("a growing integrand is averaged where its mass is, and never NaN", {
  ## exp(eta + sd z) phi(z) peaks at z = sd, where a rule centred at zero has
  ## almost no nodes: it was 0.2% low at an SD of 10, 19% low at 12, and NaN
  ## from 13.8, which draws of a variance reach with few groups
  s <- c(5, 8, 10, 12, 13.8, 16, 20)
  expect_equal(ilm_normal_expect(exp, 0.5, s), exp(0.5 + s^2 / 2),
               tolerance = 1e-12)
  ## past what double precision holds at the nodes it says so, as Inf
  v <- ilm_normal_expect(exp, c(0, 0), c(3, 30))
  expect_equal(v[1], exp(4.5), tolerance = 1e-12)
  expect_identical(v[2], Inf)
  expect_null(attributes(v))
  ## a peak far from zero at a small spread is found too: the chance of a
  ## count of 3 where the mean is exp(6) lives in the far left tail
  fine <- function(h, e0, sd) {
    z <- seq(-40, 40, by = 1e-4)
    sum(h(e0 + sd * z) * stats::dnorm(z)) * 1e-4
  }
  h3 <- function(e) stats::dpois(3, exp(e))
  expect_equal(ilm_normal_expect(h3, 6, 0.3), fine(h3, 6, 0.3),
               tolerance = 1e-8)
})

test_that("predict()'s population mean under a log link is exact at any spread", {
  set.seed(4); n <- 200
  d <- data.frame(x = stats::rnorm(n), g = factor(rep(1:20, each = 10)))
  d$y <- stats::rpois(n, exp(0.2 + 0.3 * d$x + stats::rnorm(20, 0, 0.5)[d$g]))
  f <- ilm_model(y ~ x + (1 | g), data = d, family = "poisson", verbose = FALSE)
  nd <- data.frame(x = c(-1, 0, 1))
  eta <- as.vector(cbind(1, nd$x) %*% f$beta)
  for (s in c(0.5, 3, 12)) {
    f$Sigma[["g"]][1, 1] <- s^2
    ## the argument's name changes with the `groups =` vocabulary; either
    ## spelling reaches the same average
    pm <- suppressWarnings(predict(f, nd, marginal = TRUE))
    expect_equal(as.vector(pm), exp(eta + s^2 / 2), tolerance = 1e-10)
  }
})

test_that("ilm_ranef() gives lme4's modes and condVar, labelled", {
  skip_if_not_installed("lme4")
  d <- acc_data()
  f <- ilm_model(y ~ x + (1 + t | id), data = d, family = "gaussian",
                 verbose = FALSE)
  m <- lme4::lmer(y ~ x + (1 + t | id), data = d, REML = FALSE)
  r <- ilm_ranef(f)
  expect_s3_class(r, "ilm_ranef")
  expect_equal(nrow(r), 30L)
  expect_identical(levels(r$level), levels(d$id))
  expect_identical(unique(r$factor), "id")
  expect_identical(unique(r$dim), c("(Intercept)", "t"))
  re <- lme4::ranef(m, condVar = TRUE)$id
  pv <- attr(re, "postVar")
  for (j in 1:2) {
    rj <- r[r$dim == unique(r$dim)[j], ]
    expect_identical(as.character(rj$level), rownames(re))
    expect_equal(rj$mode, re[, j], tolerance = 1e-3)
    expect_equal(rj$sd, sqrt(pv[j, j, ]), tolerance = 1e-3)
  }
  ## lme4's long names on request
  expect_named(as.data.frame(r, lme4 = TRUE),
               c("grpvar", "term", "grp", "condval", "condsd"))
  ## and through nlme's generic
  expect_identical(nlme::ranef(f), r)
})

test_that("the conditional SDs are the exact posterior's, by ML and by REML", {
  ## For a gaussian model the Laplace approximation is exact, so the SDs must
  ## be those of the posterior of the random effects given the variance
  ## components -- holding the fixed effects at their estimate under ML, and
  ## integrating them out under a flat prior under REML. Built here densely
  ## from the fit's own estimates. (An earlier draft filled the REML vector
  ## by position and was 60% off.)
  d <- acc_data()
  n <- nrow(d); nl <- nlevels(d$id); g <- as.integer(d$id)
  Zi <- outer(g, seq_len(nl), `==`) * 1
  Z <- cbind(Zi, Zi * d$t)                       # intercepts, then slopes
  X <- cbind(1, d$x)
  for (reml in c(FALSE, TRUE)) {
    f <- ilm_model(y ~ x + (1 + t | id), data = d, family = "gaussian",
                   reml = reml, verbose = FALSE)
    S <- ilm_varcorr(f)$re$id; s2 <- unname(ilm_varcorr(f)$dispersion$value)^2
    Ginv <- solve(kronecker(S[, ], diag(nl)))
    P <- if (reml)
      rbind(cbind(crossprod(X), crossprod(X, Z)),
            cbind(crossprod(Z, X), crossprod(Z) + s2 * Ginv)) / s2
    else (crossprod(Z) + s2 * Ginv) / s2
    V <- solve(P)
    want <- sqrt(diag(V))[if (reml) -(1:2) else seq_len(2 * nl)]
    r <- expect_silent(ilm_ranef(f))
    expect_equal(r$sd, want, tolerance = 1e-6, label = paste("reml =", reml))
  }
})

test_that("each mode is the fit's own value at its row, by ML and by REML", {
  ## Under REML the fixed effects are in the random block and come first, and
  ## reading the random effects by POSITION shifted every mode by p * C. The
  ## rows here are positions in the full parameter vector, so the modes must
  ## be exactly the values there -- with a random slope, where the shift shows.
  d <- acc_data()
  for (reml in c(FALSE, TRUE)) {
    f <- ilm_model(y ~ x + (1 + t | id), data = d, family = "gaussian",
                   reml = reml, verbose = FALSE)
    r <- ilm_ranef(f)
    expect_identical(r$mode, centre_by_name(f)[r$row])
    ## the fixed effects are not random effects, whatever block they sit in
    expect_false(any(names(f$obj$env$par)[r$row] == "beta"))
  }
})

test_that("a correlation over time and a smooth get their own rows", {
  d <- acc_data()
  f <- suppressMessages(ilm_model(y ~ x + (1 | id), data = d,
                                  family = "gaussian",
                                  ar = ilm_rw1(~ t | id), verbose = FALSE))
  r <- ilm_ranef(f)
  cl <- ilm_cells(f)
  w <- r[r$type == "rw1", ]
  expect_equal(nrow(w), nrow(cl))
  expect_identical(as.character(w$level), as.character(cl$group))
  expect_equal(w$time, cl$time)
  expect_identical(w$cell, seq_len(nrow(cl)))
  ## the anchors are held at zero: no row, and nothing to be uncertain about
  expect_true(all(is.na(w$row[cl$anchor])))
  expect_true(all(w$mode[cl$anchor] == 0 & w$sd[cl$anchor] == 0))
  ok <- !is.na(r$row)
  expect_identical(r$mode[ok], centre_by_name(f)[r$row[ok]])
  expect_identical(unique(w$factor), "id")
  ## a smooth's penalised coefficients
  fs <- ilm_model(y ~ s(x) + (1 | id), data = d, family = "gaussian",
                  verbose = FALSE)
  rs <- ilm_ranef(fs)
  expect_true(all(c("smooth", "re") %in% rs$type))
  expect_identical(rs$mode, centre_by_name(fs)[rs$row])
})

test_that("a multinomial term is labelled by category", {
  set.seed(4)
  d <- data.frame(id = factor(rep(1:30, each = 8)), x = stats::rnorm(240))
  d$k <- factor(sample(c("a", "b", "c"), 240, TRUE))
  f <- ilm_model(k ~ x + (1 | id), data = d, family = "multinomial",
                 verbose = FALSE)
  r <- ilm_ranef(f)
  expect_identical(unique(r$dim), c("a:(Intercept)", "b:(Intercept)"))
  expect_equal(nrow(r), 60L)
  expect_identical(r$mode, centre_by_name(f)[r$row])
})

test_that("the labels survive a refit", {
  d <- acc_data()
  f <- ilm_model(y ~ x + (1 | id), data = d, family = "gaussian", verbose = FALSE)
  g <- illume:::ilm_refit_like(f)
  expect_identical(levels(ilm_ranef(g)$level), levels(d$id))
  expect_identical(unique(ilm_ranef(g)$factor), "id")
})

test_that("ilm_varcorr() gives lme4's variance components", {
  skip_if_not_installed("lme4")
  d <- acc_data()
  f <- ilm_model(y ~ x + (1 + t | id), data = d, family = "gaussian",
                 verbose = FALSE)
  m <- lme4::lmer(y ~ x + (1 + t | id), data = d, REML = FALSE)
  v <- ilm_varcorr(f)
  expect_s3_class(v, "ilm_VarCorr")
  expect_equal(c(v$re$id), c(as.matrix(lme4::VarCorr(m)$id)), tolerance = 1e-3)
  expect_equal(unname(attr(v$re$id, "stddev")),
               unname(attr(lme4::VarCorr(m)$id, "stddev")), tolerance = 1e-3)
  ## the residual is a standard deviation, named as one
  expect_named(v$dispersion$value, "sigma")
  expect_equal(unname(v$dispersion$value), stats::sigma(m), tolerance = 1e-3)
  ## lme4's data frame: the same rows, the same columns
  a <- as.data.frame(v)
  b <- as.data.frame(lme4::VarCorr(m))
  expect_named(a, c("grp", "var1", "var2", "vcov", "sdcor"))
  expect_equal(a$vcov, b$vcov, tolerance = 1e-3)
  expect_equal(a$sdcor, b$sdcor, tolerance = 1e-3)
  expect_identical(a$grp, b$grp)
  expect_identical(nlme::VarCorr(f), v)
  expect_output(print(v), "Std.Dev.")
})

test_that("the variance components are the one transform a draw would get", {
  ## at the estimate, the same numbers from the fit and from its parameter
  ## vector run through ilm_rebuild() -- what a draw of them goes through
  d <- acc_data()
  f <- suppressMessages(ilm_model(y ~ x + (1 + t | id), data = d,
                                  family = "gaussian",
                                  ar = ilm_car1(~ t | id, verbose = FALSE),
                                  verbose = FALSE))
  a <- illume:::ilm_natural(f)
  b <- illume:::ilm_natural(f, f$opt$par)
  expect_equal(a, b, tolerance = 1e-12)
  ## and a correlation over time says what its parameters mean, in the words
  ## ilm_cells() uses
  v <- ilm_varcorr(f)
  expect_identical(v$latent$type, "car1")
  expect_equal(v$latent$rho, f$rho)
  expect_identical(v$latent$meaning, attr(ilm_cells(f), "parameterisation"))
  fr <- ilm_model(y ~ x, data = d, family = "gaussian",
                  ar = ilm_rw1(~ t | id), verbose = FALSE)
  vr <- ilm_varcorr(fr)
  expect_equal(c(vr$latent$var_per_time), c(fr$Sigma$ar))
  expect_output(print(vr), "random walk")
})
