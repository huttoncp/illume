## Denominator degrees of freedom.
##
## The reference is lmerTest with REML = FALSE, because that is the SAME
## estimator illume uses. Comparing against lmerTest's default (REML) would
## measure the difference between two fitting methods rather than the
## correctness of the approximation, and the two are worth keeping apart.

rm_data <- function(ns, nt, seed = 42) {
  set.seed(seed)
  d <- expand.grid(id = factor(seq_len(ns)), time = factor(seq_len(nt)))
  d$grp <- factor(rep(c("ctl", "trt"), each = ns / 2))[as.integer(d$id)]
  d$x <- rnorm(nrow(d))
  u <- rnorm(ns, 0, 1.2)
  d$y <- 2 + 0.6 * as.integer(d$time) + 0.8 * (d$grp == "trt") + 0.3 * d$x +
    u[as.integer(d$id)] + rnorm(nrow(d), 0, 1)
  d
}

test_that("Satterthwaite matches lmerTest fitted the same way", {
  skip_if_not_installed("lmerTest")
  for (cfg in list(list(30, 4, y ~ grp * time + x + (1 | id)),
                   list(12, 3, y ~ grp + time + (1 | id)),
                   list(20, 6, y ~ x + (1 | id)))) {
    d <- rm_data(cfg[[1]], cfg[[2]])
    f <- ilm_model(cfg[[3]], data = d, family = "gaussian", verbose = FALSE)
    m <- lmerTest::lmer(cfg[[3]], data = d, REML = FALSE)
    lt <- coef(summary(m))
    p <- length(f$beta)
    for (k in seq_len(p)) {
      l <- numeric(p); l[k] <- 1
      r <- ilm_denom_df(f, l, method = "satterthwaite")
      expect_equal(r$df, unname(lt[k, "df"]), tolerance = 1e-3)
    }
  }
})

test_that("a multi-row contrast reproduces lmerTest's F denominator", {
  skip_if_not_installed("lmerTest")
  d <- rm_data(30, 4)
  f <- ilm_model(y ~ grp + time + (1 | id), data = d, family = "gaussian",
                 verbose = FALSE)
  m <- lmerTest::lmer(y ~ grp + time + (1 | id), data = d, REML = FALSE)
  a <- anova(m, ddf = "Satterthwaite", type = 3)
  for (j in seq_along(f$term_labels)) {
    idx <- which(f$assign == j)
    L <- matrix(0, length(idx), length(f$beta))
    for (i in seq_along(idx)) L[i, idx[i]] <- 1
    r <- ilm_denom_df(f, L, method = "satterthwaite")
    expect_equal(r$df, a$DenDF[j], tolerance = 1e-3)
  }
})

test_that("with nothing integrated out the df is exact and auto finds it", {
  set.seed(1)
  d <- data.frame(x = rnorm(60), g = factor(rep(c("a", "b", "c"), 20)))
  d$y <- 0.5 * d$x + rnorm(60)
  f <- ilm_model(y ~ x + g, data = d, family = "gaussian", verbose = FALSE)
  r <- ilm_denom_df(f, c(0, 1, 0, 0))
  expect_identical(r$method, "residual")
  expect_equal(r$df, nrow(d) - length(f$beta))
  ## asking for Satterthwaite on a fit with no variance components to estimate
  ## should still land on the exact answer rather than an approximation to it
  expect_equal(ilm_denom_df(f, c(0, 1, 0, 0), method = "satterthwaite")$df,
               nrow(d) - length(f$beta))
})

test_that("asymptotic recovers the chi-square reference", {
  d <- rm_data(20, 4)
  f <- ilm_model(y ~ x + (1 | id), data = d, family = "gaussian",
                 verbose = FALSE)
  expect_identical(ilm_denom_df(f, c(0, 1), method = "asymptotic")$df, Inf)
})

test_that("auto is Satterthwaite for a gaussian mixed fit", {
  d <- rm_data(20, 4)
  f <- ilm_model(y ~ x + (1 | id), data = d, family = "gaussian",
                 verbose = FALSE)
  expect_identical(ilm_denom_df(f, c(0, 1))$method, "satterthwaite")
  expect_true(is.finite(ilm_denom_df(f, c(0, 1))$df))
})

test_that("neither approximation is offered for a non-gaussian family", {
  d <- rm_data(20, 4)
  d$yb <- as.integer(d$y > stats::median(d$y))
  f <- ilm_model(yb ~ x + (1 | id), data = d, family = "binomial",
                 verbose = FALSE)
  expect_error(ilm_denom_df(f, c(0, 1), method = "satterthwaite"),
               "LINEAR mixed models")
  expect_error(ilm_denom_df(f, c(0, 1), method = "satterthwaite"),
               "ilm_pb_lrt")
  ## auto must not quietly hand back a number it cannot stand behind
  expect_identical(ilm_denom_df(f, c(0, 1))$method, "asymptotic")
})

test_that("Kenward-Roger refuses where the covariance is not linear in theta", {
  set.seed(3)
  ncl <- 12; per <- 8
  d <- data.frame(g = factor(rep(seq_len(ncl), each = per)),
                  time = rep(seq_len(per), ncl))
  d$x <- rnorm(nrow(d))
  d$y <- 0.4 * d$x + rnorm(nrow(d))
  f <- ilm_model(y ~ x, data = d, family = "gaussian",
                 ar = ilm_ar1(d$time, d$g, verbose = FALSE), verbose = FALSE)
  expect_error(ilm_denom_df(f, c(0, 1), method = "kenward-roger"),
               "AR\\(1\\)")
  expect_error(ilm_denom_df(f, c(0, 1), method = "kenward-roger"),
               "satterthwaite")
})

test_that("a contrast of the wrong length is rejected", {
  d <- rm_data(20, 4)
  f <- ilm_model(y ~ x + (1 | id), data = d, family = "gaussian",
                 verbose = FALSE)
  expect_error(ilm_denom_df(f, c(0, 1, 1)), "column")
})

## ---- Kenward-Roger ----------------------------------------------------------
##
## The reference is pbkrtest, which computes it from the same covariance of the
## observations. An earlier version assembled it from differences of Vb(theta)
## and got the algebra wrong -- the "inflated" covariance came out smaller, and
## the only test of it asserted that its df were below Satterthwaite's, which
## the defect guaranteed. So these compare numbers.

kr_data <- function() {
  set.seed(10)
  ng <- 8; per <- 6
  d <- data.frame(id = factor(rep(seq_len(ng), each = per)),
                  t = rep(0:(per - 1), ng), x = stats::rnorm(ng * per))
  d$trt <- factor(rep(c("a", "b"), length.out = ng))[as.integer(d$id)]
  b0 <- stats::rnorm(ng, 0, 1); b1 <- stats::rnorm(ng, 0, 0.3)
  d$y <- 1 + 0.5 * d$x + 0.2 * d$t + (d$trt == "b") * 0.4 + b0[d$id] +
    b1[d$id] * d$t + stats::rnorm(nrow(d), 0, 0.7)
  d
}

test_that("Kenward-Roger agrees with pbkrtest", {
  skip_if_not_installed("pbkrtest")
  skip_if_not_installed("lme4")
  d <- kr_data()
  ## eight clusters, where the adjustment is visible: a random intercept, a
  ## correlated random slope -- which the old gate said it refused, and did
  ## not -- and an uncorrelated one
  cases <- list(
    list(f = y ~ x + t + trt + (1 | id), rs = NULL,
         lf = y ~ x + t + trt + (1 | id)),
    list(f = y ~ x + t + trt + (1 + t | id), rs = NULL,
         lf = y ~ x + t + trt + (1 + t | id)),
    list(f = y ~ x + t + trt + (1 + t | id), rs = list(id = list(d_cor = FALSE)),
         lf = y ~ x + t + trt + (1 + t || id)))
  for (cs in cases) {
    fi <- ilm_model(cs$f, data = d, family = "gaussian", reml = TRUE,
                    re_struct = cs$rs, verbose = FALSE)
    fl <- lme4::lmer(cs$lf, data = d, REML = TRUE)
    Va <- pbkrtest::vcovAdj(fl)
    lab <- deparse(cs$lf)
    for (j in 2:4) {
      l <- replace(numeric(4), j, 1)
      k <- ilm_denom_df(fi, l, method = "kenward-roger")
      expect_identical(k$method, "kenward-roger")
      expect_equal(k$df, pbkrtest::Lb_ddf(l, stats::vcov(fl), Va),
                   tolerance = 1e-4, label = lab)
      expect_equal(k$scale, 1, tolerance = 1e-12)
    }
    ## the adjusted covariance; the two fits' REML estimates differ by about
    ## 1e-5 at most, which is what the tolerance allows for
    expect_equal(unname(k$V), unname(as.matrix(Va)), tolerance = 1e-4,
                 label = lab)
    ## a two-row F test: its own df and its scale
    L <- rbind(c(0, 1, 0, 0), c(0, 0, 0, 1))
    k2 <- ilm_denom_df(fi, L, method = "kenward-roger")
    kr <- pbkrtest::KRmodcomp(fl, L)
    expect_equal(k2$df, kr$stats$ddf, tolerance = 1e-4, label = lab)
    expect_equal(k2$scale, kr$stats$F.scaling, tolerance = 1e-4, label = lab)
  }
  ## crossed grouping factors
  set.seed(11)
  dc <- expand.grid(a = factor(1:6), b = factor(1:5), rep = 1:2)
  dc$x <- stats::rnorm(nrow(dc))
  dc$y <- 0.3 * dc$x + stats::rnorm(6)[dc$a] + stats::rnorm(5, 0, 0.7)[dc$b] +
    stats::rnorm(nrow(dc))
  fi <- ilm_model(y ~ x + (1 | a) + (1 | b), data = dc, family = "gaussian",
                  reml = TRUE, verbose = FALSE)
  fl <- lme4::lmer(y ~ x + (1 | a) + (1 | b), data = dc, REML = TRUE)
  Va <- pbkrtest::vcovAdj(fl)
  k <- ilm_denom_df(fi, c(0, 1), method = "kenward-roger")
  expect_equal(k$df, pbkrtest::Lb_ddf(c(0, 1), stats::vcov(fl), Va),
               tolerance = 1e-4)
  expect_equal(unname(k$V), unname(as.matrix(Va)), tolerance = 1e-4)
})

test_that("Kenward-Roger reproduces the exact test of a balanced design", {
  ## a treatment that varies between eight equal clusters: the exact test is a
  ## two-sample t on the cluster means, with 6 df, and Kenward and Roger built
  ## their approximation to reproduce it
  set.seed(12)
  d <- data.frame(id = factor(rep(1:8, each = 5)))
  d$trt <- factor(rep(c("a", "b"), each = 20))
  d$y <- 0.5 * (d$trt == "b") + stats::rnorm(8, 0, 0.8)[d$id] + stats::rnorm(40)
  f <- ilm_model(y ~ trt + (1 | id), data = d, family = "gaussian",
                 reml = TRUE, verbose = FALSE)
  k <- ilm_denom_df(f, c(0, 1), method = "kenward-roger")
  expect_equal(k$df, 6, tolerance = 1e-8)
  m <- tapply(d$y, d$id, mean)
  g <- tapply(as.character(d$trt), d$id, `[`, 1)
  tt <- stats::t.test(m[g == "b"], m[g == "a"], var.equal = TRUE)
  expect_equal(sqrt(k$V[2, 2]), unname(tt$stderr), tolerance = 1e-6)
})

test_that("Kenward-Roger's t in ilm_trends() uses the adjusted covariance", {
  d <- kr_data()
  f <- ilm_model(y ~ x + t + trt + (1 + t | id), data = d, family = "gaussian",
                 reml = TRUE, verbose = FALSE)
  tr <- ilm_trends(f, specs = NULL, var = "x", df = "kenward-roger")
  k <- ilm_denom_df(f, c(0, 1, 0, 0), method = "kenward-roger")
  expect_equal(tr$se, sqrt(k$V[2, 2]), tolerance = 1e-10)
  expect_equal(tr$df, k$df, tolerance = 1e-10)
  expect_identical(attr(tr, "df_method"), "kenward-roger")
})

test_that("Kenward-Roger says why when it is not available", {
  d <- kr_data()
  ## by maximum likelihood: Kenward and Roger derived it for REML
  fm <- ilm_model(y ~ x + t + trt + (1 | id), data = d, family = "gaussian",
                  verbose = FALSE)
  expect_error(ilm_denom_df(fm, c(0, 1, 0, 0), method = "kenward-roger"),
               "reml = TRUE")
  expect_error(ilm_trends(fm, specs = NULL, var = "x", df = "kenward-roger"),
               "reml = TRUE")
  ## a penalised smooth's variance is a smoothing parameter
  fs <- ilm_model(y ~ s(x) + (1 | id), data = d, family = "gaussian",
                  reml = TRUE, verbose = FALSE)
  expect_error(ilm_denom_df(fs, replace(numeric(length(fs$beta)), 2, 1),
                            method = "kenward-roger"), "penalised smooth")
  ## frequency weights
  d$w <- rep(1:2, length.out = nrow(d))
  fw <- ilm_model(y ~ x + (1 | id), data = d, family = "gaussian", weights = w,
                  reml = TRUE, verbose = FALSE)
  expect_error(ilm_denom_df(fw, c(0, 1), method = "kenward-roger"),
               "frequency weights")
})
