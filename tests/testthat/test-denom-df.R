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

test_that("Kenward-Roger gives a smaller df than Satterthwaite", {
  ## KR inflates the covariance for the uncertainty in theta, so it should be
  ## the more conservative of the two. If it ever comes out larger, the
  ## adjustment has the wrong sign.
  d <- rm_data(12, 3)
  f <- ilm_model(y ~ grp + time + (1 | id), data = d, family = "gaussian",
                 verbose = FALSE)
  s <- ilm_denom_df(f, c(0, 1, 0, 0), method = "satterthwaite")$df
  k <- ilm_denom_df(f, c(0, 1, 0, 0), method = "kenward-roger")
  expect_identical(k$method, "kenward-roger")
  expect_lt(k$df, s + 1e-8)
  expect_true(is.matrix(k$V))
})
