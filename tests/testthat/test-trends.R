## Marginal slopes.
##
## The reference is emmeans::emtrends, which is the implementation people
## actually check against. The two questions this answers -- does the slope
## differ BETWEEN levels, and does it differ from zero WITHIN a level -- have
## different answers and both are tested.

trend_data <- function(seed = 11) {
  set.seed(seed)
  ns <- 24
  d <- expand.grid(id = factor(seq_len(ns)), week = 0:5)
  d$arm <- factor(rep(c("ctl", "drugA", "drugB"),
                      each = ns / 3))[as.integer(d$id)]
  d$age <- rnorm(nrow(d), 50, 8)
  sl <- c(ctl = 0.05, drugA = 0.45, drugB = 0.20)
  u <- rnorm(ns, 0, 1.0)
  d$y <- 3 + sl[as.character(d$arm)] * d$week + 0.01 * d$age +
    u[as.integer(d$id)] + rnorm(nrow(d), 0, 0.8)
  d
}

test_that("slopes and their standard errors match emtrends", {
  skip_if_not_installed("emmeans")
  d <- trend_data()
  f <- ilm_model(y ~ arm * week + age, data = d, family = "gaussian",
                 verbose = FALSE)
  l <- lm(y ~ arm * week + age, data = d)
  tr <- ilm_trends(f, "arm", var = "week")
  et <- as.data.frame(emmeans::emtrends(l, "arm", var = "week"))
  expect_equal(tr$estimate, et$week.trend, tolerance = 1e-5)
  expect_equal(tr$se, et$SE, tolerance = 1e-6)
  expect_equal(unique(tr$df), unique(et$df), tolerance = 1e-8)
})

test_that("differences between slopes match emmeans' pairs()", {
  skip_if_not_installed("emmeans")
  d <- trend_data()
  f <- ilm_model(y ~ arm * week + age, data = d, family = "gaussian",
                 verbose = FALSE)
  l <- lm(y ~ arm * week + age, data = d)
  ct <- ilm_contrast(ilm_trends(f, "arm", var = "week"), adjust = "none")
  ref <- as.data.frame(emmeans::contrast(
    emmeans::emtrends(l, "arm", var = "week"), "pairwise", adjust = "none"))
  ## emmeans orders its pairs the other way round, so compare magnitudes
  expect_equal(sort(abs(ct$estimate)), sort(abs(ref$estimate)),
               tolerance = 1e-5)
  expect_equal(sort(ct$p_adj), sort(ref$p.value), tolerance = 1e-5)
})

test_that("a mixed REML fit gets Satterthwaite df matching emtrends", {
  skip_if_not_installed("lmerTest")
  skip_if_not_installed("emmeans")
  d <- trend_data()
  f <- ilm_model(y ~ arm * week + age + (1 | id), data = d,
                 family = "gaussian", reml = TRUE, verbose = FALSE)
  m <- lmerTest::lmer(y ~ arm * week + age + (1 | id), data = d, REML = TRUE)
  tr <- ilm_trends(f, "arm", var = "week")
  et <- as.data.frame(emmeans::emtrends(m, "arm", var = "week",
                                        mode = "satterthwaite"))
  expect_equal(tr$estimate, et$week.trend, tolerance = 1e-6)
  expect_equal(tr$se, et$SE, tolerance = 1e-6)
  expect_equal(tr$df, et$df, tolerance = 1e-3)
  expect_identical(attr(tr, "df_method"), "satterthwaite")
})

test_that("the within-level test and the between-level test are different", {
  ## The point of having both. On this design the control arm's slope is not
  ## distinguishable from zero and drugA's is, while whether the two DIFFER is
  ## a separate question with its own answer.
  d <- trend_data()
  f <- ilm_model(y ~ arm * week + age, data = d, family = "gaussian",
                 verbose = FALSE)
  tr <- ilm_trends(f, "arm", var = "week")
  expect_true(all(c("statistic", "p.value", "df") %in% names(tr)))
  expect_gt(tr$p.value[tr$arm == "ctl"], 0.05)      # nothing in the control arm
  expect_lt(tr$p.value[tr$arm == "drugA"], 0.05)    # something in drugA
  ct <- ilm_contrast(tr, adjust = "none")
  expect_equal(nrow(ct), 3L)                        # 3 arms -> 3 pairs
})

test_that("a trends object is an emm object, so contrasts just work", {
  d <- trend_data()
  f <- ilm_model(y ~ arm * week, data = d, family = "gaussian",
                 verbose = FALSE)
  tr <- ilm_trends(f, "arm", var = "week")
  expect_s3_class(tr, "ilm_emm")
  expect_false(is.null(attr(tr, "L")))
  expect_equal(dim(attr(tr, "V")), c(3L, 3L))
  expect_s3_class(ilm_contrast(tr, method = "trt.vs.ctrl", ref = "ctl"),
                  "ilm_contrast")
})

test_that("the derivative is exact for a linear term whatever the step", {
  d <- trend_data()
  f <- ilm_model(y ~ arm * week, data = d, family = "gaussian",
                 verbose = FALSE)
  a <- ilm_trends(f, "arm", var = "week", delta = 1e-2)
  b <- ilm_trends(f, "arm", var = "week", delta = 1e-6)
  expect_equal(a$estimate, b$estimate, tolerance = 1e-7)
})

test_that("nonsensical requests are refused with the alternative named", {
  d <- trend_data()
  f <- ilm_model(y ~ arm * week, data = d, family = "gaussian",
                 verbose = FALSE)
  expect_error(ilm_trends(f, "arm", var = "arm"), "only defined for a continuous")
  expect_error(ilm_trends(f, "arm", var = "arm"), "ilm_emmeans")
  expect_error(ilm_trends(f, "week", var = "week"), "cannot be one of them")
  expect_error(ilm_trends(f, "arm", var = "nope"), "not a predictor")
  expect_error(ilm_trends(f, "nope", var = "week"), "not in the model")
})
