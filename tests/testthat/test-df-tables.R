## Finite degrees of freedom in a gaussian mixed model's tables: summary(),
## ilm_anova(), ilm_emmeans() and ilm_contrast() test against t and F on
## Satterthwaite's df by default, as lmerTest does (Craig's rulings D-DF1 to
## D-DF3). The reference is lmerTest fitted the same way (ML, REML = FALSE),
## so what is compared is the approximation, not two estimators.

q <- function(e) suppressMessages(suppressWarnings(e))
tab_data <- function(ns = 12, nt = 4, seed = 42) {
  set.seed(seed)
  d <- expand.grid(id = factor(seq_len(ns)), time = factor(seq_len(nt)))
  d$grp <- factor(rep(c("ctl", "trt"), each = ns / 2))[as.integer(d$id)]
  d$x <- rnorm(nrow(d))
  d$y <- 2 + 0.3 * as.integer(d$time) + 0.8 * (d$grp == "trt") + 0.3 * d$x +
    rnorm(ns, 0, 1.2)[as.integer(d$id)] + rnorm(nrow(d))
  d
}

test_that("summary()'s coefficient table is lmerTest's: df, t and p", {
  skip_if_not_installed("lmerTest")
  d <- tab_data()
  f <- q(ilm_model(reml = FALSE, y ~ grp + time + x + (1 | id), data = d, family = "gaussian",
                   verbose = FALSE))
  m <- lmerTest::lmer(y ~ grp + time + x + (1 | id), data = d, REML = FALSE)
  ct <- ilm_coef_table(f)
  lt <- coef(summary(m))
  expect_identical(names(ct), c("Estimate", "Std. Error", "df", "t value", "Pr(>|t|)"))
  expect_identical(attr(ct, "df_method"), "satterthwaite")
  expect_equal(ct$df, unname(lt[, "df"]), tolerance = 1e-3)
  expect_equal(ct[["Pr(>|t|)"]], unname(lt[, "Pr(>|t|)"]), tolerance = 1e-3)
  ## the between-cluster effect has few df, and its p is larger than the z's
  z <- ilm_coef_table(f, df = "asymptotic")
  expect_gt(ct[["Pr(>|t|)"]][2], z[["Pr(>|z|)"]][2])
  out <- capture.output(summary(f))
  expect_true(any(grepl("t tests on Satterthwaite's degrees of freedom", out, fixed = TRUE)))
})

test_that("ilm_anova() gives F on Satterthwaite's denominator df, as lmerTest", {
  skip_if_not_installed("lmerTest")
  d <- tab_data()
  f <- q(ilm_model(reml = FALSE, y ~ grp + time + x + (1 | id), data = d, family = "gaussian",
                   verbose = FALSE))
  m <- lmerTest::lmer(y ~ grp + time + x + (1 | id), data = d, REML = FALSE)
  a <- ilm_anova(f)
  la <- anova(m, type = 2)
  expect_identical(names(a), c("NumDF", "DenDF", "F value", "Pr(>F)"))
  expect_equal(a$NumDF, unname(la$NumDF))
  expect_equal(a$DenDF, unname(la$DenDF), tolerance = 1e-3)
  ## the covariance is the joint Hessian's, lmerTest's is the fixed effects'
  ## given theta: the F agree to about 1e-3
  expect_equal(a[["F value"]], unname(la[["F value"]]), tolerance = 1e-3)
  expect_equal(a[["Pr(>F)"]], unname(la[["Pr(>F)"]]), tolerance = 1e-3)
  expect_match(attr(a, "heading")[2], "Satterthwaite", fixed = TRUE)
  ## the chi-square table, as before, on request
  ch <- ilm_anova(f, statistic = "Chisq")
  expect_identical(names(ch), c("Df", "Chisq", "Pr(>Chisq)"))
})

test_that("a Type II test with a higher-order relative takes its df from the reduced fit", {
  d <- tab_data()
  f <- q(ilm_model(reml = FALSE, y ~ grp * time + (1 | id), data = d, family = "gaussian",
                   verbose = FALSE))
  a <- q(ilm_anova(f))
  expect_true(all(is.finite(a$DenDF)))
  expect_true(all(a$DenDF > 0))
  expect_equal(a$NumDF, c(1, 3, 3))
})

test_that("marginal means and their contrasts carry Satterthwaite's df, as emmeans", {
  skip_if_not_installed("lmerTest")
  skip_if_not_installed("emmeans")
  d <- tab_data()
  f <- q(ilm_model(reml = FALSE, y ~ grp + time + (1 | id), data = d, family = "gaussian",
                   verbose = FALSE))
  m <- lmerTest::lmer(y ~ grp + time + (1 | id), data = d, REML = FALSE)
  em <- ilm_emmeans(f, "grp")
  ref <- as.data.frame(emmeans::emmeans(m, "grp", lmer.df = "satterthwaite"))
  expect_equal(em$df, ref$df, tolerance = 1e-3)
  expect_equal(em$lower, ref$lower.CL, tolerance = 1e-3)
  expect_output(print(em), "Intervals are t on Satterthwaite's degrees of freedom",
                fixed = TRUE)
  ct <- ilm_contrast(em)
  refc <- as.data.frame(emmeans::contrast(emmeans::emmeans(m, "grp", lmer.df = "satterthwaite"),
                                          "revpairwise"))
  expect_identical(attr(ct, "df_method"), "satterthwaite")
  expect_equal(ct$df, refc$df, tolerance = 1e-3)
})

test_that("Kenward-Roger on request, and other families untouched", {
  d <- tab_data()
  fr <- q(ilm_model(y ~ grp + x + (1 | id), data = d, family = "gaussian",
                    reml = TRUE, verbose = FALSE))
  kr <- ilm_coef_table(fr, df = "kenward-roger")
  expect_identical(attr(kr, "df_method"), "kenward-roger")
  expect_true(all(is.finite(kr$df)))
  expect_false(isTRUE(all.equal(kr$df, ilm_coef_table(fr)$df)))
  expect_identical(names(ilm_anova(fr, df = "kenward-roger")),
                   c("NumDF", "DenDF", "F value", "Pr(>F)"))
  d$b <- rbinom(nrow(d), 1, 0.4)
  fb <- q(ilm_model(reml = FALSE, b ~ grp + (1 | id), data = d, family = "binomial", verbose = FALSE))
  expect_identical(names(ilm_coef_table(fb)), c("Estimate", "Std. Error", "z value", "Pr(>|z|)"))
  expect_identical(names(ilm_anova(fb)), c("Df", "Chisq", "Pr(>Chisq)"))
  ## and a linear model keeps its exact t and F
  fl <- q(ilm_model(reml = FALSE, y ~ grp + x, data = d, family = "gaussian", verbose = FALSE))
  expect_identical(names(ilm_coef_table(fl)), c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
})

test_that("a variance held at its boundary still gives finite df, or says why not", {
  set.seed(5)
  d <- data.frame(g = factor(rep(1:8, each = 5)), trt = rep(0:1, each = 20))
  d$y <- 1 + 0.5 * d$trt + rnorm(40)          # no cluster variance at all
  f <- q(ilm_model(reml = FALSE, y ~ trt + (1 | g), data = d, family = "gaussian", verbose = FALSE))
  ct <- q(ilm_coef_table(f))
  expect_true("df" %in% names(ct))
  expect_true(all(is.finite(ct$df) | !is.null(attr(ct, "fallback"))))
})
