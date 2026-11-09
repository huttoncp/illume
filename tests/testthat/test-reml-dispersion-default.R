# REML by default for a gaussian dispersion model with no random terms, and
# F tests on Satterthwaite's df for it under df = "auto". The Welch
# comparisons call stats::t.test() and stats::oneway.test().

q <- function(e) suppressMessages(e)
disp_data <- function(n = c(8, 12, 20), sd = c(1, 2, 4), seed = 11) {
  set.seed(seed)
  g <- factor(rep(letters[seq_along(n)], n))
  data.frame(g = g, x = rnorm(sum(n)), y = rnorm(sum(n), 0, rep(sd, n)))
}

test_that("a dispersion model with no random terms is fitted by REML by default, and says so", {
  d <- disp_data()
  expect_message(f <- ilm_model(y ~ g, data = d, family = "gaussian",
                                dispformula = ~ g, verbose = FALSE),
                 "fitted by REML")
  expect_true(isTRUE(f$reml))
  expect_true(isTRUE(f$reml_by_default))
  ## not written into the call: a refit decides again for its own model
  expect_null(f$call$reml)
  ## the same fit as asking for it
  fr <- ilm_model(y ~ g, data = d, family = "gaussian", dispformula = ~ g,
                  reml = TRUE, verbose = FALSE)
  expect_equal(coef(f, full = TRUE), coef(fr, full = TRUE), tolerance = 1e-8)
  expect_false(isTRUE(fr$reml_by_default))
})

test_that("reml = FALSE is respected, and other models keep maximum likelihood", {
  d <- disp_data()
  expect_no_message(f0 <- ilm_model(y ~ g, data = d, family = "gaussian",
                                    dispformula = ~ g, reml = FALSE, verbose = FALSE))
  expect_false(isTRUE(f0$reml))
  ## no dispersion formula
  f1 <- q(ilm_model(y ~ g, data = d, family = "gaussian", verbose = FALSE))
  expect_false(isTRUE(f1$reml))
  expect_false(isTRUE(f1$reml_by_default))
  ## a dispersion in the fitted mean
  d$k <- rpois(nrow(d), 5) + 1
  f2 <- q(ilm_model(k ~ g, data = d, family = "gaussian", dispformula = ~ mu,
                    verbose = FALSE))
  expect_false(isTRUE(f2$reml))
  ## a random term beside the dispersion formula
  d$id <- factor(rep(1:10, length.out = nrow(d)))
  f3 <- q(ilm_model(y ~ g + (1 | id), data = d, family = "gaussian",
                    dispformula = ~ g, verbose = FALSE))
  expect_false(isTRUE(f3$reml))
  ## another family
  d$cnt <- rpois(nrow(d), 3)
  f4 <- q(ilm_model(cnt ~ g, data = d, family = "poisson", verbose = FALSE))
  expect_false(isTRUE(f4$reml))
})

test_that("update() decides again for the model it describes", {
  d <- disp_data()
  d$id <- factor(rep(1:10, length.out = nrow(d)))
  f <- q(ilm_model(y ~ g, data = d, family = "gaussian", dispformula = ~ g,
                   verbose = FALSE))
  expect_true(isTRUE(q(stats::update(f, y ~ g + x))$reml))
  expect_false(isTRUE(q(stats::update(f, dispformula = NULL))$reml))
  expect_false(isTRUE(q(stats::update(f, y ~ g + (1 | id)))$reml))
  ## an explicit choice is kept
  f0 <- q(ilm_model(y ~ g, data = d, family = "gaussian", dispformula = ~ g,
                    reml = FALSE, verbose = FALSE))
  expect_false(isTRUE(q(stats::update(f0, y ~ g + x))$reml))
})

test_that("reml must be TRUE, FALSE or NULL", {
  d <- disp_data()
  expect_error(ilm_model(y ~ g, data = d, family = "gaussian", dispformula = ~ g,
                         reml = NA, verbose = FALSE), "TRUE, FALSE or NULL")
  expect_error(ilm_model(y ~ g, data = d, family = "gaussian", dispformula = ~ g,
                         reml = "yes", verbose = FALSE), "TRUE, FALSE or NULL")
})

test_that("ilm_anova()'s default gives F on Satterthwaite's df for this model", {
  d <- disp_data()
  f <- q(ilm_model(y ~ x + g, data = d, family = "gaussian", dispformula = ~ g,
                   verbose = FALSE))
  a <- as.data.frame(q(ilm_anova(f)))
  s <- as.data.frame(q(ilm_anova(f, df = "satterthwaite")))
  expect_true(all(c("NumDF", "DenDF", "F value", "Pr(>F)") %in% names(a)))
  expect_equal(a, s)
  expect_true(all(is.finite(a$DenDF)))
  ## the chi-square is still there by name
  ch <- as.data.frame(q(ilm_anova(f, statistic = "Chisq")))
  expect_true("Chisq" %in% names(ch))
  expect_equal(ch$Chisq, a$NumDF * a$`F value`, tolerance = 1e-8)
  ## and so is it on a maximum-likelihood fit, the F on its own fit
  f0 <- q(ilm_model(y ~ x + g, data = d, family = "gaussian", dispformula = ~ g,
                    reml = FALSE, verbose = FALSE))
  expect_true("DenDF" %in% names(as.data.frame(q(ilm_anova(f0)))))
})

test_that("the coefficient table gives t tests with df for this model", {
  d <- disp_data()
  f <- q(ilm_model(y ~ g, data = d, family = "gaussian", dispformula = ~ g,
                   verbose = FALSE))
  ct <- ilm_coef_table(f)
  expect_true(all(c("df", "t value", "Pr(>|t|)") %in% names(ct)))
  expect_true(all(is.finite(ct$df)))
})

test_that("with two groups the default test is Welch's t test", {
  d <- disp_data(n = c(9, 14), sd = c(1, 3))
  f <- q(ilm_model(y ~ g, data = d, family = "gaussian", dispformula = ~ g,
                   verbose = FALSE))
  a <- as.data.frame(q(ilm_anova(f)))
  w <- stats::t.test(y ~ g, data = d, var.equal = FALSE)
  expect_equal(a["g", "F value"], unname(w$statistic)^2, tolerance = 1e-5)
  expect_equal(a["g", "DenDF"], unname(w$parameter), tolerance = 1e-4)
  expect_equal(a["g", "Pr(>F)"], w$p.value, tolerance = 1e-4)
})

test_that("with three groups the F test has Welch's numerator", {
  d <- disp_data()
  f <- q(ilm_model(y ~ g, data = d, family = "gaussian", dispformula = ~ g,
                   verbose = FALSE))
  a <- as.data.frame(q(ilm_anova(f)))
  n <- tapply(d$y, d$g, length); s2 <- tapply(d$y, d$g, var); k <- length(n)
  w <- n / s2; m <- sum(w * tapply(d$y, d$g, mean)) / sum(w)
  num <- sum(w * (tapply(d$y, d$g, mean) - m)^2) / (k - 1)
  expect_equal(a["g", "F value"], num, tolerance = 1e-5)
  ## Welch's own statistic divides that by his correction
  L <- sum((1 - w / sum(w))^2 / (n - 1))
  ow <- stats::oneway.test(y ~ g, data = d, var.equal = FALSE)
  expect_equal(unname(ow$statistic), num / (1 + 2 * (k - 2) * L / (k^2 - 1)),
               tolerance = 1e-8)
})

test_that("a likelihood-ratio test on the default fit refits by maximum likelihood", {
  d <- disp_data()
  f <- q(ilm_model(y ~ g, data = d, family = "gaussian", dispformula = ~ g,
                   verbose = FALSE))
  expect_message(a <- ilm_anova(f, test = "LRT"),
                 "refitted by maximum likelihood.*Wald F test")
  f0 <- q(ilm_model(y ~ g, data = d, family = "gaussian", dispformula = ~ g,
                    reml = FALSE, verbose = FALSE))
  expect_no_message(a0 <- ilm_anova(f0, test = "LRT"))
  expect_equal(as.data.frame(a), as.data.frame(a0), tolerance = 1e-5)
})
