ctr_data <- function(n = 400L, seed = 1L) {
  set.seed(seed)
  d <- data.frame(g = factor(sample(c("a", "b", "c"), n, TRUE)),
                  h = factor(sample(c("p", "q"), n, TRUE)),
                  o = factor(sample(c("lo", "mid", "hi"), n, TRUE),
                             levels = c("lo", "mid", "hi"), ordered = TRUE),
                  x = rnorm(n))
  d$y <- 1 + 0.5 * (d$g == "b") + 0.2 * d$x + rnorm(n)
  d
}

test_that("the default coding is R's own, so a coefficient means what it does in lm", {
  d <- ctr_data()
  f <- ilm_model(y ~ g + h + x, data = d, family = "gaussian", verbose = FALSE)
  l <- stats::lm(y ~ g + h + x, data = d)
  expect_equal(names(coef(f)), names(coef(l)))
  expect_equal(unname(coef(f)), unname(coef(l)), tolerance = 1e-5)
  expect_equal(unname(unlist(f$contrasts["g"])), "contr.treatment")
  ## an ordered predictor keeps R's polynomial coding, again as lm() does
  fo <- ilm_model(y ~ o, data = d, family = "gaussian", verbose = FALSE)
  expect_equal(names(coef(fo)), names(coef(stats::lm(y ~ o, data = d))))
})

test_that("a single string sets every unordered factor at once", {
  d <- ctr_data()
  for (ct in c("sum", "contr.sum")) {
    f <- ilm_model(y ~ g + h + o + x, data = d, family = "gaussian",
                   contrasts = ct, verbose = FALSE)
    expect_true(all(c("g1", "g2", "h1") %in% names(coef(f))))
    ## the ordered factor is NOT swept up: its polynomial coding says
    ## something about spacing that sum coding would silently discard
    expect_true(all(c("o.L", "o.Q") %in% names(coef(f))))
  }
  f <- ilm_model(y ~ g + h + x, data = d, family = "gaussian",
                 contrasts = "treatment", verbose = FALSE)
  expect_true(all(c("gb", "gc", "hq") %in% names(coef(f))))
  ## a named list still works, one factor at a time, exactly as lm() takes it
  f2 <- ilm_model(y ~ g + h + x, data = d, family = "gaussian",
                  contrasts = list(g = "contr.sum"), verbose = FALSE)
  expect_true(all(c("g1", "g2", "hq") %in% names(coef(f2))))
})

test_that("the coding changes the parameterisation and not the fit", {
  d <- ctr_data()
  f1 <- ilm_model(y ~ g * h + x, data = d, family = "gaussian",
                  verbose = FALSE)
  f2 <- ilm_model(y ~ g * h + x, data = d, family = "gaussian",
                  contrasts = "sum", verbose = FALSE)
  expect_equal(as.numeric(logLik(f1)), as.numeric(logLik(f2)),
               tolerance = 1e-6)
  expect_equal(max(abs(predict(f1) - predict(f2))), 0, tolerance = 1e-4)
  ## and marginal means are invariant to it, which is the point of them
  e1 <- ilm_emmeans(f1, "g"); e2 <- ilm_emmeans(f2, "g")
  expect_equal(e1$estimate, e2$estimate, tolerance = 1e-4)
  expect_equal(e1$se, e2$se, tolerance = 1e-4)
})

test_that("a Type III main effect is not invariant to it, and the check says so", {
  d <- ctr_data()
  f1 <- ilm_model(y ~ g * h + x, data = d, family = "gaussian",
                  verbose = FALSE)
  f2 <- ilm_model(y ~ g * h + x, data = d, family = "gaussian",
                  contrasts = "sum", verbose = FALSE)
  ## treatment coding tests g at h = "p" rather than averaged over h, so the
  ## main-effect row differs -- which is why Type III refits, and why Type II
  ## is the default
  a1 <- suppressWarnings(ilm_anova(f1, type = 3, recode = FALSE))
  a2 <- suppressMessages(ilm_anova(f2, type = 3))
  gp <- function(a) a[rownames(a) == "g", ncol(a)]
  expect_false(isTRUE(all.equal(gp(a1), gp(a2))))
  ## the interaction row does not depend on the coding
  ip <- function(a) a[rownames(a) == "g:h", ncol(a)]
  expect_equal(ip(a1), ip(a2), tolerance = 1e-3)
  ## Type II needs no recoding at all, so it agrees with itself -- and it is
  ## the default, so neither of these says anything
  expect_silent(ilm_anova(f1))
  expect_equal(ilm_anova(f1, type = 2)[["Pr(>F)"]],
               ilm_anova(f2, type = 2)[["Pr(>F)"]], tolerance = 1e-3)
  ## asking for Type III on the treatment-coded fit refits it to where the
  ## sum-coded one already is
  expect_message(a3 <- ilm_anova(f1, type = 3), "refitted")
  expect_equal(a3[["Pr(>F)"]], a2[["Pr(>F)"]], tolerance = 1e-4)
  ## and with no interaction there is nothing to warn about
  f3 <- ilm_model(y ~ g + h + x, data = d, family = "gaussian",
                  verbose = FALSE)
  expect_silent(ilm_anova(f3))
})

test_that("a contrast shorthand that is not one says so", {
  d <- ctr_data(200L)
  expect_error(ilm_model(y ~ g, data = d, family = "gaussian",
                         contrasts = "nope", verbose = FALSE),
               "should be one of")
  expect_error(ilm_model(y ~ g, data = d, family = "gaussian",
                         contrasts = 1, verbose = FALSE),
               "must be NULL")
  expect_error(ilm_model(y ~ g, data = d, family = "gaussian",
                         contrasts = c("sum", "treatment"), verbose = FALSE),
               "must be NULL")
  ## a model with no unordered factor has nothing to set, and does not fail
  expect_silent(ilm_model(y ~ x, data = d, family = "gaussian",
                          contrasts = "sum", verbose = FALSE))
})
