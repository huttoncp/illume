test_that("each term is tested jointly across category dimensions", {
  # car::Anova would report df = 1 here by testing only the first category
  fit <- fit_basic()
  a <- lum_anova(fit, type = 3)
  expect_s3_class(a, "anova")
  expect_equal(rownames(a), fit$term_labels)
  ncols <- vapply(seq_along(fit$term_labels),
                  function(j) sum(!is.na(fit$assign) & fit$assign == j), 1L)
  expect_equal(a$Df, ncols * fit$C)
  expect_true(all(a$Chisq >= 0))
  expect_true(all(a[["Pr(>Chisq)"]] >= 0 & a[["Pr(>Chisq)"]] <= 1))
})

test_that("Type II and Type III agree when there are no interactions", {
  fit <- fit_basic()
  a3 <- lum_anova(fit, type = 3)
  a2 <- lum_anova(fit, type = 2)
  expect_equal(a2$Chisq, a3$Chisq)
})

test_that("Type II on an interaction model reproduces Type III on the reduced one", {
  # the strongest available correctness check for Type II
  dd <- sim_mlmm(seed = 8, n_subj = 35, per = 16)
  f_main <- lum_model(y ~ x1 + grp + (1 | subj), data = dd, family = "multinomial", verbose = FALSE)
  f_int  <- lum_model(y ~ x1 * grp + (1 | subj), data = dd, family = "multinomial", verbose = FALSE)
  a_main <- lum_anova(f_main, type = 3)
  a_int2 <- lum_anova(f_int,  type = 2)
  expect_equal(a_int2[c("x1", "grp"), "Chisq"], a_main[c("x1", "grp"), "Chisq"],
               tolerance = 1e-5)
})

test_that("likelihood-ratio tests run and are non-negative", {
  fit <- fit_basic()
  a <- lum_anova(fit, type = 3, test = "LRT")
  expect_true(all(a$Chisq >= 0 | is.na(a$Chisq)))
  expect_equal(a$Df, lum_anova(fit, type = 3)$Df)
})

test_that("a matrix-interface fit cannot be tested", {
  dd <- sim_mlmm(seed = 3)
  X <- model.matrix(~ x1, dd)
  f <- lum_fit(X, as.integer(dd$y), 3L, list(subj = dd$subj), family = "multinomial", verbose = FALSE)
  expect_error(lum_anova(f), "formula interface")
})

test_that("bad type arguments are rejected", {
  fit <- fit_basic()
  expect_error(lum_anova(fit, type = 1), "type must be")
})
