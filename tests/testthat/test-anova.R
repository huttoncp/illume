test_that("each term is tested jointly across category dimensions", {
  # car::Anova would report df = 1 here by testing only the first category
  fit <- fit_basic()
  a <- ilm_anova(fit, type = 3)
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
  a3 <- ilm_anova(fit, type = 3)
  a2 <- ilm_anova(fit, type = 2)
  expect_equal(a2$Chisq, a3$Chisq)
})

test_that("Type II on an interaction model reproduces Type III on the reduced one", {
  # the strongest available correctness check for Type II
  dd <- sim_mlmm(seed = 8, n_subj = 35, per = 16)
  f_main <- ilm_model(y ~ x1 + grp + (1 | subj), data = dd, family = "multinomial", verbose = FALSE)
  f_int  <- ilm_model(y ~ x1 * grp + (1 | subj), data = dd, family = "multinomial", verbose = FALSE)
  a_main <- ilm_anova(f_main, type = 3)
  a_int2 <- ilm_anova(f_int,  type = 2)
  expect_equal(a_int2[c("x1", "grp"), "Chisq"], a_main[c("x1", "grp"), "Chisq"],
               tolerance = 1e-5)
})

test_that("likelihood-ratio tests run and are non-negative", {
  fit <- fit_basic()
  a <- ilm_anova(fit, type = 3, test = "LRT")
  expect_true(all(a$Chisq >= 0 | is.na(a$Chisq)))
  expect_equal(a$Df, ilm_anova(fit, type = 3)$Df)
})

test_that("a matrix-interface fit cannot be tested", {
  dd <- sim_mlmm(seed = 3)
  X <- model.matrix(~ x1, dd)
  f <- ilm_fit(X, as.integer(dd$y), 3L, list(subj = dd$subj), family = "multinomial", verbose = FALSE)
  expect_error(ilm_anova(f), "formula interface")
})

test_that("bad type arguments are rejected", {
  fit <- fit_basic()
  expect_error(ilm_anova(fit, type = 1), "type must be")
})

## ---- Type III depends on the coding ----------------------------------------

test_that("Type III with interactions warns when the coding is not orthogonal", {
  set.seed(11)
  n <- 600
  d <- data.frame(id = factor(rep(1:60, each = 10)), x = stats::rnorm(n),
                  z = stats::rnorm(n),
                  g = factor(sample(c("a", "b", "c"), n, TRUE)))
  b <- stats::rnorm(60, 0, 0.5)[as.integer(d$id)]
  d$y <- 0.5 + 0.8 * d$x + 0.3 * d$z +
    c(a = 0, b = 0.9, c = -0.6)[as.character(d$g)] +
    ifelse(d$g == "b", 1, 0) * d$x + b + stats::rnorm(n, 0, 0.8)

  ftrt <- ilm_model(y ~ g * x + z + (1 | id), d, family = "gaussian",
                    verbose = FALSE)
  fsum <- ilm_model(y ~ g * x + z + (1 | id), d, family = "gaussian",
                    verbose = FALSE, contrasts = list(g = "contr.sum"))

  expect_warning(ilm_anova(ftrt, type = 3), "not Type III tests")
  expect_warning(ilm_anova(ftrt, type = 3), "contr.sum", fixed = TRUE)
  # the orthogonal coding is the one the label means, so it says nothing
  expect_silent(ilm_anova(fsum, type = 3))
  # Type II does not depend on the coding at all
  expect_silent(ilm_anova(ftrt, type = 2))

  # and the warning is earned: the main-effect row really does change
  a1 <- suppressWarnings(ilm_anova(ftrt, type = 3))
  a2 <- ilm_anova(fsum, type = 3)
  expect_gt(a2["x", "Chisq"] / a1["x", "Chisq"], 3)
  # the multi-column rows are invariant, because both codings span the same
  # subspace for those terms
  expect_equal(a1["g", "Chisq"], a2["g", "Chisq"], tolerance = 1e-6)
  expect_equal(a1["g:x", "Chisq"], a2["g:x", "Chisq"], tolerance = 1e-6)
})

test_that("a model without interactions is never warned about", {
  set.seed(12)
  d <- data.frame(id = factor(rep(1:40, each = 10)), x = stats::rnorm(400),
                  g = factor(sample(c("a", "b"), 400, TRUE)))
  d$y <- d$x + stats::rnorm(400)
  f <- ilm_model(y ~ g + x + (1 | id), d, family = "gaussian", verbose = FALSE)
  expect_silent(ilm_anova(f, type = 3))
})

test_that("interaction terms are grouped, tested and predicted as one term", {
  set.seed(13)
  n <- 600
  d <- data.frame(id = factor(rep(1:60, each = 10)), x = stats::rnorm(n),
                  z = stats::rnorm(n),
                  g = factor(sample(c("a", "b", "c"), n, TRUE)))
  b <- stats::rnorm(60, 0, 0.5)[as.integer(d$id)]
  d$y <- 0.5 + 0.8 * d$x + 0.3 * d$z + 0.7 * d$x * d$z + b +
    stats::rnorm(n, 0, 0.8)
  f <- ilm_model(y ~ x * z + g + (1 | id), d, family = "gaussian",
                 verbose = FALSE)
  expect_true("x:z" %in% f$term_labels)
  expect_equal(unname(coef(f)[["x:z"]]), 0.7, tolerance = 0.15)
  a <- suppressWarnings(ilm_anova(f, type = 3))
  expect_equal(a["g", "Df"], 2L)      # a 3-level factor is one 2-df row
  expect_equal(a["x:z", "Df"], 1L)
  # the collinearity check reports one row per term, interactions included
  cc <- ilm_check_collinearity(f)
  expect_setequal(cc$term, f$term_labels)
  expect_equal(cc$df[cc$term == "g"], 2L)
})
