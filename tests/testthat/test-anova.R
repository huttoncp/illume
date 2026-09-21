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

  # Type III on a treatment-coded fit now REFITS with sum coding and says so,
  # rather than reporting a test that is not the one the label claims
  expect_message(ilm_anova(ftrt, type = 3), "refitted")
  expect_message(ilm_anova(ftrt, type = 3), "contr.sum", fixed = TRUE)
  # recode = FALSE keeps the old behaviour: test it as coded, and warn
  expect_warning(ilm_anova(ftrt, type = 3, recode = FALSE), "not Type III tests")
  # Type II does not depend on the coding at all, and is the default
  expect_silent(ilm_anova(ftrt, type = 2))
  expect_silent(ilm_anova(ftrt))
  expect_equal(ilm_anova(ftrt)$Chisq, ilm_anova(ftrt, type = 2)$Chisq)

  # the refit lands exactly where fitting with contr.sum by hand lands
  expect_equal(suppressMessages(ilm_anova(ftrt, type = 3))$Chisq,
               suppressMessages(ilm_anova(fsum, type = 3))$Chisq,
               tolerance = 1e-4)
  # and it does not touch the fit it was given
  expect_equal(unname(unlist(ftrt$contrasts["g"])), "contr.treatment")
  expect_true("gb" %in% names(stats::coef(ftrt)))

  # The recode is earned: the main-effect row really does change. Both of
  # these leave x uncentred, so the ONLY difference between them is how g is
  # coded -- comparing against a recoded fit would move the x origin too, and
  # confound the two things this is separating.
  a1 <- suppressWarnings(ilm_anova(ftrt, type = 3, recode = FALSE))
  a2 <- suppressWarnings(ilm_anova(fsum, type = 3, recode = FALSE))
  expect_gt(a2["x", "Chisq"] / a1["x", "Chisq"], 3)
  # the multi-column rows are invariant, because both codings span the same
  # subspace for those terms
  expect_equal(a1["g", "Chisq"], a2["g", "Chisq"], tolerance = 1e-6)
  expect_equal(a1["g:x", "Chisq"], a2["g:x", "Chisq"], tolerance = 1e-6)
  # but the g row is NOT invariant to where x sits, which is the other half
  # of the same point: it is tested at x = 0 either way, and centring moves
  # that to the average
  a3 <- suppressMessages(ilm_anova(fsum, type = 3))
  expect_false(isTRUE(all.equal(a2["g", "Chisq"], a3["g", "Chisq"])))
})

test_that("an uncentred numeric in an interaction is the same problem", {
  # the main effect of x in x:z is the slope where z = 0, and only when z is
  # centred is that the average slope -- no contrast to blame, same issue
  set.seed(13); n <- 600
  d <- data.frame(x = stats::rnorm(n, mean = 5), z = stats::rnorm(n, mean = 3),
                  g = factor(sample(c("a", "b", "c"), n, TRUE)))
  d$y <- 0.5 + 0.8 * d$x + 0.3 * d$z + 0.7 * (d$x - 5) * (d$z - 3) +
    stats::rnorm(n, 0, 0.8)
  f <- ilm_model(y ~ x * z + g, data = d, family = "gaussian", verbose = FALSE)
  expect_message(ilm_anova(f, type = 3), "centred")
  # it matches fitting on centred data directly
  d2 <- d; d2$x <- d$x - mean(d$x); d2$z <- d$z - mean(d$z)
  f2 <- ilm_model(y ~ x * z + g, data = d2, family = "gaussian",
                  verbose = FALSE)
  expect_equal(suppressMessages(ilm_anova(f, type = 3))[["F value"]],
               ilm_anova(f2, type = 3)[["F value"]], tolerance = 1e-3)
  # already centred, so nothing to do and nothing to say
  expect_silent(ilm_anova(f2, type = 3))
  # and the uncorrected test really is a different number
  a_raw <- suppressWarnings(ilm_anova(f, type = 3, recode = FALSE))
  a_fix <- suppressMessages(ilm_anova(f, type = 3))
  expect_gt(abs(a_raw["x", "F value"] - a_fix["x", "F value"]), 100)
  # the interaction row does not depend on any of this
  expect_equal(a_raw["x:z", "F value"], a_fix["x:z", "F value"],
               tolerance = 1e-3)
})

test_that("a model without interactions is never warned about", {
  set.seed(12)
  d <- data.frame(id = factor(rep(1:40, each = 10)), x = stats::rnorm(400),
                  g = factor(sample(c("a", "b"), 400, TRUE)))
  d$y <- d$x + stats::rnorm(400)
  f <- ilm_model(y ~ g + x + (1 | id), d, family = "gaussian", verbose = FALSE)
  expect_silent(ilm_anova(f, type = 3))
  # with nothing in an interaction, the two types are the same test
  expect_equal(ilm_anova(f, type = 2)$Chisq, ilm_anova(f, type = 3)$Chisq,
               tolerance = 1e-8)
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
