## `groups` says which groups' random effects a prediction, a fitted value or
## a residual uses: "fitted", "typical" or "population" ("new" is taken by
## nothing yet). It replaced `marginal` and `conditional`, which still work for
## one release with a warning once per session.

grp_fit <- function() {
  set.seed(21); ng <- 40
  d <- data.frame(x = rnorm(ng * 10), g = factor(rep(seq_len(ng), each = 10)))
  d$y <- rbinom(nrow(d), 1, plogis(0.3 + 0.8 * d$x + rnorm(ng, 0, 1.2)[d$g]))
  ilm_model(y ~ x + (1 | g), data = d, family = "binomial", verbose = FALSE)
}

## the once-per-session record, cleared so each test sees its first warning
grp_reset <- function()
  rm(list = ls(ilm_deprecation_seen), envir = ilm_deprecation_seen)

test_that("the defaults are unchanged", {
  f <- grp_fit()
  b <- as.vector(f$X %*% f$beta)
  expect_equal(unname(predict(f)[, 1]), plogis(b), tolerance = 1e-12)
  expect_identical(predict(f), predict(f, groups = "typical"))
  expect_identical(ilm_fitted(f), ilm_fitted(f, groups = "fitted"))
  expect_identical(ilm_rqr(f), ilm_rqr(f, groups = "fitted"))
  expect_identical(ilm_scores(f), ilm_scores(f, groups = "fitted"))
  expect_identical(ilm_ame(f, "x"), ilm_ame(f, "x", groups = "typical"))
  ## a fitted group's own effects are not zero, so the two differ
  expect_gt(max(abs(ilm_fitted(f) - ilm_fitted(f, groups = "typical"))), 0.05)
  expect_equal(unname(ilm_fitted(f, groups = "typical")[, 1]), plogis(b),
               tolerance = 1e-12)
})

test_that("each old argument gives what its word gives, and warns once", {
  f <- grp_fit()
  grp_reset()
  expect_warning(p1 <- predict(f, marginal = TRUE),
                 "`marginal` in predict() is deprecated", fixed = TRUE)
  expect_identical(p1, predict(f, groups = "population"))
  ## once per session, whichever value
  expect_no_warning(p0 <- predict(f, marginal = FALSE))
  expect_identical(p0, predict(f, groups = "typical"))

  expect_warning(a1 <- ilm_ame(f, "x", marginal = TRUE), "ilm_ame()",
                 fixed = TRUE)
  expect_identical(a1, ilm_ame(f, "x", groups = "population"))

  expect_warning(q0 <- ilm_fitted(f, conditional = FALSE),
                 "groups = \"typical\" for conditional = FALSE", fixed = TRUE)
  expect_identical(q0, ilm_fitted(f, groups = "typical"))

  expect_warning(s0 <- ilm_scores(f, conditional = FALSE), "deprecated")
  expect_identical(s0, ilm_scores(f, groups = "typical"))

  expect_warning(u0 <- ilm_rqr(f, conditional = FALSE, seed = 3L),
                 "deprecated")
  expect_identical(u0, ilm_rqr(f, groups = "typical", seed = 3L))
})

test_that("the old argument given by position is read as before", {
  f <- grp_fit()
  grp_reset()
  ## ilm_fitted(fit, TRUE) and ilm_rqr(fit, FALSE, seed) were common
  expect_warning(q1 <- ilm_fitted(f, TRUE), "deprecated")
  expect_identical(q1, ilm_fitted(f, groups = "fitted"))
  expect_warning(u0 <- ilm_rqr(f, FALSE, 3L), "deprecated")
  expect_identical(u0, ilm_rqr(f, groups = "typical", seed = 3L))
  expect_warning(a1 <- ilm_ame(f, "x", 1e-4, TRUE), "deprecated")
  expect_identical(a1, ilm_ame(f, "x", groups = "population"))
})

test_that("the words are checked, and a word not taken says where it is", {
  f <- grp_fit()
  expect_error(ilm_fitted(f, groups = "new"), "unknown effects")
  expect_error(predict(f, groups = "new"), "a spread rather than one value")
  expect_error(ilm_ame(f, groups = "fitted"),
               "ilm_ame() takes groups = \"typical\" or \"population\"",
               fixed = TRUE)
  expect_error(ilm_fitted(f, groups = "population"),
               "predict(groups = \"population\")", fixed = TRUE)
  expect_error(ilm_rqr(f, groups = "new"), "unknown effects")
  expect_error(predict(f, groups = "everyone"),
               "must be \"fitted\", \"new\", \"typical\" or \"population\"")
  expect_error(predict(f, groups = c("typical", "fitted")), "one word")
  expect_error(predict(f, groups = "population", marginal = TRUE),
               "not both")
  expect_error(ilm_fitted(f, conditional = NA), "TRUE or FALSE")
  ## a unique abbreviation is the word, and the whole default is the default
  expect_identical(predict(f, groups = "pop"),
                   predict(f, groups = "population"))
  expect_identical(predict(f, groups = c("typical", "population", "fitted")),
                   predict(f))
})

test_that("the marginaleffects bridge averages over the groups by default", {
  f <- grp_fit()
  grp_reset()
  pop <- as.vector(predict(f, groups = "population"))
  typ <- as.vector(predict(f, groups = "typical"))
  expect_equal(get_predict.ilm_model(f)$estimate, pop, tolerance = 1e-12)
  expect_equal(get_predict.ilm_model(f, groups = "typical")$estimate, typ,
               tolerance = 1e-12)
  op <- options(ilm_model.groups = "typical")
  expect_equal(get_predict.ilm_model(f)$estimate, typ, tolerance = 1e-12)
  options(op)
  ## the old option still chooses, with its own warning
  op <- options(ilm_model.marginal = FALSE)
  expect_warning(g0 <- get_predict.ilm_model(f), "ilm_model.marginal")
  options(op)
  expect_equal(g0$estimate, typ, tolerance = 1e-12)
  expect_error(get_predict.ilm_model(f, groups = "fitted"), "get_predict()",
               fixed = TRUE)

  skip_if_not_installed("marginaleffects")
  ilm_register_marginaleffects()
  ## `groups` reaches the bridge through marginaleffects' `...`
  mt <- suppressWarnings(marginaleffects::avg_predictions(f, groups = "typical"))
  mp <- marginaleffects::avg_predictions(f)
  expect_equal(mt$estimate, mean(typ), tolerance = 1e-10)
  expect_equal(mp$estimate, mean(pop), tolerance = 1e-10)
})
