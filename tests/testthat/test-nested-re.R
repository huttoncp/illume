## Nested random effects: (1 | school / class)
##
## lme4::findbars() expands a nested bar into a grouping EXPRESSION, not a
## name: `(1 | continent/country)` becomes `(1 | country:continent)` plus
## `(1 | continent)`. illume looked the group up by name, found no column
## called "country:continent", and refused a formula lme4 accepts -- on the
## first thing anyone tries with hierarchical data. Reported from use.
##
## Two separate faults, both needed fixing: the constituent variables never
## reached the model frame (deparse() of the grouping side is not a variable
## name), and the group itself was looked up rather than evaluated.

nested_data <- function(n = 600, seed = 1) {
  set.seed(seed)
  d <- data.frame(sch = factor(sample(paste0("s", 1:12), n, TRUE)))
  d$cls <- factor(paste0(d$sch, "_", sample(1:3, n, TRUE)))
  d$x <- rnorm(n)
  u1 <- rnorm(nlevels(d$sch), 0, 0.8)
  u2 <- rnorm(nlevels(d$cls), 0, 0.5)
  d$y <- 1 + 0.5 * d$x + u1[as.integer(d$sch)] + u2[as.integer(d$cls)] +
    rnorm(n)
  d
}

test_that("a nested bar fits at all", {
  d <- nested_data()
  expect_silent(f <- ilm_model(y ~ x + (1 | sch / cls), data = d,
                               family = "gaussian",
                               verbose = FALSE))
  expect_s3_class(f, "ilm_model")
  ## two terms: the nesting expands to the inner interaction and the outer
  expect_identical(length(f$re), 2L)
  expect_true(any(grepl(":", names(f$re))))
})

test_that("it agrees with lme4 on the same nested model", {
  skip_if_not_installed("lme4")
  d <- nested_data()
  f <- ilm_model(y ~ x + (1 | sch / cls), data = d, verbose = FALSE)
  l <- lme4::lmer(y ~ x + (1 | sch / cls), data = d, REML = FALSE)
  expect_equal(unname(stats::coef(f)), unname(lme4::fixef(l)),
               tolerance = 1e-5)
  vi <- sort(vapply(f$Sigma, function(s) sqrt(s[1, 1]), 0))
  vl <- sort(as.data.frame(lme4::VarCorr(l))$sdcor[1:2])
  expect_equal(unname(vi), unname(vl), tolerance = 1e-4)
})

test_that("the explicit form gives the same fit as the slash", {
  ## (1 | a/b) IS (1 | a) + (1 | b:a); if the two disagree the expansion is
  ## being handled differently from the thing it expands to
  skip_if_not_installed("lme4")
  d <- nested_data()
  a <- ilm_model(y ~ x + (1 | sch / cls), data = d, verbose = FALSE)
  d$cls_in_sch <- factor(paste(d$cls, d$sch, sep = ":"))
  b <- ilm_model(y ~ x + (1 | sch) + (1 | cls_in_sch), data = d,
                 verbose = FALSE)
  expect_equal(unname(stats::coef(a)), unname(stats::coef(b)), tolerance = 1e-5)
})

test_that("unused interaction levels are dropped", {
  ## an interaction of two factors carries every combination as a level,
  ## including those that never occur -- each would otherwise claim a random
  ## effect that no row informs
  d <- nested_data()
  f <- ilm_model(y ~ x + (1 | sch / cls), data = d, verbose = FALSE)
  inner <- f$re[[grep(":", names(f$re))[1]]]
  expect_identical(length(unique(inner$group)), nlevels(droplevels(d$cls)))
  expect_lt(length(unique(inner$group)), nlevels(d$sch) * nlevels(d$cls))
})

test_that("three levels of nesting work", {
  set.seed(3); n <- 900
  d <- data.frame(region = factor(sample(c("n", "s"), n, TRUE)))
  d$sch <- factor(paste0(d$region, "_s", sample(1:6, n, TRUE)))
  d$cls <- factor(paste0(d$sch, "_c", sample(1:3, n, TRUE)))
  d$x <- rnorm(n)
  d$y <- 0.4 * d$x + rnorm(nlevels(d$sch), 0, .6)[as.integer(d$sch)] + rnorm(n)
  expect_silent(f <- ilm_model(y ~ x + (1 | region / sch / cls), data = d,
                               family = "gaussian",
                               verbose = FALSE))
  expect_identical(length(f$re), 3L)
})

test_that("a grouping variable that really is missing still errors clearly", {
  d <- nested_data()
  expect_error(ilm_model(y ~ x + (1 | sch / nope), data = d, verbose = FALSE),
               "could not be built")
  expect_error(ilm_model(y ~ x + (1 | sch / nope), data = d, verbose = FALSE),
               "nope")                       # names the offending variable
})

test_that("the nested fit supports what any other fit supports", {
  ## the bug was in building the term, so everything downstream needs a look
  d <- nested_data()
  f <- ilm_model(y ~ x + (1 | sch / cls), data = d, verbose = FALSE)
  expect_s3_class(ilm_anova(f), "anova")
  expect_true(is.numeric(stats::predict(f)))
  expect_identical(nrow(ilm_coef_table(f)), 2L)
  expect_silent(s <- ilm_simulate(f, nsim = 5, seed = 1))
  expect_identical(nrow(s), nrow(d))
})
