## Methods for other packages' generics have to be REGISTERED with them. An
## exported function named like a method is not found by R's method lookup,
## and for car::Anova(), performance::model_performance() and check_model()
## that is exactly what happened: each fell back to its own default -- an
## error, for check_model() -- while the docs said otherwise.

test_that("car::Anova() reaches illume's joint tests", {
  skip_if_not_installed("car")
  set.seed(1); n <- 300
  d <- data.frame(x = rnorm(n), h = factor(sample(c("u", "v", "w"), n, TRUE)))
  d$k <- factor(sample(c("a", "b", "c"), n, TRUE))
  f <- ilm_model(k ~ x + h, data = d, family = "multinomial", verbose = FALSE)
  a <- car::Anova(f, type = 3)
  ## tested across both category dimensions: car's own fallback gives 1 and 2
  expect_equal(a$Df, c(2, 4))
  expect_equal(a$Df, ilm_anova(f, type = 3)$Df)
})

test_that("performance::model_performance() reaches illume's fit indices", {
  skip_if_not_installed("performance")
  set.seed(2); n <- 200
  d <- data.frame(x = rnorm(n)); d$y <- 0.5 * d$x + rnorm(n)
  f <- ilm_model(y ~ x, data = d, family = "gaussian", verbose = FALSE)
  mp <- performance::model_performance(f)
  expect_false(is.null(mp))
  expect_true(all(c("AIC", "BIC", "R2_McFadden") %in% names(mp)))
})

test_that("performance::check_model() draws ilm_appraise()'s panels", {
  skip_if_not_installed("performance")
  set.seed(3); n <- 200
  d <- data.frame(x = rnorm(n)); d$y <- 0.5 * d$x + rnorm(n)
  f <- ilm_model(y ~ x, data = d, family = "gaussian", verbose = FALSE)
  grDevices::pdf(NULL); on.exit(grDevices::dev.off(), add = TRUE)
  ## the arguments reach ilm_appraise() too
  r <- suppressWarnings(performance::check_model(f, B = 20L))
  expect_true(all(c("rqr", "calibration", "re") %in% names(r)))
})
