# ilm_refit(): the same model, to a new response or to other data.

rf_data <- function(seed = 1) {
  set.seed(seed)
  d <- data.frame(id = factor(rep(1:12, each = 6)), t = rep(1:6, 12),
                  x = stats::rnorm(72))
  d$y <- 1 + 0.5 * d$x + stats::rnorm(12)[d$id] + stats::rnorm(72)
  d
}

test_that("a new response is refitted on the same design, as a full fit", {
  d <- rf_data()
  fit <- ilm_model(y ~ x + (1 | id), data = d, family = "gaussian",
                   verbose = FALSE)
  ystar <- ilm_simulate(fit, nsim = 2, seed = 2)
  r <- ilm_refit(fit, y = ystar[, 1])
  expect_s3_class(r, "ilm_model")
  expect_equal(coef(r), coef(illume:::ilm_refit_like(fit, y = ystar[, 1])))
  ## a one-column matrix is taken as ilm_simulate() gives it
  expect_equal(coef(ilm_refit(fit, y = ystar[, 1, drop = FALSE])), coef(r))
  ## and it is a model like any other: new rows predicted, emmeans built,
  ## the response in its model frame the new one
  nd <- d[1:3, ]
  expect_length(predict(r, newdata = nd), 3L)
  expect_s3_class(ilm_emmeans(r, "x"), "ilm_emm")
  expect_equal(r$model$y, as.numeric(ystar[, 1]))
  ## refitting its own response gives the fit back
  expect_equal(coef(ilm_refit(fit)), coef(fit), tolerance = 1e-6)
})

test_that("a categorical response is refitted from category numbers", {
  dm <- sim_mlmm(seed = 4, n_subj = 20, per = 10, J = 3)
  f <- suppressMessages(ilm_model(y ~ x1 + (1 | subj), data = dm,
                                  family = "multinomial", verbose = FALSE))
  ys <- ilm_simulate(f, nsim = 1, seed = 5)[, 1]
  r <- suppressMessages(ilm_refit(f, y = ys))
  expect_s3_class(r, "ilm_model")
  expect_identical(levels(r$model$y), levels(dm$y))
  expect_identical(as.integer(r$model$y), as.integer(ys))
  ## labels do as well as numbers
  r2 <- suppressMessages(ilm_refit(f, y = f$ylevels[ys]))
  expect_equal(coef(r2), coef(r))
})

test_that("other data rebuild the model, a correlation over time included", {
  ## one observation per cell of the AR grid, which draws warnings and a
  ## boundary message that are not what is tested here
  quiet <- function(e) suppressMessages(suppressWarnings(e))
  d <- rf_data()
  fit <- quiet(ilm_model(y ~ x + (1 | id), data = d, family = "gaussian",
                         ar = ilm_ar1(~ t | id), verbose = FALSE))
  sub <- d[as.integer(d$id) <= 8, ]
  r <- quiet(ilm_refit(fit, data = sub))
  direct <- quiet(ilm_model(y ~ x + (1 | id), data = sub, family = "gaussian",
                            ar = ilm_ar1(~ t | id), verbose = FALSE))
  expect_equal(coef(r), coef(direct))
  expect_equal(nrow(r$X), nrow(sub))
  expect_identical(nlevels(droplevels(r$model$id)), 8L)
  ## the call now names the new data
  expect_identical(deparse(r$call$data), "sub")
})

test_that("what cannot be refitted says why", {
  d <- rf_data()
  fit <- ilm_model(y ~ x + (1 | id), data = d, family = "gaussian",
                   verbose = FALSE)
  expect_error(ilm_refit(fit, data = d, y = fit$y), "but not both")
  expect_error(ilm_refit(fit, y = fit$y[-1]), "one value per row")
  expect_error(ilm_refit(fit, y = cbind(fit$y, fit$y)), "one response")
  expect_error(ilm_refit(mtcars), "must be a fitted ilm_model")
  ## a correlation over time given as vectors belongs to the old rows
  fv <- suppressWarnings(ilm_model(y ~ x, data = d, family = "gaussian",
                                   ar = ilm_ar1(d$t, d$id), verbose = FALSE))
  expect_error(ilm_refit(fv, data = d[1:36, ]), "give it by name")
})

test_that("a censored response stays censored where the data were", {
  set.seed(6); n <- 150
  d <- data.frame(x = stats::rnorm(n))
  lat <- 1 + 0.5 * d$x + stats::rnorm(n)
  d$y <- pmin(lat, 1.5)                   # censored above at 1.5
  fit <- suppressMessages(ilm_model(y ~ x, data = d, family = "gaussian",
                                    censor = ilm_censor(d$y, upper = 1.5),
                                    verbose = FALSE))
  ys <- ilm_simulate(fit, nsim = 1, seed = 7)[, 1]
  expect_true(max(ys) <= 1.5)
  r <- ilm_refit(fit, y = ys)
  expect_equal(coef(r), coef(illume:::ilm_refit_like(fit, y = ys)))
})
