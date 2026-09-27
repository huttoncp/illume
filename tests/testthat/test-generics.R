# ilm_scores(), ilm_calibration() and ilm_contrast() are S3 generics, so a
# package whose objects predict category probabilities, or hold estimates
# with a joint covariance, can give them methods.

gen_fit <- function() {
  set.seed(2); n <- 240
  d <- data.frame(x = stats::rnorm(n), g = factor(sample(c("a", "b", "c"), n, TRUE)))
  d$y <- stats::rbinom(n, 1, stats::plogis(0.3 + 0.8 * d$x))
  suppressMessages(ilm_model(y ~ x + g, data = d, family = "binomial",
                             verbose = FALSE))
}

test_that("each generic dispatches to a method another package registers", {
  on.exit({
    for (g in c("ilm_scores", "ilm_calibration", "ilm_contrast"))
      rm(list = paste0(g, ".illume_test_class"),
         envir = get(".__S3MethodsTable__.", envir = asNamespace("illume")))
  })
  obj <- structure(list(), class = "illume_test_class")
  for (g in c("ilm_scores", "ilm_calibration", "ilm_contrast")) {
    registerS3method(g, "illume_test_class",
                     function(object, ...) paste("method for", class(object)),
                     envir = asNamespace("illume"))
    expect_identical(get(g, envir = asNamespace("illume"))(obj),
                     "method for illume_test_class")
  }
})

test_that("anything without a method is refused with a sentence", {
  expect_error(ilm_scores(mtcars), "needs a fitted ilm_model")
  expect_error(ilm_calibration(mtcars), "needs a fitted ilm_model")
  expect_error(ilm_contrast(mtcars), "must be an ilm_emmeans")
})

test_that("a fitted model gets what it got before, positional arguments too", {
  f <- gen_fit()
  expect_identical(ilm_scores(f), illume:::ilm_scores.ilm_model(f))
  expect_identical(ilm_scores(f, "typical"),
                   illume:::ilm_scores.ilm_model(f, groups = "typical"))
  e <- ilm_emmeans(f, "g")
  expect_identical(ilm_contrast(e, nsim = 2000L),
                   illume:::ilm_contrast.ilm_emm(e, nsim = 2000L))
  ## B = 20 draws the warning that the envelope is unstable; the point here
  ## is only that the generic and the method agree
  cal <- suppressWarnings(ilm_calibration(f, nbins = 5L, B = 20L, seed = 3L))
  expect_identical(cal, suppressWarnings(illume:::ilm_calibration.ilm_model(
    f, nbins = 5L, B = 20L, seed = 3L)))
})
