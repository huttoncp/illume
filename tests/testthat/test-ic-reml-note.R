# AIC and BIC from a REML fit compare only with fits of the same fixed
# effects, and say so in a message: they cannot refuse, since a single
# criterion does not know what it will be compared with.

ic_data <- function() {
  set.seed(23)
  g <- factor(rep(c("a", "b", "c"), c(10, 15, 20)))
  data.frame(g = g, x = rnorm(45), y = rnorm(45, 0, c(1, 2, 4)[as.integer(g)]))
}

test_that("AIC and BIC of a REML fit carry a message, and the same value", {
  d <- ic_data()
  fr <- ilm_model(y ~ x + g, data = d, family = "gaussian", dispformula = ~ g,
                  reml = TRUE, verbose = FALSE)
  expect_message(a <- AIC(fr), "AIC from a REML fit: comparable only")
  expect_message(b <- BIC(fr), "BIC from a REML fit")
  expect_message(AIC(fr), "reml = FALSE", fixed = TRUE)
  ll <- logLik(fr)
  expect_equal(a, -2 * as.numeric(ll) + 2 * attr(ll, "df"))
  expect_equal(b, -2 * as.numeric(ll) + log(attr(ll, "nobs")) * attr(ll, "df"))
})

test_that("a maximum-likelihood fit's AIC and BIC say nothing", {
  d <- ic_data()
  fm <- ilm_model(y ~ x + g, data = d, family = "gaussian", dispformula = ~ g,
                  reml = FALSE, verbose = FALSE)
  expect_no_message(AIC(fm))
  expect_no_message(BIC(fm))
})

test_that("the default fit of a dispersion model, being REML, gets the message", {
  d <- ic_data()
  f <- suppressMessages(ilm_model(y ~ g, data = d, family = "gaussian",
                                  dispformula = ~ g, verbose = FALSE))
  expect_message(AIC(f), "REML fit")
})

test_that("model_performance() says it once, and carries it as a note", {
  skip_if_not_installed("performance")
  d <- ic_data()
  fr <- ilm_model(y ~ x + g, data = d, family = "gaussian", dispformula = ~ g,
                  reml = TRUE, verbose = FALSE)
  msgs <- character()
  mp <- withCallingHandlers(performance::model_performance(fr),
                            message = function(m) {
                              msgs <<- c(msgs, conditionMessage(m))
                              invokeRestart("muffleMessage")
                            })
  expect_length(grep("from a REML fit", msgs), 1L)
  expect_match(attr(mp, "ic_note"), "same fixed effects")
  expect_no_message(performance::model_performance(fr, verbose = FALSE))
  fm <- ilm_model(y ~ x + g, data = d, family = "gaussian", dispformula = ~ g,
                  reml = FALSE, verbose = FALSE)
  expect_null(attr(performance::model_performance(fm), "ic_note"))
})

test_that("summary() and print() of a REML fit do not repeat it", {
  d <- ic_data()
  fr <- ilm_model(y ~ x + g, data = d, family = "gaussian", dispformula = ~ g,
                  reml = TRUE, verbose = FALSE)
  expect_no_message(utils::capture.output(print(fr)))
  expect_no_message(utils::capture.output(print(summary(fr))))
})
