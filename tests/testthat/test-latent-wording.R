## The latent budget of a correlation over time is judged where the family is
## known. For a gaussian response the Laplace approximation is exact, so one
## observation per latent value is no reason to warn about it; what it does
## cost is separating the residual SD from the latent SD, which the check says.

q <- function(e) suppressMessages(suppressWarnings(e))
one_per_latent <- function(n = 160, seed = 1) {
  set.seed(seed)
  d <- data.frame(g = factor("s"), t = seq_len(n))
  d$y <- 5 + as.numeric(stats::arima.sim(list(ar = 0.7), n)) + stats::rnorm(n, 0, 0.5)
  d$k <- stats::rpois(n, exp(0.5 + 0.3 * as.numeric(scale(d$y))))
  d
}
warns_of <- function(expr) {
  w <- character(0)
  withCallingHandlers(expr, warning = function(cnd) {
    w <<- c(w, conditionMessage(cnd)); invokeRestart("muffleWarning") })
  w
}

test_that("a gaussian AR(1) at one observation per latent is not warned about the Laplace approximation", {
  d <- one_per_latent()
  w <- warns_of(f <- ilm_model(y ~ 1, data = d, family = "gaussian",
                               ar = ilm_ar1(~ t | g), verbose = FALSE))
  expect_false(any(grepl("Laplace|per latent", w)))
  ck <- f$checks[f$checks$check == "obs_per_ar_latent", ]
  expect_identical(ck$status, "OK")
  expect_match(ck$detail, "weakly separable", fixed = TRUE)
  expect_match(ck$detail, "hold pass", fixed = TRUE)
})

test_that("a count AR(1) at one observation per latent warns as its check says", {
  d <- one_per_latent()
  w <- warns_of(f <- ilm_model(k ~ 1, data = d, family = "poisson",
                               ar = ilm_ar1(~ t | g), verbose = FALSE))
  ck <- f$checks[f$checks$check == "obs_per_ar_latent", ]
  expect_identical(ck$status, "WARN")
  hit <- grep("^obs_per_ar_latent", w, value = TRUE)
  expect_length(hit, 1L)
  ## the console and fit$checks say the same thing
  expect_match(hit, ck$detail, fixed = TRUE)
  ## and a structure built quietly stays quiet
  w2 <- warns_of(ilm_model(k ~ 1, data = d, family = "poisson",
                           ar = ilm_ar1(~ t | g, verbose = FALSE), verbose = FALSE))
  expect_false(any(grepl("^obs_per_ar_latent", w2)))
})

test_that("a structure built from vectors does not warn before the family is known", {
  d <- one_per_latent()
  expect_silent(ilm_ar1(d$t, d$g))
  expect_silent(ilm_car1(d$t, d$g))
})
