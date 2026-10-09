## A resampled study keeps each row's exposure (review finding 3.4). The
## offset stayed at the fitted rows' length: every study of another size
## failed, and one of the same size refitted without its exposure (intercept
## -5.18 where the truth is -2).

offset_fit <- function() {
  set.seed(7)
  n <- 200
  d <- data.frame(x = stats::rnorm(n), expo = stats::runif(n, 1, 50))
  d$y <- stats::rpois(n, d$expo * exp(-2 + 0.3 * d$x))
  ilm_model(y ~ x + offset(log(expo)), data = d, family = "poisson", verbose = FALSE)
}

test_that("a resampled study carries its rows' offset into the draw and the refit", {
  f <- offset_fit()
  st <- illume:::ilm_refit_stub(f)
  r <- illume:::ilm_power_rows(100, NULL, seq_len(nrow(f$X)))
  sr <- illume:::ilm_power_stub(st, r$rows, r$copy, logical(0), FALSE)
  expect_length(sr$offset, 100L)
  expect_equal(sr$offset, f$offset[r$rows])
  set.seed(3)
  eta <- illume:::ilm_power_eta_stub(f, sr, f$beta)
  expect_equal(as.numeric(eta), as.numeric(sr$X %*% f$beta + sr$offset))
  y <- illume:::ilm_power_response(f, eta)
  ff <- illume:::ilm_refit_like(sr, y = y, restarts = 1L)
  expect_equal(as.numeric(ff$beta), as.numeric(f$beta), tolerance = 0.15)
})

test_that("power with an offset is computed, at the size asked for", {
  skip_on_cran()
  f <- offset_fit()
  pw <- ilm_power(f, n = c(100, 200), sims = 20, seed = 1, progress = FALSE)
  expect_true(all(pw$converged >= 0.95))
  ## the slope of 0.3 on these exposures is detected nearly every time (a
  ## glm-based simulation gives 1.00)
  expect_true(all(pw$power >= 0.9))
})
