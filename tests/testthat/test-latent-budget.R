# The latent-budget checks for families with one linear predictor: warned,
# not failed, until they are recalibrated on information per latent value.

lb_counts <- function(seed = 1) {
  set.seed(seed); G <- 20; Tn <- 15
  d <- expand.grid(t = seq_len(Tn), g = factor(seq_len(G)))
  d$x <- stats::rnorm(nrow(d))
  lat <- unlist(lapply(seq_len(G), function(i) as.numeric(
    stats::arima.sim(list(ar = 0.7), Tn, sd = 0.5 * sqrt(1 - 0.49)))))
  d$y <- stats::rpois(nrow(d), exp(log(3) + 0.3 * d$x + lat))
  d
}

test_that("a count at one observation per cell is warned, not failed", {
  ## ilm_ar1() warns about one observation per cell as it is built; the
  ## check is what is tested here
  f <- suppressWarnings(suppressMessages(ilm_model(y ~ x, data = lb_counts(), family = "poisson",
                                  ar = ilm_ar1(~ t | g), verbose = FALSE)))
  ck <- f$checks
  for (nm in c("latent_budget", "obs_per_ar_latent")) {
    expect_identical(ck$status[ck$check == nm], "WARN")
    expect_match(ck$cause[ck$check == nm], "provisional")
  }
  ## the suggestion leads with coarsening, and no longer offers s(time)
  sug <- ck$suggestion[ck$check == "obs_per_ar_latent"]
  expect_match(sug, "^coarsen")
  expect_false(grepl("s(time)", sug, fixed = TRUE))
  ## the remedies offer s(time) as a description of the past only
  r <- ilm_remedies(f)
  est <- r$remedy[r$check == "obs_per_ar_latent" & r$tier == "estimand"]
  expect_length(est, 1L)
  expect_match(est, "not to forecast")
})

test_that("a multinomial keeps the thresholds it was calibrated on", {
  d <- sim_mlmm(seed = 4, n_subj = 20, per = 3, J = 4)
  f <- suppressMessages(suppressWarnings(
    ilm_model(y ~ x1 + (1 | subj), data = d, family = "multinomial",
              verbose = FALSE)))
  st <- f$checks$status[f$checks$check == "latent_budget"]
  expect_identical(st, "FAIL")
  expect_false(grepl("provisional",
                     f$checks$cause[f$checks$check == "latent_budget"]))
})
