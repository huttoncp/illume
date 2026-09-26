# A random effect's SD that the optimiser left short of zero, on a flat
# likelihood, is held as a variance at zero is.

## Another agent's case: a random intercept beside an AR(1) that takes up each
## series' level. The intercept's variance has its optimum at zero, and the
## fit stopped at an SD of 0.0204 -- above the 1e-3 line -- with an SE of 535
## on its log scale and draws running to +/-1800.
rsd_data <- function() {
  set.seed(3e6 + 1e4 + 149)
  G <- 20L; Tn <- 24L; H <- 6L; n_t <- Tn + H
  ar1 <- function(r, s) as.numeric(stats::arima.sim(list(ar = r), n_t,
                                                    sd = s * sqrt(1 - r^2)))
  full <- expand.grid(t = seq_len(n_t), g = factor(sprintf("s%02d", seq_len(G))))
  full$x <- unlist(lapply(seq_len(G), function(k) ar1(0.8, 1)))
  full$a <- unlist(lapply(seq_len(G), function(k) ar1(0.7, 0.5)))
  b <- stats::rnorm(G, 0, 0.5)
  full$y <- stats::rnorm(nrow(full), 1 + 0.8 * full$x + b[full$g] + full$a, 0.5)
  full[full$t <= Tn, ]
}

test_that("an intercept SD stopped short of zero on a flat likelihood is held", {
  d <- rsd_data()
  ## (ilm_ar1() warns of one observation per cell as it is built, which is
  ## the design here, not what is tested)
  expect_message(f <- suppressWarnings(
    ilm_model(y ~ x + (1 | g), data = d, family = "gaussian",
              ar = ilm_ar1(~ t | g), verbose = FALSE)),
    "edge of its range")
  sd_g <- sqrt(ilm_varcorr(f)$re$g[1, 1])
  expect_gt(sd_g, 1e-3)                    # above the old line
  expect_true("g" %in% f$hessian_held)
  expect_identical(f$hessian_how, "boundary")
  ## its draws stay where it is, and every natural value is finite
  dr <- ilm_draws(f, nsim = 200, seed = 1)
  th <- dr$draws[rownames(dr$draws) == "theta", ]
  expect_true(all(th == th[1]))
  expect_true(all(is.finite(dr$natural$re$g)))
  ## the fixed effects are the model without the term's, near enough
  f0 <- suppressWarnings(ilm_model(y ~ x, data = d, family = "gaussian",
                                   ar = ilm_ar1(~ t | g), verbose = FALSE))
  expect_equal(sqrt(diag(vcov(f)))[["x"]], sqrt(diag(vcov(f0)))[["x"]],
               tolerance = 0.02)
})

test_that("a small SD the data do support is left alone", {
  ## 30 groups of 10 with a true SD of 0.3: estimated, not held
  set.seed(11)
  d <- data.frame(g = factor(rep(1:30, each = 10)), x = stats::rnorm(300))
  d$y <- 0.5 + 0.3 * d$x + stats::rnorm(30, 0, 0.3)[d$g] + stats::rnorm(300, 0, 0.5)
  f <- ilm_model(y ~ x + (1 | g), data = d, family = "gaussian", verbose = FALSE)
  expect_false("g" %in% f$hessian_held)
  expect_length(f$boundary_terms, 0L)
})
