test_that("a weighted fit equals the expanded fit exactly", {
  # frequency weights are equivalent to replicating rows, so this must hold to
  # numerical precision for coefficients AND standard errors
  dd <- sim_mlmm(seed = 11, n_subj = 25, per = 24)
  ag <- aggregate(list(w = rep(1L, nrow(dd))),
                  by = list(subj = dd$subj, grp = dd$grp, y = dd$y), FUN = sum)
  ag$x1 <- 0; dd2 <- dd; dd2$x1 <- 0     # drop the continuous term so rows pool
  f_exp <- ilm_model(y ~ grp + (1 | subj), data = dd2, family = "multinomial", verbose = FALSE)
  f_agg <- ilm_model(y ~ grp + (1 | subj), data = ag, family = "multinomial", weights = w, verbose = FALSE)
  expect_lt(max(abs(coef(f_exp) - coef(f_agg))), 1e-5)
  expect_lt(max(abs(sqrt(diag(vcov(f_exp))) - sqrt(diag(vcov(f_agg))))), 1e-5)
  expect_lt(abs(as.numeric(logLik(f_exp)) - as.numeric(logLik(f_agg))), 1e-5)
  expect_lt(nrow(ag), nrow(dd2))         # and it is genuinely smaller
})

test_that("weights that look like sampling weights are flagged", {
  dd <- sim_mlmm(seed = 12)
  dd$wn <- runif(nrow(dd), 1, 5)         # non-integer: cannot be counts
  f <- ilm_model(y ~ x1 + (1 | subj), data = dd, family = "multinomial", weights = wn, verbose = FALSE)
  wt <- f$checks[f$checks$check == "weights_type", ]
  expect_equal(wt$status, "FAIL")
  expect_match(wt$cause, "sampling weights")

  dd$wf <- rep(10, nrow(dd))             # flat scaling, not replicate counts
  f2 <- ilm_model(y ~ x1 + (1 | subj), data = dd, family = "multinomial", weights = wf, verbose = FALSE)
  expect_equal(f2$checks$status[f2$checks$check == "weights_type"], "WARN")
})

test_that("genuine replicate counts pass the check", {
  dd <- sim_mlmm(seed = 13)
  dd$wc <- sample(1:3, nrow(dd), TRUE)
  f <- ilm_model(y ~ x1 + (1 | subj), data = dd, family = "multinomial", weights = wc, verbose = FALSE)
  expect_equal(f$checks$status[f$checks$check == "weights_type"], "OK")
})

test_that("negative or mis-sized weights are refused", {
  dd <- sim_mlmm(seed = 14)
  X <- model.matrix(~ x1, dd)
  expect_error(ilm_fit(X, as.integer(dd$y), 3L, list(subj = dd$subj), family = "multinomial",
                        weights = rep(-1, nrow(dd)), verbose = FALSE),
               "non-negative")
  expect_error(ilm_fit(X, as.integer(dd$y), 3L, list(subj = dd$subj), family = "multinomial",
                        weights = c(1, 2), verbose = FALSE), "one entry per row")
})
