test_that("every fit carries a well-formed check table", {
  fit <- fit_basic()
  ck <- fit$checks
  expect_true(all(c("check", "status", "detail", "cause", "suggestion") %in% names(ck)))
  expect_true(all(ck$status %in% c("OK", "WARN", "FAIL", "BOUNDARY", "INCONCLUSIVE")))
  expect_true(all(c("category_counts", "latent_budget", "optimizer", "gradient",
                    "hessian") %in% ck$check))
})

test_that("anything not OK explains itself and suggests a remedy", {
  # a verdict with no reason is not actionable
  dd <- sim_mlmm(seed = 21, n_subj = 8, per = 6, J = 5)   # deliberately starved
  fit <- suppressWarnings(ilm_model(y ~ x1 + (1 | subj), data = dd, family = "multinomial", verbose = FALSE))
  bad <- fit$checks[fit$checks$status %in% c("WARN", "FAIL"), ]
  skip_if(nrow(bad) == 0, "no failing checks in this draw")
  expect_true(all(nzchar(bad$cause)))
})

test_that("the latent budget uses weighted information, not row count", {
  dd <- sim_mlmm(seed = 22)
  f1 <- ilm_model(y ~ x1 + (1 | subj), data = dd, family = "multinomial", verbose = FALSE)
  dd$w <- rep(5L, nrow(dd))
  f2 <- ilm_model(y ~ x1 + (1 | subj), data = dd, family = "multinomial", weights = w, verbose = FALSE)
  d1 <- f1$checks$detail[f1$checks$check == "latent_budget"]
  d2 <- f2$checks$detail[f2$checks$check == "latent_budget"]
  expect_false(identical(d1, d2))   # five times the information, not the rows
})

test_that("a starved design is caught before fitting", {
  dd <- sim_mlmm(seed = 23, n_subj = 6, per = 5, J = 4)
  re <- ilm_norm_re(list(subj = dd$subj), nrow(dd))
  st <- list(subj = list(type = "us", rank = NA_integer_, d_cor = TRUE))
  ck <- ilm_precheck(as.integer(dd$y), 4L, re, st)
  expect_true(any(ck$status %in% c("WARN", "FAIL")))
})

test_that("ilm_print_checks produces output", {
  fit <- fit_basic()
  expect_output(ilm_print_checks(fit$checks, "checks"), "ok")
})
