mar_data <- function(n = 500L, seed = 11L, on = c("outcome", "covariate", "none")) {
  on <- match.arg(on)
  set.seed(seed)
  z <- rnorm(n); x <- 0.6 * z + rnorm(n); y <- 0.5 * x + 0.3 * z + rnorm(n)
  d <- data.frame(x = x, z = z, y = y)
  p <- switch(on, outcome = plogis(1.0 * y - 0.4),
              covariate = plogis(1.2 * z - 0.4), none = rep(0.3, n))
  d$x[stats::runif(n) < p] <- NA
  d
}

test_that("check_missing reports how much is missing and where", {
  d <- mar_data()
  r <- ilm_check_missing(d, y = "y", verbose = FALSE)
  expect_s3_class(r, "ilm_missing")
  expect_equal(r$n, nrow(d))
  expect_equal(r$n_complete, sum(stats::complete.cases(d)))
  expect_equal(r$variables$n_missing[r$variables$variable == "x"],
               sum(is.na(d$x)))
  expect_equal(r$variables$n_missing[r$variables$variable == "y"], 0L)
  expect_true(nrow(r$patterns) >= 1L)

  ## nothing missing is said plainly rather than analysed
  clean <- ilm_check_missing(d[stats::complete.cases(d), ], verbose = FALSE)
  expect_equal(clean$verdict, "NONE")
  expect_null(clean$associations)
  expect_output(print(clean), "no missing values")
})

test_that("check_missing separates missingness on a covariate from on the outcome", {
  ## This is the distinction the whole function exists for: complete cases are
  ## unbiased in the first case and biased in the second, so the advice differs.
  cov <- ilm_check_missing(mar_data(on = "covariate"), y = "y", verbose = FALSE)
  expect_equal(cov$verdict, "RELATED_TO_COVARIATES")
  fl <- cov$associations[cov$associations$flag, ]
  expect_true("z" %in% fl$related_to)
  ## The MARGINAL table does flag the outcome, and correctly: missingness
  ## depends on z, the outcome depends on z, so the two are associated. That is
  ## exactly why the verdict cannot be read off this table -- the conditional
  ## test is what decides, and it clears the outcome here.
  expect_true(any(fl$is_outcome))
  expect_false(any(cov$outcome_test$flag))
  expect_lt(cov$outcome_test$effect[1], 0.1)
  expect_output(print(cov), "complete cases stay unbiased")

  out <- ilm_check_missing(mar_data(on = "outcome"), y = "y", verbose = FALSE)
  expect_equal(out$verdict, "RELATED_TO_OUTCOME")
  expect_true(any(out$outcome_test$flag))
  expect_gt(out$outcome_test$effect[1], 0.2)
  expect_output(print(out), "complete cases are biased")

  ## and missingness unrelated to anything is not talked up
  expect_equal(ilm_check_missing(mar_data(on = "none"), y = "y",
                                 verbose = FALSE)$verdict, "MCAR_NOT_REJECTED")
})

test_that("check_missing needs an effect as well as a p-value", {
  ## a trivial association at a large n is significant and not worth acting on
  set.seed(5); n <- 6000
  d <- data.frame(z = rnorm(n), y = rnorm(n))
  d$x <- rnorm(n)
  d$x[plogis(0.06 * d$z - 0.8) > runif(n)] <- NA
  loose <- ilm_check_missing(d, y = "y", min_effect = 0.001, verbose = FALSE)
  strict <- ilm_check_missing(d, y = "y", min_effect = 0.3, verbose = FALSE)
  expect_gte(sum(loose$associations$flag), sum(strict$associations$flag))
  expect_equal(strict$verdict, "MCAR_NOT_REJECTED")
})

test_that("check_missing detects a monotone pattern", {
  set.seed(3); n <- 200
  d <- data.frame(a = rnorm(n), b = rnorm(n), c = rnorm(n))
  drop <- sample(n, 60)
  d$b[drop] <- NA
  d$c[drop] <- NA; d$c[sample(setdiff(seq_len(n), drop), 30)] <- NA
  expect_true(ilm_check_missing(d, verbose = FALSE)$monotone)

  d2 <- d; d2$b[sample(n, 40)] <- NA      # now they cross
  expect_false(ilm_check_missing(d2, verbose = FALSE)$monotone)
})

test_that("imputation draws rather than fits", {
  d <- mar_data(on = "outcome")
  imp <- ilm_impute(d, m = 5, seed = 1, verbose = FALSE)
  expect_s3_class(imp, "ilm_mids")
  expect_equal(imp$m, 5L)
  expect_equal(imp$incomplete, "x")
  expect_equal(unname(imp$families[["x"]]), "gaussian")
  expect_length(imp$imputations, 5L)
  expect_true(all(vapply(imp$imputations, function(z) !anyNA(z$x), TRUE)))
  ## observed values are never altered
  obs <- !is.na(d$x)
  for (z in imp$imputations) expect_equal(z$x[obs], d$x[obs])
  ## and the fills DIFFER between imputations -- if they did not, the
  ## between-imputation variance would be zero and pooling would be pointless
  mi <- which(is.na(d$x))
  f1 <- imp$imputations[[1]]$x[mi]; f2 <- imp$imputations[[2]]$x[mi]
  expect_false(isTRUE(all.equal(f1, f2)))
  expect_gt(stats::sd(f1 - f2), 0)
})

test_that("pooling follows Rubin's rules and widens with missingness", {
  d <- mar_data(on = "outcome")
  p <- ilm_mi_pool(ilm_impute(d, m = 10, seed = 1, verbose = FALSE),
                   y ~ x + z, family = "gaussian")
  expect_s3_class(p, "ilm_pooled")
  expect_setequal(p$term, c("(Intercept)", "x", "z"))
  expect_true(all(p$se > 0))
  expect_true(all(p$lower < p$estimate & p$estimate < p$upper))
  ## fmi is a fraction, and the incomplete variable loses the most information
  expect_true(all(p$fmi >= 0 & p$fmi <= 1))
  expect_gt(p$fmi[p$term == "x"], p$fmi[p$term == "z"])
  ## Barnard-Rubin df never exceeds what the complete data could supply
  expect_true(all(p$df <= nrow(d)))
  expect_true(all(p$df > 1))
  expect_output(print(p), "Rubin")
})

test_that("multiple imputation recovers what complete cases lose", {
  ## missingness driven by the outcome: complete cases are biased, and this is
  ## the case imputation exists for
  d <- mar_data(on = "outcome")
  cc <- unname(coef(ilm_model(y ~ x + z, data = d, family = "gaussian",
                              verbose = FALSE))["x"])
  p <- ilm_mi_pool(ilm_impute(d, m = 20, seed = 1, verbose = FALSE),
                   y ~ x + z, family = "gaussian")
  mi <- p$estimate[p$term == "x"]
  expect_lt(abs(mi - 0.5), abs(cc - 0.5))      # closer to the truth
  expect_lt(p$lower[p$term == "x"], 0.5)       # and its interval covers it
  expect_gt(p$upper[p$term == "x"], 0.5)
})

test_that("single imputation is available but warns", {
  d <- mar_data()
  expect_warning(s <- ilm_impute(d, single = TRUE, seed = 1, verbose = FALSE),
                 "standard errors that are too small")
  expect_s3_class(s, "data.frame")
  expect_false(anyNA(s$x))
  expect_equal(nrow(s), nrow(d))
})

test_that("imputation handles other column types and leaves the rest alone", {
  set.seed(4); n <- 400
  d <- data.frame(num = rnorm(n), cnt = rpois(n, 3),
                  bin = factor(sample(c("no", "yes"), n, TRUE)),
                  keep = rnorm(n))
  d$num[sample(n, 60)] <- NA
  d$cnt[sample(n, 50)] <- NA
  d$bin[sample(n, 40)] <- NA
  imp <- ilm_impute(d, m = 3, seed = 1, verbose = FALSE)
  expect_setequal(imp$incomplete, c("num", "cnt", "bin"))
  expect_equal(unname(imp$families[["cnt"]]), "poisson")
  expect_equal(unname(imp$families[["bin"]]), "binomial")
  one <- imp$imputations[[1]]
  expect_false(anyNA(one$num)); expect_false(anyNA(one$cnt))
  expect_false(anyNA(one$bin))
  ## drawn counts stay non-negative integers, drawn factors stay in-level
  expect_true(all(one$cnt >= 0 & one$cnt == round(one$cnt)))
  expect_setequal(levels(one$bin), levels(d$bin))
  ## a complete column is untouched
  expect_equal(one$keep, d$keep)
})

test_that("imputation refuses or declines gracefully", {
  expect_error(ilm_impute("nope"), "must be a data frame")
  expect_error(ilm_impute(mar_data(), m = 0), "at least 1")
  ## nothing missing means nothing to do, said rather than errored
  clean <- mar_data()[stats::complete.cases(mar_data()), ]
  expect_message(r <- ilm_impute(clean, verbose = TRUE), "no missing values")
  expect_equal(r$m, 1L)
  expect_error(ilm_mi_pool(ilm_impute(mar_data(), m = 2, seed = 1,
                                      verbose = FALSE)), "`formula` is needed")
})

test_that("ilm_model says how many rows it dropped", {
  d <- mar_data()
  expect_message(f <- ilm_model(y ~ x + z, data = d, family = "gaussian"),
                 "dropped for missing values")
  expect_equal(f$n_dropped, sum(is.na(d$x)))
  ## and stays quiet when asked to
  expect_silent(ilm_model(y ~ x + z, data = d, family = "gaussian",
                          verbose = FALSE))
  ## no missing data, nothing to report
  f2 <- ilm_model(y ~ x + z, data = d[stats::complete.cases(d), ],
                  family = "gaussian", verbose = FALSE)
  expect_equal(f2$n_dropped, 0L)
})

test_that("the imputation agrees with mice", {
  skip_if_not_installed("mice")
  d <- mar_data(on = "outcome")
  a <- ilm_mi_pool(ilm_impute(d, m = 20, seed = 1, verbose = FALSE),
                   y ~ x + z, family = "gaussian")
  mi <- mice::mice(d, m = 20, printFlag = FALSE, seed = 1)
  b <- summary(mice::pool(with(mi, stats::lm(y ~ x + z))))
  for (t in c("x", "z")) {
    ia <- a[a$term == t, ]; ib <- b[b$term == t, ]
    ## different RNG and slightly different imputation models, so agreement is
    ## judged on the scale of the uncertainty rather than digit by digit
    expect_lt(abs(ia$estimate - ib$estimate), 0.75 * ia$se)
    expect_gt(ia$se / ib$std.error, 0.8)
    expect_lt(ia$se / ib$std.error, 1.25)
  }
})
