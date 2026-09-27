## Remedies from another package's diagnostics.
##
## ilm_remedies() is a generic, and ilm_remedy_table() builds the table a
## method returns, so a remedy another package names is listed, merged,
## ordered, tied to its fit and applied exactly as illume's own are -- and a
## table that illume could not refit correctly is refused when it is built,
## not when it is applied.

## counts more variable than a poisson allows: a negative binomial is the
## remedy another package's calibration check would name
nb_fit <- function() {
  set.seed(2)
  d <- data.frame(x = rnorm(300))
  d$y <- rnbinom(300, mu = exp(1 + 0.5 * d$x), size = 1.5)
  ilm_model(y ~ x, data = d, family = "poisson", verbose = FALSE)
}

test_that("ilm_remedies() dispatches, and says what it needs otherwise", {
  f <- nb_fit()
  expect_error(ilm_remedies(data.frame(a = 1)),
               "a result whose package gives ilm_remedies\\(\\) a method")
  ## a package's own result, with its own method
  diag <- structure(list(fit = f), class = "fake_calibration")
  ilm_remedies.fake_calibration <- function(object, ...)
    ilm_remedy_table(object$fit, "pit_shape", "FAIL", "structural",
                     "refit as a negative binomial",
                     args = list(list(family = "nbinom")))
  rem <- ilm_remedies(diag)
  expect_s3_class(rem, "ilm_remedies")
  expect_identical(rem$check, "pit_shape")
  expect_identical(rem$change, "family = \"nbinom\"")
})

test_that("a table from another package refits like illume's own", {
  f <- nb_fit()
  rem <- ilm_remedy_table(f, check = "pit_shape", status = "FAIL",
                          tier = "structural",
                          remedy = "the intervals are too narrow: refit as a negative binomial",
                          args = list(list(family = "nbinom")))
  expect_output(print(rem), "structural -- pit_shape (FAIL)", fixed = TRUE)
  expect_output(print(rem), "change: family = \"nbinom\"", fixed = TRUE)
  ## the report says to run the check again: it is not one the fit makes,
  ## and must not be reported as "no longer checked"
  expect_message(f2 <- ilm_apply_remedy(f, rem, 1, reason = "PIT was U-shaped"),
                 "pit_shape: FAIL before; run it again on the new fit",
                 fixed = TRUE)
  expect_identical(f2$family$name, "nbinom")
  expect_identical(f2$remedy_log$check, "pit_shape")
  expect_identical(f2$remedy_log$tier, "structural")
  expect_identical(f2$remedy_log$reason, "PIT was U-shaped")
  ## tied to its fit, as illume's own tables are
  expect_error(ilm_apply_remedy(f2, rem, 1), "different fit")
})

test_that("remedies are merged and ordered as illume's own are", {
  f <- nb_fit()
  rem <- ilm_remedy_table(f,
    check  = c("interval_width", "pit_shape", "coverage_by_horizon"),
    status = c("WARN", "FAIL", "WARN"),
    tier   = c("estimand", "structural", "structural"),
    remedy = c("add the predictor the forecasts miss",
               "refit as a negative binomial",
               "the same: a negative binomial"),
    args   = list(NULL, list(family = "nbinom"), list(family = "nbinom")))
  ## one change named by two checks is one remedy, answering both
  expect_identical(nrow(rem), 2L)
  expect_identical(rem$tier, c("structural", "estimand"))
  expect_identical(rem$check[1], "pit_shape, coverage_by_horizon")
  expect_identical(rem$status[1], "FAIL, WARN")
  ## a remedy made by hand is listed, and refused as a refit
  expect_identical(rem$change[2], "")
  expect_error(suppressMessages(ilm_apply_remedy(f, rem, 2)), "made by hand")
  expect_message(ilm_apply_remedy(f, rem, 1),
                 "coverage_by_horizon: WARN before; run it again", fixed = TRUE)
})

test_that("a formula in a remedy means what the model's own formula means", {
  f <- nb_fit()
  ## built inside another package's function, where its names would
  ## otherwise be looked up
  mk <- function() list(list(ziformula = ~ x))
  rem <- ilm_remedy_table(f, "zero_check", "WARN", "structural",
                          "add a zero part", args = mk())
  expect_identical(environment(attr(rem, "args")[["1"]]$ziformula),
                   environment(f$formula))
})

test_that("a table illume could not refit correctly is refused when built", {
  f <- nb_fit()
  one <- function(...) {
    a <- list(check = "pit_shape", status = "FAIL", tier = "structural",
              remedy = "refit as a negative binomial",
              args = list(list(family = "nbinom")))
    a[names(list(...))] <- list(...)
    do.call(ilm_remedy_table, c(list(f), a))
  }
  expect_error(ilm_remedy_table(lm(1 ~ 1), "a", "WARN", "numerical", "b"),
               "fitted ilm_model object")
  expect_error(one(status = c("FAIL", "WARN")), "one entry per remedy")
  expect_error(one(status = "OK"), "a check that is OK names no remedy")
  expect_error(one(tier = "cosmetic"), "`tier` must be one of")
  expect_error(one(remedy = ""), "no missing or empty")
  expect_error(one(check = NA_character_), "no missing or empty")
  ## illume's own checks are ilm_remedies()'s, and the report reads them
  expect_error(one(check = "hessian"), "a check of illume's own")
  expect_error(one(check = "ilm_check_zeros"), "a check of illume's own")
  expect_error(one(check = "pit, shape"), "must not contain a comma")
  ## arguments a refit can set, one list per remedy
  ## the outer list forgotten: one remedy's arguments given as the whole
  expect_error(one(args = list(family = "nbinom")), "args = list(list(",
               fixed = TRUE)
  expect_error(one(args = list(NULL, NULL)), "one element per remedy")
  expect_error(one(args = list(list("nbinom"))), "named list")
  expect_error(one(args = list(list(uncertainty = "joint"))),
               "uncertainty, which is not an argument of ilm_model")
  expect_error(one(args = list(list(data = data.frame()))),
               "data, which is not an argument")
  ## no remedies is an empty table, printed as such
  rem0 <- ilm_remedy_table(f, character(0), character(0), character(0),
                           character(0))
  expect_identical(nrow(rem0), 0L)
  expect_output(print(rem0), "No remedies")
})

test_that("c() makes one list for a fit, with no remedy or check twice", {
  f <- nb_fit()
  own <- ilm_remedies(f, dispersion = suppressMessages(
    ilm_check_dispersion(f, B = 100)))
  expect_true(any(own$change == "family = \"nbinom\""))
  ext <- ilm_remedy_table(f, c("pit_shape", "interval_width"), c("FAIL", "WARN"),
                          c("structural", "estimand"),
                          c("refit as a negative binomial",
                            "add the predictor the forecasts miss"),
                          args = list(list(family = "nbinom"), NULL))
  all <- c(own, ext)
  expect_s3_class(all, "ilm_remedies")
  expect_identical(all$id, seq_len(nrow(all)))
  ## the negative binomial once, answering both checks that named it
  nb <- which(all$change == "family = \"nbinom\"")
  expect_length(nb, 1L)
  expect_identical(all$check[nb], "ilm_check_dispersion, pit_shape")
  expect_false(anyDuplicated(ifelse(nzchar(all$change), all$change, all$remedy)) > 0L)
  ## ordered by tier
  expect_false(is.unsorted(match(all$tier, illume:::ilm_rem_tiers)))
  ## and it applies, reporting both checks as ones to run again
  msg <- capture_messages(f2 <- ilm_apply_remedy(f, all, all$id[nb]))
  expect_match(paste(msg, collapse = ""), "ilm_check_dispersion: .* run it again")
  expect_match(paste(msg, collapse = ""), "pit_shape: FAIL before; run it again")
  expect_identical(f2$family$name, "nbinom")

  ## the same list twice is the list once
  twice <- c(ext, ext)
  expect_identical(nrow(twice), nrow(ext))
  expect_identical(twice$check, ext$check)
  ## a change two lists put in different tiers takes the more cautious
  cau <- c(ext, ilm_remedy_table(f, "coverage", "WARN", "estimand",
                                 "a negative binomial", list(list(family = "nbinom"))))
  expect_identical(cau$tier[cau$change == "family = \"nbinom\""], "estimand")
})

test_that("c() refuses tables it cannot combine", {
  f <- nb_fit()
  ext <- ilm_remedy_table(f, "pit_shape", "FAIL", "structural", "negative binomial",
                          list(list(family = "nbinom")))
  g <- suppressMessages(ilm_apply_remedy(f, ext, 1))
  other <- ilm_remedy_table(g, "pit_shape", "WARN", "structural", "a zero part",
                            list(list(ziformula = ~ 1)))
  expect_error(c(ext, other), "different fits")
  expect_error(c(ext, data.frame(a = 1)), "only remedy tables")
  ## a subset of the rows keeps its tie, and its remedies are listed once
  expect_identical(nrow(c(ext, ext[1, ])), nrow(ext))
  lost <- ext; attr(lost, "fit_id") <- NULL
  expect_error(c(ext, lost), "lost what ties it to its fit")
})
