# Each smooth's effective degrees of freedom, as mgcv defines them.

edf_data <- function(seed = 2, n = 300) {
  set.seed(seed)
  d <- data.frame(x = stats::runif(n), z = stats::runif(n),
                  id = factor(sample(12, n, TRUE)))
  d$y <- sin(2 * pi * d$x) + 0.5 * d$z^2 + stats::rnorm(12, 0, 0.5)[d$id] +
    stats::rnorm(n, 0, 0.4)
  d$bin <- stats::rbinom(n, 1, stats::plogis(sin(2 * pi * d$x) +
                                               stats::rnorm(12, 0, 0.5)[d$id]))
  d$t <- stats::ave(seq_len(n), d$id, FUN = function(i) ceiling(seq_along(i) / 5))
  d
}

## mgcv's edf of a smooth: the sum over its coefficients, "s(x).1" onwards
edf_mgcv <- function(m, lab)
  sum(m$edf[startsWith(names(m$edf), paste0(lab, "."))])

test_that("a smooth's edf is mgcv's, by ML and by REML", {
  skip_if_not_installed("mgcv")
  d <- edf_data()
  for (method in c("ML", "REML")) {
    f <- ilm_model(y ~ s(x) + z + (1 | id), data = d, family = "gaussian",
                   reml = method == "REML", verbose = FALSE)
    m <- mgcv::gam(y ~ s(x) + z + s(id, bs = "re"), data = d, method = method)
    expect_equal(unname(f$edf[["s(x)"]]), edf_mgcv(m, "s(x)"), tolerance = 1e-4)
  }
  ## and for a binary response, where W is the fit's working weights
  f <- suppressMessages(ilm_model(bin ~ s(x) + (1 | id), data = d,
                                  family = "binomial", verbose = FALSE))
  m <- mgcv::gam(bin ~ s(x) + s(id, bs = "re"), data = d, family = binomial,
                 method = "ML")
  expect_equal(unname(f$edf[["s(x)"]]), edf_mgcv(m, "s(x)"), tolerance = 1e-3)
})

test_that("a smooth shrunk to its null space has the null space's edf", {
  skip_if_not_installed("mgcv")
  d <- edf_data()
  f <- ilm_model(y ~ s(x) + s(z), data = d, family = "gaussian",
                 verbose = FALSE)
  m <- mgcv::gam(y ~ s(x) + s(z), data = d, method = "ML")
  expect_named(f$edf, c("s(x)", "s(z)"))
  expect_equal(as.numeric(f$edf), c(edf_mgcv(m, "s(x)"), edf_mgcv(m, "s(z)")),
               tolerance = 1e-4)
  ## 0.5 z^2 is all but straight over (0, 1): one unpenalised column
  expect_equal(f$edf[["s(z)"]], 1, tolerance = 1e-3)
})

test_that("the edf is kept with the fit, and is within its bounds", {
  d <- edf_data()
  ## no correlation over time in y, so the AR term sits at its boundary;
  ## what matters here is that its cells are among the random effects
  f <- suppressMessages(ilm_model(y ~ s(x) + z, data = d, family = "gaussian",
                                  ar = ilm_ar1(~ t | id), verbose = FALSE))
  expect_identical(f$edf, illume:::ilm_smooth_edf(f))
  ## between the null space alone and every basis function
  nb <- ncol(f$smooths[["s(x)"]]$Xf) + f$nlk[match("s(x)", names(f$re))]
  expect_true(f$edf[["s(x)"]] > 1 && f$edf[["s(x)"]] < nb)
  ## a model without smooths has none
  f0 <- ilm_model(y ~ z + (1 | id), data = d, family = "gaussian",
                  verbose = FALSE)
  expect_null(f0$edf)
})

test_that("a multinomial smooth's edf is its total over the categories", {
  dm <- sim_mlmm(seed = 4, n_subj = 20, per = 10, J = 3)   # C = 2
  f <- suppressMessages(ilm_model(y ~ s(x1) + (1 | subj), data = dm,
                                  family = "multinomial", verbose = FALSE))
  C <- f$C
  nb <- ncol(f$smooths[["s(x1)"]]$Xf) + f$nlk[match("s(x1)", names(f$re))]
  expect_true(f$edf[["s(x1)"]] >= C - 1e-6 && f$edf[["s(x1)"]] <= nb * C)
})

test_that("summary() and ilm_varcorr() print the edf with what the SD is", {
  d <- edf_data()
  f <- ilm_model(y ~ s(x) + z + (1 | id), data = d, family = "gaussian",
                 verbose = FALSE)
  lab <- sprintf("edf %.2f of %d", f$edf[["s(x)"]],
                 as.integer(attr(f$edf, "max")[["s(x)"]]))
  out <- paste(capture.output(summary(f)), collapse = " ")
  expect_match(out, paste0("s(x)  [us]  8 basis functions, ", lab), fixed = TRUE)
  expect_match(out, "the SD is on the basis's scale", fixed = TRUE)
  v <- ilm_varcorr(f)
  expect_equal(attr(v$re[["s(x)"]], "edf"), unname(f$edf[["s(x)"]]))
  expect_identical(attr(v$re[["s(x)"]], "edf_max"), 9)
  expect_null(attr(v$re[["id"]], "edf"))
  out <- paste(capture.output(print(v)), collapse = " ")
  expect_match(out, paste0("s(x): ", lab), fixed = TRUE)
  ## the variance table itself is lme4's, unchanged
  expect_named(as.data.frame(v), c("grp", "var1", "var2", "vcov", "sdcor"))
})
