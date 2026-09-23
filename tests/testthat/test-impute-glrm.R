# Imputation through a generalized low rank model. The model itself,
# ilm_glrm(), is illumex's; filling in and drawing from it is illume's.

glrm_mixed <- function(n = 400L, seed = 1L) {
  set.seed(seed)
  f <- rnorm(n); f2 <- rnorm(n)
  d <- data.frame(a = f + rnorm(n, 0, .4), b = -f + rnorm(n, 0, .4),
                  c = f2 + rnorm(n, 0, .4))
  d$g <- factor(ifelse(stats::plogis(2.5 * f) > runif(n), "hi", "lo"))
  d
}
hide <- function(d, frac = .2, seed = 1L) {
  set.seed(seed)
  mi <- matrix(runif(nrow(d) * ncol(d)) < frac, nrow(d), ncol(d))
  dm <- d
  for (j in seq_len(ncol(d))) dm[mi[, j], j] <- NA
  list(dm = dm, mi = mi)
}

test_that("ilm_impute can use it, and then categories get imputed too", {
  d <- glrm_mixed(n = 300L, seed = 6L)
  h <- hide(d, seed = 5L)
  im <- ilm_impute(h$dm, m = 5L, method = "glrm", verbose = FALSE,
                   progress = FALSE)
  expect_s3_class(im, "ilm_mids")
  expect_equal(im$method, "glrm")
  expect_equal(im$m, 5L)
  expect_false(anyNA(im$imputations[[1]]))
  expect_s3_class(im$imputations[[1]]$g, "factor")
  ## the m datasets are DRAWS, not copies of one reconstruction: identical
  ## imputations would make the pooled variance that of a single fit
  expect_false(identical(im$imputations[[1]], im$imputations[[2]]))
  expect_gt(mean(im$imputations[[1]]$g[h$mi[, 4]] == d$g[h$mi[, 4]]), 0.45)
  ## the low-rank route leaves the factor alone, which is the gap this fills
  il <- ilm_impute(h$dm, m = 3L, method = "lowrank", verbose = FALSE,
                   progress = FALSE)
  expect_true(anyNA(il$imputations[[1]]$g))
})
