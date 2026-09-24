test_that("parameter counts match the structures they describe", {
  expect_equal(ilm_ncov(2L), 3L)          # 2x2 symmetric: 3 free entries
  expect_equal(ilm_ncov(9L), 45L)         # the J = 10 case
  expect_equal(ilm_nrr(9L, 1L), 9L)
  expect_equal(ilm_nrr(9L, 4L), 30L)
  expect_lt(ilm_nrr(9L, 4L), ilm_ncov(9L))    # reduced rank must cost less
  expect_equal(ilm_npar_d(1L), 0L)        # no within-group covariance to estimate
  expect_equal(ilm_npar_d(2L, TRUE), ilm_ncov(2L) - 1L)   # one fixed for the scale
  expect_equal(ilm_npar_d(2L, FALSE), 1L)
})

test_that("covariance factors rebuild positive-definite matrices", {
  v <- c(0.2, -0.3, 0.1)
  L <- ilm_mkL_num(v, 2L)
  S <- L %*% t(L)
  expect_true(all(eigen(S, only.values = TRUE)$values > 0))
  expect_equal(L[1, 2], 0)                       # lower triangular
  expect_equal(diag(L), exp(v[c(1, 3)]))         # exp() on the diagonal
})

test_that("diagonal structure has no off-diagonal correlation", {
  S <- ilm_mkD_num(c(0.1, -0.2), 2L); S <- S %*% t(S)
  expect_equal(S[1, 2], 0)
})

test_that("reduced-rank covariance has exactly the requested rank", {
  Lam <- ilm_mkLam_num(c(0.1, 0.2, -0.1, 0.3, 0.05), 3L, 2L)
  expect_equal(dim(Lam), c(3L, 2L))
  S <- Lam %*% t(Lam)
  ev <- eigen(S, only.values = TRUE)$values
  expect_gt(ev[2], 1e-8)              # first two dimensions are real
  expect_lt(abs(ev[3]), 1e-8)         # third is exactly absent
})

test_that("within-group factor fixes its first element at 1", {
  Ld <- ilm_mkLd_num(c(-0.2, 0.15), 2L, TRUE)
  expect_equal(Ld[1, 1], 1)           # resolves the Kronecker scale
  Ld0 <- ilm_mkLd_num(0.1, 2L, FALSE)
  expect_equal(Ld0[2, 1], 0)          # uncorrelated variant
})

test_that("ilm_msqrt handles singular matrices where chol would fail", {
  S <- matrix(c(1, 1, 1, 1), 2, 2)    # rank 1
  expect_error(chol(S))
  A <- ilm_msqrt(S)
  expect_equal(A %*% A, S, tolerance = 1e-8)
})

test_that("ilm_rec_rank suggests something the levels can support", {
  r <- ilm_rec_rank(9L, 200L)
  expect_gte(r, 1L)
  expect_lte(ilm_nrr(9L, r), 200L / 6)    # lands in the safe band
  expect_equal(ilm_rec_rank(9L, 12L), 0L) # too few levels for any rank
})

test_that("structure helpers agree with each other", {
  us <- list(type = "us"); rr <- list(type = "rr", rank = 2L)
  expect_equal(ilm_str_npar(us, 4L), ilm_ncov(4L))
  expect_equal(ilm_str_width(us, 4L), 4L)
  expect_equal(ilm_str_width(rr, 4L), 2L)   # fewer latent values, the real saving
  expect_equal(ilm_str_label(rr), "rr(2)")
})

# A reduced-rank structure is specified by hand, so its input must be checked
# properly: users type rank = 2, not rank = 2L, and a wrong rank should give a
# sentence explaining what to do rather than an internal vapply error.

test_that("a rank given as a double is accepted", {
  dd <- sim_mlmm(seed = 4, n_subj = 20, per = 10, J = 4)   # C = 3
  f <- ilm_model(y ~ x1 + (1 | subj), data = dd, family = "multinomial",
                 re_struct = list(subj = list(type = "rr", rank = 2)),
                 verbose = FALSE)
  expect_s3_class(f, "ilm_model")
  expect_identical(f$re_struct$subj$rank, 2L)
})

test_that("an unusable rank is refused with a helpful message", {
  dd <- sim_mlmm(seed = 4, n_subj = 20, per = 10, J = 4)   # C = 3
  bad <- function(s)
    tryCatch(ilm_model(y ~ x1 + (1 | subj), data = dd, family = "multinomial",
                       re_struct = list(subj = s), verbose = FALSE),
             error = conditionMessage)
  expect_match(bad(list(type = "rr")), "needs a rank")
  expect_match(bad(list(type = "rr", rank = 2.5)), "whole number")
  expect_match(bad(list(type = "rr", rank = 0)), "whole number")
  expect_match(bad(list(type = "rr", rank = 99)), "only 3 category dimensions")
  expect_match(bad(list(type = "banana")), 'one of "us", "diag" or "rr"')
})

test_that("an re_struct element may leave its type to the default", {
  ## d_cor is all an uncorrelated random slope needs outside a multinomial
  ## model, and naming it alone used to stop, asking for a type
  set.seed(3); n <- 240
  d <- data.frame(x = rnorm(n), subj = factor(sample(20, n, TRUE)))
  d$y <- 0.5 * d$x + rnorm(20)[d$subj] + rnorm(20, 0, 0.3)[d$subj] * d$x + rnorm(n)
  f <- ilm_model(y ~ x + (1 + x | subj), data = d, family = "gaussian",
                 re_struct = list(subj = list(d_cor = FALSE)), verbose = FALSE)
  expect_s3_class(f, "ilm_model")
  expect_identical(f$re_struct$subj$type, "us")
  expect_false(f$re_struct$subj$d_cor)
  ## a type that is not one of the three still stops
  expect_error(ilm_model(y ~ x + (1 | subj), data = d, family = "gaussian",
                         re_struct = list(subj = list(type = "full")),
                         verbose = FALSE),
               "must be one of")
})
