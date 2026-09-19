test_that("parameter counts match the structures they describe", {
  expect_equal(lum_ncov(2L), 3L)          # 2x2 symmetric: 3 free entries
  expect_equal(lum_ncov(9L), 45L)         # the J = 10 case
  expect_equal(lum_nrr(9L, 1L), 9L)
  expect_equal(lum_nrr(9L, 4L), 30L)
  expect_lt(lum_nrr(9L, 4L), lum_ncov(9L))    # reduced rank must cost less
  expect_equal(lum_npar_d(1L), 0L)        # no within-group covariance to estimate
  expect_equal(lum_npar_d(2L, TRUE), lum_ncov(2L) - 1L)   # one fixed for the scale
  expect_equal(lum_npar_d(2L, FALSE), 1L)
})

test_that("covariance factors rebuild positive-definite matrices", {
  v <- c(0.2, -0.3, 0.1)
  L <- lum_mkL_num(v, 2L)
  S <- L %*% t(L)
  expect_true(all(eigen(S, only.values = TRUE)$values > 0))
  expect_equal(L[1, 2], 0)                       # lower triangular
  expect_equal(diag(L), exp(v[c(1, 3)]))         # exp() on the diagonal
})

test_that("diagonal structure has no off-diagonal correlation", {
  S <- lum_mkD_num(c(0.1, -0.2), 2L); S <- S %*% t(S)
  expect_equal(S[1, 2], 0)
})

test_that("reduced-rank covariance has exactly the requested rank", {
  Lam <- lum_mkLam_num(c(0.1, 0.2, -0.1, 0.3, 0.05), 3L, 2L)
  expect_equal(dim(Lam), c(3L, 2L))
  S <- Lam %*% t(Lam)
  ev <- eigen(S, only.values = TRUE)$values
  expect_gt(ev[2], 1e-8)              # first two dimensions are real
  expect_lt(abs(ev[3]), 1e-8)         # third is exactly absent
})

test_that("within-group factor fixes its first element at 1", {
  Ld <- lum_mkLd_num(c(-0.2, 0.15), 2L, TRUE)
  expect_equal(Ld[1, 1], 1)           # resolves the Kronecker scale
  Ld0 <- lum_mkLd_num(0.1, 2L, FALSE)
  expect_equal(Ld0[2, 1], 0)          # uncorrelated variant
})

test_that("lum_msqrt handles singular matrices where chol would fail", {
  S <- matrix(c(1, 1, 1, 1), 2, 2)    # rank 1
  expect_error(chol(S))
  A <- lum_msqrt(S)
  expect_equal(A %*% A, S, tolerance = 1e-8)
})

test_that("lum_rec_rank suggests something the levels can support", {
  r <- lum_rec_rank(9L, 200L)
  expect_gte(r, 1L)
  expect_lte(lum_nrr(9L, r), 200L / 6)    # lands in the safe band
  expect_equal(lum_rec_rank(9L, 12L), 0L) # too few levels for any rank
})

test_that("structure helpers agree with each other", {
  us <- list(type = "us"); rr <- list(type = "rr", rank = 2L)
  expect_equal(lum_str_npar(us, 4L), lum_ncov(4L))
  expect_equal(lum_str_width(us, 4L), 4L)
  expect_equal(lum_str_width(rr, 4L), 2L)   # fewer latent values, the real saving
  expect_equal(lum_str_label(rr), "rr(2)")
})

# A reduced-rank structure is specified by hand, so its input must be checked
# properly: users type rank = 2, not rank = 2L, and a wrong rank should give a
# sentence explaining what to do rather than an internal vapply error.

test_that("a rank given as a double is accepted", {
  dd <- sim_mlmm(seed = 4, n_subj = 20, per = 10, J = 4)   # C = 3
  f <- lum_model(y ~ x1 + (1 | subj), data = dd, family = "multinomial",
                 re_struct = list(subj = list(type = "rr", rank = 2)),
                 verbose = FALSE)
  expect_s3_class(f, "lum_model")
  expect_identical(f$re_struct$subj$rank, 2L)
})

test_that("an unusable rank is refused with a helpful message", {
  dd <- sim_mlmm(seed = 4, n_subj = 20, per = 10, J = 4)   # C = 3
  bad <- function(s)
    tryCatch(lum_model(y ~ x1 + (1 | subj), data = dd, family = "multinomial",
                       re_struct = list(subj = s), verbose = FALSE),
             error = conditionMessage)
  expect_match(bad(list(type = "rr")), "needs a rank")
  expect_match(bad(list(type = "rr", rank = 2.5)), "whole number")
  expect_match(bad(list(type = "rr", rank = 0)), "whole number")
  expect_match(bad(list(type = "rr", rank = 99)), "only 3 category dimensions")
  expect_match(bad(list(type = "banana")), 'one of "us", "diag" or "rr"')
})
