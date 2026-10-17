# ilm_fc_ddf(): the Fai and Cornelius (1996) one-moment denominator df.
#
# Fai, A. H.-T. and Cornelius, P. L. (1996). Approximate F-tests of multiple
# degree of freedom hypotheses in generalized least squares analyses of
# unbalanced split-plot experiments. J. Statist. Comput. Simul. 54, 363-378.
#
# p. 368, eq. (4.2): E_Q = sum_m nu_m / (nu_m - 2), each nu_m > 2 (p. 369).
# p. 369: nu_GB1 = 2 E_Q / (E_Q - q).
# The paper has no numeric worked example (its Section 5 is a simulation), so
# these tests use identities that follow from the formula. Each is marked
# [OURS] where it is our derivation rather than a statement in the paper.

fc_formula <- function(nu) {
  E <- sum(ifelse(is.infinite(nu), 1, nu / (nu - 2)))
  2 * E / (E - length(nu))
}

test_that("one direction gives its own df [OURS, identity: 2E/(E - 1) = nu]", {
  expect_identical(ilm_fc_ddf(7.5), 7.5)
  expect_identical(ilm_fc_ddf(3.2), 3.2)
  # whatever its size: the paper's nu_m > 2 is needed only to combine
  # several directions (eq. 4.2), and one direction combines nothing
  expect_identical(ilm_fc_ddf(1.5), 1.5)
  expect_identical(ilm_fc_ddf(2), 2)
})

test_that("equal df give that df [OURS, identity applied q times]", {
  # identical, not just equal: the formula alone gives 12 + 5e-15, so this
  # pins the equal-df step itself
  expect_identical(ilm_fc_ddf(c(12, 12, 12)), 12)
  # within the tolerance the mean is returned: 12 + 5e-11
  expect_equal(ilm_fc_ddf(c(12, 12 + 1e-10)), 12 + 5e-11, tolerance = 1e-12)
})

test_that("unequal df follow p. 369's nu_GB1 = 2E/(E - q), eq. (4.2) for E", {
  E <- 5 / 3 + 9 / 7 + 30 / 28
  expect_equal(ilm_fc_ddf(c(5, 9, 30)), 2 * E / (E - 3), tolerance = 1e-12)
  expect_equal(ilm_fc_ddf(c(5, 9, 30)), 7.86, tolerance = 1e-3)
})

test_that("the result lies between the smallest and largest df [OURS, derived]", {
  # nu / (nu - 2) decreases in nu > 2, so E/q lies between its values at
  # max(nu) and min(nu), and so does the df that matches it.
  cases <- list(c(5, 9, 30), c(2.5, 100), c(3, 3.5, 4, 50, 1e4), c(20, 21))
  for (nu in cases) {
    d <- ilm_fc_ddf(nu)
    expect_gte(d, min(nu))
    expect_lte(d, max(nu))
    expect_equal(d, fc_formula(nu), tolerance = 1e-12)
  }
})

test_that("raising any one direction's df raises the result [OURS, derived]", {
  base <- list(c(5, 9, 30), c(3, 40), c(2.5, 6, 6, 80))
  for (nu in base) {
    for (m in seq_along(nu)) {
      up <- nu
      up[m] <- up[m] + 1
      expect_gt(ilm_fc_ddf(up), ilm_fc_ddf(nu))
    }
  }
})

test_that("with several directions, one at 2 or fewer df makes it undefined (p. 369 needs nu_m > 2)", {
  why <- "Fai and Cornelius's approximation is undefined when a direction has 2 or fewer df"
  for (nu in list(c(1.5, 30), c(2, 30), c(30, 0.5, 8))) {
    d <- ilm_fc_ddf(nu)
    expect_true(is.na(d))
    expect_identical(as.numeric(d), NA_real_)
    expect_identical(attr(d, "reason"), why)
  }
})

test_that("an infinite df contributes 1 to E [OURS, the limit of nu/(nu - 2)]", {
  E <- 1 + 10 / 8
  expect_equal(ilm_fc_ddf(c(Inf, 10)), 2 * E / (E - 2), tolerance = 1e-12)
  expect_equal(ilm_fc_ddf(c(Inf, 10)), 18, tolerance = 1e-12)
  expect_identical(ilm_fc_ddf(c(Inf, Inf)), Inf)
})

test_that("very large df neither overflow nor cancel [identity: E - q = sum 2/(nu - 2)]", {
  # nu/(nu - 2) - 1 = 2/(nu - 2), so E - q needs no subtraction of near-equal
  # numbers; written that way the result is accurate at any size
  d <- ilm_fc_ddf(c(1e17, 2e17))
  expect_true(is.finite(d))
  expect_gte(d, 1e17)
  expect_lte(d, 2e17)
  nu <- c(1e12, 2e12)
  D <- sum(2 / (nu - 2)); E <- length(nu) + D
  expect_equal(ilm_fc_ddf(nu), 2 * E / D, tolerance = 1e-12)
})

test_that("an empty or missing df is an internal error [OURS]", {
  expect_error(ilm_fc_ddf(numeric(0)), "internal")
  expect_error(ilm_fc_ddf(c(5, NA)), "internal")
  expect_error(ilm_fc_ddf(c(5, NaN)), "internal")
})
