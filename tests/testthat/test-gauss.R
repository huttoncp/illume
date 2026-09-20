# The gaussian agreement index. Two properties matter more than accuracy on
# any one distribution: it must not fire on genuinely normal data at any n,
# and it must not degenerate as n grows the way a normality test does.

test_that("normal data scores 1 at every sample size", {
  for (n in c(30, 100, 1000, 10000)) {
    set.seed(7)
    g <- ilm_gauss_check(stats::rnorm(n))
    expect_equal(g$gauss, 1, info = paste("n =", n))
    expect_equal(g$gauss_note, "", info = paste("n =", n))
  }
})

test_that("the index does not degenerate at large n", {
  # a normality test would reject here; an effect size should not
  set.seed(1)
  expect_equal(ilm_gauss_check(stats::rnorm(1e5))$gauss, 1)
})

test_that("false alarms on normal data stay rare", {
  set.seed(11)
  flagged <- mean(replicate(120, ilm_gauss_check(stats::rnorm(200))$gauss < 0.95))
  expect_lt(flagged, 0.05)
})

test_that("departures are scored low with the right reason", {
  set.seed(2)
  expect_match(ilm_gauss_check(stats::rlnorm(2000, 0, 0.6))$gauss_note,
               "right-skewed")
  expect_match(ilm_gauss_check(-stats::rlnorm(2000, 0, 0.6))$gauss_note,
               "left-skewed")
  expect_match(ilm_gauss_check(stats::rt(2000, 3))$gauss_note, "heavy-tailed")
  expect_match(ilm_gauss_check(stats::rpois(2000, 3))$gauss_note, "discrete")
  expect_match(ilm_gauss_check(stats::rexp(2000))$gauss_note, "bounded at zero")
})

test_that("a symmetric heavy-tailed sample is not called skewed", {
  # moment skewness is unstable under heavy tails and called t(3) right-skewed;
  # the quartile-based measure does not
  set.seed(3)
  expect_false(grepl("skewed", ilm_gauss_check(stats::rt(3000, 3))$gauss_note))
})

test_that("mixtures are identified as multimodal, not as some other family", {
  # this is the case a single-family fit gets confidently wrong
  set.seed(4)
  g <- ilm_gauss_check(c(stats::rnorm(1000), stats::rnorm(1000, 5)))
  expect_match(g$gauss_note, "multimodal")
  expect_lt(g$gauss, 0.6)
})

test_that("ripple is not mistaken for modes", {
  # uniform, exponential and Poisson all produced spurious peaks before
  # prominence was required
  set.seed(5)
  for (x in list(stats::runif(2000), stats::rexp(2000), stats::rpois(2000, 3)))
    expect_false(grepl("multimodal", ilm_gauss_check(x)$gauss_note))
})

test_that("small samples say so instead of guessing", {
  set.seed(6)
  g <- ilm_gauss_check(stats::rnorm(10))
  expect_true(is.na(g$gauss))
  expect_match(g$gauss_note, "n too small")
})

test_that("a constant vector is reported as constant, not as too small", {
  expect_match(ilm_gauss_check(rep(3, 100))$gauss_note, "constant")
})

test_that("ks_d is the measured quantity and does not move with cap", {
  set.seed(8)
  x <- stats::rlnorm(1000)
  a <- ilm_gauss_check(x, cap = 0.06)
  b <- ilm_gauss_check(x, cap = 0.24)
  expect_equal(a$ks_d, b$ks_d)          # the distance is a property of x
  expect_lt(a$gauss, b$gauss)           # only the index rescales
})

test_that("the gauss argument selects which columns appear", {
  d <- ilm_sim(n_id = 20)
  expect_true("gauss" %in% names(ilm_describe(d, "score", gauss = "index")))
  expect_true("ks_d" %in% names(ilm_describe(d, "score", gauss = "ks_d")))
  both <- names(ilm_describe(d, "score", gauss = "both"))
  expect_true(all(c("gauss", "ks_d") %in% both))
  none <- names(ilm_describe(d, "score", gauss = "none"))
  expect_false(any(c("gauss", "ks_d", "gauss_note") %in% none))
})
