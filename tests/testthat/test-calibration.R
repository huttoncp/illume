# Calibration curves with a simulation envelope.
#
# An in-sample calibration curve on its own is close to self-fulfilling:
# fitting with an intercept forces the mean residual to zero, so
# calibration-in-the-large is guaranteed. What the envelope tests is the shape.

sim_bin <- function(kind, seed, n = 1500) {
  set.seed(seed)
  d <- data.frame(g = factor(rep(seq_len(50), each = n / 50)),
                  x = stats::rnorm(n))
  b <- stats::rnorm(50, 0, 0.5)[as.integer(d$g)]
  eta <- 0.1 + 1.2 * d$x + b
  p <- switch(kind,
    correct = 1 / (1 + exp(-eta)),
    thresh  = 1 / (1 + exp(-(0.1 + 2.5 * (d$x > 0.5) + b))))
  d$y <- stats::rbinom(n, 1, p)
  d
}

fit_bin <- function(d)
  ilm_model(y ~ x + (1 | g), data = d, family = "binomial", verbose = FALSE)

test_that("calibration works for a binomial model", {
  r <- suppressWarnings(ilm_calibration(fit_bin(sim_bin("correct", 1)), B = 200))
  expect_length(r, 2L)                       # one entry per outcome
  cal <- r[[2]]
  expect_true(all(c("mean_p", "obs", "lo", "hi", "n_bins", "n_outside",
                    "p_outside", "status") %in% names(cal)))
  expect_true(all(cal$lo <= cal$hi, na.rm = TRUE))
  # predicted probability should increase across bins by construction
  expect_true(!is.unsorted(cal$mean_p, na.rm = TRUE))
})

test_that("a correctly specified model is not flagged as failing", {
  # The thresholds are graded against Binomial(bins, 0.05). WARN is a soft
  # signal that fires on roughly 9% of correct models by construction, so the
  # thing to assert is the absence of FAIL, not the absence of WARN.
  st <- vapply(1:4, function(s) {
    r <- suppressWarnings(ilm_calibration(fit_bin(sim_bin("correct", s)), B = 200))
    r[[2]]$status
  }, "")
  expect_false(any(st == "FAIL"))
  expect_true(sum(st == "OK") >= 3L)
})

test_that("too few replicates is warned about rather than silently trusted", {
  expect_warning(ilm_calibration(fit_bin(sim_bin("correct", 1)), B = 60),
                 "unstable envelope")
})

test_that("missing nonlinearity is detected", {
  st <- vapply(1:3, function(s) {
    r <- suppressWarnings(ilm_calibration(fit_bin(sim_bin("thresh", s)), B = 200))
    r[[2]]$status
  }, "")
  expect_true(all(st == "FAIL"))
})

test_that("the status follows the binomial tail probability", {
  r <- suppressWarnings(ilm_calibration(fit_bin(sim_bin("correct", 1)), B = 60))
  cal <- r[[2]]
  expected <- stats::pbinom(cal$n_outside - 1L, cal$n_bins, 0.05,
                            lower.tail = FALSE)
  expect_equal(cal$p_outside, round(expected, 4))
})

test_that("calibration still works for the multinomial family", {
  set.seed(3); n <- 900
  d <- data.frame(g = factor(rep(seq_len(30), each = 30)), x = stats::rnorm(n))
  b <- stats::rnorm(30, 0, 0.5)[as.integer(d$g)]
  eta <- cbind(0.3 + 0.8 * d$x + b, -0.2 + 0.4 * d$x + b)
  P <- exp(eta %*% t(stats::contr.sum(3))); P <- P / rowSums(P)
  d$y <- factor(c("a", "b", "c")[apply(P, 1, function(p) sample.int(3, 1, prob = p))])
  f <- ilm_model(y ~ x + (1 | g), data = d, family = "multinomial",
                 verbose = FALSE)
  r <- suppressWarnings(ilm_calibration(f, B = 200))
  expect_length(r, 3L)
  expect_true(all(vapply(r, function(z) is.null(z) || "status" %in% names(z), TRUE)))
})

test_that("families without a predicted probability are refused clearly", {
  set.seed(4); n <- 600
  d <- data.frame(g = factor(rep(seq_len(30), each = 20)), x = stats::rnorm(n))
  d$y <- 1 + d$x + stats::rnorm(n)
  f <- ilm_model(y ~ x + (1 | g), data = d, family = "gaussian",
                 verbose = FALSE)
  expect_error(ilm_calibration(f), "binomial and multinomial")
  expect_error(ilm_calibration(f), "ilm_rqr")
})
