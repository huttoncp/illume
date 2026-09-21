eff_data <- function(n = 600L, seed = 1L) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n), z = rnorm(n),
                  g = factor(sample(c("a", "b"), n, TRUE)))
  d$y <- rbinom(n, 1L, stats::plogis(-0.4 + 0.8 * d$x + 0.5 * (d$g == "b")))
  d
}

test_that("odds ratios and their intervals match glm", {
  d <- eff_data()
  f <- ilm_model(y ~ x + g, data = d, family = "binomial", verbose = FALSE)
  e <- ilm_effects(f)
  gl <- stats::glm(y ~ x + g, data = d, family = stats::binomial())
  ci <- exp(stats::confint.default(gl))
  expect_equal(e$odds_ratio, unname(exp(stats::coef(gl))), tolerance = 1e-6)
  expect_equal(e$lower, unname(ci[, 1]), tolerance = 1e-4)
  expect_equal(e$upper, unname(ci[, 2]), tolerance = 1e-4)
  ## the interval is built on the link scale and transformed, so it is
  ## asymmetric on the ratio scale -- which is the shape it really has
  expect_true(all((e$upper - e$odds_ratio) > (e$odds_ratio - e$lower)))
  ## and it never reaches a ratio below zero, which a symmetric one can
  expect_true(all(e$lower > 0))
})

test_that("each family gets the scale it belongs on", {
  set.seed(2); n <- 600L
  d <- data.frame(x = rnorm(n))
  d$cnt <- rpois(n, exp(0.5 + 0.4 * d$x))
  fp <- ilm_model(cnt ~ x, data = d, family = "poisson", verbose = FALSE)
  ep <- ilm_effects(fp)
  expect_true("rate_ratio" %in% names(ep))
  expect_equal(ep$rate_ratio,
               unname(exp(stats::coef(stats::glm(cnt ~ x, data = d,
                                                 family = stats::poisson())))),
               tolerance = 1e-6)

  ## an accelerated failure time model gives a TIME ratio, not a hazard ratio,
  ## and says so rather than letting the reader assume
  set.seed(3); n3 <- 500L
  d3 <- data.frame(x = rnorm(n3))
  d3$t <- rweibull(n3, 1.5, exp(1 + 0.6 * d3$x))
  f3 <- ilm_model(t ~ x, data = d3, family = "weibull", verbose = FALSE)
  e3 <- ilm_effects(f3)
  expect_true("time_ratio" %in% names(e3))
  expect_output(print(f3 <- ilm_effects(f3)), "acceleration factor")

  ## an ordinal fit gives a proportional odds ratio, and names the check
  set.seed(8); n8 <- 600L
  d8 <- data.frame(x = rnorm(n8))
  zz <- 0.8 * d8$x + stats::rlogis(n8)
  d8$y <- factor(cut(zz, c(-Inf, -0.8, 0.9, Inf),
                     labels = c("lo", "mid", "hi")), ordered = TRUE)
  f8 <- ilm_model(y ~ x, data = d8, family = "ordinal", verbose = FALSE)
  expect_output(print(ilm_effects(f8)), "ilm_check_proportional")
  expect_true("odds_ratio" %in% names(ilm_effects(f8)))
})

test_that("a gaussian fit gets no ratio and a partial R2 instead", {
  set.seed(4); n <- 600L
  d <- data.frame(x = rnorm(n), z = rnorm(n),
                  g = factor(sample(c("a", "b", "c"), n, TRUE)))
  d$y <- 1 + 0.5 * d$x + 0.2 * d$z + 0.8 * (d$g == "b") + rnorm(n)
  f <- ilm_model(y ~ x + z + g, data = d, family = "gaussian",
                 verbose = FALSE)
  e <- ilm_effects(f)
  ## no duplicate column: without a transform the effect IS the estimate
  expect_false("difference" %in% names(e))
  expect_equal(e$estimate, unname(stats::coef(f)))
  expect_output(print(e), "SHIFTS")
  ## partial R2 matches dropping the term and refitting
  pp <- attr(e, "partial")
  expect_setequal(pp$term, c("x", "z", "g"))
  full <- stats::lm(y ~ x + z + g, data = d)
  red <- stats::lm(y ~ z + g, data = d)
  r1 <- sum(stats::resid(full)^2); r0 <- sum(stats::resid(red)^2)
  expect_equal(pp$partial_r2[pp$term == "x"], (r0 - r1) / r0,
               tolerance = 1e-8)
  ## a multi-column term is one row, not one per column
  expect_equal(sum(pp$term == "g"), 1L)
})

test_that("standardisation scales numeric predictors and leaves binary alone", {
  set.seed(4); n <- 600L
  d <- data.frame(x = rnorm(n, sd = 3), g = factor(sample(c("a", "b"), n, TRUE)))
  d$y <- 1 + 0.5 * d$x + 0.8 * (d$g == "b") + rnorm(n)
  f <- ilm_model(y ~ x + g, data = d, family = "gaussian", verbose = FALSE)
  e0 <- ilm_effects(f)
  e1 <- ilm_effects(f, standardise = TRUE)
  e2 <- ilm_effects(f, standardise = "gelman")
  bx <- unname(stats::coef(f)["x"]); s <- stats::sd(d$x)
  expect_equal(e1$estimate[e1$term == "x"], bx * s, tolerance = 1e-8)
  expect_equal(e2$estimate[e2$term == "x"], bx * 2 * s, tolerance = 1e-8)
  ## the binary predictor is untouched, and the reason is printed
  expect_equal(e1$estimate[e1$term == "gb"], e0$estimate[e0$term == "gb"])
  expect_equal(attr(e1, "skipped"), "gb")
  expect_output(print(e1), "Left alone: gb", fixed = TRUE)
  ## the interval is built from the scaled coefficient, so z is unchanged
  expect_equal(e1$z, e0$z, tolerance = 1e-8)
  expect_error(ilm_effects(f, standardise = "nope"), "must be FALSE, TRUE")
})

test_that("a ratio from a mixed model is labelled conditional", {
  set.seed(5); ng <- 60L; ni <- 12L; n <- ng * ni
  d <- data.frame(id = factor(rep(seq_len(ng), each = ni)), x = rnorm(n))
  b <- rnorm(ng, 0, 1.2)
  d$y <- rbinom(n, 1L, stats::plogis(-0.3 + 0.9 * d$x + b[as.integer(d$id)]))
  f <- ilm_model(y ~ x + (1 | id), data = d, family = "binomial",
                 verbose = FALSE)
  e <- ilm_effects(f)
  expect_true(attr(e, "mixed"))
  expect_output(print(e), "CONDITIONAL, not marginal")
  expect_output(print(e), "ilm_ame")
  ## a fixed-effects fit carries no such warning
  ff <- ilm_model(y ~ x, data = d, family = "binomial", verbose = FALSE)
  expect_false(grepl("CONDITIONAL",
                     paste(utils::capture.output(print(ilm_effects(ff))),
                           collapse = " ")))
})

test_that("a robust covariance can be substituted for the model's", {
  set.seed(7); G <- 40L; ni <- 10L; n <- G * ni
  d <- data.frame(g = factor(rep(seq_len(G), each = ni)), x = rnorm(n))
  d$y <- rpois(n, exp(0.4 + 0.5 * d$x + rep(rnorm(G, 0, .7), each = ni)))
  f <- ilm_model(y ~ x, data = d, family = "poisson", verbose = FALSE)
  a <- ilm_effects(f)
  b <- ilm_effects(f, vcov = ilm_vcov_cluster(f, ~ g))
  ## the estimate does not move; only the uncertainty does
  expect_equal(a$rate_ratio, b$rate_ratio)
  expect_gt(b$upper[2] - b$lower[2], a$upper[2] - a$lower[2])
  expect_true(attr(b, "robust"))
  expect_output(print(b), "covariance supplied")
  expect_error(ilm_effects(f, vcov = matrix(1, 1, 1)), "against 2 coefficients")
})

test_that("it refuses what has no single ratio per predictor", {
  d <- eff_data(400L)
  expect_error(ilm_effects(mtcars), "must be a fitted ilm_model")
  set.seed(6); n <- 400L
  dm <- data.frame(x = rnorm(n))
  dm$y <- factor(sample(c("a", "b", "c"), n, TRUE))
  fm <- ilm_model(y ~ x, data = dm, family = "multinomial", verbose = FALSE)
  expect_error(ilm_effects(fm), "one coefficient per category dimension")
})
