# Censored responses: floors, ceilings and detection limits.
#
# A reading of "0" from an assay whose limit is 0 does not mean zero, it means
# "at or below the limit". Fitting it as a zero pulls the mean down and the
# variance in. The likelihood contribution becomes the probability of the
# interval instead of a density at a point.
#
# The fixed-effects case is pinned to survival::survreg(), which solves the
# same problem and has been solving it for thirty years.

cens_ceiling <- function(seed = 11, n = 800, up = 2.5) {
  set.seed(seed)
  d <- data.frame(x = stats::rnorm(n), z = stats::rnorm(n))
  d$ystar <- 1.0 + 1.2 * d$x - 0.7 * d$z + stats::rnorm(n, 0, 1.5)
  d$y <- pmin(d$ystar, up)
  d
}

## ---- the specification -----------------------------------------------------

test_that("codes are derived from the limits", {
  y <- c(0, 0, 1.4, 2.9, 5, 5)
  s <- ilm_censor(y, lower = 0, upper = 5)
  expect_equal(as.integer(s), c(-1L, -1L, 0L, 0L, 1L, 1L))
  expect_equal(attr(s, "lower"), 0)
  expect_equal(attr(s, "upper"), 5)
  expect_equal(as.integer(ilm_censor(y, lower = 0)), c(-1L, -1L, 0L, 0L, 0L, 0L))
  expect_equal(as.integer(ilm_censor(y, upper = 5)), c(0L, 0L, 0L, 0L, 1L, 1L))
})

test_that("a status can be given directly, and then the limits are not known", {
  s <- ilm_censor(1:5, status = c(0, 0, 1, 1, 0))
  expect_equal(as.integer(s), c(0L, 0L, 1L, 1L, 0L))
  expect_true(is.na(attr(s, "lower")))
  # with no limits the codes are held as supplied, whatever response is handed in
  expect_equal(illume:::ilm_censor_for(s, c(99, 99, 99, 99, 99)), as.integer(s))
})

test_that("codes are re-derived for a simulated response", {
  s <- ilm_censor(c(0, 1, 2, 5), lower = 0, upper = 5)
  # a replicate lands where it lands; the limits, not the original pattern,
  # decide what is censored in it
  expect_equal(illume:::ilm_censor_for(s, c(3, 0, 5, 1)), c(0L, -1L, 1L, 0L))
  expect_equal(illume:::ilm_censor_apply(s, c(-4, 3, 9)), c(0, 3, 5))
})

test_that("misspecification is named", {
  expect_error(ilm_censor(1:5), "give `lower`, `upper` or `status`")
  expect_error(ilm_censor(1:5, lower = 3, upper = 1), "must be below")
  expect_error(ilm_censor(1:5, status = c(0, 1)), "but `y` has")
  expect_error(ilm_censor(1:5, status = c(0, 1, 2, 0, 1)), "must be -1")
  expect_error(ilm_censor(1:5, status = rep(2, 5)), "1 - event", fixed = TRUE)
  expect_error(ilm_censor(letters[1:5], lower = 1), "must be numeric")
  expect_warning(ilm_censor(1:5, lower = -99), "nothing is censored")
})

test_that("printing says how much is censored", {
  out <- capture.output(print(ilm_censor(c(0, 0, 1, 2, 5), lower = 0, upper = 5)))
  expect_true(any(grepl("2 left", out)))
  expect_true(any(grepl("1 right", out)))
  expect_true(any(grepl("limits", out)))
  out2 <- capture.output(print(ilm_censor(1:4, status = c(0, 0, 1, 1))))
  expect_true(any(grepl("limits not recorded", out2)))
})

## ---- against survreg -------------------------------------------------------

test_that("a ceiling model matches survreg exactly", {
  skip_if_not_installed("survival")
  d <- cens_ceiling()
  f <- ilm_model(y ~ x + z, data = d, family = "gaussian",
                 censor = ilm_censor(d$y, upper = 2.5), verbose = FALSE)
  sv <- survival::survreg(
    survival::Surv(y, y < 2.5, type = "right") ~ x + z, data = d,
    dist = "gaussian")
  expect_equal(unname(coef(f)), unname(coef(sv)), tolerance = 1e-5)
  expect_equal(unname(f$dispersion[[1]]), sv$scale, tolerance = 1e-5)
  expect_equal(as.numeric(logLik(f)), as.numeric(logLik(sv)), tolerance = 1e-6)
  expect_equal(unname(sqrt(diag(vcov(f)))[c("x", "z")]),
               unname(sqrt(diag(vcov(sv)))[c("x", "z")]), tolerance = 1e-4)
})

test_that("a floor model matches survreg exactly", {
  skip_if_not_installed("survival")
  set.seed(12)
  n <- 800
  d <- data.frame(x = stats::rnorm(n))
  d$y <- pmax(0.4 + 1.1 * d$x + stats::rnorm(n, 0, 1.2), 0)
  f <- ilm_model(y ~ x, data = d, family = "gaussian",
                 censor = ilm_censor(d$y, lower = 0), verbose = FALSE)
  sv <- survival::survreg(survival::Surv(y, y > 0, type = "left") ~ x,
                          data = d, dist = "gaussian")
  expect_equal(unname(coef(f)), unname(coef(sv)), tolerance = 1e-5)
  expect_equal(unname(f$dispersion[[1]]), sv$scale, tolerance = 1e-5)
})

test_that("ignoring the censoring is what it looks like", {
  # the reason the feature exists: the naive fit is badly biased
  d <- cens_ceiling()
  f <- ilm_model(y ~ x + z, data = d, family = "gaussian",
                 censor = ilm_censor(d$y, upper = 2.5), verbose = FALSE)
  naive <- stats::lm(y ~ x + z, data = d)
  expect_equal(unname(coef(f)[["x"]]), 1.2, tolerance = 0.12)
  expect_lt(unname(coef(naive)[["x"]]), 1.0)     # attenuated toward zero
  expect_gt(abs(coef(naive)[["x"]] - 1.2), 3 * abs(coef(f)[["x"]] - 1.2))
})

test_that("censoring works alongside a random effect, which survreg cannot do", {
  set.seed(13)
  ng <- 60; nt <- 10
  d <- data.frame(id = factor(rep(seq_len(ng), each = nt)),
                  x = stats::rnorm(ng * nt))
  b <- stats::rnorm(ng, 0, 0.8)[as.integer(d$id)]
  d$y <- pmin(0.5 + 1.0 * d$x + b + stats::rnorm(ng * nt, 0, 1.0), 2.0)
  f <- ilm_model(y ~ x + (1 | id), data = d, family = "gaussian",
                 censor = ilm_censor(d$y, upper = 2.0), verbose = FALSE)
  g <- ilm_model(y ~ x + (1 | id), data = d, family = "gaussian",
                 verbose = FALSE)
  expect_equal(unname(coef(f)[["x"]]), 1.0, tolerance = 0.12)
  expect_equal(unname(f$dispersion[[1]]), 1.0, tolerance = 0.15)
  expect_equal(sqrt(f$Sigma[["id"]][1, 1]), 0.8, tolerance = 0.2)
  # ignoring it attenuates both the slope and the between-unit spread
  expect_lt(coef(g)[["x"]], coef(f)[["x"]])
  expect_lt(sqrt(g$Sigma[["id"]][1, 1]), sqrt(f$Sigma[["id"]][1, 1]))
})

## ---- inference and residuals ----------------------------------------------

test_that("censoring turns off the exact t and F path", {
  # exact inference rests on the response being a gaussian linear model. A
  # censored likelihood mixes densities with tail probabilities and is only
  # asymptotically normal, and the N/(N-p) correction is an OLS result.
  d <- cens_ceiling()
  f <- ilm_model(y ~ x + z, data = d, family = "gaussian",
                 censor = ilm_censor(d$y, upper = 2.5), verbose = FALSE)
  g <- ilm_model(ystar ~ x + z, data = d, family = "gaussian", verbose = FALSE)
  expect_false(isTRUE(f$exact_df))
  expect_true(isTRUE(g$exact_df))
  expect_true("z value" %in% names(ilm_coef_table(f)))
  expect_true("Chisq" %in% names(ilm_anova(f)))
})

test_that("quantile residuals are uniform for a censored fit", {
  # a censored row is known only as an interval, so its residual is drawn
  # across the probability of that interval; without that every censored row
  # lands on one quantile and the QQ plot shows a spike at the limit
  set.seed(21)
  n <- 1500
  d <- data.frame(g = factor(rep(seq_len(50), each = 30)),
                  x = stats::rnorm(n))
  b <- stats::rnorm(50, 0, 0.5)[as.integer(d$g)]
  d$y <- pmin(0.5 + 1.0 * d$x + b + stats::rnorm(n, 0, 1), 1.8)
  cs <- ilm_censor(d$y, upper = 1.8)
  f <- ilm_model(y ~ x + (1 | g), data = d, family = "gaussian", censor = cs,
                 verbose = FALSE)
  u <- ilm_rqr(f, TRUE, 1L)
  expect_gt(suppressWarnings(stats::ks.test(u, "punif")$p.value), 0.01)
  # every censored row gets its own draw, not a shared quantile
  expect_equal(length(unique(round(u[cs != 0L], 8))), sum(cs != 0L))
})

test_that("simulated replicates are censored the way the data were", {
  d <- cens_ceiling(n = 400)
  cs <- ilm_censor(d$y, upper = 2.5)
  f <- ilm_model(y ~ x + z, data = d, family = "gaussian", censor = cs,
                 verbose = FALSE)
  ys <- illume:::ilm_sim_cond(f, 25, seed = 1)
  expect_true(all(ys <= 2.5 + 1e-12))
  # the proportion censored varies between replicates, because whether a draw
  # lands beyond the limit is itself random; holding it fixed would understate
  # the spread the envelope is meant to carry
  p <- apply(ys, 2, function(v) mean(v >= 2.5))
  expect_gt(stats::sd(p), 0)
  expect_equal(mean(p), mean(d$y >= 2.5), tolerance = 0.06)
})

test_that("a family with no censored form says so", {
  set.seed(3)
  d <- data.frame(x = stats::rnorm(200))
  d$cnt <- stats::rpois(200, exp(0.4 + 0.3 * d$x))
  expect_error(
    ilm_model(cnt ~ x, data = d, family = "poisson",
              censor = ilm_censor(d$cnt, upper = max(d$cnt)), verbose = FALSE),
    "has no censored form")
  d$y <- stats::rnorm(200)
  expect_error(
    ilm_model(y ~ x, data = d, family = "gaussian",
              censor = structure(rep(0L, 5), lower = NA_real_,
                                 upper = NA_real_, class = "ilm_censor"),
              verbose = FALSE),
    "but the model matrix has 200 rows")
})

## ---- the data-side half ----------------------------------------------------

test_that("a pile-up at a limit is named before the model is fitted", {
  set.seed(1)
  y <- 0.5 + stats::rnorm(1500, 0, 1.4)
  d <- data.frame(free = y, ceiling = pmin(y, 1.2), floor = pmax(y, 0))
  r <- ilm_describe(d)
  note <- stats::setNames(r$gauss_note, r$variable)
  expect_match(note[["ceiling"]], "ceiling")
  expect_match(note[["ceiling"]], "ilm_censor", fixed = TRUE)
  expect_match(note[["floor"]], "floor")
  expect_false(grepl("ilm_censor", note[["free"]], fixed = TRUE))
})

test_that("ordinary shapes are not mistaken for limits", {
  set.seed(2)
  d <- data.frame(norm = stats::rnorm(500),
                  pois = as.numeric(stats::rpois(500, 3)),
                  lnorm = stats::rlnorm(500), unif = stats::runif(500),
                  few = as.numeric(sample(1:4, 500, TRUE)))
  r <- ilm_describe(d)
  expect_false(any(grepl("ilm_censor", r$gauss_note, fixed = TRUE)))
})

test_that("ilm_describe handles a frame, several columns, and says so otherwise", {
  # it used to hand the whole frame to the categorical branch, which failed
  # with "the condition has length > 1" for any frame of more than one column
  set.seed(4)
  d <- data.frame(a = stats::rnorm(50), b = stats::rnorm(50),
                  c = stats::rnorm(50))
  r <- ilm_describe(d)
  expect_equal(nrow(r), 3L)
  expect_equal(r$variable, c("a", "b", "c"))
  expect_equal(nrow(ilm_describe(d, c("a", "b"))), 2L)
  expect_equal(nrow(ilm_describe(d, "a")), 1L)
  expect_equal(nrow(ilm_describe(d$a)), 1L)
  expect_error(ilm_describe(data.frame(a = 1:5, s = letters[1:5])),
               "more than one kind")
  expect_error(ilm_describe(data.frame(a = 1:5, s = letters[1:5])),
               "ilm_describe_all", fixed = TRUE)
  expect_error(ilm_describe(d, "nope"), "not found in the data")
})
