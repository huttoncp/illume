## ilm_survival() follows predict()'s groups (Craig's item 213): each row's own
## group by default, a typical group, or the average of the curve over the
## groups.

q <- function(e) suppressMessages(suppressWarnings(e))

surv_data <- function(slope = FALSE, n_g = 15L, per = 12L, seed = 11) {
  set.seed(seed)
  d <- data.frame(g = factor(sprintf("g%02d", rep(seq_len(n_g), each = per))),
                  x = stats::rnorm(n_g * per))
  b <- stats::rnorm(n_g, 0, 0.6); s <- if (slope) stats::rnorm(n_g, 0, 0.3) else rep(0, n_g)
  eta <- 1 + 0.5 * d$x + b[d$g] + s[d$g] * d$x
  tt <- exp(eta + 0.7 * log(stats::rexp(nrow(d))))
  ct <- stats::rexp(nrow(d), rate = 1 / (3 * stats::median(tt)))
  d$time <- pmin(tt, ct); d$event <- as.integer(tt <= ct)
  d
}
weibull_fit <- function(d, f = time ~ x + (1 | g))
  q(ilm_model(f, data = d, family = "weibull", censor = ilm_surv(d$time, d$event),
              verbose = FALSE))
S_weib <- function(times, eta, sc) exp(-exp((log(times) - eta) / sc))
mode_of <- function(f, grp, dim = "(Intercept)") {
  r <- ilm_ranef(f)
  r$mode[r$type == "re" & as.character(r$level) == grp & r$dim == dim]
}

test_that("fitted is each row's own group's curve, built by hand from ranef()", {
  d <- surv_data(); f <- weibull_fit(d)
  nd <- data.frame(x = c(0, 0, 1), g = factor(c("g01", "g02", "g01"), levels = levels(d$g)))
  tg <- c(1, 3, 10)
  s <- ilm_survival(f, newdata = nd, times = tg, nsim = 50)
  sc <- unname(f$dispersion[[1]]); b <- coef(f)
  for (i in 1:3) {
    eta <- b[["(Intercept)"]] + b[["x"]] * nd$x[i] + mode_of(f, as.character(nd$g[i]))
    expect_equal(s$surv[s$row == i], S_weib(tg, eta, sc), tolerance = 1e-8)
  }
  ## rows that differ only in their group get different curves
  expect_false(isTRUE(all.equal(s$surv[s$row == 1], s$surv[s$row == 2])))
  ## and the intervals bracket the curve
  expect_true(all(s$lower <= s$surv + 1e-12 & s$surv <= s$upper + 1e-12))
})

test_that("typical is today's curve, with every random effect at zero", {
  d <- surv_data(); f <- weibull_fit(d)
  nd <- data.frame(x = c(-1, 1), g = factor(c("g01", "g02"), levels = levels(d$g)))
  tg <- c(1, 3, 10)
  s <- ilm_survival(f, newdata = nd, times = tg, groups = "typical")
  sc <- unname(f$dispersion[[1]]); b <- coef(f)
  for (i in 1:2)
    expect_equal(s$surv[s$row == i], S_weib(tg, b[["(Intercept)"]] + b[["x"]] * nd$x[i], sc),
                 tolerance = 1e-8)
  ## with no newdata the curve is a typical group's, and fitted is refused
  expect_equal(ilm_survival(f, times = tg)$surv,
               ilm_survival(f, times = tg, groups = "typical")$surv)
  expect_error(ilm_survival(f, times = tg, groups = "fitted"), "belongs to no group")
})

test_that("population averages the curve over the groups, as a Monte Carlo average does", {
  d <- surv_data(); f <- weibull_fit(d)
  nd <- data.frame(x = 0.5, g = factor("g01", levels = levels(d$g)))
  tg <- c(0.5, 2, 6, 15)
  s <- ilm_survival(f, newdata = nd, times = tg, groups = "population", nsim = 50)
  sc <- unname(f$dispersion[[1]]); b <- coef(f)
  sd_g <- sqrt(ilm_varcorr(f)$re$g[1, 1])
  set.seed(9); bb <- stats::rnorm(2e5, 0, sd_g)
  eta0 <- b[["(Intercept)"]] + 0.5 * b[["x"]]
  mc <- vapply(tg, function(t) mean(S_weib(t, eta0 + bb, sc)), 0)
  expect_equal(s$surv, mc, tolerance = 3e-3)
  ## the average of the curve, not the curve at zero, and not the curve at the
  ## log of the averaged mean -- the AFT path's old error
  typ <- S_weib(tg, eta0, sc)
  naive_eta <- log(mean(exp(eta0 + bb)))
  expect_gt(max(abs(s$surv - typ)), 0.01)
  expect_gt(max(abs(s$surv - S_weib(tg, naive_eta, sc))), 0.01)
})

test_that("a random slope enters the fitted curve and the population's variance", {
  d <- surv_data(slope = TRUE); f <- weibull_fit(d, time ~ x + (1 + x | g))
  nd <- data.frame(x = 2, g = factor("g03", levels = levels(d$g)))
  tg <- c(1, 5)
  s <- ilm_survival(f, newdata = nd, times = tg, nsim = 30)
  sc <- unname(f$dispersion[[1]]); b <- coef(f)
  eta <- b[["(Intercept)"]] + 2 * b[["x"]] + mode_of(f, "g03") + 2 * mode_of(f, "g03", "x")
  expect_equal(s$surv, S_weib(tg, eta, sc), tolerance = 1e-8)
  ## the population average at x = 2 integrates the intercept and the slope
  sp <- ilm_survival(f, newdata = nd, times = tg, groups = "population", nsim = 30)
  S <- as.matrix(ilm_varcorr(f)$re$g); z <- c(1, 2)
  v <- as.numeric(t(z) %*% S %*% z)
  set.seed(9); bb <- stats::rnorm(2e5, 0, sqrt(v))
  eta0 <- b[["(Intercept)"]] + 2 * b[["x"]]
  expect_equal(sp$surv, vapply(tg, function(t) mean(S_weib(t, eta0 + bb, sc)), 0),
               tolerance = 3e-3)
})

test_that("a flexible (Royston-Parmar) model follows the same groups", {
  d <- surv_data()
  f <- q(ilm_model(time ~ x + (1 | g), data = d, family = "rp",
                   censor = ilm_surv(d$time, d$event), verbose = FALSE))
  nd <- data.frame(x = c(0, 0), g = factor(c("g01", "g02"), levels = levels(d$g)))
  tg <- c(1, 4)
  sf <- ilm_survival(f, newdata = nd, times = tg, nsim = 30)
  st <- ilm_survival(f, newdata = nd, times = tg, groups = "typical")
  ## the typical curve is the same for both rows; each row's own is not
  expect_equal(st$surv[st$row == 1], st$surv[st$row == 2])
  expect_false(isTRUE(all.equal(sf$surv[sf$row == 1], sf$surv[sf$row == 2])))
  ## a row's own curve is the typical one moved by its group's effect
  rp <- f$rp; keep <- setdiff(seq_len(ncol(f$X)), rp$cols)
  base <- as.vector(ilm_rcs(log(tg), rp$knots) %*% f$beta[rp$cols, 1L])
  e1 <- sum(c(1, 0) * f$beta[keep, 1L]) + mode_of(f, "g01")
  expect_equal(sf$surv[sf$row == 1], f$family$surv(base + e1), tolerance = 1e-8)
  ## and the population's is the Monte Carlo average of the curve
  sp <- ilm_survival(f, newdata = nd[1, ], times = tg, groups = "population", nsim = 30)
  sd_g <- sqrt(ilm_varcorr(f)$re$g[1, 1])
  set.seed(9); bb <- stats::rnorm(2e5, 0, sd_g)
  e0 <- sum(c(1, 0) * f$beta[keep, 1L])
  mc <- vapply(seq_along(tg), function(k) mean(f$family$surv(base[k] + e0 + bb)), 0)
  expect_equal(sp$surv, mc, tolerance = 3e-3)
})

test_that("an unseen group or missing grouping column is refused, naming the alternatives", {
  d <- surv_data(); f <- weibull_fit(d)
  nd <- data.frame(x = 0, g = factor("new"))
  expect_error(ilm_survival(f, newdata = nd, times = 1), "typical")
  expect_error(ilm_survival(f, newdata = data.frame(x = 0), times = 1), "typical")
  ## either alternative answers for them
  expect_silent(ilm_survival(f, newdata = data.frame(x = 0), times = 1, groups = "typical"))
})

test_that("a model without random effects has one curve whatever groups says", {
  d <- surv_data()
  f <- q(ilm_model(time ~ x, data = d, family = "weibull", censor = ilm_surv(d$time, d$event),
                   verbose = FALSE))
  nd <- data.frame(x = c(-1, 1))
  expect_equal(ilm_survival(f, newdata = nd, times = c(1, 3)),
               ilm_survival(f, newdata = nd, times = c(1, 3), groups = "typical"))
})
