## ilm_survival() on a model with a smooth beside a random effect. Its curves
## for each row's own group ("fitted") and for the average over the groups
## ("population"), and a flexible-parametric model's typical curve, were built
## from the fixed design alone, so a smooth's penalised part was left out: on
## a Weibull fit the survival at t = 3 came out 8.8e-21 where the model gave
## 0.917. And "fitted" and "population" stopped when given a single time.
## Each curve here is checked against the survivor function built by hand
## from the fit's own parts.

ss_data <- function(seed = 31, n = 300) {
  set.seed(seed)
  d <- data.frame(x = stats::runif(n), h = factor(sample(letters[1:8], n, TRUE)))
  hl <- stats::rnorm(8, 0, 0.4)
  eta <- 1 + 1.2 * sin(2 * pi * d$x) + hl[as.integer(d$h)]
  tt <- exp(eta + 0.5 * log(stats::rexp(n)))
  ct <- stats::rexp(n, 1 / (2 * stats::median(tt)))
  d$time <- pmin(tt, ct); d$event <- as.integer(tt <= ct)
  d
}
ss_fit <- function(d, family, form = time ~ s(x, k = 8) + (1 | h))
  suppressMessages(suppressWarnings(
    ilm_model(form, data = d, family = family, censor = ilm_surv(d$time, d$event),
              reml = TRUE, verbose = FALSE)))

## S(t) at the rows of `nd` from the fit's parts: the fixed effects, every
## smooth's penalised part, and for "fitted" the row's level of h; for
## "population" averaged over h's normal by integrate()
ss_by_hand <- function(f, nd, times, groups) {
  ndx <- illume:::ilm_newX(f, nd)
  keep <- if (is.null(f$rp)) seq_len(ncol(f$X)) else setdiff(seq_len(ncol(f$X)), f$rp$cols)
  e <- as.vector(ndx$X %*% f$beta[keep, 1L])
  for (k in which(vapply(f$re, function(r) identical(r$kind, "basis"), TRUE)))
    e <- e + as.vector(illume:::ilm_basis_new(f, ndx, k) %*% illume:::ilm_Bhat_term(f, k))
  kh <- which(names(f$re) == "h")
  if (groups == "fitted")
    e <- e + as.vector(illume:::ilm_Bhat_term(f, kh))[match(as.character(nd$h), f$re[[kh]]$levels)]
  S <- if (is.null(f$rp)) {
    sc <- unname(f$dispersion[[1]])
    function(t, e) exp(-exp((log(t) - e) / sc))
  } else {
    function(t, e) f$family$surv(as.vector(illume:::ilm_rcs(log(t), f$rp$knots) %*%
                                             f$beta[f$rp$cols, 1L]) + e)
  }
  v <- f$Sigma[["h"]][1, 1]
  ## row by row, times within each row, as ilm_survival() returns them
  unlist(lapply(e, function(ei) vapply(times, function(t) {
    if (groups != "population") return(S(t, ei))
    stats::integrate(function(z) S(t, ei + sqrt(v) * z) * stats::dnorm(z),
                     -Inf, Inf, rel.tol = 1e-10)$value
  }, 0)))
}

test_that("a Weibull curve carries the smooth under 'fitted' and 'population'", {
  d <- ss_data(); f <- ss_fit(d, "weibull")
  nd <- data.frame(x = c(0.25, 0.5, 0.75), h = factor(c("a", "b", "a"), levels = levels(d$h)))
  tg <- c(1, 3)
  for (gr in c("typical", "fitted", "population")) {
    s <- suppressWarnings(ilm_survival(f, newdata = nd, times = tg, groups = gr, nsim = 20L))
    expect_equal(s$surv, ss_by_hand(f, nd, tg, gr), tolerance = 1e-8, label = gr)
  }
})

test_that("a flexible-parametric curve carries the smooth under every groups setting", {
  d <- ss_data(32); f <- ss_fit(d, "rp")
  nd <- data.frame(x = c(0.25, 0.5, 0.75), h = factor(c("a", "b", "a"), levels = levels(d$h)))
  tg <- c(1, 3)
  for (gr in c("typical", "fitted", "population")) {
    s <- suppressWarnings(ilm_survival(f, newdata = nd, times = tg, groups = gr, nsim = 20L))
    expect_equal(s$surv, ss_by_hand(f, nd, tg, gr), tolerance = 1e-8, label = gr)
  }
})

test_that("'fitted' and 'population' take a single time", {
  d <- ss_data(33)
  for (fam in c("weibull", "rp")) {
    f <- ss_fit(d, fam)
    nd <- data.frame(x = c(0.3, 0.6), h = factor(c("a", "c"), levels = levels(d$h)))
    for (gr in c("fitted", "population")) {
      one <- suppressWarnings(ilm_survival(f, newdata = nd, times = 2, groups = gr, nsim = 20L))
      two <- suppressWarnings(ilm_survival(f, newdata = nd, times = c(2, 4), groups = gr, nsim = 20L))
      expect_equal(nrow(one), 2L)
      expect_equal(one$surv, two$surv[two$time == 2], tolerance = 1e-12,
                   label = paste(fam, gr))
      expect_true(all(one$lower <= one$surv & one$surv <= one$upper))
    }
  }
})

test_that("without a smooth the curves are what they were", {
  ## the same by-hand construction; there is no penalised part to add
  d <- ss_data(34); f <- ss_fit(d, "weibull", time ~ x + (1 | h))
  nd <- data.frame(x = c(0.2, 0.8), h = factor("d", levels = levels(d$h)))
  for (gr in c("typical", "fitted", "population")) {
    s <- suppressWarnings(ilm_survival(f, newdata = nd, times = c(1, 3), groups = gr, nsim = 20L))
    expect_equal(s$surv, ss_by_hand(f, nd, c(1, 3), gr), tolerance = 1e-8, label = gr)
  }
})

test_that("a flexible-parametric typical curve leaves a factor-smooth's subjects out", {
  ## its subjects are a population, so a typical curve has the term at zero,
  ## as predict(groups = "typical") has it, while an ordinary smooth stays in
  d <- ss_data(35)
  f <- ss_fit(d, "rp", time ~ s(x, k = 5) + s(x, h, bs = "fs", k = 4))
  at <- function(lv) data.frame(x = c(0.25, 0.5, 0.75), h = factor(lv, levels = levels(d$h)))
  s_a <- suppressWarnings(ilm_survival(f, newdata = at("a"), times = c(1, 3), groups = "typical"))
  s_e <- suppressWarnings(ilm_survival(f, newdata = at("e"), times = c(1, 3), groups = "typical"))
  expect_equal(s_a$surv, s_e$surv)
  ## by hand: the fixed effects and s(x)'s penalised part, no fs blocks
  nd <- at("a")
  ndx <- illume:::ilm_newX(f, nd, fs_fill = TRUE)
  keep <- setdiff(seq_len(ncol(f$X)), f$rp$cols)
  e <- as.vector(ndx$X %*% f$beta[keep, 1L])
  fs <- illume:::ilm_fs_terms(f)
  expect_true(any(fs))
  for (k in which(vapply(f$re, function(r) identical(r$kind, "basis"), TRUE) & !fs))
    e <- e + as.vector(illume:::ilm_basis_new(f, ndx, k) %*% illume:::ilm_Bhat_term(f, k))
  S <- unlist(lapply(e, function(ei) f$family$surv(
    as.vector(illume:::ilm_rcs(log(c(1, 3)), f$rp$knots) %*% f$beta[f$rp$cols, 1L]) + ei)))
  expect_equal(s_a$surv, S, tolerance = 1e-10)
})

test_that("a Weibull 'fitted' curve agrees with predict(), through exported functions only", {
  ## predict(groups = "fitted", type = "link") carries the smooth and the
  ## row's group; the Weibull survivor function at that linear predictor is
  ## the curve, with nothing from the package's internals
  d <- ss_data(36); f <- ss_fit(d, "weibull")
  nd <- data.frame(x = c(0.2, 0.55, 0.85), h = factor(c("b", "c", "f"), levels = levels(d$h)))
  tg <- c(1, 2.5)
  s <- suppressWarnings(ilm_survival(f, newdata = nd, times = tg, groups = "fitted", nsim = 20L))
  eta <- as.numeric(suppressWarnings(predict(f, newdata = nd, type = "link", groups = "fitted")))
  sc <- unname(f$dispersion[[1]])
  expect_equal(s$surv, as.vector(sapply(eta, function(e) exp(-exp((log(tg) - e) / sc)))),
               tolerance = 1e-8)
})

test_that("the 'fitted' interval is the percentiles over joint draws of the coefficients, the smooth's included", {
  ## ilm_survival() draws with ilm_draws(blocks = c("beta", "bvec")); the same
  ## draws, put through the by-hand curve, give the same interval exactly. A
  ## smooth's penalised coefficients held at their estimates instead give a
  ## different, wider one at the ends (measured 14 Nov)
  d <- ss_data(37); f <- ss_fit(d, "weibull")
  nd <- data.frame(x = c(0.25, 0.75), h = factor("a", levels = levels(d$h)))
  tg <- 3
  s <- suppressWarnings(ilm_survival(f, newdata = nd, times = tg, groups = "fitted",
                                     nsim = 60L, seed = 4L))
  dr <- ilm_draws(f, nsim = 60L, seed = 4L, blocks = c("beta", "bvec"), natural = FALSE)
  w <- rownames(dr$draws)
  ndx <- illume:::ilm_newX(f, nd)
  kh <- which(names(f$re) == "h")
  ks <- which(vapply(f$re, function(r) identical(r$kind, "basis"), TRUE))
  sc <- unname(f$dispersion[[1]])
  S <- sapply(seq_len(ncol(dr$draws)), function(j) {
    b <- dr$draws[w == "beta", j]; bv <- dr$draws[w == "bvec", j]
    e <- as.vector(ndx$X %*% b)
    for (k in ks) e <- e + as.vector(illume:::ilm_basis_new(f, ndx, k) %*% illume:::ilm_Bhat_term(f, k, bv))
    e <- e + as.vector(illume:::ilm_Bhat_term(f, kh, bv))[match(as.character(nd$h), f$re[[kh]]$levels)]
    exp(-exp((log(tg) - e) / sc))
  })
  expect_equal(s$lower, apply(S, 1, stats::quantile, 0.025, names = FALSE), tolerance = 1e-10)
  expect_equal(s$upper, apply(S, 1, stats::quantile, 0.975, names = FALSE), tolerance = 1e-10)
  expect_true(all(s$lower <= s$surv & s$surv <= s$upper))
})
