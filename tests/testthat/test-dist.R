# ilm_dist(): a fit's response distribution, as functions.
#
# The load-bearing check: on a fit with nothing integrated out, the negative
# log-likelihood the optimiser minimised is the sum of the fit's own density
# over the rows. So ilm_dist()'s density, fed the fit's linear predictors and
# parameters, must give back f$opt$objective EXACTLY -- for every family, with
# censoring, weights, a zero part and a dispersion model.

## -sum(w * log d) at the fit, from the fit's own parameter vector
dist_nll <- function(f, status = NULL) {
  pe <- f$opt$par; tl <- names(pe)
  lp <- list(mu = if (f$C > 1L) f$X %*% f$beta else as.numeric(f$X %*% f$beta))
  if (!is.null(f$Zzi)) lp$zi <- as.numeric(f$Zzi %*% pe[tl == "gzi"])
  if (!is.null(f$Zd)) lp$disp <- as.numeric(f$Zd %*% pe[tl == "gamma"])
  par <- pe[tl %in% c("logdisp", "zeta_raw", "mu_pow")]
  w <- if (is.null(f$weights)) rep(1, nrow(f$X)) else f$weights
  -sum(w * ilm_dist(f)$d(f$y, lp, par, log = TRUE, status = status))
}

test_that("the density is the fit's likelihood, family by family", {
  set.seed(1); n <- 300
  d <- data.frame(x = stats::rnorm(n), g = factor(sample(c("a", "b"), n, TRUE)))
  eta <- 0.3 + 0.5 * d$x
  d$yg <- eta + stats::rnorm(n)
  d$yb <- stats::rbinom(n, 1, stats::plogis(eta))
  d$yp <- stats::rpois(n, exp(eta))
  d$yn <- stats::rnbinom(n, mu = exp(eta), size = 2)
  d$yt <- stats::rbeta(n, stats::plogis(eta) * 6, (1 - stats::plogis(eta)) * 6)
  d$yo <- factor(cut(eta + stats::rlogis(n), c(-Inf, -0.2, 0.6, 1.4, Inf),
                     labels = c("lo", "mid", "high", "top")), ordered = TRUE)
  d$ym <- factor(sample(c("u", "v", "w"), n, TRUE))
  d$wt <- sample(1:3, n, TRUE)
  cases <- list(
    gaussian = list(y ~ x + g, "gaussian"),
    binomial = list(yb ~ x, "binomial"),
    poisson = list(yp ~ x, "poisson"),
    nbinom = list(yn ~ x, "nbinom"),
    beta = list(yt ~ x, "beta"),
    ordinal = list(yo ~ x, "ordinal"),
    ordinal_probit = list(yo ~ x, "ordinal_probit"),
    multinomial = list(ym ~ x, "multinomial"))
  cases$gaussian[[1]] <- yg ~ x + g
  for (nm in names(cases)) {
    f <- ilm_model(cases[[nm]][[1]], data = d, family = cases[[nm]][[2]],
                   verbose = FALSE)
    expect_equal(dist_nll(f), f$opt$objective, tolerance = 1e-9, label = nm)
  }
  ## frequency weights
  f <- ilm_model(yp ~ x, data = d, family = "poisson", weights = wt,
                 verbose = FALSE)
  expect_equal(dist_nll(f), f$opt$objective, tolerance = 1e-9)
})

test_that("a censored row contributes its interval, as in the fit", {
  set.seed(2); n <- 400
  d <- data.frame(x = stats::rnorm(n))
  d$y <- pmax(0.5 + 0.8 * d$x + stats::rnorm(n), 0)          # a floor at zero
  f <- ilm_model(y ~ x, data = d, family = "gaussian",
                 censor = ilm_censor(d$y, lower = 0), verbose = FALSE)
  st <- illume:::ilm_censor_for(f$censor, f$y)
  expect_true(any(st < 0))
  expect_equal(dist_nll(f, status = st), f$opt$objective, tolerance = 1e-9)
  ## and right censoring in each survival family
  for (dist in c("weibull", "lognormal", "loglogistic")) {
    set.seed(3)
    s <- data.frame(x = stats::rnorm(n))
    tt <- exp(1 + 0.5 * s$x + 0.6 * switch(dist, weibull = log(stats::rexp(n)),
                                           lognormal = stats::rnorm(n),
                                           loglogistic = stats::rlogis(n)))
    ct <- stats::quantile(tt, 0.7)
    s$time <- pmin(tt, ct); s$event <- as.integer(tt <= ct)
    fs <- ilm_model(time ~ x, data = s, family = dist,
                    censor = ilm_surv(s$time, s$event), verbose = FALSE)
    st <- illume:::ilm_censor_for(fs$censor, fs$y)
    expect_equal(dist_nll(fs, status = st), fs$opt$objective, tolerance = 1e-9,
                 label = dist)
  }
})

test_that("a zero part and a dispersion model are the fit's too", {
  set.seed(4); n <- 500
  d <- data.frame(x = stats::runif(n, 0, 3), z = stats::rnorm(n),
                  g = factor(sample(c("a", "b"), n, TRUE)))
  pz <- stats::plogis(-1 + 0.8 * d$z)
  d$yz <- ifelse(stats::runif(n) < pz, 0, stats::rpois(n, exp(0.3 + 0.4 * d$x)))
  d$yh <- ifelse(stats::runif(n) < pz, 0,
                 stats::rnbinom(n, mu = exp(0.5 + 0.3 * d$x), size = 2))
  mu <- stats::plogis(-0.5 + 0.4 * d$x)
  d$yb <- ifelse(stats::runif(n) < pz, 0, stats::rbeta(n, mu * 5, (1 - mu) * 5))
  d$yd <- 1 + 0.5 * d$x + stats::rnorm(n, 0, ifelse(d$g == "b", 2, 0.5))
  fits <- list(
    zip = ilm_model(yz ~ x, data = d, family = "poisson", ziformula = ~ z,
                    verbose = FALSE),
    hnb = ilm_model(yh ~ x, data = d, family = "nbinom", ziformula = ~ z,
                    zi_type = "hurdle", verbose = FALSE),
    hbeta = ilm_model(yb ~ x, data = d, family = "beta", ziformula = ~ z,
                      zi_type = "hurdle", verbose = FALSE),
    disp = ilm_model(yd ~ x, data = d, family = "gaussian", dispformula = ~ g,
                     verbose = FALSE),
    dispmu = ilm_model(yh ~ x, data = d, family = "nbinom", dispformula = ~ mu,
                       verbose = FALSE))
  for (nm in names(fits))
    expect_equal(dist_nll(fits[[nm]]), fits[[nm]]$opt$objective,
                 tolerance = 1e-9, label = nm)
  dd <- ilm_dist(fits$zip)
  expect_identical(dd$zero, "inflated")
  expect_identical(dd$lp_names, c("mu", "zi"))
  expect_error(dd$d(1, list(mu = 0), NULL), "lp\\$zi")
})

test_that("the distribution, quantile and random functions agree", {
  set.seed(5); n <- 300
  d <- data.frame(x = stats::runif(n, 0, 3), z = stats::rnorm(n))
  d$y <- ifelse(stats::runif(n) < 0.3, 0,
                stats::rnbinom(n, mu = exp(0.5 + 0.3 * d$x), size = 2))
  f <- ilm_model(y ~ x, data = d, family = "nbinom", ziformula = ~ 1,
                 zi_type = "hurdle", verbose = FALSE)
  dd <- ilm_dist(f)
  pe <- f$opt$par
  lp <- list(mu = 1.2, zi = unname(pe[names(pe) == "gzi"]))
  par <- pe[names(pe) == "logdisp"]
  ## the density sums to one, and cumulates to the distribution function
  ys <- 0:200
  dens <- dd$d(ys, lp, par)
  expect_equal(sum(dens), 1, tolerance = 1e-8)
  expect_equal(cumsum(dens)[1:30], dd$p(0:29, lp, par), tolerance = 1e-10)
  ## the quantile inverts it
  pr <- c(0.05, 0.2, 0.5, 0.8, 0.97)
  qv <- dd$q(pr, lp, par)
  expect_true(all(dd$p(qv, lp, par) >= pr - 1e-12))
  expect_true(all(qv == 0 | dd$p(qv - 1, lp, par) < pr))
  ## the mean is the density's
  expect_equal(dd$mean(lp, par), sum(ys * dens), tolerance = 1e-8)
  ## and draws follow it
  set.seed(6)
  x <- dd$r(list(mu = rep(1.2, 40000), zi = lp$zi), par)
  expect_equal(mean(x == 0), dens[1], tolerance = 0.05)
  expect_equal(mean(x), dd$mean(lp, par), tolerance = 0.05)
  expect_true(dd$discrete)
})

test_that("an ordered response gives categories and their probabilities", {
  set.seed(7); n <- 400
  d <- data.frame(x = stats::rnorm(n))
  d$y <- factor(cut(0.8 * d$x + stats::rlogis(n), c(-Inf, -0.5, 0.5, Inf),
                    labels = c("low", "mid", "high")), ordered = TRUE)
  f <- ilm_model(y ~ x, data = d, family = "ordinal", verbose = FALSE)
  dd <- ilm_dist(f)
  pe <- f$opt$par
  par <- pe[names(pe) == "zeta_raw"]
  lp <- list(mu = as.numeric(f$X %*% f$beta))
  P <- dd$mean(lp, par)
  expect_equal(unname(P), unname(predict(f)), tolerance = 1e-10)
  expect_equal(dd$d(d$y, lp, par), P[cbind(seq_len(n), as.integer(d$y))])
  expect_equal(dd$p(2, lp, par), P[, 1] + P[, 2])
  expect_identical(dd$categories, c("low", "mid", "high"))
  expect_identical(dd$q(c(0, 1), list(mu = 0), par), c(1L, 3L))
  ## multinomial: probabilities, and no ordering
  d$m <- factor(sample(c("a", "b", "c"), n, TRUE))
  fm <- ilm_model(m ~ x, data = d, family = "multinomial", verbose = FALSE)
  dm <- ilm_dist(fm)
  expect_null(dm$p)
  expect_equal(unname(dm$mean(list(mu = fm$X %*% fm$beta), NULL)),
               unname(predict(fm)), tolerance = 1e-10)
  expect_output(print(dm), "multinomial")
})

test_that("a survival family's density, quantile and mean are exact", {
  set.seed(8); n <- 300
  s <- data.frame(x = stats::rnorm(n))
  s$time <- exp(1 + 0.5 * s$x + 0.6 * log(stats::rexp(n)))
  f <- ilm_model(time ~ x, data = s, family = "weibull", verbose = FALSE)
  dd <- ilm_dist(f)
  par <- f$opt$par[names(f$opt$par) == "logdisp"]
  sc <- exp(unname(par)); lp <- list(mu = 1.3)
  ## the mean against a numerical integral of the survivor function
  surv <- function(t) dd$p(t, lp, par, lower.tail = FALSE)
  expect_equal(dd$mean(lp, par), stats::integrate(surv, 0, Inf)$value,
               tolerance = 1e-6)
  expect_equal(dd$q(dd$p(2.5, lp, par), lp, par), 2.5, tolerance = 1e-10)
  expect_false(dd$discrete)
  expect_error(ilm_dist(mtcars), "must be a fitted ilm_model")
})
