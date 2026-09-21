zi_data <- function(n = 600L, seed = 1L, hurdle = FALSE) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n), z = rnorm(n))
  mu <- exp(0.7 + 0.5 * d$x)
  pz <- stats::plogis(-0.4 + 0.9 * d$z)
  yc <- if (hurdle) stats::qpois(exp(-mu) + runif(n) * (1 - exp(-mu)), mu)
        else stats::rpois(n, mu)
  d$y <- ifelse(runif(n) < pz, 0, yc)
  d
}
zi_fit <- function(d, type = "inflated", zf = ~ z, ...)
  ilm_model(y ~ x, data = d, family = "poisson", ziformula = zf,
            zi_type = type, verbose = FALSE, ...)

test_that("a zero part is fitted and recovers what generated it", {
  f <- zi_fit(zi_data(1500L))
  expect_s3_class(f, "ilm_model")
  expect_true(f$ok)
  expect_equal(unname(coef(f)), c(0.7, 0.5), tolerance = 0.15)
  zc <- ilm_zi_coef(f)
  expect_equal(nrow(zc), 2L)
  expect_setequal(zc$term, c("(Intercept)", "z"))
  expect_equal(zc$estimate, c(-0.4, 0.9), tolerance = 0.25)
  expect_true(all(zc$se > 0))
  expect_equal(zc$odds_ratio, exp(zc$estimate))
  ## the zero part shows up in the coefficient names the fit carries
  expect_true(any(grepl("^zi:", f$pnames)))
})

test_that("the two types are different models, and each wins on its own data", {
  di <- zi_data(2000L, seed = 5L, hurdle = FALSE)
  dh <- zi_data(2000L, seed = 5L, hurdle = TRUE)
  expect_lt(AIC(zi_fit(di, "inflated")), AIC(zi_fit(di, "hurdle")))
  expect_lt(AIC(zi_fit(dh, "hurdle")), AIC(zi_fit(dh, "inflated")))
  ## under a hurdle the probability is EVERY zero, so with a constant zero
  ## part it is exactly the observed proportion
  fh <- zi_fit(dh, "hurdle", zf = ~ 1)
  expect_equal(mean(ilm_zi_prob(fh)), mean(dh$y == 0), tolerance = 1e-4)
  ## under a mixture it is the structural share only, which is smaller
  fi <- zi_fit(dh, "inflated", zf = ~ 1)
  expect_lt(mean(ilm_zi_prob(fi)), mean(dh$y == 0))
})

test_that("the response scale carries the zeros", {
  for (tp in c("inflated", "hurdle")) {
    d <- zi_data(3000L, seed = 9L, hurdle = tp == "hurdle")
    f <- zi_fit(d, tp)
    pr <- as.numeric(predict(f))
    expect_equal(mean(pr), mean(d$y), tolerance = 0.06)
    ## and it is NOT the inverse link of the linear predictor, which is the
    ## count part on its own
    cnt <- exp(as.numeric(predict(f, type = "link")))
    expect_gt(mean(cnt), mean(pr) * 1.2)
    expect_true(all(pr >= 0))
    ## newdata goes through the zero part's own design
    nd <- data.frame(x = c(-1, 0, 1), z = c(-2, 0, 2))
    p2 <- as.numeric(predict(f, newdata = nd))
    expect_length(p2, 3L)
    ## more z means more zeros means a smaller mean, holding x fixed
    nd2 <- data.frame(x = 0, z = c(-2, 2))
    pp <- as.numeric(predict(f, newdata = nd2))
    expect_gt(pp[1], pp[2])
    expect_equal(length(ilm_zi_prob(f, newdata = nd)), 3L)
  }
})

test_that("quantile residuals stay uniform, which is what would fail silently", {
  for (tp in c("inflated", "hurdle")) {
    f <- zi_fit(zi_data(1200L, seed = 4L, hurdle = tp == "hurdle"), tp)
    u <- ilm_rqr(f, seed = 2L)
    expect_true(all(u >= 0 & u <= 1))
    expect_gt(suppressWarnings(stats::ks.test(u, "punif")$p.value), 0.01)
    ## scoring the same fit against the count part alone is what the zero
    ## branch exists to prevent, and it is not a close call
    mu <- as.numeric(f$family$linkinv(ilm_eta_hat(f, TRUE)[, 1]))
    y <- as.numeric(f$y); set.seed(2)
    lo <- stats::ppois(y - 1, mu); hi <- stats::ppois(y, mu)
    bad <- lo + runif(length(y)) * (hi - lo)
    expect_lt(suppressWarnings(stats::ks.test(bad, "punif")$p.value), 1e-6)
  }
})

test_that("simulation reproduces the zero fraction, and the checks agree", {
  for (tp in c("inflated", "hurdle")) {
    d <- zi_data(1500L, seed = 6L, hurdle = tp == "hurdle")
    f <- zi_fit(d, tp)
    ys <- ilm_simulate(f, 100L, seed = 3L)
    expect_equal(mean(ys == 0), mean(d$y == 0), tolerance = 0.03)
    expect_equal(mean(ys), mean(d$y), tolerance = 0.12)
    ## a hurdle's positives come from a count that cannot be zero, so the
    ## draws must respect that rather than being filtered after the fact
    z <- suppressMessages(ilm_check_zeros(f, B = 200L))
    expect_equal(z$status, "OK")
    expect_equal(z$expected, z$observed, tolerance = 0.06 * z$observed)
  }
  ## and the check still fires when the zero part is missing
  d <- zi_data(1500L, seed = 6L)
  f0 <- ilm_model(y ~ x, data = d, family = "poisson", verbose = FALSE)
  z0 <- suppressMessages(ilm_check_zeros(f0, B = 200L))
  expect_equal(z0$status, "FAIL")
  expect_match(z0$note, "ziformula")
})

test_that("a zero part works with a random effect and with nbinom", {
  set.seed(21); ng <- 50L; ni <- 12L; n <- ng * ni
  d <- data.frame(id = factor(rep(seq_len(ng), each = ni)),
                  x = rnorm(n), z = rnorm(n))
  b <- rnorm(ng, 0, 0.5)
  mu <- exp(0.5 + 0.5 * d$x + b[as.integer(d$id)])
  d$y <- ifelse(runif(n) < stats::plogis(-0.5 + 0.8 * d$z), 0, rpois(n, mu))
  fr <- ilm_model(y ~ x + (1 | id), data = d, family = "poisson",
                  ziformula = ~ z, verbose = FALSE)
  expect_true(fr$ok)
  expect_equal(unname(coef(fr)), c(0.5, 0.5), tolerance = 0.2)
  expect_equal(sqrt(fr$Sigma$id[1, 1]), 0.5, tolerance = 0.2)
  ## an absolute bound, not a relative one: across 150 replicates at this
  ## size the zero-part coefficients have a standard deviation of 0.13 and a
  ## bias under 0.03, so a single draw lands well outside a 30% band on -0.5
  ## without anything being wrong
  expect_lt(max(abs(ilm_zi_coef(fr)$estimate - c(-0.5, 0.8))), 0.45)

  set.seed(13); n <- 1200L
  dn <- data.frame(x = rnorm(n), z = rnorm(n))
  dn$y <- ifelse(runif(n) < stats::plogis(-0.4 + 0.9 * dn$z), 0,
                 rnbinom(n, size = 2, mu = exp(0.7 + 0.5 * dn$x)))
  fn <- ilm_model(y ~ x, data = dn, family = "nbinom", ziformula = ~ z,
                  verbose = FALSE)
  expect_true(fn$ok)
  expect_equal(unname(coef(fn)), c(0.7, 0.5), tolerance = 0.2)
  expect_equal(unname(fn$dispersion), 2, tolerance = 0.8)
})

test_that("a zero part is refused where a zero is not a special value", {
  set.seed(2); d <- data.frame(x = rnorm(200), z = rnorm(200))
  d$y <- rnorm(200)
  expect_error(
    ilm_model(y ~ x, data = d, family = "gaussian", ziformula = ~ z,
              verbose = FALSE),
    "applies to counts")
  d$y <- rpois(200, 2)
  expect_error(ilm_model(y ~ x, data = d, family = "poisson",
                         ziformula = ~ nope, verbose = FALSE),
               "not in the data")
  expect_error(ilm_model(y ~ x, data = d, family = "poisson",
                         ziformula = "z", verbose = FALSE),
               "one-sided formula")
  ## a bar is refused rather than silently dropped
  d$g <- factor(rep(1:20, 10))
  expect_error(ilm_model(y ~ x, data = d, family = "poisson",
                         ziformula = ~ (1 | g), verbose = FALSE),
               "fixed effects only")
  ## and the accessors say what to do when there is no zero part
  f <- ilm_model(y ~ x, data = d, family = "poisson", verbose = FALSE)
  expect_error(ilm_zi_coef(f), "no zero part")
  expect_error(ilm_zi_prob(f), "no zero part")
})

test_that("the fit matches pscl, which is the reference implementation", {
  skip_if_not_installed("pscl")
  d <- zi_data(800L, seed = 11L)
  f <- zi_fit(d, "inflated")
  m <- pscl::zeroinfl(y ~ x | z, data = d, dist = "poisson")
  expect_equal(unname(coef(f)), unname(coef(m)[1:2]), tolerance = 1e-4)
  expect_equal(unname(ilm_zi_coef(f)$estimate), unname(coef(m)[3:4]),
               tolerance = 1e-4)
  expect_equal(unname(ilm_zi_coef(f)$se),
               unname(sqrt(diag(vcov(m)))[3:4]), tolerance = 1e-4)
  expect_equal(as.numeric(logLik(f)), as.numeric(logLik(m)), tolerance = 1e-5)

  dh <- zi_data(800L, seed = 12L, hurdle = TRUE)
  fh <- zi_fit(dh, "hurdle")
  mh <- pscl::hurdle(y ~ x | z, data = dh, dist = "poisson",
                     zero.dist = "binomial")
  expect_equal(unname(coef(fh)), unname(coef(mh)[1:2]), tolerance = 1e-4)
  ## pscl models P(y > 0) where illume models P(y = 0), so the zero part is
  ## the same numbers with the opposite sign
  expect_equal(unname(ilm_zi_coef(fh)$estimate), -unname(coef(mh)[3:4]),
               tolerance = 1e-4)
  expect_equal(as.numeric(logLik(fh)), as.numeric(logLik(mh)), tolerance = 1e-5)

  dn <- zi_data(800L, seed = 13L)
  dn$y <- ifelse(dn$y == 0, 0, rnbinom(sum(dn$y != 0), size = 2, mu = 3))
  fn <- ilm_model(y ~ x, data = dn, family = "nbinom", ziformula = ~ z,
                  verbose = FALSE)
  mn <- pscl::zeroinfl(y ~ x | z, data = dn, dist = "negbin")
  expect_equal(as.numeric(logLik(fn)), as.numeric(logLik(mn)), tolerance = 1e-4)
})
