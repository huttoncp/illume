med_data <- function(n = 3000L, seed = 1L, theta = 0, lam = 0) {
  set.seed(seed)
  u <- rnorm(n)
  d <- data.frame(x = rbinom(n, 1L, .5), c = rnorm(n), u = u)
  d$m <- 0.3 + 0.7 * d$x + 0.2 * d$c + lam * u + rnorm(n)
  d$y <- 1 + 0.4 * d$x + 0.6 * d$m + theta * d$x * d$m + 0.1 * d$c +
    lam * u + rnorm(n)
  d
}
med_fit <- function(d, interaction = FALSE) {
  fm <- ilm_model(m ~ x + c, data = d, family = "gaussian", verbose = FALSE)
  fy <- if (interaction)
    ilm_model(y ~ x * m + c, data = d, family = "gaussian", verbose = FALSE)
  else ilm_model(y ~ x + m + c, data = d, family = "gaussian", verbose = FALSE)
  list(m = fm, y = fy)
}

test_that("with a linear outcome and no interaction, ACME is exactly a*b", {
  ## the Baron-Kenny product is the special case the counterfactual definition
  ## has to reproduce, and reproducing it is the test that the counterfactual
  ## machinery is assembled correctly
  d <- med_data()
  f <- med_fit(d)
  md <- ilm_mediate(f$m, f$y, "x", "m", sims = 4000L, progress = FALSE)
  expect_s3_class(md, "ilm_mediate")
  a <- unname(stats::coef(f$m)["x"])
  b <- unname(stats::coef(f$y)["m"])
  cp <- unname(stats::coef(f$y)["x"])
  get <- function(e) md$estimate[md$effect == e]
  expect_equal(get("ACME (control)"), a * b, tolerance = 0.01)
  expect_equal(get("ACME (treated)"), a * b, tolerance = 0.01)
  expect_equal(get("ADE (control)"), cp, tolerance = 0.02)
  expect_equal(get("Total effect"), a * b + cp, tolerance = 0.02)
  ## with no interaction the two arms agree by construction, and it says so
  expect_equal(get("ACME (control)"), get("ACME (treated)"))
  expect_output(print(md), "no x:m interaction")
  ## the decomposition adds up, whichever pairing is taken
  expect_equal(get("ACME (control)") + get("ADE (treated)"),
               get("Total effect"), tolerance = 1e-8)
  expect_true(all(md$lower < md$estimate & md$estimate < md$upper))
})

test_that("with an interaction the two ACMEs differ, and match the closed form", {
  ## a single product cannot be two numbers, which is why the product is not
  ## the general answer
  set.seed(2); n <- 4000L
  d <- data.frame(x = rbinom(n, 1L, .5))
  d$m <- 0.3 + 0.7 * d$x + rnorm(n)
  d$y <- 1 + 0.4 * d$x + 0.6 * d$m + 0.5 * d$x * d$m + rnorm(n)
  fm <- ilm_model(m ~ x, data = d, family = "gaussian", verbose = FALSE)
  fy <- ilm_model(y ~ x * m, data = d, family = "gaussian", verbose = FALSE)
  md <- ilm_mediate(fm, fy, "x", "m", sims = 4000L, progress = FALSE)
  a <- unname(stats::coef(fm)["x"]); b <- unname(stats::coef(fy)["m"])
  th <- unname(stats::coef(fy)["x:m"])
  ## ACME(t) = a * (b + theta * t)
  expect_equal(md$estimate[1], a * b, tolerance = 0.01)
  expect_equal(md$estimate[2], a * (b + th), tolerance = 0.01)
  expect_gt(md$estimate[2] / md$estimate[1], 1.5)
  expect_true(attr(md, "interaction"))
  expect_false(grepl("no x:m interaction",
                     paste(utils::capture.output(print(md)), collapse = " ")))
})

test_that("a binary outcome gets an effect on the probability scale", {
  set.seed(3); n <- 4000L
  d <- data.frame(x = rbinom(n, 1L, .5))
  d$m <- 0.2 + 0.8 * d$x + rnorm(n)
  d$y <- rbinom(n, 1L, stats::plogis(-0.5 + 0.5 * d$x + 0.9 * d$m))
  fm <- ilm_model(m ~ x, data = d, family = "gaussian", verbose = FALSE)
  fy <- ilm_model(y ~ x + m, data = d, family = "binomial", verbose = FALSE)
  md <- ilm_mediate(fm, fy, "x", "m", sims = 2000L, progress = FALSE)
  ## the total effect is a risk difference, and matches the observed one
  obs <- mean(d$y[d$x == 1L]) - mean(d$y[d$x == 0L])
  expect_equal(md$estimate[md$effect == "Total effect"], obs, tolerance = 0.03)
  expect_true(all(abs(md$estimate[1:5]) < 1))
  ## the product of coefficients is on the logit scale and is a different
  ## quantity entirely -- it is not even close
  prod <- unname(stats::coef(fm)["x"]) * unname(stats::coef(fy)["m"])
  expect_gt(prod / md$estimate[1], 3)
})

test_that("the counterfactual mediator is drawn, not fixed at its mean", {
  ## Y(t, M(t')) averages over the mediator's own variation. Fixing it at the
  ## mean is right only when the outcome is linear in it.
  d <- med_data(1500L, seed = 5L)
  f <- med_fit(d)
  a <- ilm_mediate(f$m, f$y, "x", "m", sims = 500L, seed = 1L,
                   progress = FALSE)
  b <- ilm_mediate(f$m, f$y, "x", "m", sims = 500L, seed = 2L,
                   progress = FALSE)
  ## different seeds give different draws, so the answer moves a little
  expect_false(identical(a$estimate, b$estimate))
  expect_equal(a$estimate[1], b$estimate[1], tolerance = 0.05)
  ## and the draws are kept, so the interval is a percentile of them
  dr <- attr(a, "draws")
  expect_equal(dim(dr), c(500L, 5L))
  expect_equal(a$estimate[1], mean(dr[, "acme0"]))
})

test_that("the sensitivity correction recovers a hidden confounder", {
  ## generate data WITH a confounder of known strength, hide it, and check
  ## that the correction at the true strength returns the true effect
  skip_on_cran()
  for (lam in c(0.4, 0.8)) {
    obs <- cor_ <- numeric(10L)
    for (i in seq_len(10L)) {
      set.seed(3000L + i); n <- 4000L
      u <- rnorm(n)
      d <- data.frame(x = rbinom(n, 1L, .5))
      d$m <- 0.3 + 0.7 * d$x + lam * u + rnorm(n)
      d$y <- 1 + 0.4 * d$x + 0.6 * d$m + lam * u + rnorm(n)
      fm <- ilm_model(m ~ x, data = d, family = "gaussian", verbose = FALSE)
      fy <- ilm_model(y ~ x + m, data = d, family = "gaussian",
                      verbose = FALSE)
      md <- ilm_mediate(fm, fy, "x", "m", sims = 200L, progress = FALSE)
      obs[i] <- md$estimate[1]
      cor_[i] <- ilm_mediate_sens(md, lambda_m = lam, lambda_y = lam)$acme
    }
    ## The observed effect is inflated and the correction removes it. Ten
    ## replicates leave a percent or two of Monte Carlo error on the mean; a
    ## 30-replicate run returned 0.4205, 0.4188 and 0.4167 against a true
    ## 0.4200 at confounder strengths of 0.4, 0.8 and 1.2, and the predicted
    ## bias in the mediator-outcome coefficient matched the realised one to
    ## three decimals at every strength.
    expect_gt(mean(obs), 0.42 + 0.05)
    expect_equal(mean(cor_), 0.7 * 0.6, tolerance = 0.05)
  }
})

test_that("the sensitivity grid reports where the effect would vanish", {
  d <- med_data(1000L, seed = 7L)
  f <- med_fit(d)
  md <- ilm_mediate(f$m, f$y, "x", "m", sims = 500L, progress = FALSE)
  ss <- ilm_mediate_sens(md)
  expect_s3_class(ss, "ilm_mediate_sens")
  expect_setequal(names(ss), c("lambda_m", "lambda_y", "acme"))
  ## no confounding leaves the observed effect alone
  z <- ss$acme[ss$lambda_m == 0 & ss$lambda_y == 0]
  expect_equal(z, attr(ss, "observed"), tolerance = 1e-10)
  ## stronger confounding shrinks it monotonically
  sub <- ss[ss$lambda_y == max(ss$lambda_y), ]
  expect_false(is.unsorted(rev(sub$acme)))
  ## and the reported threshold is where it hits zero
  thr <- attr(ss, "threshold")
  chk <- attr(ss, "observed") - attr(ss, "a") * thr / attr(ss, "s2m")
  expect_equal(chk, 0, tolerance = 1e-10)
  expect_output(print(ss), "reaches zero")
  ## it is not offered where it was not derived
  db <- med_data(800L, seed = 8L)
  db$yb <- rbinom(nrow(db), 1L, stats::plogis(db$y - mean(db$y)))
  fyb <- ilm_model(yb ~ x + m + c, data = db, family = "binomial",
                   verbose = FALSE)
  fmb <- ilm_model(m ~ x + c, data = db, family = "gaussian", verbose = FALSE)
  mb <- ilm_mediate(fmb, fyb, "x", "m", sims = 200L, progress = FALSE)
  expect_error(ilm_mediate_sens(mb), "exact for linear models")
})

test_that("it refuses models that do not describe the same decomposition", {
  d <- med_data(800L, seed = 9L)
  f <- med_fit(d)
  expect_error(ilm_mediate(f$m, f$y, "zz", "m", sims = 10L),
               "must be a predictor in BOTH")
  expect_error(ilm_mediate(f$m, f$y, "x", "zz", sims = 10L),
               "not a predictor in the outcome model")
  ## the mediator has to be the mediator model's response
  fw <- ilm_model(c ~ x, data = d, family = "gaussian", verbose = FALSE)
  expect_error(ilm_mediate(fw, f$y, "x", "m", sims = 10L),
               "not the response of the mediator model")
  expect_error(ilm_mediate(mtcars, f$y, "x", "m"), "must be a fitted ilm_model")
  ## different row counts mean different units
  d2 <- d[-1, ]
  f2 <- ilm_model(y ~ x + m + c, data = d2, family = "gaussian",
                  verbose = FALSE)
  expect_error(ilm_mediate(f$m, f2, "x", "m", sims = 10L), "same units")
  ## a random effect makes the counterfactual conditional rather than
  ## population-level
  d$g <- factor(rep(1:40, length.out = nrow(d)))
  fr <- ilm_model(m ~ x + c + (1 | g), data = d, family = "gaussian",
                  verbose = FALSE)
  expect_error(ilm_mediate(fr, f$y, "x", "m", sims = 10L), "random effects")
})
