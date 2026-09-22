## Every refit -- a reduced model for a test, a bootstrap replicate, a
## simulation envelope, a consistency check -- has to be the SAME model as the
## fit it comes from. Each of those paths used to write out its own list of
## what "the same model" means, and each list went out of date: without the
## censoring the likelihood-ratio test rejected a true null every time, and
## without the zero part it put a predictor with no effect at chi-square 68.8
## where the right answer was 0.02. They are now all built from
## ilm_refit_stub(), and these tests hold that in place.

## The invariant: refitting a model on its OWN response, through the stub,
## gives back its own likelihood. A part of the model the stub forgot would
## make it a different model, and the likelihood would move.
test_that("a refit through the stub is the same model, whatever it contains", {
  set.seed(21); n <- 300
  d <- data.frame(x = rnorm(n), g = factor(sample(c("north", "central", "south"), n, TRUE)),
                  id = factor(sample(15, n, TRUE)))
  re <- rnorm(15, 0, .6)[d$id]
  d$y <- 1 + 0.5 * d$x + re + rnorm(n, 0, ifelse(d$g == "south", 2, 1))
  d$cnt <- ifelse(runif(n) < 0.3, 0, rpois(n, exp(0.4 + 0.3 * d$x)))
  d$yc <- pmin(d$y, 2)                                   # a ceiling at 2
  fits <- list(
    zero_inflated = ilm_model(cnt ~ x, data = d, family = "poisson",
                              ziformula = ~ 1, verbose = FALSE),
    hurdle        = ilm_model(cnt ~ x, data = d, family = "poisson",
                              ziformula = ~ x, zi_type = "hurdle",
                              verbose = FALSE),
    dispersion    = ilm_model(y ~ x, data = d, dispformula = ~ g,
                              verbose = FALSE),
    censored      = ilm_model(yc ~ x, data = d,
                              censor = ilm_censor(d$yc, upper = 2),
                              verbose = FALSE),
    reml          = ilm_model(y ~ x + (1 | id), data = d, reml = TRUE,
                              verbose = FALSE))
  for (nm in names(fits)) {
    f <- fits[[nm]]
    r <- ilm_refit_like(ilm_refit_stub(f), restarts = 1L)
    expect_equal(-r$opt$objective, -f$opt$objective, tolerance = 1e-6,
                 label = paste("refit likelihood,", nm))
  }
})

test_that("a likelihood-ratio test of a zero-inflated model keeps its zero part", {
  set.seed(11); n <- 500
  d <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  d$y <- ifelse(runif(n) < 0.35, 0, rpois(n, exp(0.5 + 0.4 * d$x1)))
  full <- ilm_model(y ~ x1 + x2, data = d, family = "poisson",
                    ziformula = ~ 1, verbose = FALSE)
  red <- ilm_model(y ~ x1, data = d, family = "poisson",
                   ziformula = ~ 1, verbose = FALSE)
  by_hand <- 2 * (as.numeric(logLik(full)) - as.numeric(logLik(red)))
  a <- suppressMessages(ilm_anova(full, test = "LRT"))
  expect_equal(a["x2", 2], by_hand, tolerance = 1e-4)
  expect_gt(a["x2", 3], 0.05)                     # x2 has no effect
  ## the parametric bootstrap's observed statistic is the same test
  pb <- ilm_pb_lrt(full, "x2", B = 5, verbose = FALSE)
  expect_equal(pb$LR, by_hand, tolerance = 1e-4)
})

test_that("McFadden's R-squared is not formed from restricted likelihoods", {
  skip_if_not_installed("performance")
  set.seed(3); n <- 200
  d <- data.frame(x = rnorm(n), id = factor(sample(10, n, TRUE)))
  d$y <- 0.5 * d$x + rnorm(10)[d$id] + rnorm(n)
  ml <- ilm_model(y ~ x + (1 | id), data = d, verbose = FALSE)
  rl <- ilm_model(y ~ x + (1 | id), data = d, reml = TRUE, verbose = FALSE)
  expect_true(is.finite(model_performance.ilm_model(ml, verbose = FALSE)$R2_McFadden))
  pr <- model_performance.ilm_model(rl, verbose = FALSE)
  expect_true(is.na(pr$R2_McFadden))
  expect_match(attr(pr, "r2_note"), "reml = FALSE", fixed = TRUE)
})

test_that("fit indices work for every family, scoring only what can be scored", {
  ## model_performance() used to fail on any fit whose outcome is not a set of
  ## categories, inside the scoring rules, on "subscript out of bounds"
  set.seed(4); n <- 200
  d <- data.frame(x = rnorm(n))
  d$y <- 0.5 * d$x + rnorm(n)
  d$b <- rbinom(n, 1, stats::plogis(0.8 * d$x))
  pg <- model_performance.ilm_model(
    ilm_model(y ~ x, data = d, verbose = FALSE), verbose = FALSE)
  expect_true(is.na(pg$Brier))
  expect_match(attr(pg, "score_note"), "NA for this gaussian fit", fixed = TRUE)
  expect_error(ilm_scores(ilm_model(y ~ x, data = d, verbose = FALSE)),
               "categorical outcome")
  ## a binary outcome is two categories, and is scored as such
  fb <- ilm_model(b ~ x, data = d, family = "binomial", verbose = FALSE)
  s <- ilm_scores(fb)
  p <- ilm_fitted(fb, TRUE)[, 1]
  expect_equal(unname(s["brier"]), mean(2 * (d$b - p)^2))
  expect_equal(unname(s["log_score"]),
               mean(-log(ifelse(d$b == 1, p, 1 - p))))
})
