## Every likelihood a result is built from goes through logLik() (review
## findings 2.2, 2.13, 2.19). The optimiser's objective leaves out the latent
## values' (q/2) log(2 pi) and includes any boundary penalty, so it is not
## the log-likelihood; comparisons that used it compared unlike things.

glmm_data <- function() {
  set.seed(1)
  ng <- 40
  d <- data.frame(id = factor(rep(1:ng, each = 6)))
  d$x <- stats::rnorm(nrow(d))
  d$y <- stats::rbinom(nrow(d), 1, stats::plogis(-0.3 + 0.8 * d$x + stats::rnorm(ng)[d$id]))
  d
}

test_that("McFadden's R-squared of a mixed model is built from logLik(), as lme4's", {
  skip_if_not_installed("lme4")
  d <- glmm_data()
  f <- ilm_model(y ~ x + (1 | id), data = d, family = "binomial", verbose = FALSE)
  f0 <- ilm_model(y ~ 1 + (1 | id), data = d, family = "binomial", verbose = FALSE)
  expect_equal(illume:::ilm_null_ll(f), as.numeric(logLik(f0)), tolerance = 1e-6)
  r2 <- illume:::model_performance.ilm_model(f)$R2_McFadden
  g <- lme4::glmer(y ~ x + (1 | id), data = d, family = stats::binomial)
  g0 <- lme4::glmer(y ~ 1 + (1 | id), data = d, family = stats::binomial)
  expect_equal(r2, 1 - as.numeric(logLik(g)) / as.numeric(logLik(g0)), tolerance = 1e-4)
  expect_gt(r2, 0)
})

test_that("print() shows logLik(), not the objective", {
  d <- glmm_data()
  f <- ilm_model(y ~ x + (1 | id), data = d, family = "binomial", verbose = FALSE)
  out <- capture.output(print(f))
  expect_true(any(grepl(sprintf("logLik %.2f", as.numeric(logLik(f))), out, fixed = TRUE)))
  expect_false(any(grepl(sprintf("logLik %.2f", -f$opt$objective), out, fixed = TRUE)))
})

test_that("a likelihood-ratio test under boundary = \"avoid\" leaves the penalty out", {
  set.seed(8)
  ng <- 15
  d <- data.frame(id = factor(rep(1:ng, each = 8)), x = stats::rnorm(120), z = stats::rnorm(120))
  d$y <- stats::rbinom(120, 1, stats::plogis(-0.2 + 0.7 * d$x + stats::rnorm(ng, 0, 0.8)[d$id]))
  f1 <- ilm_model(y ~ x + z + (1 | id), data = d, family = "binomial", boundary = "avoid",
                  verbose = FALSE)
  f0 <- ilm_model(y ~ z + (1 | id), data = d, family = "binomial", boundary = "avoid",
                  verbose = FALSE)
  a <- suppressMessages(ilm_anova(f1, test = "LRT"))
  expect_equal(a["x", "Chisq"], 2 * (as.numeric(logLik(f1)) - as.numeric(logLik(f0))),
               tolerance = 1e-5)
})
