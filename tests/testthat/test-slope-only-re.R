## A slope-only random effect, (0 + x | g), is a slope on new rows too
## (review finding 3.1). A one-column term was read as an intercept wherever
## the bar was evaluated afresh: new-data predictions, population averages,
## ilm_matrices() and power simulation.

slope_only_data <- function() {
  set.seed(25)
  ng <- 25; m <- 10
  d <- data.frame(id = factor(sprintf("s%02d", rep(1:ng, each = m))))
  d$x <- stats::rnorm(nrow(d), 1, 1)
  b1 <- stats::rnorm(ng, 0, 0.8)
  d$y <- 1 + (0.5 + b1[d$id]) * d$x + stats::rnorm(nrow(d), 0, 0.5)
  d$yc <- stats::rpois(nrow(d), exp(0.2 + (0.3 + b1[d$id]) * d$x))
  d
}

test_that("new-data predictions of a slope-only term match lme4 and the fitted rows", {
  skip_if_not_installed("lme4")
  d <- slope_only_data()
  f <- ilm_model(y ~ x + (0 + x | id), data = d, family = "gaussian", verbose = FALSE)
  l <- lme4::lmer(y ~ x + (0 + x | id), data = d)
  nd <- data.frame(id = c("s01", "s02", "s01"), x = c(3, -2, 0))
  expect_equal(as.numeric(predict(f, newdata = nd)), as.numeric(predict(l, newdata = nd)),
               tolerance = 1e-5)
  expect_equal(as.numeric(predict(f, newdata = d[1:3, ])), as.numeric(predict(f))[1:3],
               tolerance = 1e-10)
})

test_that("a slope-only term's population average and design are a slope's", {
  d <- slope_only_data()
  fp <- ilm_model(yc ~ x + (0 + x | id), data = d, family = "poisson", verbose = FALSE)
  s2 <- unlist(fp$Sigma)
  ndp <- data.frame(x = c(0, 2))
  pp <- as.numeric(predict(fp, newdata = ndp, groups = "population"))
  ## E[exp(eta + b x)] for b ~ N(0, s2) is exp(eta + s2 x^2 / 2)
  expect_equal(pp, as.numeric(exp(fp$beta[1] + fp$beta[2] * ndp$x + s2 * ndp$x^2 / 2)),
               tolerance = 1e-6)
  mm <- ilm_matrices(fp, data.frame(id = "s01", x = 3))
  expect_equal(as.numeric(mm$re$id$Z), 3)
})

test_that("a power scaffold with a slope-only term simulates a slope, not an intercept", {
  skip_if_not_installed("lme4")
  s <- suppressWarnings(ilm_scaffold(y ~ time + (0 + time | id),
         design = list(time = c(0, 1, 2, 3)), within = "time", n_unit = 300,
         coefs = c("(Intercept)" = 10, time = 0.2), sd = 1, re_sd = list(id = 0.8),
         verbose = FALSE))
  g <- illume:::ilm_scaffold_grid(s$scaffold$design, 300, s$scaffold$within, s$scaffold$group)
  set.seed(1)
  g$y <- illume:::ilm_power_response(s, illume:::ilm_power_eta_grid(s, g, s$beta))
  l <- suppressMessages(lme4::lmer(y ~ time + (1 | id) + (0 + time | id), data = g))
  vc <- as.data.frame(lme4::VarCorr(l))
  expect_gt(vc$sdcor[vc$grp == "id.1"], 0.6)       # the slope, assumed 0.8
  expect_lt(vc$sdcor[vc$grp == "id"], 0.3)         # no intercept was assumed
})
