## A term built from a variable -- log(x), a Fourier basis -- sits in the model
## frame as a column of its own, without the variable underneath. Everything
## that builds new rows from the fit's own rows (marginal means and slopes,
## average effects, scenarios) could not rebuild the term: "object 'x' not
## found", or, for a time column called t or time, R's own t() and time() in
## its place. The fit now keeps the variable.

dt_data <- function(seed = 1) {
  set.seed(seed); n <- 240
  d <- data.frame(t = rep(1:24, 10), g = factor(rep(1:10, each = 24)),
                  grp = factor(sample(c("a", "b", "c"), n, TRUE)),
                  x = stats::runif(n, 1, 10), x1 = stats::rnorm(n))
  d$y <- 0.5 * log(d$x) + 0.3 * (d$grp == "b") + d$x1 +
    sin(2 * pi * d$t / 12) + stats::rnorm(10, 0, 0.3)[d$g] + stats::rnorm(n)
  d
}

test_that("the fit keeps the variables its transformed terms are built from", {
  d <- dt_data(); d$x[3] <- NA          # a dropped row: the rows must line up
  f <- suppressMessages(ilm_model(y ~ log(x) + grp + (1 | g), data = d,
                                  family = "gaussian", verbose = FALSE))
  full <- ilm_data(f)
  expect_identical(names(f$data_extra), "x")
  expect_identical(nrow(full), nrow(f$model))
  expect_equal(log(full$x), full[["log(x)"]])
  ## a fit made before they were kept reads its frame, as it always did
  f$data_extra <- NULL
  expect_identical(ilm_data(f), f$model)
})

test_that("marginal means through log(x) are emmeans'", {
  skip_if_not_installed("emmeans")
  d <- dt_data()
  f <- ilm_model(y ~ log(x) + grp, data = d, family = "gaussian",
                 verbose = FALSE)
  m <- stats::lm(y ~ log(x) + grp, data = d)
  e1 <- as.data.frame(ilm_emmeans(f, "grp"))
  e2 <- as.data.frame(emmeans::emmeans(m, "grp"))
  ## the fit is an optimiser's, so it agrees with lm() to its tolerance
  expect_equal(e1$estimate, e2$emmean, tolerance = 1e-5)
  expect_equal(e1$se, e2$SE, tolerance = 1e-5)
})

test_that("slopes through log(x) are the exact derivative", {
  d <- dt_data()
  f <- ilm_model(y ~ log(x) * grp, data = d, family = "gaussian",
                 verbose = FALSE)
  b <- stats::coef(f)
  xbar <- mean(d$x)
  exact <- c(b[["log(x)"]], b[["log(x)"]] + b[["log(x):grpb"]],
             b[["log(x)"]] + b[["log(x):grpc"]]) / xbar
  tr <- as.data.frame(ilm_trends(f, "grp", var = "x"))
  expect_equal(tr$estimate, exact, tolerance = 1e-6)
})

test_that("average effects and scenarios run, and log(x) gets no silent zero", {
  d <- dt_data()
  d$yb <- stats::rbinom(nrow(d), 1, stats::plogis(0.8 * log(d$x) - 1 +
                                                    0.5 * (d$grp == "b")))
  f <- ilm_model(yb ~ log(x) + grp, data = d, family = "binomial",
                 verbose = FALSE)
  a <- ilm_ame(f)
  ## grp is a variable; the column "log(x)" is not, and changing it would move
  ## nothing, since the prediction rebuilds log(x) from x
  expect_setequal(a$term, "grp")
  b <- stats::coef(f)
  base <- b[[1]] + b[["log(x)"]] * log(d$x)
  expect_equal(a$estimate, c(mean(stats::plogis(base + b[["grpb"]]) -
                                    stats::plogis(base)),
                             mean(stats::plogis(base + b[["grpc"]]) -
                                    stats::plogis(base))),
               tolerance = 1e-8)
  ## a scenario on x moves log(x) with it
  s <- ilm_scenario(f, x = c(2, 8), progress = FALSE)
  nd2 <- d; nd2$x <- 2; nd8 <- d; nd8$x <- 8
  expect_equal(s$estimate, c(mean(stats::predict(f, nd2)),
                             mean(stats::predict(f, nd8))), tolerance = 1e-12)
  ## and what it holds is named as variables, not as the term's column
  expect_identical(attr(s, "held"), "grp")
})

test_that("average effects agree with marginaleffects through get_data()", {
  skip_if_not_installed("marginaleffects")
  d <- dt_data()
  d$yb <- stats::rbinom(nrow(d), 1, stats::plogis(0.8 * log(d$x) - 1 +
                                                    0.5 * (d$grp == "b")))
  f <- ilm_model(yb ~ log(x) + grp, data = d, family = "binomial",
                 verbose = FALSE)
  ilm_register_marginaleffects()
  mc <- as.data.frame(marginaleffects::avg_comparisons(f, variables = "grp"))
  expect_equal(ilm_ame(f)$estimate, mc$estimate, tolerance = 1e-6)
})

test_that("a Fourier term in t or time is rebuilt, not taken for t() or time()", {
  d <- dt_data()
  for (nm in c("t", "time", "mo")) {
    d2 <- d; d2[[nm]] <- d$t
    fo <- stats::as.formula(sprintf(
      "y ~ x1 + grp + ilm_fourier(%s, 12, K = 2) + (1 | g)", nm))
    f <- suppressMessages(ilm_model(fo, data = d2, family = "gaussian",
                                    verbose = FALSE))
    e <- as.data.frame(ilm_emmeans(f, "grp"))
    ## each mean is the typical group's prediction at the grid: x1 and the
    ## time at their means, the Fourier basis rebuilt there
    nd <- data.frame(x1 = mean(d$x1),
                     grp = factor(levels(d$grp), levels = levels(d$grp)),
                     g = d$g[1])
    nd[[nm]] <- mean(d$t)
    expect_equal(e$estimate,
                 as.vector(stats::predict(f, nd, type = "link")),
                 tolerance = 1e-10)
    expect_s3_class(ilm_trends(f, "grp", var = "x1"), "ilm_trends")
    expect_setequal(ilm_ame(f)$term, c("x1", "grp"))
  }
})
