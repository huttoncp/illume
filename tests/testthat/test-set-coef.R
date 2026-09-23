## set_coef() has to rebuild EVERYTHING the parameters determine
##
## get_coef() deliberately hands back the variance parameters beside the fixed
## effects, so a caller who sets one reasonably expects it to take effect.
## It used to rebuild only the mean and the covariance structure: setting a
## dispersion changed the stored vector and nothing else, and ilm_simulate()
## went on drawing from the fitted residual SD. Nothing errored. The draws
## simply had the wrong spread.
##
## Reported from outside the package -- which is, again, where every real
## defect in illume has come from, never from this suite.
##
## Each test below sets one parameter to a value chosen so the correct answer
## is known in closed form, because "it changed" is not the same as "it changed
## to the right thing": a rebuilt CAR(1) correlation was wrong by 7% and still
## looked like a perfectly ordinary correlation.

pv <- function(fit) names(fit$opt$par)     # TYPE labels, not the pretty names

test_that("a dispersion set through set_coef() reaches ilm_simulate()", {
  set.seed(1); n <- 400
  d <- data.frame(x = rnorm(n), g = factor(rep(1:40, each = 10)))
  u <- rnorm(40, 0, 0.5)
  d$y <- 1 + 0.5 * d$x + u[as.integer(d$g)] + rnorm(n, 0, 1)
  f <- ilm_model(y ~ x + (1 | g), data = d, family = "gaussian", verbose = FALSE)

  p <- get_coef.ilm_model(f)
  p[pv(f) == "logdisp"] <- log(4)
  g <- set_coef.ilm_model(f, p)

  expect_equal(unname(g$dispersion), 4, tolerance = 1e-8)
  ## the point of the whole thing: the draws have to have the new spread
  ys <- ilm_simulate(g, nsim = 100, seed = 3)
  expect_equal(mean(apply(ys, 2, stats::sd)), sqrt(16 + 0.25 + 0.25),
               tolerance = 0.05)
  ## and the fit's own value must be far enough away that passing is not luck
  expect_lt(unname(f$dispersion), 2)
})

test_that("the parameters that already worked still work", {
  set.seed(1); n <- 400
  d <- data.frame(x = rnorm(n), g = factor(rep(1:40, each = 10)))
  u <- rnorm(40, 0, 0.5)
  d$y <- 1 + 0.5 * d$x + u[as.integer(d$g)] + rnorm(n, 0, 1)
  f <- ilm_model(y ~ x + (1 | g), data = d, family = "gaussian", verbose = FALSE)

  p <- get_coef.ilm_model(f); p[pv(f) == "theta"] <- log(3)
  expect_equal(sqrt(set_coef.ilm_model(f, p)$Sigma[[1]][1, 1]), 3,
               tolerance = 1e-8)
  p <- get_coef.ilm_model(f); p[2] <- 5
  expect_equal(unname(set_coef.ilm_model(f, p)$beta[2, 1]), 5, tolerance = 1e-8)
})

test_that("setting the parameters back gives the fit back, unchanged", {
  ## marginaleffects sets the parameters back after every nudge, so a round
  ## trip that drifts would poison every later derivative
  set.seed(1); n <- 300
  d <- data.frame(x = rnorm(n), g = factor(rep(1:30, each = 10)))
  d$y <- 1 + 0.5 * d$x + rnorm(30, 0, 0.5)[as.integer(d$g)] + rnorm(n)
  f <- ilm_model(y ~ x + (1 | g), data = d, family = "gaussian", verbose = FALSE)
  g <- set_coef.ilm_model(f, get_coef.ilm_model(f))
  expect_equal(g$dispersion, f$dispersion)
  expect_equal(g$beta, f$beta)
  expect_equal(g$Sigma, f$Sigma)
})

test_that("an exact fit keeps reporting the dispersion on the n - p scale", {
  ## ilm_model() reports sigma on the unbiased scale wherever it can, so it
  ## agrees with lm(). A rebuilt object that used the ML scale would disagree
  ## with the object it came from.
  set.seed(8); de <- data.frame(x = rnorm(200))
  de$y <- 1 + de$x + rnorm(200, 0, 2)
  f <- ilm_model(y ~ x, data = de, family = "gaussian", verbose = FALSE)
  expect_true(isTRUE(f$exact_df))
  p <- get_coef.ilm_model(f); p[pv(f) == "logdisp"] <- log(5)
  expect_equal(unname(set_coef.ilm_model(f, p)$dispersion),
               5 * sqrt(200 / 198), tolerance = 1e-8)
  expect_equal(unname(set_coef.ilm_model(f, get_coef.ilm_model(f))$dispersion),
               unname(f$dispersion), tolerance = 1e-10)
})

test_that("ordinal thresholds are rebuilt, not copied", {
  ## zeta is built from increments so it stays increasing; it is a transform
  ## of the parameters, which is exactly why storing the vector is not enough
  set.seed(2); m <- 500
  d <- data.frame(x = rnorm(m))
  d$y <- factor(cut(0.8 * d$x + stats::rlogis(m), c(-Inf, -0.5, 0.8, Inf),
                    labels = c("a", "b", "c")), ordered = TRUE)
  f <- ilm_model(y ~ x, data = d, family = "ordinal", verbose = FALSE)
  p <- get_coef.ilm_model(f)
  p[pv(f) == "zeta_raw"] <- c(-2, log(1.5))          # cutpoints at -2 and -0.5
  g <- set_coef.ilm_model(f, p)
  expect_equal(unname(g$zeta), c(-2, -0.5), tolerance = 1e-8)
  expect_named(g$zeta, names(f$zeta))
})

test_that("a zero-inflation parameter is rebuilt", {
  set.seed(3); k <- 600
  d <- data.frame(x = rnorm(k))
  d$y <- ifelse(runif(k) < 0.3, 0, rpois(k, exp(1 + 0.4 * d$x)))
  f <- ilm_model(y ~ x, data = d, family = "poisson", ziformula = ~ 1,
                 verbose = FALSE)
  p <- get_coef.ilm_model(f); p[pv(f) == "gzi"] <- -1.5
  expect_equal(unname(set_coef.ilm_model(f, p)$zi_gamma[1]), -1.5,
               tolerance = 1e-8)
})

test_that("a dispersion MODEL is rebuilt row by row", {
  set.seed(5); q <- 500
  d <- data.frame(x = rnorm(q), grp = factor(sample(c("lo", "hi"), q, TRUE)))
  d$y <- 1 + d$x + rnorm(q, 0, ifelse(d$grp == "hi", 3, 1))
  f <- ilm_model(y ~ x, data = d, family = "gaussian", dispformula = ~ grp,
                 verbose = FALSE)
  p <- get_coef.ilm_model(f)
  ## raising the dispersion INTERCEPT by log(2) doubles every row, so it
  ## doubles the median the object reports
  p[which(pv(f) == "gamma")[1]] <- f$disp_gamma[1] + log(2)
  g <- set_coef.ilm_model(f, p)
  expect_equal(unname(g$dispersion), unname(f$dispersion) * 2, tolerance = 1e-6)
  expect_equal(g$disp_gamma[1], f$disp_gamma[1] + log(2), tolerance = 1e-8)
})

test_that("CAR(1) is rebuilt with the CAR(1) transform, not AR(1)'s", {
  ## AR(1) stores atanh(rho); CAR(1) stores log(range) and rho is
  ## exp(-1/range). Both come back inside (-1, 1), so using the wrong one is
  ## invisible in the value: here tanh() of the raw parameter gives 0.70 where
  ## the fit reported 0.66.
  set.seed(11); ni <- 60; nt <- 8
  d <- expand.grid(t = 1:nt, id = factor(seq_len(ni)))
  d$x <- rnorm(nrow(d))
  d$y <- 1 + 0.5 * d$x + unlist(lapply(seq_len(ni), function(i)
    as.numeric(stats::arima.sim(list(ar = 0.6), nt))))
  f <- ilm_model(y ~ x, data = d, family = "gaussian",
                 ar = suppressWarnings(ilm_car1(d$t, d$id)), verbose = FALSE)
  g <- set_coef.ilm_model(f, get_coef.ilm_model(f))
  expect_equal(g$rho, f$rho, tolerance = 1e-10)
  expect_equal(g$ar_range, f$ar_range, tolerance = 1e-10)
  ## the wrong transform really would have been wrong, by more than tolerance
  rr <- unname(g$opt$par[g$pnames == "ar:rho_raw"])
  expect_gt(abs(tanh(rr) - f$rho), 1e-3)
})

test_that("AR(1) is still rebuilt with tanh()", {
  set.seed(11); ni <- 60; nt <- 8
  d <- expand.grid(t = 1:nt, id = factor(seq_len(ni)))
  d$x <- rnorm(nrow(d))
  d$y <- 1 + 0.5 * d$x + unlist(lapply(seq_len(ni), function(i)
    as.numeric(stats::arima.sim(list(ar = 0.6), nt))))
  f <- ilm_model(y ~ x, data = d, family = "gaussian",
                 ar = suppressWarnings(ilm_ar1(d$t, d$id)), verbose = FALSE)
  g <- set_coef.ilm_model(f, get_coef.ilm_model(f))
  expect_equal(g$rho, f$rho, tolerance = 1e-10)
})

test_that("set_coef takes the fixed effects alone, and says what it wants otherwise", {
  set.seed(1); d <- data.frame(x = rnorm(100))
  d$y <- 1 + 0.5 * d$x + rnorm(100)
  f <- ilm_model(y ~ x, data = d, family = "gaussian", verbose = FALSE)
  g <- set_coef.ilm_model(f, c(2, -1))
  expect_equal(unname(coef(g)), c(2, -1))
  ## everything else stays where the fit put it
  expect_equal(unname(coef(g, full = TRUE)[-(1:2)]),
               unname(coef(f, full = TRUE)[-(1:2)]))
  expect_error(set_coef.ilm_model(f, 1:7), "coef(model, full = TRUE)",
               fixed = TRUE)
})
