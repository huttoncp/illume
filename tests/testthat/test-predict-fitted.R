## predict(groups = "fitted"): each row with its own group's estimated effects,
## as lme4's predict() gives by default. A row with no estimated effect -- a
## new group, a time off the fitted cells -- is an error, never a silent zero.

pf_data <- function(seed = 5) {
  set.seed(seed); ng <- 30; ni <- 8
  d <- data.frame(g = factor(rep(sprintf("g%02d", seq_len(ng)), each = ni)),
                  x = stats::rnorm(ng * ni))
  b0 <- stats::rnorm(ng, 0, 0.8); b1 <- stats::rnorm(ng, 0, 0.4)
  d$y <- stats::rbinom(nrow(d), 1, stats::plogis(-0.2 + 0.7 * d$x + b0[d$g] +
                                                   b1[d$g] * d$x))
  d$yg <- 1 + 0.5 * d$x + b0[d$g] + b1[d$g] * d$x +
    stats::rnorm(nrow(d), 0, 0.7)
  d
}

## the fixed part plus a group's own modes, by hand from ilm_ranef()
pf_by_hand <- function(f, nd) {
  r <- ilm_ranef(f); b <- stats::coef(f)
  vapply(seq_len(nrow(nd)), function(i) {
    ri <- r[r$level == nd$g[i], ]
    b[[1]] + b[[2]] * nd$x[i] + ri$mode[ri$dim == "(Intercept)"] +
      ri$mode[ri$dim == "x"] * nd$x[i]
  }, 0)
}

test_that("the fitted rows' prediction is ilm_fitted()'s, and new rows' by hand", {
  d <- pf_data()
  f <- ilm_model(y ~ x + (1 + x | g), data = d, family = "binomial",
                 verbose = FALSE)
  expect_equal(as.vector(predict(f, groups = "fitted")),
               as.vector(ilm_fitted(f)), tolerance = 1e-12)
  ## new rows of fitted groups, with the random slope at the new x
  nd <- data.frame(g = c("g03", "g17", "g03"), x = c(-1, 0.5, 2))
  eta <- pf_by_hand(f, nd)
  expect_equal(as.vector(predict(f, nd, groups = "fitted")),
               stats::plogis(eta), tolerance = 1e-10)
  expect_equal(as.vector(predict(f, nd, groups = "fitted", type = "link")),
               eta, tolerance = 1e-10)
  ## and the rows of one group differ from a typical group's by its effect
  expect_false(isTRUE(all.equal(predict(f, nd, groups = "fitted"),
                                predict(f, nd, groups = "typical"))))
})

test_that("a row with no estimated effect is an error, not a zero", {
  d <- pf_data()
  f <- ilm_model(y ~ x + (1 | g), data = d, family = "binomial",
                 verbose = FALSE)
  nd <- data.frame(g = c("zz", "g01"), x = 0)
  expect_error(predict(f, nd, groups = "fitted"),
               "(zz) was not in the fit", fixed = TRUE)
  expect_error(predict(f, nd, groups = "fitted"),
               "groups = \"typical\" or groups = \"population\"", fixed = TRUE)
  expect_error(predict(f, data.frame(x = 0), groups = "fitted"),
               "grouping variable")
  ## ilm_ame() does not take the word, and says where it is answered
  expect_error(ilm_ame(f, groups = "fitted"), "predict(groups = \"fitted\")",
               fixed = TRUE)
})

test_that("the intervals carry each group's own uncertainty, exactly", {
  ## gaussian, so a prediction is linear in the parameters and its variance
  ## is a quadratic form in the joint covariance
  d <- pf_data()
  f <- ilm_model(yg ~ x + (1 + x | g), data = d, family = "gaussian",
                 verbose = FALSE)
  nd <- data.frame(g = c("g05", "g21"), x = c(-0.5, 1.5))
  expect_equal(as.vector(predict(f, nd, groups = "fitted")),
               pf_by_hand(f, nd), tolerance = 1e-10)
  S <- solve(as.matrix(illume:::ilm_joint_prec(f)))
  mp <- ilm_draws(f, nsim = 1, seed = 1, natural = FALSE)$map
  v <- vapply(seq_len(nrow(nd)), function(i) {
    a <- numeric(nrow(mp))
    a[mp$block == "beta"] <- c(1, nd$x[i])
    a[mp$block == "bvec" & mp$level == nd$g[i] & mp$dim == "(Intercept)"] <- 1
    a[mp$block == "bvec" & mp$level == nd$g[i] & mp$dim == "x"] <- nd$x[i]
    drop(a %*% S %*% a)
  }, 0)
  pu <- predict(f, nd, groups = "fitted", se.fit = TRUE, nsim = 4000)
  expect_true(pu$joint)
  expect_equal(as.vector(pu$se.fit), sqrt(v), tolerance = 0.05)
  expect_true(all(pu$lower < pu$fit & pu$fit < pu$upper))
})

test_that("a correlation over time contributes the fitted value of each cell", {
  set.seed(2)
  d <- data.frame(id = factor(rep(c("a", "b", "c"), each = 10)),
                  t = rep(1:10, 3))
  d$y <- stats::rnorm(30, 0, 0.5) + rep(cumsum(stats::rnorm(10, 0, 0.5)), 3) +
    rep(stats::rnorm(3), each = 10)
  f <- ilm_model(y ~ 1, data = d, family = "gaussian",
                 ar = ilm_rw1(~ t | id), verbose = FALSE)
  fv <- as.vector(ilm_fitted(f))
  expect_equal(as.vector(predict(f, groups = "fitted")), fv,
               tolerance = 1e-12)
  ## rows at fitted cells, in another order
  nd <- data.frame(t = c(7, 2, 10), id = c("b", "a", "c"))
  at <- match(paste(nd$id, nd$t), paste(d$id, d$t))
  expect_equal(as.vector(predict(f, nd, groups = "fitted")), fv[at],
               tolerance = 1e-12)
  ## past a group's last time the latent value would be a forecast
  expect_error(predict(f, data.frame(t = 12, id = "a"), groups = "fitted"),
               "forecast")
  expect_error(predict(f, data.frame(t = 3, id = "z"), groups = "fitted"),
               "a group the fit has not seen")
  ## a term built from vectors cannot place new rows
  fv2 <- ilm_model(y ~ 1, data = d, family = "gaussian",
                   ar = ilm_rw1(d$t, d$id), verbose = FALSE)
  expect_error(predict(fv2, nd, groups = "fitted"), "built from vectors")
  expect_equal(as.vector(predict(fv2, groups = "fitted")),
               as.vector(ilm_fitted(fv2)), tolerance = 1e-12)
})

test_that("an AR(1) grid cell without observations has its own fitted value", {
  set.seed(4)
  d <- data.frame(id = factor(rep(c("a", "b"), each = 36)),
                  t = rep(rep(1:12, each = 3), 2))
  d$y <- stats::rnorm(72) + rep(stats::arima.sim(list(ar = 0.6), 24), each = 3)
  d <- d[!(d$id == "a" & d$t == 6), ]            # a gap on a's grid
  f <- ilm_model(y ~ 1, data = d, family = "gaussian",
                 ar = ilm_ar1(~ t | id), verbose = FALSE)
  cl <- ilm_cells(f)
  k <- which(cl$group == "a" & cl$time == 6)
  expect_identical(cl$n_obs[k], 0L)
  Ba <- illume:::ilm_Bar_hat(f)
  expect_equal(as.vector(predict(f, data.frame(t = 6, id = "a"),
                                 groups = "fitted", type = "link")),
               unname(stats::coef(f)[1] + Ba[k, 1]), tolerance = 1e-12)
})

test_that("REML and a multinomial outcome read their effects by name", {
  d <- pf_data()
  fr <- ilm_model(yg ~ x + (1 + x | g), data = d, family = "gaussian",
                  reml = TRUE, verbose = FALSE)
  expect_equal(as.vector(predict(fr, groups = "fitted")),
               as.vector(ilm_fitted(fr)), tolerance = 1e-12)
  nd <- data.frame(g = c("g05", "g21"), x = c(-0.5, 1.5))
  expect_equal(as.vector(predict(fr, nd, groups = "fitted")),
               pf_by_hand(fr, nd), tolerance = 1e-10)
  set.seed(9)
  d$m <- factor(sample(c("p", "q", "r"), nrow(d), TRUE))
  fm <- suppressMessages(ilm_model(m ~ x + (1 | g), data = d,
                                   family = "multinomial", verbose = FALSE))
  expect_equal(unname(predict(fm, groups = "fitted")), unname(ilm_fitted(fm)),
               tolerance = 1e-12)
})
