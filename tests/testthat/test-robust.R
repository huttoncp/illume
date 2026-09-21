rob_data <- function(G = 40L, ni = 10L, seed = 1L) {
  set.seed(seed)
  n <- G * ni
  d <- data.frame(g = factor(rep(seq_len(G), each = ni)), x = rnorm(n),
                  z = rnorm(n))
  d$y <- 1 + 0.5 * d$x - 0.3 * d$z + rep(rnorm(G, 0, 1.5), each = ni) + rnorm(n)
  d
}

test_that("the score contributions are the ones sandwich builds", {
  skip_if_not_installed("sandwich")
  d <- rob_data()
  f <- ilm_model(y ~ x + z, data = d, family = "gaussian", verbose = FALSE)
  l <- stats::lm(y ~ x + z, data = d)
  ## illume's carry the 1/sigma^2 that lm's estfun leaves out; the sandwich
  ## cancels it, so the comparison puts it back
  s2 <- unname(f$dispersion)^2
  expect_equal(max(abs(ilm_estfun(f) * s2 - sandwich::estfun(l))), 0,
               tolerance = 1e-5)
})

test_that("with one row per cluster these are HC0, HC1 and HC2", {
  skip_if_not_installed("sandwich")
  d <- rob_data()
  f <- ilm_model(y ~ x + z, data = d, family = "gaussian", verbose = FALSE)
  l <- stats::lm(y ~ x + z, data = d)
  for (tp in c("CR0", "CR1", "CR2")) {
    mine <- sqrt(diag(ilm_vcov_cluster(f, NULL, type = tp)))
    theirs <- sqrt(diag(sandwich::vcovHC(l, type = sub("CR", "HC", tp))))
    expect_equal(unname(mine), unname(theirs), tolerance = 1e-7, info = tp)
  }
  ## CR2 reducing exactly to HC2 is the check on the leverage adjustment:
  ## with 1x1 blocks, (I - H_gg)^(-1/2) is 1 / sqrt(1 - h_ii)
})

test_that("clustered CR0 and CR1 agree with sandwich::vcovCL", {
  skip_if_not_installed("sandwich")
  d <- rob_data()
  f <- ilm_model(y ~ x + z, data = d, family = "gaussian", verbose = FALSE)
  l <- stats::lm(y ~ x + z, data = d)
  for (tp in c("CR0", "CR1")) {
    mine <- sqrt(diag(ilm_vcov_cluster(f, ~ g, type = tp)))
    theirs <- sqrt(diag(sandwich::vcovCL(l, cluster = d$g,
                                         type = sub("CR", "HC", tp),
                                         cadjust = (tp == "CR1"))))
    expect_equal(unname(mine), unname(theirs), tolerance = 1e-6, info = tp)
  }
  ## and a count model, where there is no dispersion to cancel
  set.seed(3); G <- 40L; ni <- 10L; n <- G * ni
  dp <- data.frame(g = factor(rep(seq_len(G), each = ni)), x = rnorm(n))
  dp$y <- rpois(n, exp(0.5 + 0.4 * dp$x + rep(rnorm(G, 0, 0.6), each = ni)))
  fp <- ilm_model(y ~ x, data = dp, family = "poisson", verbose = FALSE)
  gp <- stats::glm(y ~ x, data = dp, family = stats::poisson())
  expect_equal(
    unname(sqrt(diag(ilm_vcov_cluster(fp, ~ g, type = "CR0")))),
    unname(sqrt(diag(sandwich::vcovCL(gp, cluster = dp$g, type = "HC0",
                                      cadjust = FALSE)))),
    tolerance = 1e-6)
})

test_that("the ordinal score matches a numerical derivative", {
  set.seed(5); n <- 500L
  d <- data.frame(x = rnorm(n), w = rnorm(n))
  z <- 0.8 * d$x - 0.4 * d$w + stats::rlogis(n)
  d$y <- factor(cut(z, c(-Inf, -0.8, 0.9, Inf),
                    labels = c("lo", "mid", "hi")), ordered = TRUE)
  f <- ilm_model(y ~ x + w, data = d, family = "ordinal", verbose = FALSE)
  S <- ilm_estfun(f)
  ## d(log P)/d(eta) by central differences on the row probabilities
  eta <- as.numeric(ilm_eta_hat(f, FALSE)[, 1])
  yi <- as.integer(f$y); h <- 1e-5
  lp <- function(e) {
    P <- ilm_ord_probs(e, f$zeta, f$family$pfun)
    log(P[cbind(seq_along(e), yi)])
  }
  u_num <- (lp(eta + h) - lp(eta - h)) / (2 * h)
  expect_equal(unname(S[, "x"]), unname(u_num * f$X[, "x"]), tolerance = 1e-5)
  expect_equal(unname(S[, "w"]), unname(u_num * f$X[, "w"]), tolerance = 1e-5)
  ## and the score sums to zero at the maximum, which is what being there means
  expect_equal(unname(colSums(S)), c(0, 0), tolerance = 1e-3)
})

test_that("robust standard errors are larger when the clustering is ignored", {
  d <- rob_data()
  f <- ilm_model(y ~ x + z, data = d, family = "gaussian", verbose = FALSE)
  r <- ilm_robust(f, ~ g)
  expect_s3_class(r, "ilm_robust")
  expect_equal(nrow(r), 3L)
  expect_equal(attr(r, "n_clusters"), 40L)
  ## the cluster effect is on the intercept here, so that is the one that
  ## the model-based standard error badly understates
  i <- which(r$term == "(Intercept)")
  expect_gt(r$se[i] / r$se_model[i], 2)
  ## the estimates themselves do not move -- only the uncertainty does
  expect_equal(r$estimate, unname(stats::coef(f)))
  ## Bell-McCaffrey degrees of freedom sit below G - 1
  expect_true(all(r$df < 40))
  expect_true(all(r$df > 1))
  ## and the reference can be asked for explicitly
  expect_true(all(is.infinite(ilm_robust(f, ~ g, df = "normal")$df)))
  expect_true(all(ilm_robust(f, ~ g, df = "G-1")$df == 39))
  ## a wider reference gives a wider interval
  w <- function(x) mean(x$upper - x$lower)
  expect_gt(w(ilm_robust(f, ~ g, df = "bm")), w(ilm_robust(f, ~ g, df = "normal")))
})

test_that("CR2 and the Bell-McCaffrey reference actually cover", {
  ## the claim the defaults rest on, at a size where CR0 plainly does not
  skip_on_cran()
  G <- 15L; ni <- 20L
  hit0 <- hit2 <- logical(150L)
  for (i in seq_len(150L)) {
    set.seed(6000 + i)
    xg <- rnorm(G)
    g <- factor(rep(seq_len(G), each = ni)); n <- G * ni
    dd <- data.frame(g = g, x = xg[as.integer(g)] + rnorm(n, 0, 0.3))
    dd$y <- 1 + 0.5 * dd$x + rep(rnorm(G, 0, 1.2), each = ni) + rnorm(n)
    ff <- try(ilm_model(y ~ x, data = dd, family = "gaussian",
                        verbose = FALSE), silent = TRUE)
    if (inherits(ff, "try-error")) next
    a <- ilm_robust(ff, ~ g, type = "CR0", df = "normal")
    b <- ilm_robust(ff, ~ g, type = "CR2", df = "bm")
    j <- which(a$term == "x")
    hit0[i] <- a$lower[j] < 0.5 && 0.5 < a$upper[j]
    hit2[i] <- b$lower[j] < 0.5 && 0.5 < b$upper[j]
  }
  ## 150 replicates put a standard error of about 0.024 on each of these, so
  ## the assertion is the comparison rather than a threshold either side of
  ## one of them: CR2 with a Bell-McCaffrey reference lands nearer the
  ## nominal 0.95 than CR0 with a normal one does. In a 600-replicate run the
  ## two were 0.798 and 0.950 at G = 10, and 0.932 and 0.952 at G = 40.
  expect_lt(abs(mean(hit2) - 0.95), abs(mean(hit0) - 0.95))
  expect_lt(mean(hit0), 0.94)
  expect_gt(mean(hit2), 0.90)
})

test_that("it refuses the cases where the score does not decompose", {
  d <- rob_data()
  fr <- ilm_model(y ~ x + (1 | g), data = d, family = "gaussian",
                  verbose = FALSE)
  expect_error(ilm_estfun(fr), "random effects")
  expect_error(ilm_robust(fr, ~ g), "random effects")

  f <- ilm_model(y ~ x + z, data = d, family = "gaussian", verbose = FALSE)
  expect_error(ilm_vcov_cluster(f, ~ nope), "neither in the model frame")
  expect_error(ilm_vcov_cluster(f, rep(1, 5)), "but the model was fitted")
  expect_error(ilm_vcov_cluster(f, rep(1, nrow(d))), "at least two clusters")
  expect_error(ilm_vcov_cluster(mtcars), "must be a fitted ilm_model")

  ## a zero part and a multinomial each give a row more than one linear
  ## predictor, and the score for one is not separable from the other
  set.seed(4); n <- 300L
  dz <- data.frame(x = rnorm(n), g = factor(rep(1:30, each = 10)))
  dz$y <- ifelse(runif(n) < 0.3, 0, rpois(n, exp(0.5 + 0.3 * dz$x)))
  fz <- ilm_model(y ~ x, data = dz, family = "poisson", ziformula = ~ 1,
                  verbose = FALSE)
  expect_error(ilm_estfun(fz), "zero part")
})

test_that("a cluster found outside the model frame lines up after dropped rows", {
  d <- rob_data(G = 30L)
  d$x[c(3, 17, 55)] <- NA                 # three rows the fit will drop
  f <- ilm_model(y ~ x + z, data = d, family = "gaussian", verbose = FALSE)
  expect_equal(nrow(f$X), nrow(d) - 3L)
  ## g is not a predictor, so it comes from the original data -- and the same
  ## three rows have to come out of it, which is the step that otherwise
  ## recycles a vector three too long without complaint
  V <- ilm_vcov_cluster(f, ~ g)
  expect_equal(sum(attr(V, "cluster_sizes")), nrow(d) - 3L)
  expect_equal(attr(V, "n_clusters"), 30L)
  ## passing the full-length vector instead is refused rather than recycled
  expect_error(ilm_vcov_cluster(f, d$g), "but the model was fitted")
})
