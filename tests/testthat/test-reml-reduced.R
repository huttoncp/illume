## Under REML the fixed effects are integrated out, so a count or yes/no
## model's only outer parameters are its variances. With every one of them at
## its boundary the fit is the plain model: all are held, the check line says
## so (WARN, never a clean pass), and group-level predictions say they are
## conditional on zero variance. Before this, such fits fell through unheld --
## one with no intervals, one graded usable whose region intervals covered a
## third of the regions, and a binomial (1 | area) at zero graded clean.

q <- function(e) suppressMessages(suppressWarnings(e))
hess_line <- function(f) f$checks[f$checks$check == "hessian", ]

test_that("a binomial REML fit whose only variance is at zero is the plain model", {
  set.seed(10); A <- 30; n <- 600
  d <- data.frame(area = factor(sample(A, n, TRUE)), x = rnorm(n))
  d$y <- rbinom(n, 1, plogis(-0.3 + 0.5 * d$x))      # no area variance at all
  f <- q(ilm_model(y ~ x + (1 | area), data = d, family = "binomial", reml = TRUE, verbose = FALSE))
  expect_lt(exp(f$opt$par[names(f$opt$par) == "theta"]), 1e-3)
  expect_identical(f$hessian_how, "reduced")
  expect_identical(f$hessian_held, "area")
  h <- hess_line(f)
  expect_identical(h$status, "WARN")
  expect_match(h$detail, "reduces to the plain model", fixed = TRUE)
  ## the plain model's coefficients and standard errors
  g <- stats::glm(y ~ x, data = d, family = stats::binomial())
  expect_equal(unname(coef(f)), unname(coef(g)), tolerance = 1e-6)
  expect_equal(unname(sqrt(diag(vcov(f)))), unname(sqrt(diag(vcov(g)))), tolerance = 1e-5)
  ## maximum likelihood on the same data holds the term the usual way
  f0 <- q(ilm_model(y ~ x + (1 | area), data = d, family = "binomial", verbose = FALSE))
  expect_identical(f0$hessian_how, "boundary")
})

test_that("a gaussian REML fit keeps its residual SD, and holds the usual way", {
  set.seed(10); A <- 30; n <- 600
  d <- data.frame(area = factor(sample(A, n, TRUE)), x = rnorm(n))
  d$y <- 1 + 0.5 * d$x + rnorm(n)
  f <- q(ilm_model(y ~ x + (1 | area), data = d, family = "gaussian", reml = TRUE, verbose = FALSE))
  expect_identical(f$hessian_how, "boundary")
  expect_identical(hess_line(f)$status, "BOUNDARY")
})

## The BYM study's count arm (studies/scripts/spatial_bym.R): a 6 x 6 rook
## lattice, cell 5 (Ebar 5, total SD 0.6, structured share 0.2), replicates
## 229 and 421, where REML took both variances to zero.
bym_data <- function(rep) {
  k <- 6L; R <- 36L
  xy <- expand.grid(i = seq_len(k), j = seq_len(k))
  lev <- sprintf("r%03d", seq_len(R))
  nb <- lapply(seq_len(R), function(r) which(abs(xy$i - xy$i[r]) + abs(xy$j - xy$j[r]) == 1))
  names(nb) <- lev
  W <- matrix(0, R, R); for (r in seq_len(R)) W[r, nb[[r]]] <- 1
  e <- eigen(diag(rowSums(W)) - W, symmetric = TRUE); pos <- e$values > 1e-8
  sc <- sqrt(mean(rowSums(sweep(e$vectors[, pos]^2, 2, e$values[pos], "/"))))
  set.seed(7919L * 5L + rep)
  phi <- drop(e$vectors[, pos] %*% (stats::rnorm(sum(pos)) / sqrt(e$values[pos]))) / sc * sqrt(0.2) * 0.6
  theta <- stats::rnorm(R, 0, sqrt(0.8) * 0.6)
  x <- stats::rnorm(R); E <- stats::rgamma(R, shape = 2, rate = 2 / 5)
  reg <- factor(lev, levels = lev)
  list(d = data.frame(region = reg, x = x, E = E,
                      y = stats::rpois(R, E * exp(-0.2 + 0.3 * x + phi + theta))),
       nd = data.frame(region = reg, x = x, E = 1), nb = nb)
}

test_that("both ways a BYM REML fit fell through are now held, and say so", {
  skip_if_not_installed("mgcv")
  for (rep in c(229L, 421L)) {                 # 229: Hessian passed; 421: it failed
    s <- bym_data(rep); nb <- s$nb
    fo <- y ~ x + offset(log(E)) + s(region, bs = "mrf", xt = list(nb = nb)) + (1 | region)
    f <- q(ilm_model(fo, data = s$d, family = "poisson", reml = TRUE, verbose = FALSE))
    expect_identical(f$hessian_how, "reduced")
    expect_setequal(f$hessian_held, c("s(region)", "region"))
    expect_identical(hess_line(f)$status, "WARN")
    expect_true(all(is.finite(sqrt(diag(vcov(f))))))
    ## group-level intervals exist, and are not printed silently
    expect_warning(p <- stats::predict(f, newdata = s$nd, groups = "fitted",
                                       interval = "confidence", nsim = 100, seed = 1),
                   "conditional on those variances being zero")
    expect_true(is.list(p) && all(is.finite(p$lower)))
    ## maximum likelihood holds both terms the usual way
    f0 <- q(ilm_model(fo, data = s$d, family = "poisson", verbose = FALSE))
    expect_identical(f0$hessian_how, "boundary")
  }
})

test_that("no fit with every variance at its boundary grades as a clean pass", {
  ## Poisson counts with no area variance, 8 data sets, REML: every fit whose
  ## area SD is at zero is reduced and carries a WARN; at least one is
  at_zero <- 0L
  for (sd in 1:8) {
    set.seed(sd); A <- 25; n <- 500
    d <- data.frame(area = factor(sample(A, n, TRUE)), x = rnorm(n))
    d$y <- rpois(n, exp(0.5 + 0.3 * d$x))
    f <- q(ilm_model(y ~ x + (1 | area), data = d, family = "poisson", reml = TRUE, verbose = FALSE))
    if (exp(f$opt$par[names(f$opt$par) == "theta"]) < 1e-3) {
      at_zero <- at_zero + 1L
      expect_identical(f$hessian_how, "reduced")
      expect_false(hess_line(f)$status %in% c("OK", "BOUNDARY"))
    }
  }
  expect_gt(at_zero, 0L)
})
