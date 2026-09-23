## boundary = "avoid": the boundary-avoiding penalty of Chung et al. (2013,
## 2015), half the log-determinant of each grouping term's covariance added to
## the log-likelihood. Checked against blme while it was built -- a gamma(2) on
## an SD agreed to 3e-7 (gaussian) and 5e-6 (binomial, Laplace), a Wishart(4,
## scale = Inf) on a 2 x 2 covariance to 5e-4 (binomial). Here the check needs
## nothing outside the package: for a gaussian model the likelihood is exact,
## so the penalised estimate can be found by brute force and compared.

bnd_data <- function(seed) {
  set.seed(seed); n <- 300
  d <- data.frame(x = rnorm(n),
                  g = factor(sample(c("north", "central", "south"), n, TRUE)),
                  site = factor(sample(30, n, TRUE)))
  eta <- cbind(0.4 * d$x, -0.3 * d$x + 0.5 * (d$g == "south"))
  P <- exp(cbind(eta, 0)); P <- P / rowSums(P)
  d$y <- factor(apply(P, 1, function(p) sample(c("lo", "mid", "hi"), 1, prob = p)),
                levels = c("lo", "mid", "hi"))
  d
}

## the exact marginal log-likelihood of a gaussian model with random
## coefficients Z_g on each group, covariance Sig, residual SD sigma
gauss_ll <- function(d, form_z, beta, Sig, sigma) {
  X <- stats::model.matrix(~ t, d); Z <- stats::model.matrix(form_z, d)
  out <- 0
  for (gg in levels(d$g)) {
    i <- d$g == gg; Zi <- Z[i, , drop = FALSE]
    R <- chol(Zi %*% Sig %*% t(Zi) + sigma^2 * diag(sum(i)))
    e <- backsolve(R, d$y[i] - X[i, , drop = FALSE] %*% beta, transpose = TRUE)
    out <- out - 0.5 * (sum(i) * log(2 * pi) + 2 * sum(log(diag(R))) + sum(e^2))
  }
  out
}

test_that("the estimate is the maximum of the likelihood plus half log|Sigma|", {
  ## a random slope whose unpenalised fit is at a correlation of +1, so the
  ## penalty is doing real work, in both of its parts: Sigma_cat and Sigma_d
  set.seed(11); m <- 15; ni <- 8
  d <- data.frame(g = factor(rep(seq_len(m), each = ni)),
                  t = rep(seq(-1, 1, length.out = ni), m))
  u0 <- rnorm(m, 0, 0.6); u1 <- rnorm(m, 0, 0.05)
  d$y <- 2 + 0.4 * d$t + u0[as.integer(d$g)] + u1[as.integer(d$g)] * d$t +
    rnorm(nrow(d), 0, 0.8)
  f <- ilm_model(y ~ t + (1 + t | g), data = d, family = "gaussian",
                 boundary = "avoid", verbose = FALSE)
  V <- f$Sigma[[1]][1, 1] * f$Sigma_d[["g"]]
  negpl <- function(p) {
    L <- matrix(c(exp(p[3]), p[4], 0, exp(p[5])), 2, 2); S <- L %*% t(L)
    v <- tryCatch(gauss_ll(d, ~ t, p[1:2], S, exp(p[6])), error = function(e) -1e10)
    -(v + 0.5 * as.numeric(determinant(S, logarithm = TRUE)$modulus))
  }
  ## started well away from the estimate: a simplex to get near, then BFGS
  o <- stats::optim(c(2, 0, log(0.5), 0, log(0.5), 0), negpl,
                    control = list(maxit = 5000, reltol = 1e-10))
  o <- stats::optim(o$par, negpl, method = "BFGS",
                    control = list(reltol = 1e-12, maxit = 2000))
  L <- matrix(c(exp(o$par[3]), o$par[4], 0, exp(o$par[5])), 2, 2); S <- L %*% t(L)
  expect_lt(max(abs(coef(f) - o$par[1:2])), 1e-4)
  expect_lt(max(abs(sqrt(diag(V)) - sqrt(diag(S)))), 1e-3)
  expect_lt(abs(stats::cov2cor(V)[1, 2] - stats::cov2cor(S)[1, 2]), 1e-3)
  expect_lt(abs(unname(f$dispersion) - exp(o$par[6])), 1e-4)
  ## and the fit is inside its range, where the unpenalised one is not
  expect_lt(abs(stats::cov2cor(V)[1, 2]), 0.9)
  expect_length(f$boundary_terms, 0L)
})

test_that("logLik() is the likelihood of the data, not the penalised objective", {
  d <- bnd_data(4)                               # the site variance is zero
  fh <- suppressMessages(ilm_model(y ~ x + g + (1 | site), data = d,
                                   family = "multinomial", verbose = FALSE))
  fa <- ilm_model(y ~ x + g + (1 | site), data = d, family = "multinomial",
                  boundary = "avoid", verbose = FALSE)
  expect_identical(fa$boundary, "avoid")
  expect_true(is.finite(fa$re_penalty) && fa$re_penalty != 0)
  ## the unpenalised objective of the same model, evaluated at the penalised
  ## estimate, is what logLik() must report
  const <- as.numeric(logLik(fh)) + as.numeric(fh$opt$objective)
  expect_equal(as.numeric(logLik(fa)), const - as.numeric(fh$obj$fn(fa$opt$par)),
               tolerance = 1e-8)
  ## which cannot beat the maximum likelihood
  expect_lte(as.numeric(logLik(fa)), as.numeric(logLik(fh)) + 1e-6)
  ## away from the boundary, and usable as an ordinary fit
  expect_equal(fa$hessian_how, "tmb")
  expect_true(fa$ok)
  expect_true(all(sqrt(diag(fa$Sigma$site)) > 1e-3))
})

test_that("a boundary fit names the penalty with its cost, and a penalised fit does not", {
  d <- bnd_data(4)
  expect_message(ilm_model(y ~ x + g + (1 | site), data = d,
                           family = "multinomial", verbose = FALSE),
                 "boundary = \"avoid\"", fixed = TRUE)
  expect_message(ilm_model(y ~ x + g + (1 | site), data = d,
                           family = "multinomial", verbose = FALSE),
                 "further from zero", fixed = TRUE)
  expect_no_message(ilm_model(y ~ x + g + (1 | site), data = d,
                              family = "multinomial", boundary = "avoid",
                              verbose = FALSE),
                    message = "edge of its range")
  ## the verdict under summary(): the remedy is offered only where it is one
  held <- paste(utils::capture.output(ilm_trust_note("site")), collapse = "\n")
  done <- paste(utils::capture.output(ilm_trust_note("site", avoided = TRUE)),
                collapse = "\n")
  expect_match(held, "boundary = \"avoid\"", fixed = TRUE)
  expect_match(held, "variance comes out larger", fixed = TRUE)
  expect_no_match(done, "boundary = \"avoid\"", fixed = TRUE)
})

test_that("refits keep the penalty", {
  ## a refit without it is a different estimator: every bootstrap, simulation
  ## envelope and reduced model would have compared a penalised fit with an
  ## unpenalised one
  d <- bnd_data(4)
  fa <- ilm_model(y ~ x + g + (1 | site), data = d, family = "multinomial",
                  boundary = "avoid", verbose = FALSE)
  expect_identical(ilm_refit_stub(fa)$boundary, "avoid")
  fr <- ilm_refit_like(fa)
  expect_identical(fr$boundary, "avoid")
  expect_equal(fr$opt$objective, fa$opt$objective, tolerance = 1e-6)
})
