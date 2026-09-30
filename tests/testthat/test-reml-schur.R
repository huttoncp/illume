## Under REML the coefficients' covariance is the Schur complement of the
## random effects' block of the joint precision. With the variances near
## zero that block's diagonal runs to about 1e20 beside the data's O(1), and
## a dense solve() stopped the fit as "computationally singular". It is now
## solved on unit-diagonal copies, which is exact.

q <- function(e) suppressMessages(suppressWarnings(e))

naive <- function(Hbb, Hbu, Huu) solve(Hbb - Hbu %*% solve(Huu, t(Hbu)))

test_that("the scaled solve is the plain one wherever the plain one works", {
  set.seed(1)
  for (k in 1:5) {
    nb <- 3; nu <- 12
    A <- crossprod(matrix(rnorm((nb + nu) * 40), 40)) + diag(nb + nu)
    Hbb <- A[1:nb, 1:nb]; Hbu <- A[1:nb, nb + 1:nu]; Huu <- A[nb + 1:nu, nb + 1:nu]
    v1 <- illume:::ilm_schur_vb(Hbb, Hbu, Huu); v0 <- naive(Hbb, Hbu, Huu)
    expect_lt(max(abs(v1 - v0) / abs(v0)), 1e-10)
  }
  expect_lt(max(abs(illume:::ilm_schur_vb(A) - solve(A)) / abs(solve(A))), 1e-10)
})

test_that("a well-conditioned REML fit's covariance does not move", {
  set.seed(1)
  d <- data.frame(g = factor(rep(1:12, each = 10)), x = rnorm(120))
  d$k <- rpois(120, exp(0.5 + 0.2 * d$x + rnorm(12, 0, 0.5)[d$g]))
  f <- q(ilm_model(k ~ x + (1 | g), data = d, family = "poisson", reml = TRUE, verbose = FALSE))
  jp <- f$sdr$jointPrecision; rn <- rownames(jp)
  ib <- which(rn == "beta"); iu <- which(rn == "bvec")
  v0 <- naive(as.matrix(jp[ib, ib]), as.matrix(jp[ib, iu]), as.matrix(jp[iu, iu]))
  v1 <- unname(as.matrix(stats::vcov(f)))
  expect_lt(max(abs(v1 - unname(v0)) / abs(v0)), 1e-10)
})

## The BYM study's count arm, cell 1: a 6 x 6 rook lattice, expected count 5,
## total SD 0.3, structured share 0.2 (studies/scripts/spatial_bym.R)
bym_cell1 <- function(rep) {
  k <- 6L; R <- 36L
  xy <- expand.grid(i = seq_len(k), j = seq_len(k))
  lev <- sprintf("r%03d", seq_len(R))
  nb <- lapply(seq_len(R), function(r) which(abs(xy$i - xy$i[r]) + abs(xy$j - xy$j[r]) == 1))
  names(nb) <- lev
  W <- matrix(0, R, R); for (r in seq_len(R)) W[r, nb[[r]]] <- 1
  e <- eigen(diag(rowSums(W)) - W, symmetric = TRUE); pos <- e$values > 1e-8
  sc <- sqrt(mean(rowSums(sweep(e$vectors[, pos]^2, 2, e$values[pos], "/"))))
  set.seed(7919L * 1L + rep)
  phi <- drop(e$vectors[, pos] %*% (stats::rnorm(sum(pos)) / sqrt(e$values[pos]))) / sc * sqrt(0.2) * 0.3
  theta <- stats::rnorm(R, 0, sqrt(0.8) * 0.3)
  x <- stats::rnorm(R); E <- stats::rgamma(R, shape = 2, rate = 2 / 5)
  list(d = data.frame(region = factor(lev, levels = lev), x = x, E = E,
                      y = stats::rpois(R, E * exp(-0.2 + 0.3 * x + phi + theta))), nb = nb)
}

test_that("REML fits whose variances approach zero now fit", {
  skip_if_not_installed("mgcv")
  for (rep in c(20L, 102L)) {
    s <- bym_cell1(rep); nb <- s$nb
    fo <- y ~ x + offset(log(E)) + s(region, bs = "mrf", xt = list(nb = nb)) + (1 | region)
    f <- NULL
    expect_no_error(f <- q(ilm_model(fo, data = s$d, family = "poisson", reml = TRUE, verbose = FALSE)))
    expect_true(all(is.finite(sqrt(diag(stats::vcov(f))))))
  }
})
