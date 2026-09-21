glrm_num <- function(n = 300L, p = 6L, k = 2L, seed = 1L, sd_e = 0.4) {
  set.seed(seed)
  F <- matrix(rnorm(n * k), n, k); L <- matrix(rnorm(k * p), k, p)
  as.data.frame(F %*% L + matrix(rnorm(n * p, 0, sd_e), n, p))
}
glrm_mixed <- function(n = 400L, seed = 1L) {
  set.seed(seed)
  f <- rnorm(n); f2 <- rnorm(n)
  d <- data.frame(a = f + rnorm(n, 0, .4), b = -f + rnorm(n, 0, .4),
                  c = f2 + rnorm(n, 0, .4))
  d$g <- factor(ifelse(stats::plogis(2.5 * f) > runif(n), "hi", "lo"))
  d
}
hide <- function(d, frac = .2, seed = 1L) {
  set.seed(seed)
  mi <- matrix(runif(nrow(d) * ncol(d)) < frac, nrow(d), ncol(d))
  dm <- d
  for (j in seq_len(ncol(d))) dm[mi[, j], j] <- NA
  list(dm = dm, mi = mi)
}

test_that("quadratic loss on numeric data IS principal components", {
  ## the identity the whole thing rests on: change the loss and it generalises
  ## PCA, keep it quadratic and it must reproduce PCA exactly
  d <- glrm_num()
  n <- nrow(d)
  g <- ilm_glrm(d, rank = 2L, lambda = 0, maxit = 5000L, tol = 1e-13,
                progress = FALSE)
  Z <- scale(as.matrix(d))
  sv <- svd(Z, nu = 2L, nv = 2L)
  rec <- sv$u[, 1:2] %*% diag(sv$d[1:2]) %*% t(sv$v[, 1:2])
  U <- as.matrix(g$scores) %*% g$archetypes + rep(g$offset, each = n)
  expect_lt(max(abs(U - rec)), 1e-8)
  ## and the subspaces coincide, whatever basis each chose for it
  a <- qr.Q(qr(t(g$archetypes))); b <- qr.Q(qr(sv$v[, 1:2]))
  expect_lt(max(acos(pmin(svd(t(a) %*% b)$d, 1))), 1e-6)
})

test_that("a category comes back as a category", {
  d <- glrm_mixed()
  g <- ilm_glrm(d, rank = 1L, progress = FALSE)
  expect_s3_class(g, "ilm_glrm")
  expect_equal(g$loss[["g"]], "logistic")
  expect_equal(unname(g$loss[c("a", "b", "c")]), rep("quadratic", 3L))
  expect_s3_class(g$fitted$g, "factor")
  expect_setequal(levels(g$fitted$g), levels(d$g))
  expect_true(is.numeric(g$fitted$a))
  expect_equal(nrow(g$scores), nrow(d))
  expect_equal(ncol(g$scores), 1L)
  ## three levels get a multinomial block, one column per level
  d3 <- d; d3$g <- factor(sample(c("x", "y", "z"), nrow(d), TRUE))
  g3 <- ilm_glrm(d3, rank = 2L, progress = FALSE)
  expect_equal(g3$loss[["g"]], "multinomial")
  expect_equal(length(g3$encoding$blocks$g$cols), 3L)
  expect_setequal(levels(g3$fitted$g), c("x", "y", "z"))
})

test_that("the ridge penalty is chosen, and it matters", {
  ## a row with few observed cells has as many scores as observations and fits
  ## them exactly, so too little penalty sends everything else in that row
  ## wherever the algebra points
  d <- glrm_num(n = 400L, p = 4L, k = 2L, seed = 7L)
  h <- hide(d, seed = 2L)
  err <- function(lam) {
    g <- ilm_glrm(h$dm, rank = 2L, lambda = lam, maxit = 2000L,
                  progress = FALSE)
    sqrt(mean((as.matrix(g$fitted)[h$mi] - as.matrix(d)[h$mi])^2))
  }
  e_small <- err(0.01); e_mid <- err(2); e_big <- err(50)
  expect_lt(e_mid, e_small)
  expect_lt(e_mid, e_big)
  ## chosen automatically, it lands in the good region
  auto <- ilm_glrm(h$dm, rank = 2L, maxit = 2000L, progress = FALSE)
  expect_gte(auto$lambda, 0.5)
  expect_lte(auto$lambda, 10)
  e_auto <- sqrt(mean((as.matrix(auto$fitted)[h$mi] - as.matrix(d)[h$mi])^2))
  expect_lt(e_auto, e_small)
  ## and it beats an iterative SVD on the same held-out cells
  M <- scale(as.matrix(h$dm)); ok <- !is.na(M); M[!ok] <- 0
  ctr <- attr(M, "scaled:center"); scl <- attr(M, "scaled:scale")
  for (it in 1:200) {
    s <- svd(M, nu = 2L, nv = 2L)
    R <- s$u %*% diag(s$d[1:2]) %*% t(s$v); M[!ok] <- R[!ok]
  }
  Rb <- sweep(sweep(R, 2L, scl, "*"), 2L, ctr, "+")
  expect_lt(e_auto, sqrt(mean((Rb[h$mi] - as.matrix(d)[h$mi])^2)))
})

test_that("the optimiser keeps going after a rejected step", {
  ## a first version tested the objective at the top of each pass, so the pass
  ## after a rejected step saw the same point and read no-change as
  ## convergence; the fit stopped after five iterations
  d <- glrm_mixed(seed = 3L)
  g <- ilm_glrm(d, rank = 2L, maxit = 2000L, progress = FALSE)
  expect_gt(g$iterations, 20L)
  ## the objective never rises, because a step that raised it was not taken
  expect_true(all(diff(g$objective) <= 1e-9))
  expect_true(g$converged)
})

test_that("missing cells are skipped by the loss and filled by the fit", {
  d <- glrm_mixed(seed = 4L)
  h <- hide(d, seed = 3L)
  g <- ilm_glrm(h$dm, rank = 1L, progress = FALSE)
  expect_false(anyNA(g$fitted))
  expect_equal(nrow(g$fitted), nrow(d))
  ## the hidden categories are recovered well above chance
  expect_gt(mean(g$fitted$g[h$mi[, 4]] == d$g[h$mi[, 4]]), 0.6)
  ## a row whose category was never observed still gets a valid subscript in
  ## the multinomial gradient rather than an NA one
  d3 <- d; d3$g <- factor(sample(c("x", "y", "z"), nrow(d), TRUE))
  h3 <- hide(d3, frac = 0.3, seed = 4L)
  expect_error(ilm_glrm(h3$dm, rank = 2L, progress = FALSE), NA)
})

test_that("ilm_reduce takes it as a method and hands back the same shape", {
  d <- glrm_mixed(seed = 5L)
  r1 <- ilm_reduce(d, ndim = 2L)
  r2 <- ilm_reduce(d, ndim = 2L, method = "glrm", progress = FALSE)
  expect_s3_class(r2, "ilm_reduce")
  expect_equal(r2$method, "glrm")
  expect_setequal(names(r1), names(r2))
  expect_equal(nrow(r2$ind_coord), nrow(d))
  expect_equal(r2$ndim, 2L)
  ## the share of variation explained is against the TOTAL, so it is not 100%
  ## by construction
  expect_lt(max(r2$eig$cum_pct_var), 100)
  expect_false(is.unsorted(rev(r2$eig$pct_var)))
  ## the two routes find the same leading direction
  expect_gt(abs(stats::cor(r1$ind_coord[["dim1"]], r2$ind_coord$dim1)), 0.9)
  ## and each variable's contribution sums to one within a dimension
  agg <- tapply(r2$var_contrib$sqload, r2$var_contrib$dim, sum)
  expect_equal(unname(as.numeric(agg)), rep(1, 2L), tolerance = 1e-8)
  ## clustering and profiling do not care which produced the coordinates
  cl <- ilm_cluster(r2, k_max = 4L, B = 20L, seed = 1L)
  expect_s3_class(cl, "ilm_cluster")
  expect_equal(nrow(cl$ind_cluster), nrow(d))
  pr <- ilm_profile(d, ndim = 2L, method = "glrm", k_max = 4L, B = 20L,
                    seed = 1L, progress = FALSE)
  expect_equal(pr$reduce$method, "glrm")
})

test_that("ilm_impute can use it, and then categories get imputed too", {
  d <- glrm_mixed(n = 300L, seed = 6L)
  h <- hide(d, seed = 5L)
  im <- ilm_impute(h$dm, m = 5L, method = "glrm", verbose = FALSE,
                   progress = FALSE)
  expect_s3_class(im, "ilm_mids")
  expect_equal(im$method, "glrm")
  expect_equal(im$m, 5L)
  expect_false(anyNA(im$imputations[[1]]))
  expect_s3_class(im$imputations[[1]]$g, "factor")
  ## the m datasets are DRAWS, not copies of one reconstruction: identical
  ## imputations would make the pooled variance that of a single fit
  expect_false(identical(im$imputations[[1]], im$imputations[[2]]))
  expect_gt(mean(im$imputations[[1]]$g[h$mi[, 4]] == d$g[h$mi[, 4]]), 0.45)
  ## the low-rank route leaves the factor alone, which is the gap this fills
  il <- ilm_impute(h$dm, m = 3L, method = "lowrank", verbose = FALSE,
                   progress = FALSE)
  expect_true(anyNA(il$imputations[[1]]$g))
})

test_that("it says what it cannot do", {
  d <- glrm_mixed(n = 200L, seed = 8L)
  expect_error(ilm_glrm(d, rank = 99L), "low-rank model needs a rank below")
  expect_error(ilm_glrm(d, loss = c(a = "nope")), "unknown loss")
  expect_error(ilm_glrm(d, loss = c(zz = "quadratic")), "not being used")
  expect_error(ilm_glrm(d, lambda = -1), "non-negative")
  expect_error(ilm_glrm(d, weights = c(1, 2)), "one non-negative")
  expect_error(ilm_glrm("nope"), "must be a data frame")
  ## a column with no loss to give it is dropped with a note
  dd <- d; dd$when <- Sys.Date() + seq_len(nrow(d))
  expect_message(ilm_glrm(dd, rank = 1L, progress = FALSE), "no loss to give")
})
