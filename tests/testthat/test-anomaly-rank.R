# The anomaly scan's rank comes from parallel analysis, not from the
# cross-validation illume's low-rank imputation uses. The scan is illumex's,
# so its selector is reached with :::; the comparison needs both packages.

## data with a k-dimensional shared structure, plus rows pushed OFF it. The
## push is orthogonal to the structure, so no single column is extreme -- which
## is the case this function exists for and the one ilm_outliers() cannot see.
anom_data <- function(n = 400L, p = 8L, k = 2L, seed = 1L, contam = 0.03,
                      sd_e = 0.4, shared = FALSE) {
  set.seed(seed)
  F <- matrix(rnorm(n * k), n, k)
  L <- matrix(rnorm(k * p), k, p)
  X <- F %*% L + matrix(rnorm(n * p, 0, sd_e), n, p)
  bad <- if (contam > 0) sample(n, max(1L, round(contam * n))) else integer(0)
  if (length(bad)) {
    P <- t(L) %*% solve(L %*% t(L)) %*% L
    if (shared) {
      v <- rnorm(p); v <- as.numeric(v - P %*% v); v <- v / sqrt(sum(v^2))
      X[bad, ] <- X[bad, ] +
        outer(3.2 * runif(length(bad), .8, 1.2), v * sqrt(p) * sd_e)
    } else for (i in bad) {
      v <- rnorm(p); v <- as.numeric(v - P %*% v)
      X[i, ] <- X[i, ] + 3.2 * v / sqrt(sum(v^2)) * sqrt(p) * sd_e
    }
  }
  list(d = as.data.frame(X), bad = sort(bad))
}

test_that("the anomaly rank is not the imputation's cross-validated rank", {
  ## the CV selector optimises held-out CELL prediction and chose 6 or 7 on a
  ## rank-2 structure, which spans the directions the anomalies depart along
  g <- anom_data(seed = 7L, contam = 0)
  Z <- scale(as.matrix(g$d))
  expect_equal(illumex:::ilm_anom_rank(Z), 2L)
  cv <- as.integer(ilm_lowrank_ncp(as.matrix(g$d),
                                   matrix(FALSE, nrow(g$d), ncol(g$d)),
                                   ncp_max = 7L, seed = 1L))
  expect_gt(cv, 3L)
})
