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

test_that("it finds rows that no single column flags", {
  g <- anom_data(seed = 3L)
  r <- suppressMessages(ilm_anomaly(g$d, seed = 1L, progress = FALSE))
  expect_s3_class(r, "ilm_anomaly")
  expect_equal(nrow(r), nrow(g$d))
  expect_setequal(names(r), c("row", "score", "p", "p_adj", "flag", "driver"))
  expect_true(all(r$p > 0 & r$p <= 1))
  ## the planted rows are found
  expect_gt(length(intersect(r$row[r$flag], g$bad)) / length(g$bad), 0.8)
  ## and almost nothing else is
  expect_lt(length(setdiff(r$row[r$flag], g$bad)) /
              (nrow(g$d) - length(g$bad)), 0.02)
  ## and the univariate scan does not see them, which is the point. A push
  ## orthogonal to the shared structure can still leave one column on the
  ## large side by chance, so the claim is about what the column-at-a-time
  ## method CATCHES rather than about every value being unremarkable.
  uni <- suppressMessages(ilm_outliers_all(g$d))
  urow <- if (is.data.frame(uni) && "row" %in% names(uni))
    unique(uni$row) else integer(0)
  expect_lt(length(intersect(urow, g$bad)) / length(g$bad), 0.2)
  ## sorted by score, so head() is the shortlist
  expect_false(is.unsorted(rev(r$score)))
})

test_that("it almost never flags anything when there is nothing to flag", {
  ## Over 100 clean datasets of 400 rows the rate is 11 in 40,000, with 95 of
  ## them producing none -- so the assertion is on the RATE across several
  ## datasets, not on every one being empty. One dataset in twenty has a row.
  fl <- 0L; tot <- 0L
  for (s in 1:6) {
    g <- anom_data(seed = 50L + s, contam = 0)
    r <- suppressMessages(ilm_anomaly(g$d, seed = s, progress = FALSE))
    fl <- fl + sum(r$flag); tot <- tot + nrow(r)
  }
  expect_lt(fl / tot, 0.005)
})

test_that("the rank comes from parallel analysis, not cross-validation", {
  ## the CV selector optimises held-out CELL prediction and chose 6 or 7 on a
  ## rank-2 structure, which spans the directions the anomalies depart along
  g <- anom_data(seed = 7L, contam = 0)
  Z <- scale(as.matrix(g$d))
  expect_equal(ilm_anom_rank(Z), 2L)
  cv <- as.integer(ilm_lowrank_ncp(as.matrix(g$d),
                                   matrix(FALSE, nrow(g$d), ncol(g$d)),
                                   ncp_max = 7L, seed = 1L))
  expect_gt(cv, 3L)
  ## and it finds the right rank for a rank-3 structure too
  g3 <- anom_data(k = 3L, p = 10L, seed = 8L, contam = 0)
  expect_equal(ilm_anom_rank(scale(as.matrix(g3$d))), 3L)
  ## the fitted rank is reported and can be overridden
  r <- suppressMessages(ilm_anomaly(g$d, seed = 1L, progress = FALSE))
  expect_equal(attr(r, "rank"), 2L)
  expect_equal(attr(suppressMessages(
    ilm_anomaly(g$d, rank = 4L, seed = 1L, progress = FALSE)), "rank"), 4L)
})

test_that("trimming earns its place when the anomalies share a direction", {
  ## they can form a component between them, and a fit that includes them
  ## reconstructs them well
  d0 <- d1 <- numeric(6)
  for (i in 1:6) {
    g <- anom_data(n = 500L, seed = 300L + i, contam = 0.10, shared = TRUE)
    r0 <- suppressMessages(ilm_anomaly(g$d, trim = 0, seed = i,
                                       progress = FALSE))
    r1 <- suppressMessages(ilm_anomaly(g$d, trim = 0.25, seed = i,
                                       progress = FALSE))
    d0[i] <- length(intersect(r0$row[r0$flag], g$bad)) / length(g$bad)
    d1[i] <- length(intersect(r1$row[r1$flag], g$bad)) / length(g$bad)
  }
  expect_gt(mean(d1), mean(d0))
  ## and it costs nothing when they do not share one
  g <- anom_data(seed = 9L)
  a <- suppressMessages(ilm_anomaly(g$d, trim = 0, seed = 1L, progress = FALSE))
  b <- suppressMessages(ilm_anomaly(g$d, trim = 0.25, seed = 1L,
                                    progress = FALSE))
  expect_gt(stats::cor(a$score[order(a$row)], b$score[order(b$row)]), 0.95)
  expect_error(ilm_anomaly(g$d, trim = 0.7), "in \\[0, 0.5\\)")
})

test_that("the driver names a column that is actually responsible", {
  set.seed(4); n <- 400L
  f <- rnorm(n)
  d <- data.frame(a = f + rnorm(n, 0, .3), b = 2 * f + rnorm(n, 0, .3),
                  c = -f + rnorm(n, 0, .3), e = f + rnorm(n, 0, .3))
  ## a = 1.2 implies the shared factor is about 1.2, so b should be near 2.4
  ## and is -2.4 instead; c and e are set consistent with it. Every value is
  ## well inside its own column's range.
  d[1, ] <- c(1.2, -2.4, -1.2, 1.2)
  r <- suppressMessages(ilm_anomaly(d, seed = 1L, progress = FALSE))
  i <- which(r$row == 1L)
  expect_true(r$flag[i])
  expect_equal(i, 1L)                       # the worst-fitting row
  expect_equal(r$driver[i], "b")
  expect_lt(max(abs(scale(as.matrix(d))[1, ])), 3)
})

test_that("B has to be large enough for a flag to be reachable", {
  g <- anom_data(seed = 11L)
  expect_warning(ilm_anomaly(g$d, B = 9L, progress = FALSE), "cannot be flagged")
  expect_silent(suppressMessages(
    ilm_anomaly(g$d, B = 39L, seed = 1L, progress = FALSE)))
  ## a lone anomaly among many rows needs the larger B: the smallest
  ## Benjamini-Hochberg value is about 1 / B
  set.seed(1); n <- 300L
  f <- rnorm(n)
  d <- data.frame(a = f + rnorm(n, 0, .3), b = 2 * f + rnorm(n, 0, .3),
                  c = -f + rnorm(n, 0, .3))
  d[1, ] <- c(1.2, -2.4, 1.2)
  expect_equal(sum(suppressWarnings(
    ilm_anomaly(d, B = 19L, progress = FALSE))$flag), 0L)
  r <- suppressMessages(ilm_anomaly(d, B = 39L, progress = FALSE))
  expect_equal(sum(r$flag), 1L)
  expect_equal(r$row[1], 1L)
})

test_that("it refuses what it cannot do, and says what to use instead", {
  g <- anom_data(seed = 12L)
  expect_error(ilm_anomaly(g$d[, 1:2]), "at least 3 numeric columns")
  expect_error(ilm_anomaly(mtcars[1:10, ]), "complete rows")
  expect_error(ilm_anomaly("nope"), "must be a data frame")
  ## a categorical column is dropped with a note naming the alternative
  dm <- g$d; dm$grp <- factor(sample(c("a", "b"), nrow(dm), TRUE))
  expect_message(ilm_anomaly(dm, B = 39L, seed = 1L, progress = FALSE),
                 "ilm_reduce")
  r <- suppressMessages(ilm_anomaly(dm, B = 39L, seed = 1L, progress = FALSE))
  expect_equal(attr(r, "dropped"), "grp")
  expect_false("grp" %in% attr(r, "columns"))
})

test_that("rows with missing values are dropped and the index still points home", {
  g <- anom_data(seed = 13L)
  d <- g$d
  d[c(5, 50, 120), 2] <- NA
  r <- suppressMessages(ilm_anomaly(d, seed = 1L, progress = FALSE))
  expect_equal(nrow(r), nrow(d) - 3L)
  expect_false(any(c(5L, 50L, 120L) %in% r$row))
  ## `row` indexes the ORIGINAL data, so a flagged row can be looked up
  expect_true(all(r$row %in% seq_len(nrow(d))))
  expect_equal(nrow(d[r$row[r$flag], , drop = FALSE]), sum(r$flag))
})
