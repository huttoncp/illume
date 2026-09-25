# ilm_matrices(): the designs of a fitted model for new rows.

mat_data <- function(seed = 3) {
  set.seed(seed)
  d <- data.frame(id = factor(rep(sprintf("s%02d", 1:15), each = 8)),
                  t = rep(0:7, 15))
  d$x <- stats::rnorm(nrow(d))
  d$y <- 1 + 0.5 * d$x + stats::rnorm(15)[d$id] +
    stats::rnorm(15, 0, 0.3)[d$id] * d$t + stats::rnorm(nrow(d))
  d
}

test_that("the designs are the fit's, and groups are matched by label", {
  d <- mat_data()
  f <- ilm_model(y ~ x + (1 + t | id), data = d, family = "gaussian",
                 verbose = FALSE)
  nd <- data.frame(x = c(0.5, -1, 2), t = c(2, 9, 4),
                   id = c("s03", "zz", "s15"))
  m <- ilm_matrices(f, nd)
  expect_equal(matrix(m$X, nrow(nd)), cbind(1, nd$x))
  expect_identical(colnames(m$X), colnames(f$X))
  expect_equal(matrix(m$re$id$Z, nrow(nd)), cbind(1, nd$t))
  expect_identical(m$re$id$group, c(3L, NA, 15L))
  expect_identical(m$re$id$new_group, c(FALSE, TRUE, FALSE))
  expect_identical(m$re$id$level, nd$id)
  expect_identical(m$re$id$factor, "id")
  ## the typical group's prediction is X beta, as predict() gives it
  expect_equal(as.numeric(m$X %*% f$beta),
               unname(predict(f, newdata = nd, type = "link")[, 1]),
               tolerance = 1e-12)
  ## and a fitted group's is X beta + Z b, with b its modes from ilm_ranef()
  r <- ilm_ranef(f)
  B <- matrix(r$mode, ncol = 2)                   # levels by (intercept, t)
  i <- 1L; g <- m$re$id$group[i]
  expect_equal(sum(m$X[i, ] * f$beta) + sum(m$re$id$Z[i, ] * B[g, ]),
               sum(nd$x[i] * f$beta[2], f$beta[1]) + B[g, 1] + nd$t[i] * B[g, 2],
               tolerance = 1e-12)
})

test_that("a smooth's basis rebuilds predict()'s prediction exactly", {
  d <- mat_data()
  f <- ilm_model(y ~ s(x) + (1 | id), data = d, family = "gaussian",
                 verbose = FALSE)
  nd <- data.frame(x = c(-1.5, 0, 0.7, 1.9), id = "s01")
  m <- ilm_matrices(f, nd)
  expect_length(m$smooth, 1L)
  r <- ilm_ranef(f)
  bs <- r$mode[r$type == "smooth"]
  eta <- as.numeric(m$X %*% f$beta + m$smooth[[1]] %*% bs)
  expect_equal(eta, unname(predict(f, newdata = nd, type = "link")[, 1]),
               tolerance = 1e-10)
})

test_that("new rows are placed among a walk's cells", {
  d <- mat_data()
  d2 <- d[!(d$id == "s01" & d$t %in% c(3, 4)), ]      # a gap in one group
  f <- ilm_model(y ~ x, data = d2, family = "gaussian",
                 ar = ilm_rw1(~ t | id), verbose = FALSE)
  cl <- ilm_cells(f)
  nd <- data.frame(x = 0, t = c(2, 3.5, 7, 10, -1, 5),
                   id = c("s01", "s01", "s01", "s01", "s01", "new"))
  a <- ilm_matrices(f, nd)$ar
  s01 <- which(cl$group == "s01")                     # times 0, 1, 2, 5, 6, 7
  expect_identical(a$cell, c(s01[3], NA, s01[6], NA, NA, NA))
  ## between the cells either side of the gap
  expect_identical(c(a$prev_cell[2], a$next_cell[2]), c(s01[3], s01[4]))
  expect_equal(c(a$dt_prev[2], a$dt_next[2]), c(1.5, 1.5))
  ## past the last cell: the forecast case
  expect_identical(c(a$prev_cell[4], a$next_cell[4]), c(s01[6], NA))
  expect_equal(a$dt_prev[4], 3)
  ## before the first
  expect_identical(c(a$prev_cell[5], a$next_cell[5]), c(NA, s01[1]))
  expect_equal(a$dt_next[5], 1)
  ## a group the fit has not seen
  expect_true(a$new_group[6])
  expect_true(all(is.na(unlist(a[6, c("cell", "prev_cell", "next_cell")]))))
  expect_identical(a$time, nd$t)
})

test_that("an AR(1) grid places on its steps and refuses between them", {
  d <- mat_data()
  f <- suppressWarnings(ilm_model(y ~ x, data = d, family = "gaussian",
                                  ar = ilm_ar1(~ t | id), verbose = FALSE))
  cl <- ilm_cells(f)
  nd <- data.frame(x = 0, t = c(3, 9), id = "s02")
  a <- ilm_matrices(f, nd)$ar
  s02 <- which(cl$group == "s02")
  expect_identical(a$cell[1], s02[4])
  ## two steps past the grid's end
  expect_identical(c(a$prev_cell[2], a$next_cell[2]), c(s02[8], NA))
  expect_equal(a$dt_prev[2], 2)
  expect_error(ilm_matrices(f, data.frame(x = 0, t = 2.5, id = "s02")),
               "between the AR\\(1\\) grid's steps")
})

test_that("a structure built from vectors needs the times and groups given", {
  d <- mat_data()
  f <- ilm_model(y ~ x, data = d, family = "gaussian",
                 ar = ilm_rw1(d$t, d$id), verbose = FALSE)
  nd <- data.frame(x = 0, when = 4, who = "s02")
  expect_error(ilm_matrices(f, nd), "pass `time` and `group`")
  a <- ilm_matrices(f, nd, time = nd$when, group = nd$who)$ar
  expect_false(is.na(a$cell))
})

test_that("the zero part's and the dispersion model's designs come too", {
  set.seed(8); n <- 400
  d <- data.frame(x = stats::runif(n, 0, 3), z = stats::rnorm(n),
                  g = factor(sample(c("a", "b"), n, TRUE)))
  d$y <- ifelse(stats::runif(n) < stats::plogis(-1 + 0.6 * d$z), 0,
                stats::rnbinom(n, mu = exp(0.5 + 0.3 * d$x), size = 2))
  f <- ilm_model(y ~ x, data = d, family = "nbinom", ziformula = ~ z,
                 dispformula = ~ g, verbose = FALSE)
  nd <- data.frame(x = c(1, 2), z = c(0, 1), g = c("b", "a"))
  m <- ilm_matrices(f, nd)
  expect_identical(colnames(m$zi), paste0("zi:", colnames(f$Zzi)))
  expect_equal(matrix(m$zi, nrow(nd)), cbind(1, nd$z))
  expect_identical(colnames(m$disp), paste0("disp:", colnames(f$Zd)))
  expect_equal(unname(m$disp[, "disp:gb"]), c(1, 0))
  ## each column finds its coefficient by name, in coef() and in the draws
  b <- stats::coef(f, full = TRUE)
  expect_true(all(c(colnames(m$zi), colnames(m$disp)) %in% names(b)))
  mp <- ilm_draws(f, nsim = 2, seed = 1)$map
  expect_true(all(c(colnames(m$zi), colnames(m$disp)) %in% mp$term))
  ## and the zero part's linear predictor for the new rows is the fit's own
  eta_zi <- as.vector(m$zi %*% b[colnames(m$zi)])
  expect_equal(eta_zi, as.vector(cbind(1, nd$z) %*% f$zi_gamma),
               tolerance = 1e-12)
})

test_that("each Z's columns are named by dimension, a lone intercept too", {
  ## unnamed for a random intercept alone, named with a slope: code reading
  ## columns by name found no intercept column and dropped the group effects
  d <- mat_data()
  f1 <- ilm_model(y ~ x + (1 | id), data = d, family = "gaussian",
                  verbose = FALSE)
  f2 <- ilm_model(y ~ x + (1 + t | id), data = d, family = "gaussian",
                  verbose = FALSE)
  nd <- data.frame(x = 0, t = 3, id = "s02")
  z1 <- ilm_matrices(f1, nd)$re$id$Z
  z2 <- ilm_matrices(f2, nd)$re$id$Z
  expect_identical(colnames(z1), "(Intercept)")
  expect_identical(colnames(z2), c("(Intercept)", "t"))
  ## the same names ilm_ranef() and the draws map give the dimensions
  expect_setequal(colnames(z1), unique(ilm_ranef(f1)$dim))
  expect_setequal(colnames(z2), unique(ilm_ranef(f2)$dim))
  mp <- ilm_draws(f2, nsim = 1, seed = 1, natural = FALSE)$map
  expect_setequal(colnames(z2), unique(mp$dim[mp$block == "bvec"]))
})
