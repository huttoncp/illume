# ilm_draws(): joint draws of every parameter.

drw_data <- function(seed = 3) {
  set.seed(seed)
  d <- data.frame(id = factor(rep(sprintf("s%02d", 1:15), each = 8)),
                  t = rep(0:7, 15))
  d$x <- stats::rnorm(nrow(d))
  d$y <- 1 + 0.5 * d$x + stats::rnorm(15)[d$id] +
    stats::rnorm(15, 0, 0.3)[d$id] * d$t + stats::rnorm(nrow(d))
  d
}

test_that("the draws are centred on the fit, row for row", {
  ## the check agreed with the package the draws were written for: every
  ## random effect's mode is the draws' centre at its row -- by ML and by
  ## REML, with a random slope, where reading by position would shift it
  d <- drw_data()
  for (reml in c(FALSE, TRUE)) {
    f <- ilm_model(y ~ x + (1 + t | id), data = d, family = "gaussian",
                   reml = reml, verbose = FALSE)
    dr <- ilm_draws(f, nsim = 5, seed = 1)
    r <- ilm_ranef(f)
    expect_identical(r$mode, unname(dr$mode[r$row]))
    ## the map names every row, and the fixed effects by their coefficients
    b <- dr$map[dr$map$block == "beta", ]
    expect_identical(b$term, names(coef(f)))
    rb <- dr$map[dr$map$block == "bvec", ]
    expect_identical(rb$level, as.character(r$level))
    expect_identical(rb$dim, r$dim)
  }
})

test_that("the fixed effects' draws have the fit's covariance", {
  d <- drw_data()
  f <- ilm_model(y ~ x + (1 | id), data = d, family = "gaussian",
                 verbose = FALSE)
  dr <- ilm_draws(f, nsim = 4000, seed = 2, natural = FALSE)
  B <- dr$draws[dr$map$block == "beta", , drop = FALSE]
  V <- vcov(f)
  expect_equal(unname(stats::cov(t(B))), unname(V), tolerance = 0.1)
  se <- sqrt(diag(V))
  expect_true(all(abs(rowMeans(B) - coef(f)) < 4 * se / sqrt(4000)))
  ## the same seed gives the same draws
  expect_identical(ilm_draws(f, nsim = 20, seed = 7, natural = FALSE)$draws,
                   ilm_draws(f, nsim = 20, seed = 7, natural = FALSE)$draws)
})

test_that("given = 'theta' holds the variance parameters", {
  d <- drw_data()
  f <- ilm_model(y ~ x + (1 + t | id), data = d, family = "gaussian",
                 verbose = FALSE)
  dr <- ilm_draws(f, nsim = 50, seed = 3, given = "theta")
  th <- dr$draws[dr$map$block == "theta", , drop = FALSE]
  expect_true(all(th == dr$mode[dr$map$block == "theta"]))
  expect_gt(stats::sd(dr$draws[dr$map$block == "beta", ][1, ]), 0)
  ## and so the variance components do not move
  expect_equal(stats::sd(dr$natural$re$id[1, 1, ]), 0)
  expect_identical(dr$given, "theta")
})

test_that("each draw's variance components are the fit's transform of it", {
  d <- drw_data()
  f <- ilm_model(y ~ x + (1 + t | id), data = d, family = "gaussian",
                 ar = ilm_rw1(~ t | id), verbose = FALSE)
  dr <- ilm_draws(f, nsim = 30, seed = 4)
  expect_identical(dim(dr$natural$re$id), c(2L, 2L, 30L))
  expect_length(dr$natural$latent$var_per_time[1, 1, ], 30L)
  expect_identical(dr$natural$latent$meaning,
                   attr(ilm_cells(f), "parameterisation"))
  ## at the estimate the transform is exactly ilm_varcorr()'s
  full <- illume:::ilm_full_par(f)
  at <- illume:::ilm_natural(f, illume:::ilm_coef_order(f, full))
  v <- ilm_varcorr(f)
  expect_equal(at$re$id[, ], v$re$id[, ], tolerance = 1e-12)
  expect_equal(at$latent$var_per_time[, ], v$latent$var_per_time[, ],
               tolerance = 1e-12)
  expect_equal(at$dispersion$value, v$dispersion$value, tolerance = 1e-12)
  ## the walk's cells are mapped to ilm_cells()' rows
  ba <- dr$map[dr$map$block == "B_ar", ]
  cl <- ilm_cells(f)
  expect_identical(ba$cell, which(!cl$anchor))
})

test_that("a boundary's held direction is held in every draw", {
  set.seed(4); n <- 300
  d <- data.frame(x = stats::rnorm(n),
                  g = factor(sample(c("north", "central", "south"), n, TRUE)),
                  site = factor(sample(30, n, TRUE)))
  eta <- cbind(0.4 * d$x, -0.3 * d$x + 0.5 * (d$g == "south"))
  P <- exp(cbind(eta, 0)); P <- P / rowSums(P)
  d$y <- factor(apply(P, 1, function(p) sample(c("lo", "mid", "hi"), 1, prob = p)),
                levels = c("lo", "mid", "hi"))
  f <- suppressMessages(ilm_model(y ~ x + g + (1 | site), data = d,
                                  family = "multinomial", verbose = FALSE))
  expect_identical(f$hessian_how, "boundary")
  expect_false(is.null(f$hessian_dirs))
  dr <- ilm_draws(f, nsim = 200, seed = 5, natural = FALSE)
  expect_equal(dr$held$n, ncol(f$hessian_dirs))
  expect_identical(dr$held$terms, "site")
  ## along each held direction the draws do not move
  fp <- setdiff(seq_along(dr$mode), f$obj$env$random)
  dev <- dr$draws[fp, , drop = FALSE] - dr$mode[fp]
  expect_lt(max(abs(crossprod(f$hessian_dirs, dev))), 1e-8)
  ## and everything else still does
  expect_gt(stats::sd(dr$draws[dr$map$block == "beta", ][1, ]), 0)
})

test_that("a fit without the precision stored forms it, and blocks subset", {
  d <- drw_data()
  f0 <- ilm_model(y ~ x + (1 | id), data = d, family = "gaussian",
                  joint = FALSE, verbose = FALSE)
  f1 <- ilm_model(y ~ x + (1 | id), data = d, family = "gaussian",
                  joint = TRUE, verbose = FALSE)
  expect_null(f0$jointPrecision)
  expect_equal(as.matrix(illume:::ilm_joint_prec(f0)),
               as.matrix(f1$jointPrecision), tolerance = 1e-6)
  dr <- ilm_draws(f0, nsim = 10, seed = 6, blocks = c("beta", "bvec"))
  expect_identical(unique(rownames(dr$draws)), c("beta", "bvec"))
  expect_identical(nrow(dr$map), nrow(dr$draws))
  expect_error(ilm_draws(f0, blocks = "nope"), "no block called")
  expect_output(print(dr), "joint draws")
})
