## The gradient check reports the gradient along a held direction but does
## not judge it: a term held at its boundary estimate sits where the
## likelihood is flat or running off to a limit, and its estimate is not
## interpreted. The rest of the gradient says whether the fit converged.

test_that("the held part of the gradient is split off and reported", {
  ## the forecasting package's 38,400-row fit: logdisp held (k at 3.7e8),
  ## its gradient 1.41e-2, lchol_ar's -6.42e-3, the betas at most 2.8e-4
  gv <- c(beta = 2.8e-4, beta = -1.2e-5, theta = -2.4e-5, lchol_ar = -6.42e-3, logdisp = 1.41e-2)
  D <- diag(1, length(gv))[, 5, drop = FALSE]
  gs <- illume:::ilm_grad_judged(gv, D)
  expect_equal(gs$held, 1.41e-2)
  expect_equal(gs$judged, 6.42e-3)          # WARN, not FAIL, on today's lines
  ## nothing held: all of it is judged
  gs0 <- illume:::ilm_grad_judged(gv, NULL)
  expect_equal(gs0$judged, 1.41e-2)
  expect_true(is.na(gs0$held))
  ## a held direction that is not a single parameter: a rotated flat direction
  u <- c(0, 0, 1, -1, 0) / sqrt(2)
  gv2 <- c(0, 0, 0.02, -0.02, 1e-4)
  gs2 <- illume:::ilm_grad_judged(gv2, matrix(u, ncol = 1))
  expect_equal(gs2$held, 0.02, tolerance = 1e-12)
  expect_equal(gs2$judged, 1e-4, tolerance = 1e-12)
  ## a gradient that is not finite is judged whole, as far from a stationary
  ## point as any
  expect_identical(illume:::ilm_grad_judged(c(NA, 1e-4), D[1:2, , drop = FALSE])$judged, NA_real_)
})

test_that("a fit whose held dispersion carries the largest gradient is judged on the rest", {
  ## Poisson counts fitted as a negative binomial with a random walk: k runs
  ## off and is held, and here the largest component of the gradient is
  ## along it
  set.seed(28)
  G <- 8; Tn <- 24
  d <- expand.grid(t = seq_len(Tn), g = factor(seq_len(G)))
  lvl <- rnorm(G, 3, 1.5)
  rw <- unlist(lapply(seq_len(G), function(i) cumsum(rnorm(Tn, 0, 0.1))))
  d$s12 <- sin(2 * pi * d$t / 12); d$c12 <- cos(2 * pi * d$t / 12)
  d$y <- rpois(nrow(d), exp(lvl[d$g] + 0.3 * d$s12 + rw))
  f <- suppressMessages(suppressWarnings(ilm_model(y ~ s12 + c12 + (1 | g), data = d,
         family = "nbinom", ar = ilm_rw1(~ t | g), verbose = FALSE)))
  expect_true("dispersion" %in% f$hessian_held)
  gr <- as.numeric(f$obj$gr(f$opt$par)); pn <- names(f$opt$par)
  gd <- abs(gr[pn == "logdisp"]); go <- max(abs(gr[pn != "logdisp"]))
  expect_gt(gd, go)                          # the held part is the largest
  r <- f$checks[f$checks$check == "gradient", ]
  expect_match(r$detail, sprintf("max |gradient| = %.2e", go), fixed = TRUE)
  expect_match(r$detail, sprintf("along the held dispersion: %.2e, not judged", gd), fixed = TRUE)
})

test_that("a fit with nothing held is judged on its whole gradient, as before", {
  set.seed(1)
  d <- data.frame(g = factor(rep(1:20, each = 10)), x = rnorm(200))
  d$y <- 1 + 0.5 * d$x + rnorm(20)[d$g] + rnorm(200)
  f <- ilm_model(y ~ x + (1 | g), data = d, family = "gaussian", verbose = FALSE)
  r <- f$checks[f$checks$check == "gradient", ]
  expect_identical(r$status, "OK")
  expect_false(grepl("held", r$detail))
  expect_match(r$detail, sprintf("max |gradient| = %.2e", max(abs(f$obj$gr(f$opt$par)))), fixed = TRUE)
})
