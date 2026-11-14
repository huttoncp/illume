## predict() on a factor-smooth interaction, s(x, id, bs = "fs"). For this
## smooth mgcv::smooth2random() returns no trans.U, only trans.D and `rind`,
## and predict() stopped on every fs fit, with or without new data. The new
## data's design is now (PredictMat(sm, data) %*% diag(trans.D))[, order(rind)],
## an identity found by measurement with mgcv 1.9.1 and checked here first.
##
## Compared with mgcv: mgcv::PredictMat() and mgcv::smooth2random() for the
## identity, and mgcv::gam() with method = "REML" and its predict() method,
## mgcv::predict.gam(), for the predictions at new x (mgcv 1.9.1 when
## written).

fsp_data <- function(seed, n_id = 30, n_per = 10, unbal = FALSE) {
  set.seed(seed)
  sizes <- if (unbal) sample(3:18, n_id, replace = TRUE) else rep(n_per, n_id)
  id <- factor(rep(seq_len(n_id), times = sizes))
  x <- stats::runif(length(id))
  a <- stats::rnorm(n_id, 0, 0.6)[id]; b <- stats::rnorm(n_id, 0, 0.8)[id]
  w <- stats::rnorm(n_id, 0, 0.3)[id]
  y <- sin(2 * pi * x) + a + b * (x - 0.5) + w * cos(3 * pi * x) +
    stats::rnorm(length(x), 0, 0.4)
  data.frame(y = y, x = x, id = id)
}

fsp_cases <- function() list(
  "30 levels, k = 5"              = list(d = fsp_data(1), form = y ~ s(x, id, bs = "fs", k = 5)),
  "5 levels, k = 4"               = list(d = fsp_data(8, 5, 12), form = y ~ s(x, id, bs = "fs", k = 4)),
  "k = 3"                         = list(d = fsp_data(3), form = y ~ s(x, id, bs = "fs", k = 3)),
  "beside s(x), k = 8, m = 1"     = list(d = fsp_data(4), form = y ~ s(x) + s(x, id, bs = "fs", k = 8, m = 1)),
  "unbalanced, k = 6, cr basis"   = list(d = fsp_data(5, unbal = TRUE),
                                         form = y ~ s(x, id, bs = "fs", k = 6, xt = list(bs = "cr"))))

fsp_fit <- function(z) suppressMessages(suppressWarnings(
  ilm_model(z$form, data = z$d, family = "gaussian", reml = TRUE, verbose = FALSE)))
fsp_sob <- function(f) {
  i <- vapply(f$smooths, function(s) inherits(s$sm, "fs.interaction"), TRUE)
  f$smooths[[which(i)[1L]]]
}

test_that("the design rebuilt from the fitting rows is the fit's own (the identity)", {
  for (z in fsp_cases()) {
    f <- fsp_fit(z); sob <- fsp_sob(f)
    expect_null(sob$re$trans.U)
    rebuilt <- do.call(cbind, ilm_smooth_design(sob, f$model)$Xr)
    expect_lt(max(abs(rebuilt - do.call(cbind, sob$rand))), 1e-10)
  }
})

test_that("predict() on an fs fit, without new data and with its own rows, is the fit", {
  for (z in fsp_cases()) {
    f <- fsp_fit(z)
    fv <- as.numeric(ilm_fitted(f))
    expect_equal(as.numeric(predict(f)), fv, tolerance = 1e-8)
    expect_equal(as.numeric(predict(f, newdata = z$d)), fv, tolerance = 1e-8)
  }
})

test_that("predict() at new x for fitted subjects agrees with mgcv::predict.gam()", {
  skip_on_cran()
  for (z in fsp_cases()) {
    f <- fsp_fit(z)
    g <- suppressWarnings(mgcv::gam(z$form, data = z$d, method = "REML"))
    lv <- levels(z$d$id)
    nd <- expand.grid(x = seq(0.03, 0.97, length.out = 9),
                      id = lv[unique(c(1L, 2L, length(lv)))])
    nd$id <- factor(as.character(nd$id), levels = lv)
    pg <- as.numeric(stats::predict(g, newdata = nd))
    pf <- as.numeric(predict(f, newdata = nd))
    ## the two fits agree to their optimisers' tolerance, not exactly
    expect_lt(max(abs(pf - pg)), 1e-3 * diff(range(pg)))
  }
})
