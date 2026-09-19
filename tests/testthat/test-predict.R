test_that("predict returns probabilities for every row and category", {
  fit <- fit_basic()
  P <- predict(fit)
  expect_equal(dim(P), c(nrow(fit$X), fit$J))
  expect_equal(colnames(P), fit$ylevels)
  expect_true(all(abs(rowSums(P) - 1) < 1e-9))
  expect_true(all(P >= 0 & P <= 1))
})

test_that("link predictions sum to zero across categories", {
  fit <- fit_basic()
  L <- predict(fit, type = "link")
  expect_true(all(abs(rowSums(L)) < 1e-9))   # sum-to-zero coding
})

test_that("type = class returns the most likely category", {
  fit <- fit_basic()
  cl <- predict(fit, type = "class")
  expect_s3_class(cl, "factor")
  expect_equal(levels(cl), fit$ylevels)
  P <- predict(fit)
  expect_equal(as.character(cl), fit$ylevels[max.col(P, ties.method = "first")])
})

test_that("new data reproduces the fitted design exactly", {
  fit <- fit_basic()
  nd <- ilm_newX(fit, fit$model)
  expect_equal(unname(nd$X), unname(fit$X))
})

test_that("predictions at new data have the right shape", {
  fit <- fit_basic()
  grid <- data.frame(x1 = c(-1, 0, 1), grp = factor("a", levels = levels(fit$model$grp)))
  P <- predict(fit, newdata = grid)
  expect_equal(dim(P), c(3L, fit$J))
  expect_true(all(abs(rowSums(P) - 1) < 1e-9))
})

test_that("population-averaged differs from conditional and stays valid", {
  fit <- fit_basic()
  pc <- predict(fit, type = "response", marginal = FALSE)
  pm <- predict(fit, type = "response", marginal = TRUE, ndraw = 60)
  expect_true(all(abs(rowSums(pm) - 1) < 1e-9))
  # softmax(E[eta]) != E[softmax(eta)]: they must not coincide
  expect_gt(max(abs(pc - pm)), 1e-4)
})

test_that("predictions are reproducible across calls", {
  # marginaleffects computes derivatives numerically, so identical inputs must
  # give identical outputs or the derivative is noise
  fit <- fit_basic()
  a <- predict(fit, type = "response", marginal = TRUE, ndraw = 40)
  b <- predict(fit, type = "response", marginal = TRUE, ndraw = 40)
  expect_equal(a, b)
})

test_that("intervals stay inside [0, 1] and narrow with the level", {
  fit <- fit_basic()
  p95 <- predict(fit, se.fit = TRUE, level = 0.95, nsim = 40)
  expect_true(all(p95$lower >= 0 & p95$upper <= 1))
  expect_true(all(p95$se.fit > 0))
  p80 <- predict(fit, se.fit = TRUE, level = 0.80, nsim = 40)
  expect_true(all((p80$upper - p80$lower) <= (p95$upper - p95$lower) + 1e-8))
})

test_that("marginal predictions are refused on the link scale", {
  fit <- fit_basic()
  expect_error(predict(fit, type = "link", marginal = TRUE), "response scale")
})

test_that("fitted probabilities agree with predict", {
  fit <- fit_basic()
  expect_equal(unname(ilm_fitted(fit, conditional = FALSE)),
               unname(predict(fit, type = "response")))
})
