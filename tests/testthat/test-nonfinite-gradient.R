## A gradient that is not finite gives a fit whose optimizer check FAILs, not
## an error. Seen on CI (macOS arm64) for a beta whose phi ran to its limit;
## not reproducible on every machine, so the non-finite gradient is injected.

## a quadratic whose optimum at (3, 1) lies where the gradient is NaN (x > 2)
toy <- function(bad = function(p, n) p[1] > 2) {
  n <- 0
  list(fn = function(p) sum((p - c(3, 1))^2),
       gr = function(p) {
         n <<- n + 1
         g <- 2 * (p - c(3, 1))
         if (bad(p, n)) g[] <- NaN
         g
       })
}
ctl <- list(iter.max = 3000, eval.max = 3000)

test_that("a gradient that breaks down is stopped at the best finite point", {
  obj <- toy()
  st <- c(a = 0, b = 0)
  ## what nlminb itself does
  expect_error(stats::nlminb(st, obj$fn, obj$gr, control = ctl), "NA/NaN gradient")
  o <- ilm_nlminb(st, obj, ctl)
  expect_true(isTRUE(o$nonfinite))
  expect_identical(o$convergence, 1L)
  expect_named(o$par, c("a", "b"))
  ## a point the gradient was finite at, with its own objective
  expect_lte(o$par[["a"]], 2)
  expect_equal(o$objective, obj$fn(o$par))
  expect_match(o$message, "gradient was not finite")
})

test_that("a restart from the best finite point is kept when it finishes", {
  ## the gradient breaks down once only, on its fourth evaluation
  obj <- toy(function(p, n) n == 4L)
  o <- ilm_nlminb(c(a = 0, b = 0), obj, ctl)
  expect_null(o$nonfinite)
  expect_identical(o$convergence, 0L)
  expect_equal(unname(o$par), c(3, 1), tolerance = 1e-6)
})

test_that("with no finite point at all, the error is nlminb's", {
  obj <- toy(function(p, n) TRUE)
  expect_error(ilm_nlminb(c(a = 0, b = 0), obj, ctl), "NA/NaN gradient")
})

test_that("a fit whose gradient breaks down returns, and its checks say so", {
  set.seed(4)
  d <- data.frame(g = factor(rep(1:15, each = 8)), x = stats::rnorm(120))
  d$y <- 1 + 0.5 * d$x + stats::rnorm(15, 0, 0.7)[d$g] + stats::rnorm(120)
  real <- ilm_nlminb
  ## every optimisation's gradient breaks down after its third evaluation
  local_mocked_bindings(ilm_nlminb = function(start, obj, ctl) {
    n <- 0; gr0 <- obj$gr; o2 <- obj
    o2$gr <- function(p) { n <<- n + 1; g <- gr0(p); if (n > 3) g[] <- NaN; g }
    real(start, o2, ctl)
  })
  f <- suppressWarnings(suppressMessages(
    ilm_model(y ~ x + (1 | g), data = d, verbose = FALSE)))
  expect_true(isTRUE(f$opt$nonfinite))
  ck <- f$checks[f$checks$check == "optimizer", ]
  expect_identical(ck$status, "FAIL")
  expect_match(ck$cause, "arithmetic broke down")
  expect_match(ck$detail, "gradient was not finite")
  ## the remedy is one the package makes -- listed once for every check it
  ## answers -- and the draws warn
  r <- ilm_remedies(f)
  i <- grep("optimizer", r$check, fixed = TRUE)
  expect_length(i, 1L)
  expect_match(r$change[i], "restarts")
  expect_warning(ilm_draws(f, nsim = 20, seed = 1))
})
