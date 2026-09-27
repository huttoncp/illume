# A seeded function puts the user's random stream back as it leaves.

rng_fit <- function() {
  set.seed(3); n <- 240
  d <- data.frame(id = factor(rep(1:20, each = 12)), x = stats::rnorm(n),
                  g = factor(sample(c("north", "central", "south"), n, TRUE)))
  d$y <- stats::rpois(n, exp(0.3 + 0.4 * d$x + stats::rnorm(20, 0, 0.5)[d$id]))
  ilm_model(y ~ x + g + (1 | id), data = d, family = "poisson", verbose = FALSE)
}

## the stream after a call, against the stream after the same user code
## without it
untouched <- function(call) {
  set.seed(99); stats::runif(1); a <- .Random.seed
  set.seed(99); stats::runif(1); force(call()); b <- .Random.seed
  identical(a, b)
}

test_that("every seeded function starts by arranging to restore the stream", {
  ## the list is every export with a `seed` argument and every method with
  ## one, so a new seeded function cannot be added without it
  ns <- asNamespace("illume")
  ex <- getNamespaceExports("illume")
  ex <- ex[!startsWith(ex, "iml_")]
  meth <- ls(ns, all.names = TRUE)
  meth <- meth[grepl("[.]ilm_", meth)]
  cand <- unique(c(ex, meth))
  seeded <- Filter(function(f) {
    fn <- get(f, envir = ns)
    is.function(fn) && "seed" %in% names(formals(fn)) &&
      !any(grepl("UseMethod", deparse(body(fn)), fixed = TRUE))
  }, cand)
  expect_gt(length(seeded), 25L)
  first <- vapply(seeded, function(f) {
    b <- body(get(f, envir = ns))
    deparse(if (is.call(b) && identical(b[[1L]], as.name("{"))) b[[2L]] else b)[1L]
  }, "")
  bad <- seeded[!startsWith(first, "ilm_rng_restore(seed")]
  expect_identical(bad, character(0),
                   label = "seeded functions that do not restore the stream")
})

test_that("the user's stream is untouched and the results reproducible", {
  f <- rng_fit()
  nd <- f$model[1:5, ]
  expect_true(untouched(function() predict(f, newdata = nd, groups = "population")))
  expect_true(untouched(function() ilm_simulate(f, nsim = 2, seed = 4)))
  expect_true(untouched(function() ilm_draws(f, nsim = 20, seed = 5)))
  expect_true(untouched(function() ilm_rqr(f)))
  expect_true(untouched(function() ilm_contrast(ilm_emmeans(f, "g"), nsim = 2000L)))
  expect_true(untouched(function() suppressWarnings(ilm_check_zeros(f, B = 20))))
  expect_true(untouched(function() ilm_scenario(f, x = c(-1, 1))))
  ## and a seeded result does not depend on the user's stream
  set.seed(1); r1 <- ilm_rqr(f)
  set.seed(2); r2 <- ilm_rqr(f)
  expect_identical(r1, r2)
})

test_that("an absent stream stays absent, and seed = NULL draws from it", {
  f <- rng_fit()
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
    rm(".Random.seed", envir = globalenv())
  invisible(ilm_simulate(f, nsim = 1, seed = 6))
  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
  ## with no seed the function is ordinary R code: it moves the stream on
  set.seed(7); a <- ilm_simulate(f, nsim = 1); b <- ilm_simulate(f, nsim = 1)
  expect_false(identical(a, b))
})

test_that("the stream is put back even when the function stops", {
  set.seed(8); stats::runif(1); a <- .Random.seed
  set.seed(8); stats::runif(1)
  expect_error(ilm_surv(time = c(1, 2, 3)), "`event` is required")
  expect_identical(.Random.seed, a)
})
