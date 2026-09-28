## What a DAG licenses, and what it does not.
##
## A valid adjustment set licenses causal language for the exposure's effect
## and for nothing else. The covariates are in the model to close back-door
## paths; their coefficients are not their effects on the outcome -- here `a`
## has no arrow into `y` at all, yet carries a coefficient of about 0.3 -- and
## reading them as effects is the "Table 2 fallacy". The interpretation once
## said "a affects y".

## a -> x, a -> b -> y, x -> y: {a} and {b} each close the back door
two_sets <- function(seed = 2, n = 300) {
  set.seed(seed)
  a <- rnorm(n); b <- 0.7 * a + rnorm(n); x <- 0.6 * a + rnorm(n)
  d <- data.frame(a = a, b = b, x = x, y = 0.5 * x + 0.6 * b + rnorm(n))
  g <- ilm_dag("dag { x [exposure] ; y [outcome] ; a -> x ; a -> b -> y ; x -> y }")
  ilm_dag_model(g, d, verbose = FALSE)
}

test_that("only the exposure is said as an effect; a covariate is not", {
  m <- two_sets()
  expect_length(m$sets, 2L)
  it <- ilm_interpret(m, ame = FALSE)
  eff <- it$sections$effects
  txt <- paste(unlist(it$sections), collapse = " ")
  expect_match(txt, "x affects y", fixed = TRUE)
  ## no effect sentence for the covariate, and no causal claim about it
  expect_false(any(startsWith(eff, "a:")))
  expect_false(grepl("a affects y", txt, fixed = TRUE))
  expect_false(grepl("per unit of a", txt, fixed = TRUE))
  ## it is named as what it is
  expect_match(txt, "Adjusted for a, to close the back-door paths the graph identifies.",
               fixed = TRUE)
  expect_match(txt, "Its coefficient is not its effect on y and is not interpreted here.",
               fixed = TRUE)
  ## the same when causal language is refused: still the exposure only
  it0 <- ilm_interpret(m, causal = FALSE, ame = FALSE)
  expect_false(any(startsWith(it0$sections$effects, "a:")))
  expect_match(paste(it0$sections$effects, collapse = " "), "x is associated with y",
               fixed = TRUE)
})

test_that("the prose says which adjustment set it shows", {
  m <- two_sets()
  mod <- paste(ilm_interpret(m, ame = FALSE)$sections$model, collapse = " ")
  expect_match(mod, "identified by adjusting for a -- the first of 2 minimal sufficient sets",
               fixed = TRUE)
  expect_match(mod, "every set's estimate is under How far to trust it", fixed = TRUE)
  ## and every set's estimate is there
  cav <- paste(ilm_interpret(m, ame = FALSE)$sections$caveats, collapse = " ")
  expect_match(cav, "2 different adjustment sets identify this effect", fixed = TRUE)

  ## with one set, it is simply the set
  set.seed(3); n <- 300
  z <- rnorm(n); w <- rnorm(n); x <- 0.5 * z + rnorm(n)
  d <- data.frame(x = x, z = z, w = w, y = 0.4 * x + 0.6 * z + 0.3 * w + rnorm(n))
  g <- ilm_dag("dag { x [exposure] ; y [outcome] ; z -> x -> y ; z -> y ; w -> y }")
  m1 <- ilm_dag_model(g, d, verbose = FALSE)
  expect_length(m1$sets, 1L)
  it1 <- ilm_interpret(m1, ame = FALSE)
  mod1 <- paste(it1$sections$model, collapse = " ")
  expect_match(mod1, "a minimal sufficient set under the supplied causal graph.", fixed = TRUE)
  expect_false(grepl("first of", mod1, fixed = TRUE))
  expect_false(any(startsWith(it1$sections$effects, "z:")))
})

test_that("an exposure in an interaction keeps its interaction's sentence", {
  ## the exposure's terms are every term it appears in
  set.seed(4); n <- 400
  z <- rnorm(n); x <- 0.5 * z + rnorm(n)
  d <- data.frame(x = x, z = z, y = 0.4 * x + 0.6 * z + 0.2 * x * z + rnorm(n))
  g <- ilm_dag("dag { x [exposure] ; y [outcome] ; z -> x -> y ; z -> y }")
  m <- ilm_dag_model(g, d, verbose = FALSE)
  ## the DAG fits no interaction itself; put one in its place
  m$fits[[1]] <- ilm_model(y ~ x * z, data = d, family = "gaussian",
                           verbose = FALSE)
  eff <- ilm_interpret(m, ame = FALSE)$sections$effects
  expect_true(any(grepl("^x:z", eff)))
  expect_false(any(startsWith(eff, "z:")))
})
