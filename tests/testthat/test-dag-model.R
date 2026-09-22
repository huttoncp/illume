## The workflow's contract, in order of how much it would cost to get wrong:
## the mean structure comes from the graph and is never searched; a
## non-identified effect is reported rather than estimated; every admissible
## set is fitted; and the narration matches what was actually done.

make_conf <- function(n = 600, bx = 0.4, seed = 1) {
  set.seed(seed)
  z <- rnorm(n); x <- 0.5 * z + rnorm(n)
  data.frame(x = x, y = bx * x + 0.6 * z + rnorm(n), z = z)
}
g_conf <- function()
  ilm_dag("dag { x [exposure] ; y [outcome] ; z -> x -> y ; z -> y }")

test_that("the graph's adjustment set becomes the mean structure", {
  d <- make_conf()
  m <- ilm_dag_model(g_conf(), d, verbose = FALSE)
  expect_s3_class(m, "ilm_dag_model")
  expect_true(m$identified)
  expect_equal(m$exposure, "x")
  expect_equal(m$outcome, "y")
  expect_equal(m$family, "gaussian")
  expect_length(m$sets, 1L)
  expect_equal(m$sets[[1]], "z")

  ## the fitted formula is exposure + adjustment set, and nothing else
  tl <- attr(stats::terms(m$fits[[1]]), "term.labels")
  expect_setequal(tl, c("x", "z"))

  ## adjusting recovers the truth where not adjusting does not
  expect_equal(nrow(m$effects), 1L)
  expect_lt(m$effects$lower, 0.4); expect_gt(m$effects$upper, 0.4)
  naive <- unname(coef(stats::lm(y ~ x, d))[2])
  expect_gt(naive, m$effects$upper)
})

test_that("a mediator is never adjusted for a total effect", {
  set.seed(3); n <- 500
  x <- rnorm(n); mm <- 0.7 * x + rnorm(n); y <- 0.6 * mm + rnorm(n)
  d <- data.frame(x = x, m = mm, y = y)
  g <- ilm_dag("dag { x [exposure] ; y [outcome] ; x -> m -> y }")
  m <- ilm_dag_model(g, d, verbose = FALSE)
  expect_setequal(attr(stats::terms(m$fits[[1]]), "term.labels"), "x")
  ## the total effect is 0.7 * 0.6 = 0.42; conditioning on m would give ~0
  expect_gt(m$effects$estimate, 0.25)
})

test_that("a non-identified effect is reported, not estimated", {
  d <- make_conf()
  g <- ilm_dag("dag { x [exposure] ; y [outcome] ; u [latent]
                      u -> x ; u -> y ; x -> y }")
  m <- ilm_dag_model(g, d[c("x", "y")], verbose = FALSE)
  expect_false(m$identified)
  expect_length(m$fits, 0L)
  expect_null(m$effects)
  expect_true(any(grepl("not identifiable", m$steps)))
  expect_output(print(m), "NOT IDENTIFIABLE")

  ## the same graph is identified once the confounder is measured
  g2 <- ilm_dag("dag { x [exposure] ; y [outcome] ; z -> x ; z -> y ; x -> y }")
  expect_true(ilm_dag_model(g2, make_conf(), verbose = FALSE)$identified)
})

test_that("every admissible set is fitted and compared", {
  ## x <- a <- u -> b -> y: the back door can be blocked at a, u or b
  set.seed(2); n <- 900
  u <- rnorm(n); a <- 0.9 * u + rnorm(n); b <- 0.9 * u + rnorm(n)
  x <- 0.8 * a + rnorm(n); y <- 0.5 * x + 0.8 * b + rnorm(n)
  d <- data.frame(x = x, y = y, a = a, b = b, u = u)
  g <- ilm_dag("dag { x [exposure] ; y [outcome]
                      u -> a -> x ; u -> b -> y ; x -> y }")
  m <- ilm_dag_model(g, d, verbose = FALSE)
  expect_length(m$sets, 3L)
  expect_equal(nrow(m$effects), 3L)
  expect_setequal(unlist(m$sets), c("a", "b", "u"))
  ## they estimate the same quantity, so every interval should cover the truth
  expect_true(all(m$effects$lower < 0.5 & m$effects$upper > 0.5))
  expect_true(any(grepl("spread across sets", m$steps)))
  ## max_sets caps the work
  expect_length(ilm_dag_model(g, d, max_sets = 2L, verbose = FALSE)$sets, 2L)
})

test_that("the family is inferred from the outcome and can be overridden", {
  set.seed(5); n <- 1200
  z <- rnorm(n); x <- 0.5 * z + rnorm(n)
  db <- data.frame(x = x, z = z, y = rbinom(n, 1, plogis(0.8 * x + 0.6 * z)))
  mb <- ilm_dag_model(g_conf(), db, verbose = FALSE)
  expect_equal(mb$family, "binomial")
  expect_true(any(grepl("family inferred as binomial", mb$steps)))
  expect_lt(mb$effects$lower, 0.8); expect_gt(mb$effects$upper, 0.8)

  dp <- data.frame(x = x, z = z, y = rpois(n, exp(0.3 * x + 0.2 * z)))
  expect_equal(ilm_dag_model(g_conf(), dp, verbose = FALSE)$family, "poisson")

  ## and an explicit family wins
  expect_equal(ilm_dag_model(g_conf(), dp, family = "nbinom",
                             verbose = FALSE)$family, "nbinom")
})

test_that("clustering is taken only from outside the graph", {
  set.seed(4); ncl <- 30; per <- 20; N <- ncl * per
  site <- factor(rep(seq_len(ncl), each = per))
  ru <- rnorm(ncl, 0, 0.8)[as.integer(site)]
  z <- rnorm(N); x <- 0.5 * z + rnorm(N)
  d <- data.frame(x = x, z = z, site = site,
                  y = 0.4 * x + 0.6 * z + ru + rnorm(N))
  m <- ilm_dag_model(g_conf(), d, verbose = FALSE)
  expect_equal(m$random, "site")
  expect_true(any(grepl("\\(1 \\| site\\)", m$steps)))
  expect_lt(m$effects$lower, 0.4); expect_gt(m$effects$upper, 0.4)

  ## cluster = character(0) fits none
  expect_length(ilm_dag_model(g_conf(), d, cluster = character(),
                              verbose = FALSE)$random, 0L)

  ## a variable the graph DOES mention is never taken as a random effect, even
  ## when it looks like a grouping factor: it has a causal role, and a
  ## confounder belongs in the mean structure
  d2 <- d
  d2$z <- factor(rep(seq_len(20L), each = N / 20L))
  g2 <- ilm_dag("dag { x [exposure] ; y [outcome] ; z -> x -> y ; z -> y }")
  m2 <- ilm_dag_model(g2, d2, verbose = FALSE)
  expect_false("z" %in% m2$random)
  expect_true("z" %in% attr(stats::terms(m2$fits[[1]]), "term.labels"))
})

test_that("the workflow narrates what it did", {
  d <- make_conf()
  expect_message(ilm_dag_model(g_conf(), d), "DAG-guided analysis")
  expect_message(ilm_dag_model(g_conf(), d), "identification")
  expect_message(ilm_dag_model(g_conf(), d), "minimal adjustment set")
  expect_silent(ilm_dag_model(g_conf(), d, verbose = FALSE))
  ## every step claimed is recorded on the object too
  m <- ilm_dag_model(g_conf(), d, verbose = FALSE)
  expect_type(m$steps, "character")
  expect_gt(length(m$steps), 2L)
})

test_that("the workflow refuses what it cannot do", {
  d <- make_conf()
  expect_error(ilm_dag_model(ilm_dag("a -> b"), d, verbose = FALSE),
               "one `exposure` and one")
  expect_error(ilm_dag_model(g_conf(), d[c("y", "z")], verbose = FALSE),
               "exposure `x` is not a column")
  expect_error(ilm_dag_model(g_conf(), "not a frame", verbose = FALSE),
               "must be a data frame")
  ## a multi-level factor outcome is a multinomial, which the package fits
  dd <- d; dd$y <- factor(sample(letters[1:4], nrow(d), TRUE))
  expect_equal(ilm_dag_model(g_conf(), dd, verbose = FALSE)$family,
               "multinomial")
  ## but a response it genuinely cannot place asks rather than guessing
  dc <- d; dc$y <- 1
  expect_error(ilm_dag_model(g_conf(), dc, verbose = FALSE),
               "cannot tell which family")
  dt <- d; dt$y <- as.Date("2024-01-01") + seq_len(nrow(d))
  expect_error(ilm_dag_model(g_conf(), dt, verbose = FALSE),
               "cannot tell which family")
})

test_that("dag_test runs inside the workflow and can be switched off", {
  set.seed(9); n <- 700
  z <- rnorm(n); x <- 0.5 * z + rnorm(n)
  d <- data.frame(x = x, z = z, w = rnorm(n),
                  y = 0.4 * x + 0.6 * z + rnorm(n))
  g <- ilm_dag("dag { x [exposure] ; y [outcome] ; z -> x -> y ; z -> y ; w }")
  m <- ilm_dag_model(g, d, verbose = FALSE)
  expect_false(is.null(m$dag_test))
  expect_equal(attr(m$dag_test, "verdict"), "OK")
  expect_true(any(grepl("graph vs data: OK", m$steps)))
  ## w is unconnected, so the claims about it are the ones tested
  expect_true(all(c("w") %in% c(m$dag_test$x, m$dag_test$y)))

  expect_null(ilm_dag_model(g, d, test_dag = FALSE, verbose = FALSE)$dag_test)
})
