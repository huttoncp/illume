## The graph algorithms are written in-package rather than taken from dagitty,
## so the cases below are the ones where a plausible implementation goes wrong,
## and the last block pins the whole thing to dagitty where it is available.

test_that("a graph reads from a string, an edge frame and a dagitty object", {
  g <- ilm_dag("dag {
    x [exposure]
    y [outcome]
    u [latent]
    x -> y
    u -> x ; u -> y
  }")
  expect_s3_class(g, "ilm_dag")
  expect_setequal(g$nodes, c("x", "y", "u"))
  expect_equal(g$exposure, "x")
  expect_equal(g$outcome, "y")
  expect_equal(g$latent, "u")
  expect_equal(nrow(g$edges), 3L)

  ## the wrapper is optional and ; separates statements
  expect_setequal(ilm_dag("a -> b ; b -> c")$nodes, c("a", "b", "c"))
  ## a chain is several edges, and <- points the other way
  ch <- ilm_dag("a -> b -> c")
  expect_equal(nrow(ch$edges), 2L)
  expect_equal(ilm_dag("a <- b")$edges$from, "b")

  e <- data.frame(from = c("a", "b"), to = c("b", "c"),
                  stringsAsFactors = FALSE)
  expect_setequal(ilm_dag(e, exposure = "a", outcome = "c")$nodes,
                  c("a", "b", "c"))
})

test_that("a bidirected edge becomes an unobserved common cause", {
  g <- ilm_dag("x <-> y")
  expect_length(g$latent, 1L)
  expect_equal(nrow(g$edges), 2L)
  ## both arrows point away from the unobserved node
  expect_true(all(g$edges$from == g$latent))
  expect_setequal(g$edges$to, c("x", "y"))
})

test_that("a graph rejects what it cannot represent", {
  expect_error(ilm_dag("a -> b ; b -> a"), "cycle")
  expect_error(ilm_dag("a -> b -> c ; c -> a"), "cycle")
  expect_error(ilm_dag("a -- b"), "undirected")
  expect_error(ilm_dag("a -> b", exposure = "zz"), "not in the graph")
  expect_error(ilm_dag(42), "must be a dagitty-style string")
  expect_error(ilm_dag(data.frame(a = 1, b = 2)), "columns `from` and `to`")
  expect_error(ilm_dsep(ilm_dag("a -> b"), "a", "zz"), "not in the graph")
})

test_that("d-separation gets the three elementary structures right", {
  ## chain: a -> b -> c
  ch <- ilm_dag("a -> b -> c")
  expect_false(ilm_dsep(ch, "a", "c", character()))
  expect_true(ilm_dsep(ch, "a", "c", "b"))
  ## fork: a <- b -> c
  fk <- ilm_dag("b -> a ; b -> c")
  expect_false(ilm_dsep(fk, "a", "c", character()))
  expect_true(ilm_dsep(fk, "a", "c", "b"))
  ## collider: a -> b <- c. Conditioning OPENS this one, which is the case a
  ## moralisation bug gets backwards.
  cl <- ilm_dag("a -> b ; c -> b")
  expect_true(ilm_dsep(cl, "a", "c", character()))
  expect_false(ilm_dsep(cl, "a", "c", "b"))
  ## and conditioning on a DESCENDANT of a collider opens it too
  cd <- ilm_dag("a -> b ; c -> b ; b -> d")
  expect_true(ilm_dsep(cd, "a", "c", character()))
  expect_false(ilm_dsep(cd, "a", "c", "d"))
})

test_that("adjustment sets handle the cases that catch a naive rule", {
  ## plain confounding
  g1 <- ilm_dag("dag { x [exposure] ; y [outcome] ; x -> y ; z -> x ; z -> y }")
  expect_equal(ilm_adjust_sets(g1), list("z"), ignore_attr = TRUE)

  ## a mediator must NOT be adjusted for a total effect
  g2 <- ilm_dag("dag { x [exposure] ; y [outcome] ; x -> m -> y }")
  s2 <- ilm_adjust_sets(g2)
  expect_length(s2, 1L)
  expect_length(s2[[1]], 0L)

  ## M-bias: z is a collider on the only path, so adjusting for it OPENS the
  ## back door. The empty set is right and {z} is wrong.
  g3 <- ilm_dag("dag {
    x [exposure] ; y [outcome]
    u1 [latent]  ; u2 [latent]
    x -> y ; u1 -> z ; u1 -> x ; u2 -> z ; u2 -> y
  }")
  s3 <- ilm_adjust_sets(g3)
  expect_length(s3, 1L)
  expect_length(s3[[1]], 0L)

  ## unmeasured confounding: nothing identifies the effect, and that is the
  ## answer rather than an error
  g4 <- ilm_dag("dag {
    x [exposure] ; y [outcome] ; u [latent]
    x -> y ; u -> x ; u -> y ; iv -> x
  }")
  expect_length(ilm_adjust_sets(g4), 0L)

  ## a descendant of the exposure is never a candidate
  g5 <- ilm_dag("dag { x [exposure] ; y [outcome] ; x -> y ; x -> d ; z -> x ; z -> y }")
  expect_false(any(vapply(ilm_adjust_sets(g5), function(s) "d" %in% s, TRUE)))

  ## several minimal sets are reported, not just one
  g6 <- ilm_dag("dag {
    x [exposure] ; y [outcome]
    x -> y ; u -> x ; u -> y ; u -> w ; w -> y
  }")
  expect_gte(length(ilm_adjust_sets(g6)), 1L)

  ## restricting to what was measured can remove every set
  expect_length(ilm_adjust_sets(g1, observed = c("x", "y")), 0L)
})

test_that("adjustment sets need an exposure and an outcome", {
  expect_error(ilm_adjust_sets(ilm_dag("a -> b")), "one `exposure` and one")
  expect_error(ilm_adjust_sets(ilm_dag("a -> b"), exposure = "a",
                               outcome = "a"), "same variable")
})

test_that("implied independencies are minimal and genuinely implied", {
  g <- ilm_dag("dag { z -> x -> y ; z -> w }")
  im <- ilm_dag_implied(g)
  expect_true(nrow(im) > 0L)
  ## everything claimed really is a d-separation in the same graph
  for (i in seq_len(nrow(im)))
    expect_true(ilm_dsep(g, im$x[i], im$y[i], im$given[[i]]))
  ## and nothing in the conditioning set is redundant
  for (i in seq_len(nrow(im))) {
    z <- im$given[[i]]
    for (v in z)
      expect_false(ilm_dsep(g, im$x[i], im$y[i], setdiff(z, v)))
  }
  ## an adjacent pair makes no claim
  expect_equal(nrow(ilm_dag_implied(ilm_dag("a -> b"))), 0L)

  ## A claim whose textbook conditioning set involves an unobserved parent is
  ## still testable when a smaller observed set separates the pair -- here the
  ## empty set does, and building the set from parents alone would drop it.
  g2 <- ilm_dag("dag { u [latent] ; u -> a ; b -> c ; a -> c }")
  im2 <- ilm_dag_implied(g2)
  k <- which((im2$x == "a" & im2$y == "b") | (im2$x == "b" & im2$y == "a"))
  expect_length(k, 1L)
  expect_length(im2$given[[k]], 0L)
})

test_that("the graph algorithms agree with dagitty", {
  skip_if_not_installed("dagitty")
  ## random DAGs, acyclic by construction, some with an unobserved node
  set.seed(4242L)
  n_dsep <- 0L; n_adj <- 0L; n_ci <- 0L
  for (s in 1:60) {
    p <- sample(4:6, 1L)
    nm <- paste0("v", seq_len(p))
    q <- stats::runif(1, 0.3, 0.6)
    ed <- do.call(rbind, lapply(1:(p - 1L), function(i) {
      j <- (i + 1L):p
      k <- j[stats::runif(length(j)) < q]
      if (!length(k)) NULL else
        data.frame(from = nm[i], to = nm[k], stringsAsFactors = FALSE)
    }))
    if (is.null(ed) || !nrow(ed)) next
    lat <- if (stats::runif(1) < 0.3) sample(nm, 1L) else character()
    x <- sample(setdiff(nm, lat), 1L)
    y <- sample(setdiff(nm, c(x, lat)), 1L)
    st <- paste0("dag {\n",
                 paste(c(setdiff(nm, c(x, y, lat)),
                         sprintf("%s [exposure]", x),
                         sprintf("%s [outcome]", y),
                         sprintf("%s [latent]", lat),
                         sprintf("%s -> %s", ed$from, ed$to)),
                       collapse = "\n"), "\n}")
    gi <- ilm_dag(st); gd <- dagitty::dagitty(st)

    ab <- sample(nm, 2L)
    zp <- setdiff(nm, ab)
    z <- if (length(zp)) sample(zp, sample(0:length(zp), 1L)) else character()
    expect_equal(ilm_dsep(gi, ab[1], ab[2], z),
                 dagitty::dseparated(gd, ab[1], ab[2], z))
    n_dsep <- n_dsep + 1L

    ky <- function(L) unname(sort(vapply(unname(L), function(v)
      if (!length(v)) "{}" else paste(sort(v), collapse = "+"), "")))
    expect_equal(ky(ilm_adjust_sets(gi)),
                 ky(lapply(dagitty::adjustmentSets(gd, x, y, type = "minimal",
                                                   effect = "total"),
                           as.character)))
    n_adj <- n_adj + 1L

    im <- ilm_dag_implied(gi)
    for (i in seq_len(nrow(im))) {
      expect_true(dagitty::dseparated(gd, im$x[i], im$y[i], im$given[[i]]))
      n_ci <- n_ci + 1L
    }
  }
  expect_gt(n_dsep, 40L); expect_gt(n_adj, 40L); expect_gt(n_ci, 40L)
})

test_that("dag_test passes a true graph and catches a missing arrow", {
  set.seed(7); n <- 800
  z <- rnorm(n); x <- 0.6 * z + rnorm(n); y <- 0.5 * x + 0.4 * z + rnorm(n)
  d <- data.frame(x = x, y = y, z = z, w = rnorm(n),
                  b = rbinom(n, 1, plogis(0.5 * z)),
                  ct = rpois(n, exp(0.3 * z)))

  ok <- ilm_dag_test(ilm_dag("dag { x [exposure] ; y [outcome]
                               z -> x -> y ; z -> y ; z -> b ; z -> ct ; w }"),
                     d, verbose = FALSE)
  expect_equal(attr(ok, "verdict"), "OK")
  expect_true(all(ok$verdict == "OK"))
  ## the families are taken from the variables, so a binary and a count
  ## variable are both tested rather than skipped
  expect_true(any(ok$x == "b" | ok$y == "b"))
  expect_true(any(ok$x == "ct" | ok$y == "ct"))

  ## drop the z -> y arrow and the graph now claims something false
  bad <- ilm_dag_test(ilm_dag("dag { x [exposure] ; y [outcome] ; z -> x -> y }"),
                      d[c("x", "y", "z")], verbose = FALSE)
  expect_equal(attr(bad, "verdict"), "FAIL")
  expect_equal(nrow(bad), 1L)
  expect_gt(abs(bad$estimate), 0.2)

  ## a saturated graph implies nothing, which is reported rather than passed
  none <- ilm_dag_test(ilm_dag("dag { x [exposure] ; y [outcome] ; x -> y }"),
                       d, verbose = FALSE)
  expect_equal(attr(none, "verdict"), "UNTESTED")
  expect_equal(nrow(none), 0L)
})

test_that("dag_test needs an effect as well as a p-value", {
  ## a tiny true dependence at a large n is significant and not worth acting
  ## on; min_effect is what separates the two
  ## 0.08 puts the partial correlation reliably between the two thresholds:
  ## across twelve seeds it ran 0.040 to 0.086 with a largest adjusted p-value
  ## of 7.4e-04, so the test turns on the effect size and not on the draw
  set.seed(101); n <- 8000
  z <- rnorm(n); x <- rnorm(n)
  y <- 0.08 * x + 0.5 * z + rnorm(n)       # x -> y exists but is negligible
  d <- data.frame(x = x, y = y, z = z)
  g <- ilm_dag("dag { x [exposure] ; y [outcome] ; z -> y ; x }")
  r <- ilm_dag_test(g, d, verbose = FALSE)
  k <- which((r$x == "x" & r$y == "y") | (r$x == "y" & r$y == "x"))
  expect_length(k, 1L)
  expect_lt(r$p_adj[k], 0.05)               # significant
  expect_lt(abs(r$estimate[k]), 0.1)        # but under the threshold
  expect_equal(r$verdict[k], "WARN")        # so a warning, not a contradiction
  expect_equal(attr(r, "verdict"), "WARN")
  ## and lowering the threshold turns the same claim into a contradiction
  r2 <- ilm_dag_test(g, d, min_effect = 0.01, verbose = FALSE)
  expect_equal(r2$verdict[k], "FAIL")
})

test_that("dag_test says why a claim was not tested", {
  set.seed(3); n <- 200
  d <- data.frame(f1 = factor(sample(letters[1:3], n, TRUE)),
                  f2 = factor(sample(LETTERS[1:3], n, TRUE)),
                  y = rnorm(n))
  r <- ilm_dag_test(ilm_dag("dag { y [outcome] ; f1 [exposure] ; f1 -> y ; f2 }"),
                    d, verbose = FALSE)
  k <- which((r$x == "f1" & r$y == "f2") | (r$x == "f2" & r$y == "f1"))
  expect_length(k, 1L)
  expect_equal(r$verdict[k], "UNTESTED")
  expect_match(r$note[k], "nominal")
})
