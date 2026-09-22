## Moderation search
##
## The measurements this design rests on, six candidates, n = 500, 150 reps,
## every test the JOINT test of the interaction block:
##
##                            naive   holm   BH     bonf   honest split
##   null (want 0.05)         0.293   0.027  0.027  0.027  0.053
##   moderation by numeric    0.807   0.547  0.553  0.547  0.320
##   moderation by 3-level    0.707   0.380  0.380  0.380  0.200
##
## Naive pick-the-best is the thing being prevented. Holm is the default
## because splitting costs about half the power and the gap WIDENS with df --
## which is why the three-level factor is in these tests, not just a numeric.

mod_data <- function(n = 400, seed = 1, eff = 0.9) {
  set.seed(seed)
  d <- data.frame(tx = rbinom(n, 1, .5), age = rnorm(n), score = rnorm(n),
                  site = factor(sample(c("a", "b", "c"), n, TRUE)))
  d$y <- .4 * d$tx + .3 * d$age + eff * d$tx * (d$site == "c") + rnorm(n)
  d
}

test_that("it finds a planted moderator and not the noise", {
  d <- mod_data()
  f <- ilm_model(y ~ tx + age + score + site, data = d, verbose = FALSE)
  m <- ilm_moderation(f, x = "tx", progress = FALSE)
  expect_s3_class(m, "ilm_moderation")
  expect_identical(m$moderator[1], "site")
  expect_lt(m$p_adj[1], 0.05)
  ## and the two irrelevant ones survive adjustment
  expect_true(all(m$p_adj[m$moderator %in% c("age", "score")] > 0.05))
  ## the test is the JOINT one, so a three-level factor costs 2 df
  expect_identical(m$df[m$moderator == "site"], 2)
  expect_identical(m$df[m$moderator == "age"], 1)
})

test_that("the exposure must be named, because it cannot be guessed", {
  ## x:m is the same term whichever is called the moderator; only the
  ## interpretation distinguishes them, so there is no defensible default
  d <- mod_data()
  f <- ilm_model(y ~ tx + age + site, data = d, verbose = FALSE)
  expect_error(ilm_moderation(f), "must name the treatment")
  expect_error(ilm_moderation(f), "tx, age, site")     # lists the choices
  expect_error(ilm_moderation(f, x = "nope"), "not a fixed effect")
})

test_that("an a priori interaction is excluded, not charged for the search", {
  ## a moderator already in the model was specified in advance and owes no
  ## multiplicity penalty; charging it one would make a pre-registered
  ## hypothesis weaker for having been tested beside exploratory ones
  d <- mod_data()
  f <- ilm_model(y ~ tx * site + age + score, data = d, verbose = FALSE)
  expect_message(m <- ilm_moderation(f, x = "tx", progress = FALSE),
                 "a priori")
  expect_false("site" %in% m$moderator)
  expect_identical(sort(m$moderator), c("age", "score"))
  expect_identical(attr(m, "apriori"), "site")
  out <- utils::capture.output(print(m))
  expect_true(any(grepl("ilm_anova", out)))
})

test_that("the adjustment is p.adjust and any of its methods works", {
  d <- mod_data()
  f <- ilm_model(y ~ tx + age + score + site, data = d, verbose = FALSE)
  raw <- ilm_moderation(f, x = "tx", adjust = "none", progress = FALSE)
  expect_equal(raw$p, raw$p_adj)
  for (a in c("holm", "BH", "bonferroni", "BY", "hochberg")) {
    m <- ilm_moderation(f, x = "tx", adjust = a, progress = FALSE)
    expect_identical(attr(m, "adjust"), a)
    ## an adjusted p is never smaller than the raw one
    expect_true(all(m$p_adj >= m$p - 1e-12))
  }
  expect_error(ilm_moderation(f, x = "tx", adjust = "nonsense"),
               "p.adjust.methods")
})

test_that("adjustment actually changes the answer on the null", {
  ## the whole point: unadjusted, the best of several candidates is not a test
  set.seed(7); n <- 400
  d <- as.data.frame(matrix(rnorm(n * 5), n, 5))
  names(d) <- paste0("v", 1:5)
  d$tx <- rbinom(n, 1, .5)
  d$y <- .4 * d$tx + rnorm(n)                 # NO moderation at all
  f <- ilm_model(y ~ tx + v1 + v2 + v3 + v4 + v5, data = d, verbose = FALSE)
  m <- ilm_moderation(f, x = "tx", progress = FALSE)
  expect_true(all(m$p_adj >= m$p))
  expect_gt(max(m$p_adj / pmax(m$p, 1e-12)), 1)   # something was adjusted
})

test_that("the honest split path runs and reports itself as one", {
  skip_on_cran()
  d <- mod_data(600, seed = 3)
  f <- ilm_model(y ~ tx + age + score + site, data = d, verbose = FALSE)
  m <- ilm_moderation(f, x = "tx", split = TRUE, seed = 2, progress = FALSE)
  expect_true(attr(m, "split"))
  expect_identical(nrow(m), 1L)                # one winner, tested once
  expect_identical(sum(attr(m, "n_split")), nrow(d))
  out <- utils::capture.output(print(m))
  expect_true(any(grepl("honest split", out)))
})

test_that("the refit keeps the model's structure, not just its formula", {
  ## rebuilding ilm_model(f, data = ) by hand drops ziformula, dispformula,
  ## ar, weights, contrasts and reml silently -- a zero-inflated model would
  ## be tested without its zero part and nothing would say so
  set.seed(4); n <- 400
  d <- data.frame(tx = rbinom(n, 1, .5), g = factor(sample(letters[1:3], n, TRUE)))
  d$y <- rpois(n, exp(.6 + .3 * d$tx))
  f <- ilm_model(y ~ tx + g, data = d, family = "poisson", ziformula = ~ 1,
                 verbose = FALSE)
  m <- ilm_moderation(f, x = "tx", progress = FALSE)
  kept <- attr(m, "fits")[[1]]
  expect_false(is.null(kept))
  expect_identical(kept$family$name, "poisson")
  expect_false(is.null(kept$zi_formula))       # the zero part survived
})

test_that("ilm_plot_moderation draws each shape", {
  skip_if_not_installed("tinyplot")
  grDevices::pdf(NULL); on.exit(grDevices::dev.off(), add = TRUE)
  d <- mod_data()
  f <- ilm_model(y ~ tx + age + score + site, data = d, verbose = FALSE)
  m <- ilm_moderation(f, x = "tx", progress = FALSE)

  ## binary treatment stored as a NUMBER still gets both arms: ilm_emmeans()
  ## would otherwise hold it at its mean and draw one curve where two are meant
  r <- ilm_plot_moderation(m)
  expect_identical(nrow(r), 6L)                # 2 arms x 3 sites
  expect_identical(nlevels(factor(r$by)), 2L)

  ## a continuous moderator is evaluated at its quartiles
  r2 <- ilm_plot_moderation(m, moderator = "age")
  expect_identical(nrow(r2), 6L)

  ## a continuous exposure is drawn as SLOPES, which is the quantity differing
  d2 <- d; d2$y <- .3 * d2$age + .8 * d2$age * (d2$site == "c") + rnorm(nrow(d2))
  f2 <- ilm_model(y ~ age + tx + site, data = d2, verbose = FALSE)
  m2 <- ilm_moderation(f2, x = "age", progress = FALSE)
  r3 <- ilm_plot_moderation(m2)
  expect_identical(nrow(r3), 3L)               # one slope per site
  expect_silent(ilm_plot_moderation(m, pch = "filled circle"))
})

test_that("realistic level names do not break the legend", {
  ## This caught a real bug the rest of the file missed. The call is built
  ## with do.call(), which puts the EVALUATED `by` vector into it, and
  ## tinyplot deparses its arguments to title a legend -- so the width
  ## computation throws "invalid graphics state" once the level names are long
  ## enough. Levels of "a"/"b"/"c" pass; "north"/"central"/"south" did not.
  skip_if_not_installed("tinyplot")
  grDevices::pdf(NULL); on.exit(grDevices::dev.off(), add = TRUE)
  set.seed(2026); n <- 300
  d <- data.frame(
    arm = factor(sample(c("control", "treatment"), n, TRUE)),
    site = factor(sample(c("north", "central", "south"), n, TRUE)),
    baseline = rnorm(n, 50, 8))
  tr <- d$arm == "treatment"
  d$outcome <- 50 + 3 * tr + 5 * tr * (d$site == "south") + rnorm(n, 0, 6)
  f <- ilm_model(outcome ~ arm + site + baseline, data = d, verbose = FALSE)
  m <- ilm_moderation(f, x = "arm", progress = FALSE)
  expect_silent(ilm_plot_moderation(m))
  expect_silent(ilm_plot_moderation(m, moderator = "baseline"))
})

test_that("a refit that was not retained says so rather than failing oddly", {
  skip_if_not_installed("tinyplot")
  grDevices::pdf(NULL); on.exit(grDevices::dev.off(), add = TRUE)
  d <- mod_data()
  f <- ilm_model(y ~ tx + age + score + site, data = d, verbose = FALSE)
  m <- ilm_moderation(f, x = "tx", n_keep = 1L, progress = FALSE)
  expect_error(ilm_plot_moderation(m, moderator = m$moderator[3]),
               "not retained")
  expect_error(ilm_plot_moderation(m, moderator = "nope"), "not among")
})

test_that("a plain model with an interaction can be drawn directly", {
  skip_if_not_installed("tinyplot")
  grDevices::pdf(NULL); on.exit(grDevices::dev.off(), add = TRUE)
  d <- mod_data()
  f <- ilm_model(y ~ tx * site + age, data = d, verbose = FALSE)
  expect_message(r <- ilm_plot_moderation(f), "drawing")
  expect_identical(nrow(r), 6L)
  expect_silent(r2 <- ilm_plot_moderation(f, exposure = "tx",
                                          moderator = "site"))
  expect_error(ilm_plot_moderation(
    ilm_model(y ~ tx + age, data = d, verbose = FALSE)), "no interaction")
})
