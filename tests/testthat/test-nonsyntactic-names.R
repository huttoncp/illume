## Column names that need backticks -- `my y`, `x 1`, `group id` -- must give
## exactly the fit the same data gives under ordinary names. Every one of
## these failed before 0.0.7.9000: all.vars() and a bare deparse() drop the
## backticks, mgcv::interpret.gam() rebuilds formulas from text, and mgcv's
## get.var() parses a smooth's variable names as code.
##
## The ORACLE is the same data with syntactic names: renaming a column cannot
## change a likelihood, so anything but identity is a defect.

nsn_data <- function(n = 240L, seed = 1L) {
  set.seed(seed)
  d <- data.frame(x1 = rnorm(n),
                  g = factor(sample(c("north", "central", "south"), n, TRUE)),
                  id = factor(sample(12, n, TRUE)), w = rpois(n, 2) + 1)
  re <- rnorm(12, 0, .7)[d$id]
  d$y <- 1 + 0.5 * d$x1 + (d$g == "south") + re + rnorm(n)
  d$cnt <- rpois(n, exp(0.2 + 0.3 * d$x1 + re / 2))
  eta <- cbind(0.3 * d$x1 + re, -0.4 * d$x1)
  P <- exp(cbind(eta, 0)); P <- P / rowSums(P)
  d$cat <- factor(apply(P, 1, function(p)
    sample(c("lo", "mid", "hi"), 1, prob = p)))
  b <- d
  names(b) <- c("x 1", "grp a", "group id", "w t", "my y", "the count",
                "out come")
  list(d = d, b = b)
}

nsn_same <- function(f1, f2) {
  expect_equal(as.numeric(logLik(f2)), as.numeric(logLik(f1)), tolerance = 1e-8)
  expect_equal(unname(coef(f2)), unname(coef(f1)), tolerance = 1e-8)
  expect_equal(unname(sqrt(diag(stats::vcov(f2)))),
               unname(sqrt(diag(stats::vcov(f1)))), tolerance = 1e-6)
}

test_that("the names that need backticks are quoted, and only those", {
  expect_equal(ilm_bq(c("x", "x 1", "if", "x.1", "_y")),
               c("x", "`x 1`", "`if`", "x.1", "`_y`"))
  expect_equal(ilm_unbq(c("`x 1`", "x", "log(`x 1`)")),
               c("x 1", "x", "log(`x 1`)"))
  expect_equal(ilm_as_label(c("x 1", "`x 1`", "g", "nope"), c("`x 1`", "g")),
               c("`x 1`", "`x 1`", "g", "nope"))
  ## a long expression stays on one line; a bare deparse() splits it at 60
  ## characters and pasting the pieces keeps only the first
  long <- quote(1 + a_long_variable_name + another_long_variable_name +
                  yet_another_long_variable_name)
  expect_length(ilm_term_text(long), 1L)
  expect_equal(ilm_mf_name(quote(`x 1`)), "x 1")
  expect_equal(ilm_mf_name(quote(log(`x 1`))), "log(`x 1`)")
})

test_that("every formula feature fits identically under backticked names", {
  z <- nsn_data(); d <- z$d; b <- z$b
  q <- function(...) suppressMessages(ilm_model(..., verbose = FALSE))
  nsn_same(q(y ~ x1 + g + (1 | id), data = d),
           q(`my y` ~ `x 1` + `grp a` + (1 | `group id`), data = b))
  nsn_same(q(y ~ x1 + (1 + x1 | id), data = d),
           q(`my y` ~ `x 1` + (1 + `x 1` | `group id`), data = b))
  nsn_same(q(y ~ x1 + (1 | g / id), data = d),
           q(`my y` ~ `x 1` + (1 | `grp a` / `group id`), data = b))
  nsn_same(q(cat ~ x1 + (1 | id), data = d, family = "multinomial"),
           q(`out come` ~ `x 1` + (1 | `group id`), data = b,
             family = "multinomial"))
  nsn_same(q(y ~ s(x1) + (1 | id), data = d),
           q(`my y` ~ s(`x 1`) + (1 | `group id`), data = b))
  nsn_same(q(y ~ x1, data = d, dispformula = ~ g),
           q(`my y` ~ `x 1`, data = b, dispformula = ~ `grp a`))
  nsn_same(q(cnt ~ x1, data = d, family = "poisson", ziformula = ~ x1),
           q(`the count` ~ `x 1`, data = b, family = "poisson",
             ziformula = ~ `x 1`))
  ## weights are evaluated in the data, so they cannot go through q()
  nsn_same(ilm_model(y ~ x1, data = d, weights = w, verbose = FALSE),
           ilm_model(`my y` ~ `x 1`, data = b, weights = `w t`,
                     verbose = FALSE))
  ## the grouping factor keeps its plain name
  f <- q(`my y` ~ `x 1` + (1 | `group id`), data = b)
  expect_true("group id" %in% names(f$re))
})

test_that("prediction on new data works for a smooth and a random slope", {
  z <- nsn_data(); d <- z$d; b <- z$b
  q <- function(...) suppressMessages(ilm_model(..., verbose = FALSE))
  s1 <- q(y ~ s(x1) + (1 | id), data = d)
  s2 <- q(`my y` ~ s(`x 1`) + (1 | `group id`), data = b)
  expect_equal(unname(predict(s2, newdata = b[1:15, ])),
               unname(predict(s1, newdata = d[1:15, ])), tolerance = 1e-8)
  ## the random slope's design is rebuilt on the new rows; it used to fail to
  ## parse and be dropped from the draws in silence
  r2 <- q(`my y` ~ `x 1` + (1 + `x 1` | `group id`), data = b)
  k <- which(vapply(r2$re, function(e) !identical(e$kind, "basis"), TRUE))[1]
  Z <- ilm_re_design(r2, k, b[1:15, ], 2L)
  expect_equal(dim(Z), c(15L, 2L))
})

test_that("a smooth's `by` variable reaches the model frame", {
  ## it used to be left out, and mgcv stopped with "Can't find by variable"
  ## unless the variable also appeared somewhere else in the formula
  set.seed(7); n <- 300
  d <- data.frame(x1 = runif(n), z = rnorm(n))
  d$y <- sin(2 * pi * d$x1) * d$z + rnorm(n, 0, .3)
  f <- ilm_model(y ~ s(x1, by = z), data = d, verbose = FALSE)
  g <- mgcv::gam(y ~ s(x1, by = z), data = d, method = "ML")
  expect_gt(stats::cor(predict(f), stats::fitted(g)), 0.9999)
})

test_that("functions that take names take them the way a person writes them", {
  z <- nsn_data(); d <- z$d; b <- z$b
  f1 <- ilm_model(y ~ x1 + g, data = d, verbose = FALSE)
  f2 <- ilm_model(`my y` ~ `x 1` + `grp a`, data = b, verbose = FALSE)
  m1 <- ilm_moderation(f1, x = "x1", progress = FALSE)
  m2 <- ilm_moderation(f2, x = "x 1", progress = FALSE)
  expect_equal(m2$p, m1$p)
  expect_equal(m2$moderator, "grp a")
  expect_equal(attr(m2, "x"), "x 1")
  grDevices::pdf(NULL); on.exit(grDevices::dev.off(), add = TRUE)
  expect_silent(ilm_plot_moderation(m2))
  expect_equal(ilm_term_cols(f2, "x 1"), ilm_term_cols(f1, "x1"))
  expect_equal(ilm_term_cols(f2, "`x 1`"), ilm_term_cols(f1, "x1"))
  g1 <- ilm_model(y ~ x1 + g + (1 | id), data = d, verbose = FALSE)
  g2 <- ilm_model(`my y` ~ `x 1` + `grp a` + (1 | `group id`), data = b,
                  verbose = FALSE)
  p1 <- ilm_pb_lrt(g1, "g", B = 5, verbose = FALSE)
  p2 <- ilm_pb_lrt(g2, "grp a", B = 5, verbose = FALSE)
  expect_equal(p2$LR, p1$LR, tolerance = 1e-8)
})

test_that("functions that build formulas from column names quote them", {
  ## instrumental variables
  set.seed(2); n <- 300
  iv <- data.frame(z = rnorm(n), wc = rnorm(n)); u <- rnorm(n)
  iv$x <- 0.8 * iv$z + 0.3 * iv$wc + u + rnorm(n)
  iv$y <- 1 + 0.5 * iv$x + 0.2 * iv$wc + u + rnorm(n)
  ivb <- iv; names(ivb) <- c("z inst", "w cov", "x 1", "my y")
  expect_equal(
    unname(ilm_iv(`my y` ~ `x 1` + `w cov` | `z inst` + `w cov`,
                  data = ivb)$coefficients),
    unname(ilm_iv(y ~ x + wc | z + wc, data = iv)$coefficients))

  ## factorial and repeated measures: effect labels come back plain, so an
  ## `observed` factor is still recognised in them
  set.seed(5)
  ad <- expand.grid(subj = factor(1:16), time = factor(c("t1", "t2", "t3")))
  ad$dose <- factor(ifelse(as.integer(ad$subj) <= 8, "low", "high"))
  ad$score <- 10 + (ad$time == "t3") + 0.8 * (ad$dose == "high") +
    rnorm(16)[ad$subj] + rnorm(nrow(ad))
  adb <- ad; names(adb) <- c("subj id", "time point", "dose group", "the score")
  a1 <- ilm_aov_ez("subj", "score", ad, between = "dose", within = "time",
                   observed = "dose", verbose = FALSE)
  a2 <- ilm_aov_ez("subj id", "the score", adb, between = "dose group",
                   within = "time point", observed = "dose group",
                   verbose = FALSE)
  expect_equal(a2$anova$effect,
               c("dose group", "time point", "dose group:time point"))
  expect_equal(a2$anova$F_value, a1$anova$F_value)
  expect_equal(a2$anova$ges, a1$anova$ges)

  ## regression discontinuity and the missing-data check
  set.seed(4)
  rd <- data.frame(r = runif(400, -1, 1), cv = rnorm(400))
  rd$y <- 0.5 * rd$r + 0.8 * (rd$r >= 0) + 0.2 * rd$cv + rnorm(400, 0, .5)
  rdb <- rd; names(rdb) <- c("run var", "cov 1", "my y")
  expect_equal(
    ilm_rdd(rdb, "my y", "run var", covariates = "cov 1", verbose = FALSE)$jump,
    ilm_rdd(rd, "y", "r", covariates = "cv", verbose = FALSE)$jump)
  set.seed(6)
  md <- data.frame(y = rnorm(200), x1 = rnorm(200),
                   g = factor(sample(c("a", "b"), 200, TRUE)))
  md$x1[sample(200, 30)] <- NA
  mdb <- md; names(mdb) <- c("my y", "x 1", "grp a")
  expect_equal(
    ilm_check_missing(mdb, y = "my y", verbose = FALSE)$associations$p_value,
    ilm_check_missing(md, y = "y", verbose = FALSE)$associations$p_value)
})
