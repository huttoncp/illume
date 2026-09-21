## Factorial and repeated-measures ANOVA.
##
## The reference is afex, which is what this is for: someone who has always
## written aov_ez() should get the same numbers. Agreement is to machine
## precision, so these are equality checks.

aov_rm <- function(ns = 30, nt = 4, seed = 42) {
  set.seed(seed)
  d <- expand.grid(id = factor(seq_len(ns)), time = factor(seq_len(nt)))
  d$grp <- factor(rep(c("ctl", "trt"), each = ns / 2))[as.integer(d$id)]
  u <- rnorm(ns, 0, 1.2)
  d$y <- 2 + 0.6 * as.integer(d$time) + 0.8 * (d$grp == "trt") +
    u[as.integer(d$id)] + rnorm(nrow(d))
  d
}

expect_matches_afex <- function(a, ref, tol = 1e-6) {
  A <- a$anova
  R <- ref$anova_table
  R$effect <- rownames(R)
  for (e in A$effect) {
    j <- match(e, R$effect)
    expect_false(is.na(j), info = paste("effect missing from afex:", e))
    i <- which(A$effect == e)
    expect_equal(A$F_value[i], R$F[j], tolerance = tol, info = e)
    expect_equal(A$ges[i], R$ges[j], tolerance = tol, info = e)
    expect_equal(A$p_value[i], R$`Pr(>F)`[j], tolerance = tol, info = e)
    expect_equal(A$num_df[i], R$`num Df`[j], tolerance = tol, info = e)
    expect_equal(A$den_df[i], R$`den Df`[j], tolerance = tol, info = e)
  }
}

test_that("a between x within design matches afex", {
  skip_if_not_installed("afex")
  d <- aov_rm()
  a <- ilm_aov_ez("id", "y", d, between = "grp", within = "time",
                  verbose = FALSE, posthoc = FALSE, correction = "none")
  r <- afex::aov_ez("id", "y", d, between = "grp", within = "time",
                    anova_table = list(es = "ges", correction = "none"))
  expect_matches_afex(a, r)
})

test_that("the Greenhouse-Geisser correction matches afex", {
  skip_if_not_installed("afex")
  d <- aov_rm()
  a <- ilm_aov_ez("id", "y", d, between = "grp", within = "time",
                  verbose = FALSE, posthoc = FALSE, correction = "GG")
  r <- afex::aov_ez("id", "y", d, between = "grp", within = "time",
                    anova_table = list(es = "ges", correction = "GG"))
  expect_matches_afex(a, r)
  ## the corrected degrees of freedom are fractional, which is the point
  expect_true(any(a$anova$num_df %% 1 != 0))
})

test_that("two within factors match afex, including the ges denominator", {
  ## This is where the subject stratum has to be put back: with no
  ## between-participants effect its error never reaches the table, and every
  ## ges came out more than twice too large.
  skip_if_not_installed("afex")
  set.seed(7); ns <- 20
  d <- expand.grid(id = factor(seq_len(ns)), A = factor(1:3), B = factor(1:2))
  u <- rnorm(ns, 0, 1)
  d$y <- 1 + 0.5 * as.integer(d$A) - 0.4 * as.integer(d$B) +
    0.3 * as.integer(d$A) * as.integer(d$B) + u[as.integer(d$id)] +
    rnorm(nrow(d))
  a <- ilm_aov_ez("id", "y", d, within = c("A", "B"), verbose = FALSE,
                  posthoc = FALSE, correction = "none")
  r <- afex::aov_ez("id", "y", d, within = c("A", "B"),
                    anova_table = list(es = "ges", correction = "none"))
  expect_matches_afex(a, r)
})

test_that("a purely between-participants factorial matches afex", {
  skip_if_not_installed("afex")
  set.seed(3)
  d <- data.frame(id = factor(1:120),
                  g1 = factor(rep(c("a", "b", "c"), each = 40)),
                  g2 = factor(rep(c("x", "y"), 60)))
  d$y <- as.integer(d$g1) * 0.5 + as.integer(d$g2) * 0.3 + rnorm(120)
  a <- ilm_aov_ez("id", "y", d, between = c("g1", "g2"), verbose = FALSE,
                  posthoc = FALSE)
  r <- afex::aov_ez("id", "y", d, between = c("g1", "g2"),
                    anova_table = list(es = "ges"))
  expect_matches_afex(a, r)
})

test_that("a between-participants covariate matches afex", {
  skip_if_not_installed("afex")
  set.seed(5); ns <- 36
  d <- expand.grid(id = factor(seq_len(ns)), time = factor(1:3))
  d$grp <- factor(rep(c("ctl", "trt"), each = ns / 2))[as.integer(d$id)]
  age <- rnorm(ns, 50, 8); d$age <- age[as.integer(d$id)]
  u <- rnorm(ns, 0, 1)
  d$y <- 2 + 0.5 * as.integer(d$time) + 0.7 * (d$grp == "trt") + 0.05 * d$age +
    u[as.integer(d$id)] + rnorm(nrow(d))
  a <- ilm_aov_ez("id", "y", d, between = "grp", within = "time",
                  covariate = "age", verbose = FALSE, posthoc = FALSE,
                  correction = "none")
  r <- afex::aov_ez("id", "y", d, between = "grp", within = "time",
                    covariate = "age", factorize = FALSE, observed = "age",
                    anova_table = list(es = "ges", correction = "none"))
  expect_matches_afex(a, r, tol = 1e-5)
})

test_that("sphericity matches car's arithmetic", {
  skip_if_not_installed("afex")
  d <- aov_rm()
  a <- ilm_aov_ez("id", "y", d, between = "grp", within = "time",
                  verbose = FALSE, posthoc = FALSE, correction = "none")
  r <- afex::aov_ez("id", "y", d, between = "grp", within = "time",
                    anova_table = list(es = "ges", correction = "none"))
  sp <- suppressWarnings(summary(r$Anova))
  mt <- sp$sphericity.tests
  for (e in rownames(mt)) {
    z <- a$sphericity[[e]]
    expect_false(is.null(z), info = e)
    expect_equal(z$W, unname(mt[e, 1]), tolerance = 1e-6, info = e)
    expect_equal(z$p_mauchly, unname(mt[e, 2]), tolerance = 1e-6, info = e)
    expect_equal(z$gg, unname(sp$pval.adjustments[e, "GG eps"]),
                 tolerance = 1e-6, info = e)
  }
})

test_that("a time-varying covariate is refused with the remedy named", {
  d <- aov_rm()
  d$stress <- rnorm(nrow(d))                 # varies within participant
  expect_error(ilm_aov_ez("id", "y", d, between = "grp", within = "time",
                          covariate = "stress", verbose = FALSE),
               "time-varying covariate")
  expect_error(ilm_aov_ez("id", "y", d, between = "grp", within = "time",
                          covariate = "stress", verbose = FALSE),
               "engine = \"mixed\"")
  ## one that IS constant within participant is accepted
  base <- rnorm(30); d$base <- base[as.integer(d$id)]
  expect_s3_class(ilm_aov_ez("id", "y", d, between = "grp", within = "time",
                             covariate = "base", verbose = FALSE,
                             posthoc = FALSE), "ilm_aov_ez")
})

test_that("incomplete participants are dropped, counted, and the remedy named", {
  d <- aov_rm()
  d2 <- d[!(d$id == "1" & d$time == "2"), ]
  a <- ilm_aov_ez("id", "y", d2, between = "grp", within = "time",
                  verbose = FALSE, posthoc = FALSE)
  expect_equal(a$dropped, 1L)
  expect_equal(a$nsub, 29L)
  expect_true(any(grepl("engine = \"mixed\"", a$notes)))
})

test_that("repeats within a cell are averaged and reported", {
  d <- aov_rm()
  a <- ilm_aov_ez("id", "y", rbind(d, d), between = "grp", within = "time",
                  verbose = FALSE, posthoc = FALSE)
  expect_true(any(grepl("averaged within cells", a$notes)))
  expect_equal(a$nsub, 30L)
})

test_that("the specification translates into the right ilm_model formula", {
  d <- aov_rm()
  a <- ilm_aov_ez("id", "y", d, between = "grp", within = "time",
                  verbose = FALSE, posthoc = FALSE)
  expect_equal(deparse(a$formula), "y ~ grp * time + (1 | id)")
  expect_s3_class(a$fit, "ilm_model")
  expect_true(isTRUE(a$fit$reml))            # the design was fixed in advance
})

test_that("an interaction is followed up with simple effects, not every cell", {
  set.seed(2024); ns <- 36
  d <- expand.grid(id = factor(seq_len(ns)), session = factor(1:3))
  d$arm <- factor(rep(c("p", "l", "h"), each = ns / 3))[as.integer(d$id)]
  eff <- c(p = 0, l = 0.6, h = 1.4)
  u <- rnorm(ns, 0, 1.1)
  d$y <- 50 + eff[as.character(d$arm)] * as.integer(d$session) * 2 +
    u[as.integer(d$id)] + rnorm(nrow(d), 0, 1.2)
  a <- ilm_aov_ez("id", "y", d, between = "arm", within = "session",
                  verbose = FALSE)
  ph <- a$posthoc[["arm:session"]]
  expect_s3_class(ph, "ilm_simple_effects")
  ## both directions, each holding one factor fixed
  expect_length(ph, 2L)
  expect_true(all(vapply(ph, length, 1L) == 3L))
  ## every family is 3 comparisons, not the 36 that all-cells would give
  expect_equal(nrow(ph[[1]][[1]]), 3L)
})

test_that("a two-level factor is not followed up, because the F already said it", {
  d <- aov_rm(ns = 40, nt = 2)
  a <- ilm_aov_ez("id", "y", d, between = "grp", within = "time",
                  verbose = FALSE)
  expect_false("grp" %in% names(a$posthoc))
  expect_false("time" %in% names(a$posthoc))
})

test_that("the specification is checked before anything is computed", {
  d <- aov_rm()
  expect_error(ilm_aov_ez("id", "y", d, between = "time", within = "time",
                          verbose = FALSE), "both between- and within")
  expect_error(ilm_aov_ez("id", "y", d, verbose = FALSE), "at least one")
  expect_error(ilm_aov_ez("id", "nope", d, within = "time", verbose = FALSE),
               "not in the data")
  d$ycat <- factor(d$y > 0)
  expect_error(ilm_aov_ez("id", "ycat", d, within = "time", verbose = FALSE),
               "numeric outcome")
})
