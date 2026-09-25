## A model built from assumptions instead of data
##
## The thing worth testing is not that ilm_scaffold() returns an object -- it
## is that the object really carries the numbers that were asked for, and that
## every route into it agrees. A scaffold whose parameters drifted would give
## a power curve for a study nobody is planning, and nothing about the output
## would look wrong.
##
## So the checks here come back through independent paths: coef() reads the
## parameters, ilm_emmeans() recomputes the cell means from them, ilm_simulate()
## draws data whose spread can be measured, and ilm_power() is compared against
## the closed form for a two-sample t test.

test_that("cell means become the coefficients they imply", {
  s <- ilm_scaffold(y ~ arm, design = list(arm = c("control", "treatment")),
                    n_unit = 120, cells = c(control = 12, treatment = 14.5),
                    sd = 4, verbose = FALSE)
  expect_s3_class(s, "ilm_scaffold")
  expect_s3_class(s, "ilm_model")
  expect_equal(unname(stats::coef(s)), c(12, 2.5), tolerance = 1e-8)
  expect_equal(unname(s$dispersion), 4, tolerance = 1e-8)
  expect_identical(nrow(s$model), 120L)
  ## the allocation is balanced, not merely random
  expect_identical(unname(as.integer(table(s$model$arm))), c(60L, 60L))
})

test_that("the two ways of saying it give the same model", {
  a <- ilm_scaffold(y ~ arm, design = list(arm = c("control", "treatment")),
                    n_unit = 120, cells = c(control = 12, treatment = 14.5),
                    sd = 4, verbose = FALSE)
  b <- ilm_scaffold(y ~ arm, design = list(arm = c("control", "treatment")),
                    n_unit = 120, sd = 4, verbose = FALSE,
                    coefs = c("(Intercept)" = 12, "armtreatment" = 2.5))
  expect_equal(stats::coef(a), stats::coef(b), tolerance = 1e-10)
  expect_equal(a$dispersion, b$dispersion, tolerance = 1e-10)
})

test_that("emmeans recovers the cell means that were assumed", {
  ## a completely independent route back: the solve went one way, this comes
  ## the other, so a sign or contrast error could not pass both
  s <- ilm_scaffold(y ~ arm * time + (1 | id),
                    design = list(arm = c("control", "treatment"),
                                  time = c("pre", "post")),
                    within = "time", n_unit = 60,
                    cells = c(control.pre = 12, control.post = 12.2,
                              treatment.pre = 12, treatment.post = 14.5),
                    sd = 4, icc = 0.5, verbose = FALSE)
  em <- as.data.frame(ilm_emmeans(s, c("arm", "time")))
  key <- paste(as.character(em$arm), as.character(em$time), sep = ".")
  want <- c(control.pre = 12, control.post = 12.2,
            treatment.pre = 12, treatment.post = 14.5)
  expect_equal(em$estimate, unname(want[key]), tolerance = 1e-6)
})

test_that("a repeated-measures grid is a real repeated-measures grid", {
  s <- ilm_scaffold(y ~ arm * time + (1 | id),
                    design = list(arm = c("control", "treatment"),
                                  time = c("pre", "post")),
                    within = "time", n_unit = 60,
                    cells = c(control.pre = 12, control.post = 12.2,
                              treatment.pre = 12, treatment.post = 14.5),
                    sd = 4, icc = 0.5, verbose = FALSE)
  expect_identical(nrow(s$model), 120L)
  ## a between-unit variable allocated to ROWS would put a participant in both
  ## arms, which is the mistake this design is built to avoid
  expect_true(all(tapply(as.character(s$model$arm), s$model$id,
                         function(z) length(unique(z))) == 1L))
  expect_true(all(tapply(as.character(s$model$time), s$model$id,
                         function(z) length(unique(z))) == 2L))
  ## icc = 0.5 with sd = 4 means a random-intercept sd of 4 as well
  expect_equal(sqrt(s$Sigma[[1]][1, 1]), 4, tolerance = 1e-6)
  expect_equal(unname(s$dispersion), 4, tolerance = 1e-8)
})

test_that("cell means the formula cannot produce are refused", {
  ## a least-squares solve ALWAYS returns something, so the failure mode here
  ## is a scaffold for a different study than the one being planned
  dd <- list(arm = c("control", "treatment"), time = c("pre", "post"))
  crossed <- c(control.pre = 12, control.post = 12.2,
               treatment.pre = 12, treatment.post = 14.5)
  expect_error(ilm_scaffold(y ~ arm + time, dd, 80, cells = crossed, sd = 4,
                            verbose = FALSE), "cannot produce the cell means")
  ## and it names the remedy
  expect_error(ilm_scaffold(y ~ arm + time, dd, 80, cells = crossed, sd = 4,
                            verbose = FALSE), "arm \\* time")
  ## the same means under the formula that CAN produce them are accepted
  expect_silent(ilm_scaffold(y ~ arm * time, dd, 80, cells = crossed, sd = 4,
                             verbose = FALSE))
  ## and additive means under the additive formula are accepted
  expect_silent(ilm_scaffold(y ~ arm + time, dd, 80, sd = 4, verbose = FALSE,
                 cells = c(control.pre = 12, control.post = 12.2,
                           treatment.pre = 14.5, treatment.post = 14.7)))
})

test_that("random slopes and their correlation are imposed exactly", {
  s <- ilm_scaffold(y ~ arm * time + (1 + time | id),
                    design = list(arm = c("control", "treatment"),
                                  time = 0:3),
                    n_unit = 60, sd = 1, re_sd = list(id = c(2, 1.2)),
                    re_cor = list(id = -0.4), verbose = FALSE,
                    coefs = c("(Intercept)" = 5, "armtreatment" = 0,
                              "time" = 0.3, "armtreatment:time" = 0.5))
  V <- s$Sigma[[1]][1, 1] * s$Sigma_d[["id"]]
  expect_equal(sqrt(V[1, 1]), 2, tolerance = 1e-6)
  expect_equal(sqrt(V[2, 2]), 1.2, tolerance = 1e-6)
  expect_equal(V[1, 2] / sqrt(V[1, 1] * V[2, 2]), -0.4, tolerance = 1e-6)
  ## `time` is in the bar, so it is within-unit without being named
  expect_true("time" %in% s$scaffold$within)
})

test_that("simulated data really has the assumed spread", {
  s <- ilm_scaffold(y ~ arm, design = list(arm = c("control", "treatment")),
                    n_unit = 120, cells = c(control = 12, treatment = 14.5),
                    sd = 4, verbose = FALSE)
  ys <- ilm_simulate(s, nsim = 100, seed = 4)
  ## residual sd 4, plus a two-group mean difference of 2.5 split evenly
  expect_equal(mean(apply(ys, 2, stats::sd)), sqrt(16 + (2.5 / 2)^2),
               tolerance = 0.04)
})

test_that("the power curve matches the closed form for a two-sample t", {
  skip_on_cran()
  p <- ilm_power_design(y ~ arm,
                        design = list(arm = c("control", "treatment")),
                        n_unit = c(100, 200), sd = 4, term = "arm",
                        cells = c(control = 12, treatment = 14.5),
                        sims = 300, seed = 3, verbose = FALSE)
  expect_identical(p$n_unit, c(100L, 200L))
  for (i in seq_len(nrow(p))) {
    want <- stats::power.t.test(n = p$n_unit[i] / 2, delta = 2.5,
                                sd = 4)$power
    expect_equal(p$power[i], want, tolerance = 0.07)
  }
})

test_that("the specification is checked, with the remedy named", {
  d <- list(arm = c("a", "b"))
  expect_error(ilm_scaffold(y ~ arm, d, 100, sd = 4), "exactly one")
  expect_error(ilm_scaffold(y ~ arm, d, 100, cells = c(a = 1, b = 2)), "`sd`")
  expect_error(ilm_scaffold(y ~ arm + (1 | id), d, 40, sd = 1,
                            cells = c(a = 1, b = 2)), "`re_sd` or `icc`")
  expect_error(ilm_scaffold(y ~ arm, d, 100, family = "poisson", sd = 4,
                            cells = c(a = 1, b = 2)), "no separate dispersion")
  expect_error(ilm_scaffold(y ~ arm, d, 2, cells = c(a = 1, b = 2), sd = 1),
               "not a study")
  expect_error(ilm_scaffold(y ~ arm, d, 100, sd = 4,
                            coefs = c("(Intercept)" = 1, nonsense = 2)),
               "does not have")
  ## a cell name that does not match says which names it wanted
  expect_error(ilm_scaffold(y ~ arm, d, 100, sd = 4, cells = c(a = 1, z = 2)),
               "Expected names")
  ## the grouping factor is counted by n_unit, not described in the design
  expect_error(ilm_scaffold(y ~ arm + (1 | id),
                            c(d, list(id = as.character(1:20))), 40, sd = 1,
                            icc = 0.5, cells = c(a = 1, b = 2)),
               "built from `n_unit`")
})

test_that("a non-gaussian scaffold works on the response scale", {
  s <- ilm_scaffold(y ~ arm, design = list(arm = c("control", "treatment")),
                    n_unit = 200, family = "poisson",
                    cells = c(control = 4, treatment = 6), verbose = FALSE)
  ## cells are means; the coefficients are their logs
  expect_equal(unname(stats::coef(s)), c(log(4), log(6) - log(4)),
               tolerance = 1e-8)
  ## a mean that cannot go on the log scale is refused, not silently clamped
  expect_error(ilm_scaffold(y ~ arm, list(arm = c("a", "b")), 100,
                            family = "poisson", cells = c(a = 0, b = 2),
                            verbose = FALSE), "log scale")
})

test_that("between-unit factors are balanced by cell, not only by margin", {
  ## allocated one factor at a time, 18 participants in a 3 x 3 design came
  ## out with an empty cell more often than not
  set.seed(1)
  for (n in c(18L, 20L, 25L)) {
    g <- ilm_scaffold_grid(list(a = c("a1", "a2", "a3"), b = c("b1", "b2", "b3")),
                           n, NULL, NULL)
    tab <- table(g$a, g$b)
    expect_lte(max(tab) - min(tab), 1L)
    expect_lte(diff(range(table(g$a))), 1L)
    expect_lte(diff(range(table(g$b))), 1L)
  }
})

test_that("a multinomial scaffold reproduces the probabilities it was given", {
  cells <- rbind(control = c(none = 0.5, some = 0.3, full = 0.2),
                 treatment = c(none = 0.35, some = 0.35, full = 0.3))
  s <- ilm_scaffold(y ~ arm, design = list(arm = c("control", "treatment")),
                    n_unit = 200, family = "multinomial", cells = cells,
                    verbose = FALSE)
  expect_identical(s$scaffold$categories, c("none", "some", "full"))
  nd <- data.frame(arm = factor(c("control", "treatment"),
                                levels = c("control", "treatment")))
  expect_equal(unname(as.matrix(predict(s, newdata = nd, type = "response"))),
               unname(cells), tolerance = 1e-8)
  ## sum-to-zero: a category's intercept is its log-probability less the
  ## average of all three
  lp <- log(cells["control", ])
  expect_equal(unname(nlme::fixef(s)["(Intercept)", ]), unname((lp - mean(lp))[1:2]),
               tolerance = 1e-8)
  ## the same study stated as coefficients, in either layout
  s2 <- ilm_scaffold(y ~ arm, design = list(arm = c("control", "treatment")),
                     n_unit = 200, family = "multinomial",
                     categories = c("none", "some", "full"),
                     coefs = coef(s), verbose = FALSE)
  expect_equal(coef(s2), coef(s))
  s3 <- ilm_scaffold(y ~ arm, design = list(arm = c("control", "treatment")),
                     n_unit = 200, family = "multinomial",
                     categories = c("none", "some", "full"),
                     coefs = nlme::fixef(s), verbose = FALSE)
  expect_equal(coef(s3), coef(s))
  expect_output(print(s), "minus the sum of the others")
})

test_that("a multinomial random intercept is a shift in every category", {
  s <- ilm_scaffold(y ~ arm + (1 | id),
                    design = list(arm = c("control", "treatment"), time = 1:4),
                    within = "time", n_unit = 60, family = "multinomial",
                    categories = c("a", "b", "c"),
                    coefs = matrix(0, 2, 2,
                                   dimnames = list(c("(Intercept)", "armtreatment"),
                                                   c("a", "b"))),
                    re_sd = list(id = 0.8), verbose = FALSE)
  ## exchangeable shifts, centred: 0.8^2 (I - 1/3)
  expect_equal(unname(s$Sigma[[1]]), 0.64 * (diag(2) - 1 / 3), tolerance = 1e-8)
})

test_that("an ordinal scaffold takes thresholds, or probabilities that are proportional odds", {
  s <- ilm_scaffold(y ~ arm, design = list(arm = c("control", "treatment")),
                    n_unit = 150, family = "ordinal",
                    categories = c("lo", "mid", "hi"),
                    coefs = c(armtreatment = 0.7), thresholds = c(-0.5, 1),
                    verbose = FALSE)
  expect_equal(unname(s$zeta), c(-0.5, 1), tolerance = 1e-8)
  expect_equal(unname(coef(s)), 0.7, tolerance = 1e-8)
  pc <- function(eta) { cp <- stats::plogis(c(-0.5, 1) - eta)
    c(lo = cp[1], mid = cp[2] - cp[1], hi = 1 - cp[2]) }
  cells <- rbind(control = pc(0), treatment = pc(0.7))
  s2 <- ilm_scaffold(y ~ arm, design = list(arm = c("control", "treatment")),
                     n_unit = 150, family = "ordinal", cells = cells,
                     verbose = FALSE)
  expect_equal(unname(coef(s2)), 0.7, tolerance = 1e-8)
  expect_equal(unname(s2$zeta), c(-0.5, 1), tolerance = 1e-8)
  ## probabilities no proportional-odds model can give are refused, and the
  ## model that can give them is named
  bad <- cells; bad["treatment", ] <- c(0.2, 0.6, 0.2)
  expect_error(ilm_scaffold(y ~ arm, design = list(arm = c("control", "treatment")),
                            n_unit = 150, family = "ordinal", cells = bad,
                            verbose = FALSE),
               "not ones proportional odds can produce.*multinomial")
  expect_error(ilm_scaffold(y ~ arm, design = list(arm = c("control", "treatment")),
                            n_unit = 150, family = "ordinal",
                            categories = c("lo", "mid", "hi"),
                            coefs = c(armtreatment = 0.7), verbose = FALSE),
               "needs `thresholds`")
})

test_that("a categorical outcome's specification is checked", {
  d <- list(arm = c("a", "b"))
  cells <- rbind(a = c(x = 0.5, y = 0.3, z = 0.2), b = c(x = 0.4, y = 0.4, z = 0.2))
  expect_error(ilm_scaffold(y ~ arm, d, 100, family = "multinomial",
                            coefs = c("x:(Intercept)" = 0), verbose = FALSE),
               "needs `categories`")
  expect_error(ilm_scaffold(y ~ arm, d, 100, family = "multinomial",
                            categories = c("x", "y"), cells = cells[, 1:2],
                            verbose = FALSE), "at least 3 categories")
  bad <- cells; bad[1, 1] <- 0.6
  expect_error(ilm_scaffold(y ~ arm, d, 100, family = "multinomial",
                            cells = bad, verbose = FALSE), "sum to 1")
  expect_error(ilm_scaffold(y ~ arm, d, 100, family = "multinomial",
                            categories = c("x", "y", "z"),
                            coefs = c("z:armb" = 1), verbose = FALSE),
               "minus the sum of the others")
  expect_error(ilm_scaffold(y ~ arm, d, 100, sd = 1, cells = c(a = 1, b = 2),
                            categories = c("x", "y", "z"), verbose = FALSE),
               "multinomial or ordinal")
})

test_that("an ICC for a binary outcome is on the latent scale", {
  s <- ilm_scaffold(y ~ arm + (1 | id),
                    design = list(arm = c("control", "treatment"), time = 1:4),
                    within = "time", n_unit = 60, family = "binomial",
                    cells = c(control = 0.3, treatment = 0.45), icc = 0.1,
                    verbose = FALSE)
  expect_equal(sqrt(s$Sigma[[1]][1, 1]), sqrt(pi^2 / 3 * 0.1 / 0.9),
               tolerance = 1e-6)
  ## a count has no latent scale to take one from, and says what to give
  expect_error(ilm_scaffold(y ~ arm + (1 | id),
                            design = list(arm = c("control", "treatment"),
                                          time = 1:4),
                            within = "time", n_unit = 60, family = "poisson",
                            cells = c(control = 3, treatment = 4), icc = 0.1,
                            verbose = FALSE), "Give `re_sd`")
})

test_that("power for a planned multinomial study rises with its size", {
  skip_on_cran()
  cells <- rbind(control = c(none = 0.5, some = 0.3, full = 0.2),
                 treatment = c(none = 0.3, some = 0.35, full = 0.35))
  p <- ilm_power_design(y ~ arm, design = list(arm = c("control", "treatment")),
                        n_unit = c(60, 240), family = "multinomial",
                        cells = cells, term = "arm", sims = 60, seed = 2,
                        verbose = FALSE, progress = FALSE)
  expect_identical(attr(p, "test"), "Wald chi-square")
  expect_true(isTRUE(attr(p, "redrawn")))
  expect_true(all(p$converged > 0.9))
  expect_gt(p$power[2], p$power[1])
  expect_output(print(p), "fresh draw of the planned design")
})

test_that("a planned small factorial is always estimable", {
  skip_on_cran()
  cf <- c("(Intercept)" = 0, aa2 = 0.5, aa3 = 0.5, bb2 = 0, bb3 = 0,
          "aa2:bb2" = 0, "aa3:bb2" = 0, "aa2:bb3" = 0, "aa3:bb3" = 0)
  p <- ilm_power_design(y ~ a * b, n_unit = c(18, 36), coefs = cf, sd = 1,
                        design = list(a = c("a1", "a2", "a3"),
                                      b = c("b1", "b2", "b3")),
                        term = "a", sims = 40, verbose = FALSE,
                        progress = FALSE)
  expect_true(all(p$converged == 1))
  expect_identical(attr(p, "test"), "F")
})

test_that("a scaffold can plan for a REML analysis", {
  s <- ilm_scaffold(y ~ arm + (1 | id),
                    design = list(arm = c("control", "treatment"), time = 1:3),
                    within = "time", n_unit = 30, sd = 1, icc = 0.3,
                    cells = c(control = 0, treatment = 0.5), reml = TRUE,
                    verbose = FALSE)
  expect_true(isTRUE(s$reml))
  ## the assumptions are imposed whatever the estimator
  expect_equal(unname(coef(s)), c(0, 0.5), tolerance = 1e-8)
  ## and the simulated studies are refitted the same way
  st <- ilm_refit_stub(s)
  expect_true(st$reml)
})
