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
