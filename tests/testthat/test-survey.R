## a stratified two-stage sample: weights differ by stratum, as they do when
## some group is deliberately oversampled, and observations cluster within PSU
svy_data <- function(nstr = 4L, npsu = 15L, nper = 12L, seed = 1L,
                     w = c(80, 40, 15, 5), sd_u = 0.7) {
  set.seed(seed)
  n <- nstr * npsu * nper
  d <- data.frame(st = factor(rep(seq_len(nstr), each = npsu * nper)),
                  psu = factor(rep(seq_len(nstr * npsu), each = nper)),
                  x = rnorm(n))
  d$w <- w[as.integer(d$st)]
  d$y <- 1 + 0.5 * d$x +
    rep(rnorm(nstr * npsu, 0, sd_u), each = nper) + rnorm(n)
  d
}

test_that("a design records the sample and its degrees of freedom", {
  d <- svy_data()
  des <- ilm_design(d, weights = ~ w, ids = ~ psu, strata = ~ st)
  expect_s3_class(des, "ilm_design")
  expect_equal(des$n, nrow(d))
  expect_equal(des$n_psu, 60L)
  expect_equal(des$n_strata, 4L)
  ## PSUs less strata, which is the number that sets the reference
  ## distribution -- not the 720 rows
  expect_equal(des$df, 56L)
  expect_output(print(des), "56")
  ## unclustered: every row is its own sampling unit
  du <- ilm_design(d, weights = ~ w)
  expect_equal(du$n_psu, nrow(d))
  expect_equal(du$df, nrow(d) - 1L)
})

test_that("it agrees with survey::svyglm", {
  skip_if_not_installed("survey")
  d <- svy_data()
  des <- ilm_design(d, weights = ~ w, ids = ~ psu, strata = ~ st)
  f <- ilm_model(y ~ x, data = d, family = "gaussian", design = des,
                 verbose = FALSE)
  sc <- ilm_svy_coef(f)
  sd_ <- survey::svydesign(ids = ~ psu, strata = ~ st, weights = ~ w,
                           data = d, nest = TRUE)
  m <- survey::svyglm(y ~ x, design = sd_, family = stats::gaussian())
  expect_equal(unname(coef(f)), unname(coef(m)), tolerance = 1e-6)
  expect_equal(sc$se, unname(sqrt(diag(vcov(m)))), tolerance = 1e-7)
  expect_equal(attr(sc, "df"), as.integer(survey::degf(sd_)))

  ## unclustered and unstratified, where the design reduces to weighting alone
  set.seed(4); n2 <- 800L
  d2 <- data.frame(x = rnorm(n2), w = runif(n2, 1, 20))
  d2$y <- 1 + 0.4 * d2$x + rnorm(n2)
  des2 <- ilm_design(d2, weights = ~ w)
  f2 <- ilm_model(y ~ x, data = d2, family = "gaussian", design = des2,
                  verbose = FALSE)
  mm <- survey::svyglm(y ~ x, family = stats::gaussian(),
                       design = survey::svydesign(ids = ~ 1, weights = ~ w,
                                                  data = d2))
  expect_equal(ilm_svy_coef(f2)$se, unname(sqrt(diag(vcov(mm)))),
               tolerance = 1e-7)
})

test_that("it agrees with survey for binomial and poisson too", {
  skip_if_not_installed("survey")
  for (fam in c("binomial", "poisson")) {
    set.seed(6); n <- 900L
    d <- data.frame(x = rnorm(n), psu = factor(rep(1:45, each = 20)),
                    st = factor(rep(1:3, each = 300)))
    d$w <- c(50, 20, 8)[as.integer(d$st)]
    u <- rep(rnorm(45, 0, .5), each = 20)
    d$y <- if (fam == "binomial") rbinom(n, 1, stats::plogis(.2 + .7 * d$x + u))
           else rpois(n, exp(.4 + .5 * d$x + u))
    de <- ilm_design(d, weights = ~ w, ids = ~ psu, strata = ~ st)
    f <- ilm_model(y ~ x, data = d, family = fam, design = de, verbose = FALSE)
    sd_ <- survey::svydesign(ids = ~ psu, strata = ~ st, weights = ~ w,
                             data = d, nest = TRUE)
    m <- survey::svyglm(y ~ x, design = sd_,
                        family = if (fam == "binomial")
                          stats::quasibinomial() else stats::quasipoisson())
    expect_equal(unname(coef(f)), unname(coef(m)), tolerance = 1e-5,
                 info = fam)
    expect_equal(ilm_svy_coef(f)$se, unname(sqrt(diag(vcov(m)))),
                 tolerance = 1e-6, info = fam)
  }
})

test_that("the model-based standard errors are the thing being replaced", {
  d <- svy_data()
  des <- ilm_design(d, weights = ~ w, ids = ~ psu, strata = ~ st)
  f <- ilm_model(y ~ x, data = d, family = "gaussian", design = des,
                 verbose = FALSE)
  sc <- ilm_svy_coef(f)
  ## a weighted likelihood believes it saw sum(w) observations, so its errors
  ## are far too small; clustering makes the gap larger still
  expect_gt(min(sc$se / sc$se_model), sqrt(sum(d$w) / nrow(d)) * 0.9)
  expect_gt(max(sc$se / sc$se_model), 10)
  ## and the coefficients are NOT the problem -- they are the design-consistent
  ## ones, which is what makes the mistake easy to miss
  expect_equal(sc$estimate[sc$term == "x"], 0.5, tolerance = 0.2)
  ## the fit carries its design, so no second argument is needed
  expect_equal(ilm_svy_coef(f)$se, ilm_svy_coef(f, des)$se)
  ## passing weights that look like sampling weights to a plain fit is warned
  ## about, and the warning now names the remedy
  pre <- ilm_model(y ~ x, data = d, family = "gaussian", weights = w,
                   verbose = FALSE)
  ck <- pre$checks[pre$checks$check == "weights_type", , drop = FALSE]
  ## these weights are whole numbers, so being integers proves nothing; what
  ## gives them away is that they total a population rather than a count of
  ## replicates
  expect_equal(ck$status, "WARN")
  expect_match(ck$cause, "population rather than a count", fixed = TRUE)
  expect_match(ck$suggestion, "ilm_design", fixed = TRUE)
  ## in a FIXED-effects fit the coefficients are the design-consistent ones
  expect_match(ck$cause, "COEFFICIENTS are fine in a fixed-effects fit",
               fixed = TRUE)
  ## in a MIXED one they are not, because the weight multiplies the
  ## conditional likelihood inside the Laplace integral
  pm <- ilm_model(y ~ x + (1 | psu), data = d, family = "gaussian",
                  weights = w, verbose = FALSE)
  ckm <- pm$checks[pm$checks$check == "weights_type", , drop = FALSE]
  expect_match(ckm$cause, "under-shrunk", fixed = TRUE)
  ## grouped binomial trial counts are genuinely large and are not flagged
  set.seed(9); nb <- 200L
  db <- data.frame(x = rnorm(nb), trials = rep(60, nb))
  db$y <- rbinom(nb, 60, stats::plogis(0.2 + 0.5 * db$x)) / 60
  fb <- ilm_model(y ~ x, data = db, family = "binomial", weights = trials,
                  verbose = FALSE)
  expect_equal(fb$checks$status[fb$checks$check == "weights_type"], "OK")
})

test_that("a stratum with one sampling unit is handled and reported", {
  d <- svy_data(nstr = 3L, npsu = 6L)
  ## give one stratum a single PSU
  d$st <- factor(ifelse(d$psu == "1", "lonely", as.character(d$st)))
  d$w <- 20
  des <- ilm_design(d, weights = ~ w, ids = ~ psu, strata = ~ st)
  f <- ilm_model(y ~ x, data = d, family = "gaussian", design = des,
                 verbose = FALSE)
  expect_error(ilm_svy_vcov(f, lonely = "fail"), "single primary sampling unit")
  v_adj <- ilm_svy_vcov(f, lonely = "adjust")
  v_cer <- ilm_svy_vcov(f, lonely = "certainty")
  expect_equal(attr(v_adj, "lonely"), 1L)
  ## treating it as a certainty unit contributes nothing, so it is the smaller
  expect_lte(sum(diag(v_cer)), sum(diag(v_adj)))
  expect_output(print(ilm_svy_coef(f)), "single sampling unit")
})

test_that("the design and the fit have to describe the same rows", {
  d <- svy_data()
  des <- ilm_design(d, weights = ~ w, ids = ~ psu, strata = ~ st)
  ## weights cannot come from two places at once
  expect_error(ilm_model(y ~ x, data = d, family = "gaussian", design = des,
                         weights = w, verbose = FALSE), "not both")
  expect_error(ilm_model(y ~ x, data = d, family = "gaussian",
                         design = mtcars, verbose = FALSE), "ilm_design")
  expect_error(ilm_design(d, weights = ~ nope), "not in the data")
  d0 <- d; d0$w[1] <- 0
  expect_error(ilm_design(d0, weights = ~ w), "must be positive")
  ## PSU labels that restart inside each stratum are not nested
  d2 <- d
  d2$psu <- factor(rep(rep(1:15, each = 12), 4))
  expect_error(ilm_design(d2, weights = ~ w, ids = ~ psu, strata = ~ st),
               "more than one stratum")

  ## rows dropped as missing come out of the design too, or the weights line
  ## up against the wrong observations
  dm <- d; dm$x[c(5, 200, 700)] <- NA
  desm <- ilm_design(dm, weights = ~ w, ids = ~ psu, strata = ~ st)
  fm <- ilm_model(y ~ x, data = dm, family = "gaussian", design = desm,
                  verbose = FALSE)
  expect_equal(nrow(fm$X), nrow(d) - 3L)
  expect_equal(fm$design$n, nrow(d) - 3L)
  expect_silent(ilm_svy_coef(fm))
})

test_that("a random effect is refused, since the score is not a row sum", {
  d <- svy_data()
  des <- ilm_design(d, weights = ~ w, ids = ~ psu, strata = ~ st)
  fr <- ilm_model(y ~ x + (1 | psu), data = d, family = "gaussian",
                  design = des, verbose = FALSE)
  expect_error(ilm_svy_coef(fr), "random effects")
})
