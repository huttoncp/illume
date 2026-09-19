# Families other than the multinomial.  The random-effect machinery is shared,
# so these check that the family-specific pieces -- likelihood, response
# validation, dispersion, link -- behave, not that mixed models work at all.

sim_glmm <- function(seed = 1, n = 900, ng = 30, sd_g = 0.6) {
  set.seed(seed)
  dd <- data.frame(g = factor(sample(ng, n, TRUE)), x = stats::rnorm(n))
  b <- stats::rnorm(ng, 0, sd_g)[as.integer(dd$g)]
  lp <- 0.5 + 0.7 * dd$x + b
  dd$y_gauss  <- lp + stats::rnorm(n, 0, 0.5)
  dd$y_binom  <- stats::rbinom(n, 1, 1 / (1 + exp(-lp)))
  dd$y_pois   <- stats::rpois(n, exp(lp))
  dd$y_nbinom <- stats::rnbinom(n, size = 2, mu = exp(lp))
  dd
}

test_that("every family is described consistently", {
  for (fm in c("gaussian", "binomial", "poisson", "nbinom", "multinomial")) {
    f <- ilm_family(fm)
    expect_equal(f$name, fm)
    expect_true(is.function(f$nll))
    expect_true(is.function(f$linkinv))
    expect_true(is.function(f$C_of))
    expect_length(f$disp_names, f$n_disp)
  }
  expect_equal(ilm_family("gaussian")$n_disp, 1L)   # residual SD
  expect_equal(ilm_family("nbinom")$n_disp, 1L)     # overdispersion
  expect_equal(ilm_family("poisson")$n_disp, 0L)
  expect_equal(ilm_family("binomial")$n_disp, 0L)
  expect_error(ilm_family("wibble"))
})

test_that("only the multinomial uses more than one linear predictor", {
  expect_equal(ilm_family("multinomial")$C_of(5L), 4L)
  for (fm in c("gaussian", "binomial", "poisson", "nbinom"))
    expect_equal(ilm_family(fm)$C_of(5L), 1L)
})

test_that("gaussian recovers its parameters and reports a residual SD", {
  dd <- sim_glmm(2)
  f <- ilm_model(y_gauss ~ x + (1 | g), data = dd, family = "gaussian",
                 verbose = FALSE)
  expect_equal(f$opt$convergence, 0L)
  expect_equal(unname(coef(f)[["x"]]), 0.7, tolerance = 0.15)
  expect_named(f$dispersion, "log_sigma")
  expect_equal(unname(f$dispersion), 0.5, tolerance = 0.12)
  expect_equal(sqrt(f$Sigma$g[1, 1]), 0.6, tolerance = 0.25)
})

test_that("poisson and binomial fit without a dispersion parameter", {
  dd <- sim_glmm(3)
  fp <- ilm_model(y_pois ~ x + (1 | g), data = dd, family = "poisson", verbose = FALSE)
  fb <- ilm_model(y_binom ~ x + (1 | g), data = dd, family = "binomial", verbose = FALSE)
  expect_null(fp$dispersion)
  expect_null(fb$dispersion)
  expect_equal(unname(coef(fp)[["x"]]), 0.7, tolerance = 0.2)
  expect_equal(unname(coef(fb)[["x"]]), 0.7, tolerance = 0.35)  # binary is noisy
})

test_that("negative binomial recovers its overdispersion parameter", {
  dd <- sim_glmm(4, n = 1500)
  f <- ilm_model(y_nbinom ~ x + (1 | g), data = dd, family = "nbinom", verbose = FALSE)
  expect_named(f$dispersion, "log_k")
  expect_equal(unname(f$dispersion), 2, tolerance = 1.0)   # k is hard to pin down
  expect_gt(unname(f$dispersion), 0)
})

test_that("coefficient names carry no category prefix for univariate families", {
  dd <- sim_glmm(5)
  f <- ilm_model(y_gauss ~ x + (1 | g), data = dd, family = "gaussian", verbose = FALSE)
  expect_equal(names(coef(f)), c("(Intercept)", "x"))
  expect_false(any(grepl(":", names(coef(f)))))
})

test_that("predictions are one column on the response scale", {
  dd <- sim_glmm(6)
  fp <- ilm_model(y_pois ~ x + (1 | g), data = dd, family = "poisson", verbose = FALSE)
  P <- predict(fp)
  expect_equal(ncol(P), 1L)
  expect_true(all(P > 0))                       # log link, so strictly positive
  L <- predict(fp, type = "link")
  expect_equal(ncol(L), 1L)
  expect_equal(unname(P[, 1]), exp(unname(L[, 1])), tolerance = 1e-8)
  expect_error(predict(fp, type = "class"), "multinomial")
})

test_that("the response is validated against the family", {
  dd <- sim_glmm(7)
  # positive but fractional: isolates the whole-number check from the sign check
  dd$frac <- abs(dd$y_gauss) + 0.5
  expect_error(ilm_model(frac ~ x + (1 | g), data = dd, family = "poisson",
                         verbose = FALSE), "whole-number")
  dd$neg <- -1L
  expect_error(ilm_model(neg ~ x + (1 | g), data = dd, family = "poisson",
                         verbose = FALSE), "non-negative")
  dd$big <- 5
  expect_error(ilm_model(big ~ x + (1 | g), data = dd, family = "binomial",
                         verbose = FALSE), "0/1")
})

test_that("a 3-level factor is refused by the univariate families", {
  dd <- sim_glmm(8)
  dd$f3 <- factor(sample(c("a", "b", "c"), nrow(dd), TRUE))
  expect_error(ilm_model(f3 ~ x + (1 | g), data = dd, family = "binomial",
                         verbose = FALSE), "multinomial")
})

test_that("sparse-category checks apply only to the multinomial", {
  dd <- sim_glmm(9)
  f <- ilm_model(y_gauss ~ x + (1 | g), data = dd, family = "gaussian",
                 verbose = FALSE)
  expect_false("category_counts" %in% f$checks$check)
  expect_true("latent_budget" %in% f$checks$check)
})

test_that("summary and anova adapt to the family", {
  dd <- sim_glmm(10)
  f <- ilm_model(y_gauss ~ x + (1 | g), data = dd, family = "gaussian",
                 verbose = FALSE)
  out <- capture.output(print(summary(f)))
  expect_true(any(grepl("Family: gaussian", out)))
  expect_true(any(grepl("Dispersion", out)))
  expect_false(any(grepl("SUM-TO-ZERO", out)))   # meaningless with one dimension
  a <- ilm_anova(f, type = 3)
  expect_equal(a$Df, 1L)                          # one column, one dimension
})

test_that("the binomial likelihood does not overflow at extreme eta", {
  # logspace_add is exact where log(1 + exp(eta)) would return Inf
  expect_true(is.finite(RTMB::logspace_add(0, 800)))
  expect_equal(RTMB::logspace_add(0, 800), 800, tolerance = 1e-8)
  expect_false(is.finite(log(1 + exp(800))))      # the naive form
})

# A two-category outcome arrives in whatever type the user's data happens to
# hold: 0/1, a factor, TRUE/FALSE, or a character column read from a csv.  All
# four mean the same thing and must give the same fit.

sim_bin2 <- function(seed = 7, n = 600) {
  set.seed(seed)
  dd <- data.frame(x = stats::rnorm(n))
  dd$num <- stats::rbinom(n, 1, 1 / (1 + exp(-(-0.3 + 0.9 * dd$x))))
  dd$fac <- factor(c("no", "yes")[dd$num + 1], levels = c("no", "yes"))
  dd$lgl <- dd$num == 1
  dd$chr <- as.character(dd$fac)
  dd
}

test_that("0/1, factor, logical and character responses agree", {
  dd <- sim_bin2()
  fits <- lapply(c("num", "fac", "lgl", "chr"), function(v)
    ilm_model(stats::as.formula(paste(v, "~ x")), data = dd,
              family = "binomial", verbose = FALSE))
  for (f in fits[-1]) expect_equal(coef(f), coef(fits[[1]]), tolerance = 1e-8)
  # the second level is the modelled outcome, as in glm()
  expect_equal(unname(coef(fits[[2]])), unname(coef(fits[[1]])), tolerance = 1e-8)
})

test_that("binomial matches glm() on the same data", {
  dd <- sim_bin2(8)
  f <- ilm_model(fac ~ x, data = dd, family = "binomial", verbose = FALSE)
  g <- stats::glm(fac ~ x, data = dd, family = stats::binomial)
  expect_equal(unname(coef(f)), unname(coef(g)), tolerance = 1e-5)
  expect_equal(unname(sqrt(diag(vcov(f)))), unname(sqrt(diag(vcov(g)))),
               tolerance = 1e-4)
})

test_that("a factor response is refused by families that cannot use one", {
  dd <- sim_bin2(9)
  bad <- function(fam)
    tryCatch(ilm_model(fac ~ x, data = dd, family = fam, verbose = FALSE),
             error = conditionMessage)
  expect_match(bad("gaussian"), "needs a numeric response")
  expect_match(bad("gaussian"), 'family = "binomial"')
  expect_match(bad("poisson"), "needs a numeric response")
})

test_that("binomial refuses more than two categories, multinomial fewer than three", {
  dd <- sim_bin2(10)
  dd$f3 <- factor(sample(c("a", "b", "c"), nrow(dd), TRUE))
  expect_error(ilm_model(f3 ~ x, data = dd, family = "binomial", verbose = FALSE),
               "2 categories")
  expect_error(ilm_model(fac ~ x, data = dd, family = "multinomial", verbose = FALSE),
               "at least 3 outcome categories")
})

# A binomial coefficient describes the probability of ONE of the two categories,
# and nothing in the coefficient names says which -- unlike a multinomial fit,
# where every name carries its category.  The summary has to say so.

test_that("summary names the category a binomial model is modelling", {
  dd <- sim_bin2(11)
  hdr <- function(v) {
    f <- ilm_model(stats::as.formula(paste(v, "~ x")), data = dd,
                   family = "binomial", verbose = FALSE)
    paste(capture.output(print(summary(f))), collapse = " ")
  }
  # the second level is modelled, so "yes" is the outcome and "no" the reference
  expect_match(hdr("fac"), "modelling P(fac = 'yes'), with 'no' as the reference",
               fixed = TRUE)
  expect_match(hdr("chr"), "modelling P(chr = 'yes')", fixed = TRUE)
  # a logical is not a string, so its levels are printed unquoted
  expect_match(hdr("lgl"), "modelling P(lgl = TRUE), with FALSE as the reference",
               fixed = TRUE)
  # a bare 0/1 response has no level names to report
  expect_match(hdr("num"), "modelling P(num = 1)", fixed = TRUE)
})

test_that("a proportion response is described as one", {
  set.seed(12)
  dd <- data.frame(x = stats::rnorm(60), k = rep(20, 60))
  dd$p <- stats::rbinom(60, 20, 0.4) / 20
  f <- ilm_model(p ~ x, data = dd, family = "binomial", weights = k,
                 verbose = FALSE)
  out <- paste(capture.output(print(summary(f))), collapse = " ")
  expect_match(out, "response is a proportion", fixed = TRUE)
})

test_that("other families do not gain the binomial line", {
  dd <- sim_bin2(13)
  dd$z <- stats::rnorm(nrow(dd))
  f <- ilm_model(z ~ x, data = dd, family = "gaussian", verbose = FALSE)
  out <- paste(capture.output(print(summary(f))), collapse = " ")
  expect_false(grepl("modelling P(", out, fixed = TRUE))
})
