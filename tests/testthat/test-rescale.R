## The fit rescales its fixed-effect columns in extreme units and converts
## back exactly (ilm_rescale.R). Income in dollars made a zero-inflated
## negative binomial fail every time; in thousands it fitted. Now the units do
## not matter: the answers are those of rescaling by hand, in every piece of
## the fit. A column whose SD is inside ilm_rescale_band is left as it is, so
## the hand-rescaled comparison divides income by its own SD: the fit on
## dollars then runs on exactly the columns the hand fit gives it.

inc_data <- function(seed = 8, n = 500) {
  set.seed(seed)
  d <- data.frame(g = factor(sample(25, n, TRUE)), income = round(stats::rnorm(n, 33000, 9000)),
                  age = stats::rnorm(n, 45, 12), x = stats::runif(n))
  eta <- 0.4 + 0.00003 * (d$income - 33000) + 0.01 * (d$age - 45) +
    stats::rnorm(25, 0, 0.4)[d$g] + sin(3 * d$x)
  d$y <- stats::rpois(n, exp(eta))
  d$inc_k <- d$income / stats::sd(d$income)          # by hand, to unit SD
  d$e <- stats::runif(n, 0.5, 2)
  d
}
q <- function(e) suppressWarnings(suppressMessages(e))
## the rescaling variant for one expression, put back afterwards whatever
## happens (on.exit in a function, where it is certain to run)
with_rescale <- function(mode, expr) {
  op <- options(illume.rescale = mode); on.exit(options(op))
  expr
}
## coefficients and SEs of a dollars fit against the hand-rescaled fit, on
## the dollars scale, as the largest relative difference
rel_gap <- function(a, b, d, k = "income", kk = "inc_k") {
  u <- stats::sd(d$income)
  ca <- coef(a); cb <- coef(b); cb[kk] <- cb[kk] / u
  sa <- sqrt(diag(vcov(a))); sb <- sqrt(diag(vcov(b))); sb[match(kk, names(cb))] <- sb[match(kk, names(cb))] / u
  max(abs(ca - unname(cb)) / abs(cb), abs(sa - unname(sb)) / sb)
}

test_that("the dollars fit is the hand-rescaled fit", {
  d <- inc_data()
  a <- q(ilm_model(y ~ income + age + (1 | g), data = d, family = "poisson", verbose = FALSE))
  b <- q(ilm_model(y ~ inc_k + age + (1 | g), data = d, family = "poisson", verbose = FALSE))
  expect_true(a$rescale$x$any)                       # income is rescaled...
  expect_null(b$rescale)                             # ...and nothing in b is
  expect_lt(rel_gap(a, b, d), 1e-6)
  expect_equal(as.numeric(logLik(a)), as.numeric(logLik(b)), tolerance = 1e-8)
  expect_equal(as.numeric(predict(a, newdata = d[1:5, ])),
               as.numeric(predict(b, newdata = d[1:5, ])), tolerance = 1e-7)
})

test_that("a model whose columns are all well scaled is fitted exactly as without rescaling", {
  ## E3 of the pre-registered study, by construction: SDs of 1, 12 and 0.29
  ## are inside the band, so "scale" leaves every column as it is
  d <- inc_data()
  a <- q(ilm_model(y ~ inc_k + age + x + (1 | g), data = d, family = "poisson", verbose = FALSE))
  b <- with_rescale("none", q(ilm_model(y ~ inc_k + age + x + (1 | g), data = d, family = "poisson",
                                        verbose = FALSE)))
  expect_identical(a$opt$par, b$opt$par)
  expect_identical(a$opt$objective, b$opt$objective)
  expect_identical(vcov(a), vcov(b))
  ## the band's edges: a column is rescaled only outside [0.01, 100]
  sc <- illume:::ilm_col_scales(cbind(1, a = c(0, 50, 100, 150, 200), b = c(0, 0.004, 0.008, 0.012, 0.016),
                                      c = c(0, 300, 600, 900, 1200)))
  expect_identical(sc$s[1:2], c(1, 1))
  expect_true(all(sc$s[3:4] != 1))
})

test_that("a zero-inflated negative binomial in dollars now fits", {
  set.seed(1); n <- 500
  d <- data.frame(g = factor(sample(25, n, TRUE)), income = round(stats::rnorm(n, 33000, 9000)),
                  age = stats::rnorm(n, 45, 12))
  eta <- 0.4 + 0.00003 * (d$income - 33000) + 0.01 * (d$age - 45) + stats::rnorm(25, 0, 0.4)[d$g]
  pz <- stats::plogis(-1 + 0.00004 * (d$income - 33000))
  d$y <- ifelse(stats::runif(n) < pz, 0, stats::rnbinom(n, mu = exp(eta), size = 2))
  f <- q(ilm_model(y ~ income + age + (1 | g), data = d, family = "nbinom",
                   ziformula = ~ income, verbose = FALSE))
  expect_true(f$ok)
  expect_true(all(is.finite(sqrt(diag(vcov(f))))))
  ## and without the rescaling it did not
  f0 <- with_rescale("none", q(ilm_model(y ~ income + age + (1 | g), data = d, family = "nbinom",
                                         ziformula = ~ income, verbose = FALSE)))
  expect_false(f0$ok)
})

test_that("a smooth beside a covariate in dollars is unchanged, fitted values and edf", {
  d <- inc_data(n = 400)
  a <- q(ilm_model(y ~ income + s(x) + (1 | g), data = d, family = "poisson", verbose = FALSE))
  b <- q(ilm_model(y ~ inc_k + s(x) + (1 | g), data = d, family = "poisson", verbose = FALSE))
  expect_equal(unname(coef(a)[["income"]] * stats::sd(d$income)), unname(coef(b)[["inc_k"]]), tolerance = 1e-6)
  expect_equal(unname(unlist(a$edf)), unname(unlist(b$edf)), tolerance = 1e-5)
  expect_equal(as.numeric(predict(a, newdata = d[1:8, ])),
               as.numeric(predict(b, newdata = d[1:8, ])), tolerance = 1e-6)
})

test_that("an offset is not rescaled, and a REML fit's logLik does not move", {
  d <- inc_data(n = 300)
  a <- q(ilm_model(y ~ income + offset(log(e)) + (1 | g), data = d, family = "poisson", verbose = FALSE))
  b <- q(ilm_model(y ~ inc_k + offset(log(e)) + (1 | g), data = d, family = "poisson", verbose = FALSE))
  expect_lt(rel_gap(a, b, d), 1e-6)
  d$z <- 2 + 0.0001 * d$income + stats::rnorm(25)[d$g] + stats::rnorm(nrow(d))
  ra <- q(ilm_model(z ~ income + (1 | g), data = d, family = "gaussian", reml = TRUE, verbose = FALSE))
  rb <- q(ilm_model(z ~ inc_k + (1 | g), data = d, family = "gaussian", reml = TRUE, verbose = FALSE))
  r0 <- with_rescale("none", q(ilm_model(z ~ income + (1 | g), data = d, family = "gaussian",
                                         reml = TRUE, verbose = FALSE)))
  ## the restricted likelihood is in the user's units, whatever the fit ran on
  expect_equal(as.numeric(logLik(ra)), as.numeric(logLik(r0)), tolerance = 1e-6)
  expect_equal(unname(coef(ra)[["income"]] * stats::sd(d$income)), unname(coef(rb)[["inc_k"]]), tolerance = 1e-6)
})

test_that("separation is caught the same with and without the rescaling", {
  d <- data.frame(x = seq(-2, 2, length.out = 60)); d$y <- as.integer(d$x > 0)
  d$big <- d$x * 1e4
  for (mode in c("none", "scale")) {
    f <- with_rescale(mode, q(ilm_model(y ~ big, data = d, family = "binomial", verbose = FALSE)))
    expect_identical(f$checks$status[f$checks$check == "separation"], "FAIL")
    expect_true("big" %in% f$separation$flat$coef)
  }
})

test_that("a random slope on a covariate in dollars has rescaling named first", {
  fit <- list(bars = list(quote(1 + income | g)),
              model = data.frame(income = stats::rnorm(50, 33000, 9000), g = factor(1:50)))
  r <- illume:::ilm_rem_rescale(fit)
  expect_length(r, 1L)
  expect_match(r[[1]]$remedy, "rescale income by hand", fixed = TRUE)
  expect_match(r[[1]]$remedy, "divide it by 1e+04", fixed = TRUE)
  fit$model$income <- stats::rnorm(50)
  expect_length(illume:::ilm_rem_rescale(fit), 0L)
})

test_that("a term held at its boundary keeps finite SEs and its Satterthwaite df", {
  ## the held term's rows of the covariance are NA, and mapping the whole
  ## matrix spread that NA to every entry: Satterthwaite then fell back to z
  set.seed(7); n <- 200                              # as test-fixed-only's sim_lm(7)
  d <- data.frame(x = stats::rnorm(n), z = factor(sample(c("a", "b", "c"), n, TRUE)))
  d$y <- 2 + 1.5 * d$x - 0.8 * (d$z == "b") + stats::rnorm(n, 0, 1.2)
  d$g <- factor(sample(20, n, TRUE))                 # no group variance at all
  d$x <- d$x * 1e4                                   # and x in large units
  f <- q(ilm_model(y ~ x + (1 | g), data = d, family = "gaussian", verbose = FALSE))
  expect_true("g" %in% f$hessian_held)
  expect_true(all(is.finite(sqrt(diag(vcov(f))))))
  ct <- ilm_coef_table(f)
  expect_identical(attr(ct, "df_method"), "satterthwaite")
  expect_true(all(is.finite(ct$df)))
})
