ord_data <- function(n = 800L, seed = 1L, link = "logit", b = 0.8) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n), g = factor(sample(c("a", "b"), n, TRUE)))
  e <- switch(link, logit = stats::rlogis(n), probit = rnorm(n),
              cloglog = log(-log(runif(n))))
  z <- b * d$x - 0.5 * (d$g == "b") + e
  d$y <- factor(cut(z, c(-Inf, -0.8, 0.9, Inf),
                    labels = c("lo", "mid", "hi")), ordered = TRUE)
  d
}
ord_fam <- c(logit = "ordinal", probit = "ordinal_probit",
             cloglog = "ordinal_cloglog")

test_that("a cumulative link model fits, with thresholds and no intercept", {
  f <- ilm_model(y ~ x + g, data = ord_data(), family = "ordinal",
                 verbose = FALSE)
  expect_s3_class(f, "ilm_model")
  expect_true(f$ok)
  expect_true(isTRUE(f$ordinal))
  expect_equal(f$J, 3L)
  ## the thresholds are the intercepts, so the design carries none
  expect_false("(Intercept)" %in% colnames(f$X))
  expect_equal(ncol(f$X), 2L)
  ## J - 1 of them, increasing, and named for the cuts they sit between
  expect_length(f$zeta, 2L)
  expect_true(all(diff(f$zeta) > 0))
  expect_equal(names(f$zeta), c("lo|mid", "mid|hi"))
  th <- ilm_thresholds(f)
  expect_equal(nrow(th), 2L)
  expect_true(all(th$se > 0))
  expect_true(all(th$lower < th$estimate & th$estimate < th$upper))
  ## and they recover the cuts they were generated from
  expect_equal(unname(f$zeta), c(-0.8, 0.9), tolerance = 0.25)
  expect_equal(unname(coef(f)), c(0.8, -0.5), tolerance = 0.3)
})

test_that("the ordering is structural, not something the optimiser respects", {
  ## many random starts, and no fit can produce an out-of-order threshold
  for (s in 1:8) {
    f <- ilm_model(y ~ x, data = ord_data(200L, seed = s), family = "ordinal",
                   verbose = FALSE)
    expect_true(all(diff(f$zeta) > 0), info = paste("seed", s))
  }
})

test_that("all three links agree with ordinal::clm", {
  skip_if_not_installed("ordinal")
  for (lk in names(ord_fam)) {
    d <- ord_data(link = lk)
    f <- ilm_model(y ~ x + g, data = d, family = ord_fam[[lk]],
                   verbose = FALSE)
    m <- ordinal::clm(y ~ x + g, data = d, link = lk)
    cf <- coef(m); nc <- length(cf) - 2L
    se <- sqrt(diag(vcov(m)))
    expect_equal(unname(coef(f)), unname(cf[-(1:nc)]), tolerance = 1e-4,
                 info = lk)
    expect_equal(unname(f$zeta), unname(cf[1:nc]), tolerance = 1e-4, info = lk)
    expect_equal(unname(sqrt(diag(vcov(f)))), unname(se[-(1:nc)]),
                 tolerance = 1e-4, info = lk)
    expect_equal(unname(f$zeta_se), unname(se[1:nc]), tolerance = 1e-4,
                 info = lk)
    expect_equal(as.numeric(logLik(f)), as.numeric(logLik(m)), tolerance = 1e-6,
                 info = lk)
  }
})

test_that("a mixed cumulative link model agrees with ordinal::clmm", {
  skip_if_not_installed("ordinal")
  set.seed(7); ng <- 60L; ni <- 12L; n <- ng * ni
  d <- data.frame(id = factor(rep(seq_len(ng), each = ni)), x = rnorm(n))
  b <- rnorm(ng, 0, 0.8)
  z <- 0.8 * d$x + b[as.integer(d$id)] + stats::rlogis(n)
  d$y <- factor(cut(z, c(-Inf, -0.8, 0.9, Inf),
                    labels = c("lo", "mid", "hi")), ordered = TRUE)
  f <- ilm_model(y ~ x + (1 | id), data = d, family = "ordinal",
                 verbose = FALSE)
  m <- ordinal::clmm(y ~ x + (1 | id), data = d)
  expect_equal(unname(coef(f)), unname(m$beta), tolerance = 1e-3)
  expect_equal(unname(f$zeta), unname(m$alpha), tolerance = 1e-3)
  expect_equal(sqrt(f$Sigma$id[1, 1]), unname(m$ST$id[1, 1]), tolerance = 1e-3)
  expect_equal(as.numeric(logLik(f)), as.numeric(logLik(m)), tolerance = 1e-4)
})

test_that("predictions are category probabilities", {
  d <- ord_data()
  f <- ilm_model(y ~ x + g, data = d, family = "ordinal", verbose = FALSE)
  P <- predict(f)
  expect_equal(dim(P), c(nrow(d), 3L))
  expect_equal(colnames(P), c("lo", "mid", "hi"))
  expect_equal(unname(rowSums(P)), rep(1, nrow(d)), tolerance = 1e-10)
  expect_true(all(P >= 0))
  ## they track the observed shares
  expect_equal(unname(colMeans(P)),
               unname(as.numeric(prop.table(table(d$y)))), tolerance = 0.02)
  ## the link scale is the latent linear predictor, one column
  expect_equal(ncol(predict(f, type = "link")), 1L)
  ## class gives an ORDERED factor, so it compares with the response
  cls <- predict(f, type = "class")
  expect_s3_class(cls, "ordered")
  expect_equal(levels(cls), levels(d$y))
  ## new data works, and a bigger x means more probability up the scale
  nd <- data.frame(x = c(-2, 2), g = factor("a", levels = c("a", "b")))
  P2 <- predict(f, newdata = nd)
  expect_equal(dim(P2), c(2L, 3L))
  expect_gt(P2[2, "hi"], P2[1, "hi"])
  expect_lt(P2[2, "lo"], P2[1, "lo"])

  skip_if_not_installed("ordinal")
  m <- ordinal::clm(y ~ x + g, data = d)
  Pm <- predict(m, newdata = subset(d, select = -y), type = "prob")$fit
  expect_equal(max(abs(P - Pm)), 0, tolerance = 1e-5)
})

test_that("residuals and simulation know about the ordering", {
  d <- ord_data(1000L)
  f <- ilm_model(y ~ x + g, data = d, family = "ordinal", verbose = FALSE)
  u <- ilm_rqr(f, seed = 1L)
  expect_true(all(u > 0 & u < 1))
  expect_gt(suppressWarnings(stats::ks.test(u, "punif")$p.value), 0.01)
  ys <- ilm_simulate(f, 50L, seed = 2L)
  expect_true(all(ys %in% 1:3))
  expect_equal(unname(as.numeric(prop.table(table(ys)))),
               unname(as.numeric(prop.table(table(d$y)))), tolerance = 0.03)
})

test_that("marginal means work on the latent scale and refuse the response one", {
  set.seed(3); n <- 700L
  d <- data.frame(x = rnorm(n), g = factor(sample(c("a", "b", "c"), n, TRUE)))
  z <- 0.8 * d$x - 0.5 * (d$g == "b") + 0.9 * (d$g == "c") + stats::rlogis(n)
  d$y <- factor(cut(z, c(-Inf, -0.8, 0.9, Inf),
                    labels = c("lo", "mid", "hi")), ordered = TRUE)
  f <- ilm_model(y ~ x + g, data = d, family = "ordinal", verbose = FALSE)
  e <- ilm_emmeans(f, "g")
  expect_equal(nrow(e), 3L)
  ## the latent scale has no intercept, so only differences mean anything --
  ## and those differences ARE the coefficients
  cc <- ilm_contrast(e)
  expect_equal(cc$estimate[cc$contrast == "b - a"], unname(coef(f)["gb"]),
               tolerance = 1e-6)
  expect_equal(cc$estimate[cc$contrast == "c - a"], unname(coef(f)["gc"]),
               tolerance = 1e-6)
  ## on the response scale, a probability for every category in every group
  er <- ilm_emmeans(f, "g", type = "response")
  expect_equal(nrow(er), 9L)
  expect_equal(as.numeric(tapply(er$estimate, er$g, sum)), c(1, 1, 1),
               tolerance = 1e-10)
  ## and an average marginal effect still works -- one per category, since
  ## an ordered outcome has a probability for each. It used to report only
  ## the last column, unlabelled.
  a <- ilm_ame(f, "x")
  expect_equal(nrow(a), 3L)
  expect_equal(a$category, c("lo", "mid", "hi"))
  expect_equal(sum(a$estimate), 0, tolerance = 1e-8)   # probabilities sum to 1
  expect_true(a$estimate[a$category == "hi"] > 0)
  expect_true(a$estimate[a$category == "lo"] < 0)
})

test_that("an ordinal family refuses what it cannot model", {
  d <- ord_data()
  ## two categories is a binomial model
  d2 <- d; d2$y <- factor(ifelse(d$y == "lo", "lo", "hi"),
                          levels = c("lo", "hi"), ordered = TRUE)
  expect_error(ilm_model(y ~ x, data = d2, family = "ordinal", verbose = FALSE),
               "at least 3")
  ## no predictors leaves only the thresholds
  expect_error(ilm_model(y ~ 1, data = d, family = "ordinal", verbose = FALSE),
               "at least one predictor")
  ## an unordered factor is used in its stored order, and says so
  d3 <- d; d3$y <- factor(as.character(d$y))
  expect_message(ilm_model(y ~ x, data = d3, family = "ordinal",
                           verbose = FALSE), "unordered factor")
  ## and a B too small for the strongest verdict says so rather than
  ## quietly capping every violation at WARN
  fo <- ilm_model(y ~ x, data = d, family = "ordinal", verbose = FALSE)
  expect_warning(ilm_check_proportional(fo, B = 19L, progress = FALSE),
                 "unreachable")
  ## the accessors explain themselves on a model without thresholds
  fg <- ilm_model(x ~ 1, data = d, family = "gaussian", verbose = FALSE)
  expect_error(ilm_thresholds(fg), "no thresholds")
  expect_error(ilm_check_proportional(fg), "no thresholds")
})

test_that("the proportional-odds check separates the two cases", {
  ## proportional odds holds: one latent variable, one set of cuts
  f_ok <- ilm_model(y ~ x + g, data = ord_data(800L, seed = 11L),
                    family = "ordinal", verbose = FALSE)
  ## B must be large enough that the strongest verdict is reachable at all:
  ## the p-value floor is 1 / (B + 1) against a FAIL threshold of alpha / 5
  r_ok <- suppressMessages(ilm_check_proportional(f_ok, B = 199L,
                                                  progress = FALSE))
  expect_s3_class(r_ok, "ilm_prop_check")
  expect_equal(nrow(r_ok), 2L)
  expect_setequal(r_ok$term, c("x", "g"))
  expect_true(all(r_ok$p >= 0 & r_ok$p <= 1))
  expect_equal(r_ok$status[r_ok$term == "x"], "OK")

  ## and it does not: x acts on the first cut only, built cut by cut
  set.seed(4); n <- 800L
  d <- data.frame(x = rnorm(n), g = factor(sample(c("a", "b"), n, TRUE)))
  p1 <- stats::plogis(0.8 + 1.4 * d$x - 0.5 * (d$g == "b"))
  p2 <- pmin(stats::plogis(-0.9 - 0.5 * (d$g == "b")), p1)
  u <- runif(n); yi <- 1L + (u < p1) + (u < p2)
  d$y <- factor(c("lo", "mid", "hi")[yi], levels = c("lo", "mid", "hi"),
                ordered = TRUE)
  f_bad <- ilm_model(y ~ x + g, data = d, family = "ordinal", verbose = FALSE)
  r_bad <- suppressMessages(ilm_check_proportional(f_bad, B = 199L,
                                                   progress = FALSE))
  expect_equal(r_bad$status[r_bad$term == "x"], "FAIL")
  expect_gt(r_bad$max_gap[r_bad$term == "x"],
            r_ok$max_gap[r_ok$term == "x"])
  ## the message names the remedy, and the remedy is in this package
  expect_message(ilm_check_proportional(f_bad, B = 199L, progress = FALSE),
                 "multinomial")
  ## which does win on AIC here, and lose when the assumption holds
  expect_lt(AIC(ilm_model(y ~ x + g, data = d, family = "multinomial",
                          verbose = FALSE)), AIC(f_bad))
})

test_that("the pre-fit checks size an ordinal fit as the fit is sized", {
  ## an ordered response has J categories but one linear predictor. The
  ## checks took J - 1 dimensions, the multinomial's, and so counted this
  ## four-category fit's 40 random intercepts as 120 latent values: a WARN
  ## where it has ten observations to each
  set.seed(4); n <- 400
  d <- data.frame(x = rnorm(n), g = factor(rep(1:40, each = 10)))
  lat <- 0.8 * d$x + rnorm(40, 0, 0.7)[d$g] + stats::rlogis(n)
  d$y <- cut(lat, c(-Inf, -1, 0.5, 2, Inf), labels = letters[1:4],
             ordered_result = TRUE)
  f <- ilm_model(y ~ x + (1 | g), data = d, family = "ordinal",
                 verbose = FALSE)
  ck <- f$checks
  expect_match(ck$detail[ck$check == "re_levels[g]"],
               "40 levels for 1 covariance parameters (40.0 per parameter); us, C = 1",
               fixed = TRUE)
  expect_match(ck$detail[ck$check == "latent_budget"],
               "(400 observations, 40 latent values: g 40)", fixed = TRUE)
  expect_identical(ck$status[ck$check == "latent_budget"], "OK")
})
