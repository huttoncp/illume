emm_data <- function(n = 300L, seed = 1L, unbalanced = FALSE) {
  set.seed(seed)
  p <- if (unbalanced) c(.6, .25, .15) else c(1, 1, 1) / 3
  d <- data.frame(g = factor(sample(c("a", "b", "c"), n, TRUE, prob = p)),
                  h = factor(sample(c("lo", "hi"), n, TRUE)), x = rnorm(n))
  d$y <- 1 + 0.5 * (d$g == "b") - 0.3 * (d$g == "c") +
    0.4 * (d$h == "hi") + 0.2 * d$x + rnorm(n)
  d
}
emm_fit <- function(d)
  ilm_model(y ~ g + h + x, data = d, family = "gaussian", verbose = FALSE)

## a design where g and h are associated, so cell weights and marginal
## weights have to disagree
emm_assoc <- function(n = 400L, seed = 7L) {
  set.seed(seed)
  h <- factor(sample(c("lo", "hi"), n, TRUE))
  g <- factor(ifelse(h == "hi", sample(c("a", "b", "c"), n, TRUE, c(.7, .2, .1)),
                                sample(c("a", "b", "c"), n, TRUE, c(.1, .2, .7))))
  d <- data.frame(g = g, h = h, x = rnorm(n))
  d$y <- 1 + 0.9 * (d$g == "b") - 1.4 * (d$g == "c") +
    0.4 * (d$h == "hi") + 0.2 * d$x + rnorm(n)
  d
}

test_that("marginal means average over what they should", {
  d <- emm_data(); f <- emm_fit(d)
  e <- ilm_emmeans(f, "g")
  expect_s3_class(e, "ilm_emm")
  expect_equal(nrow(e), 3L)
  expect_setequal(e$g, c("a", "b", "c"))
  expect_true(all(e$se > 0))
  expect_true(all(e$lower < e$estimate & e$estimate < e$upper))
  ## the differences between them are the model's group differences, since the
  ## other predictors are held at common values
  expect_equal(e$estimate[e$g == "b"] - e$estimate[e$g == "a"],
               unname(coef(f)["gb"]), tolerance = 1e-6)
  ## two-way grids give a row per combination
  expect_equal(nrow(ilm_emmeans(f, c("g", "h"))), 6L)
})

test_that("equal and proportional weights differ in level but not in contrast", {
  d <- emm_data(unbalanced = TRUE); f <- emm_fit(d)
  eq <- ilm_emmeans(f, "h", weights = "equal")
  pr <- ilm_emmeans(f, "h", weights = "proportional")
  ## averaging over an unbalanced g moves the level of each mean
  expect_false(isTRUE(all.equal(eq$estimate, pr$estimate)))
  ## but NOT the difference between them: both give every level of h the same
  ## mix of g, so the h difference is the h coefficient either way
  expect_equal(diff(eq$estimate), diff(pr$estimate), tolerance = 1e-10)
  expect_equal(abs(diff(eq$estimate)), abs(unname(coef(f)["hlo"])),
               tolerance = 1e-6)
  expect_equal(attr(eq, "weights"), "equal")
  expect_equal(attr(pr, "weights"), "proportional")
  ## with a balanced design the levels agree too
  b <- emm_fit(emm_data(unbalanced = FALSE))
  expect_equal(ilm_emmeans(b, "h", weights = "equal")$estimate,
               ilm_emmeans(b, "h", weights = "proportional")$estimate,
               tolerance = 0.02)
})

test_that("cell weights give each group its own mix, and say so", {
  d <- emm_assoc(); f <- emm_fit(d)
  eq <- ilm_emmeans(f, "h", weights = "equal")
  ce <- ilm_emmeans(f, "h", weights = "cells")
  ## the h difference under equal weights is the coefficient; under cell
  ## weights it is that plus the difference in g composition, which here is
  ## large enough to swamp it
  expect_equal(abs(diff(eq$estimate)), abs(unname(coef(f)["hlo"])),
               tolerance = 1e-6)
  expect_gt(abs(diff(ce$estimate) - diff(eq$estimate)), 0.5)
  ## and contrasting them warns rather than passing it off as adjusted
  expect_warning(ilm_contrast(ce), "own mix")
  expect_silent(ilm_contrast(eq))
})

test_that("contrasts are the differences, with a joint adjustment", {
  d <- emm_data(); f <- emm_fit(d)
  e <- ilm_emmeans(f, "g")
  cc <- ilm_contrast(e)
  expect_s3_class(cc, "ilm_contrast")
  expect_equal(nrow(cc), 3L)                    # 3 groups -> 3 pairs
  ## every contrast is a difference of two marginal means
  bm <- e$estimate[e$g == "b"] - e$estimate[e$g == "a"]
  expect_equal(cc$estimate[cc$contrast == "b - a"], bm, tolerance = 1e-8)
  ## adjustment widens, and max_t sits between none and bonferroni
  w <- function(a) mean(ilm_contrast(e, adjust = a)$upper -
                        ilm_contrast(e, adjust = a)$lower)
  expect_lt(w("none"), w("max_t"))
  expect_lt(w("max_t"), w("bonferroni"))
  ## a single contrast has nothing to adjust
  one <- ilm_contrast(ilm_emmeans(f, "h"))
  expect_equal(nrow(one), 1L)
  expect_equal(one$adjust, "none")
})

test_that("the contrast schemes do what they say", {
  d <- emm_data(); f <- emm_fit(d)
  e <- ilm_emmeans(f, "g")
  tv <- ilm_contrast(e, method = "trt.vs.ctrl")
  expect_equal(nrow(tv), 2L)
  expect_true(all(grepl("- a$", tv$contrast)))
  tv2 <- ilm_contrast(e, method = "trt.vs.ctrl", ref = "c")
  expect_true(all(grepl("- c$", tv2$contrast)))
  expect_error(ilm_contrast(e, method = "trt.vs.ctrl", ref = "zz"),
               "not one of")
  po <- ilm_contrast(e, method = "poly")
  expect_true(any(grepl("^\\.L$|linear", po$contrast)))
})

test_that("emmeans and contrasts refuse what they cannot do", {
  d <- emm_data(); f <- emm_fit(d)
  expect_error(ilm_emmeans(f, "nope"), "not in the model")
  expect_error(ilm_emmeans(1, "g"), "must be a fitted ilm_model")
  expect_error(ilm_contrast(mtcars), "must be an ilm_emmeans")
  ## a smooth has no marginal mean in this sense, and says so
  skip_if_not_installed("mgcv")
  set.seed(2); dd <- data.frame(z = runif(200, -3, 3),
                                g = factor(sample(c("a", "b"), 200, TRUE)))
  dd$y <- sin(dd$z) + rnorm(200)
  fs <- ilm_model(y ~ g + t2(z), data = dd, family = "gaussian",
                  verbose = FALSE)
  expect_error(ilm_emmeans(fs, "g"), "not a plain variable")
})

test_that("all three weightings agree with emmeans", {
  skip_if_not_installed("emmeans")
  d <- emm_assoc(); f <- emm_fit(d)
  l <- stats::lm(y ~ g + h + x, data = d)
  for (w in c("equal", "proportional", "cells")) {
    mine <- ilm_emmeans(f, "h", weights = w)
    theirs <- as.data.frame(emmeans::emmeans(l, ~ h, weights = w))
    expect_equal(mine$estimate, theirs$emmean, tolerance = 1e-4,
                 info = paste("weights =", w))
    expect_equal(mine$se, theirs$SE, tolerance = 1e-4,
                 info = paste("weights =", w))
  }
})

test_that("contrasts agree with emmeans", {
  skip_if_not_installed("emmeans")
  d <- emm_data(n = 300L); f <- emm_fit(d)
  l <- stats::lm(y ~ g + h + x, data = d)
  eg <- emmeans::emmeans(l, ~ g)

  mine <- ilm_emmeans(f, "g")
  expect_equal(mine$estimate, as.data.frame(eg)$emmean, tolerance = 1e-4)
  expect_equal(mine$se, as.data.frame(eg)$SE, tolerance = 1e-5)

  ## unadjusted p-values are the same quantity computed the same way
  m0 <- ilm_contrast(mine, adjust = "none")
  t0 <- as.data.frame(emmeans::contrast(eg, "pairwise", adjust = "none"))
  expect_equal(abs(m0$estimate), abs(t0$estimate), tolerance = 1e-4)
  expect_equal(m0$se, t0$SE, tolerance = 1e-5)
  expect_equal(m0$p_adj, t0$p.value, tolerance = 1e-4)

  ## bonferroni likewise
  mb <- ilm_contrast(mine, adjust = "bonferroni")
  tb <- as.data.frame(emmeans::contrast(eg, "pairwise", adjust = "bonferroni"))
  expect_equal(mb$p_adj, tb$p.value, tolerance = 1e-4)

  ## max_t is a DIFFERENT adjustment from Tukey -- the maximum of a
  ## multivariate t against the studentized range -- so they agree closely
  ## without being identical
  mt <- ilm_contrast(mine, adjust = "max_t")
  tt <- as.data.frame(emmeans::contrast(eg, "pairwise", adjust = "tukey"))
  expect_lt(max(abs(mt$p_adj - tt$p.value)), 0.05)
})

test_that("a link-scale model reports on the link scale unless asked", {
  set.seed(3); n <- 600
  d <- data.frame(g = factor(sample(c("a", "b"), n, TRUE)))
  d$y <- rbinom(n, 1, plogis(-0.5 + 1.2 * (d$g == "b")))
  f <- ilm_model(y ~ g, data = d, family = "binomial", verbose = FALSE)
  lk <- ilm_emmeans(f, "g")
  rs <- ilm_emmeans(f, "g", type = "response")
  expect_equal(attr(lk, "type"), "link")
  ## the response scale is the inverse logit of the link scale
  expect_equal(rs$estimate, plogis(lk$estimate), tolerance = 1e-8)
  expect_true(all(rs$estimate > 0 & rs$estimate < 1))
  ## and the difference on the link scale is the coefficient
  expect_equal(diff(lk$estimate), unname(coef(f)["gb"]), tolerance = 1e-6)
})

test_that("a multinomial fit gets a mean for every category, matching emmeans", {
  skip_if_not_installed("emmeans")
  skip_if_not_installed("nnet")
  set.seed(3); n <- 600
  d <- data.frame(g = factor(sample(c("a", "b", "c"), n, TRUE, prob = c(.5, .3, .2))),
                  h = factor(sample(c("u", "v"), n, TRUE)), x = rnorm(n))
  eta <- cbind(0, 0.3 + 0.6 * d$x + 0.5 * (d$g == "b") - 0.3 * (d$h == "v"),
               -0.2 - 0.4 * d$x + 0.8 * (d$g == "c"))
  P <- exp(eta) / rowSums(exp(eta))
  d$y <- factor(apply(P, 1, function(p) sample(c("lo", "mid", "hi"), 1, prob = p)),
                levels = c("lo", "mid", "hi"))
  f <- ilm_model(y ~ g + h + x, data = d, family = "multinomial", verbose = FALSE)
  m <- nnet::multinom(y ~ g + h + x, data = d, trace = FALSE, reltol = 1e-12,
                      maxit = 1000)
  ## the coefficients differ, because nnet codes against a baseline category;
  ## the probabilities and the centred log-odds do not
  for (w in c("equal", "proportional", "cells")) for (ty in c("response", "link")) {
    mine <- ilm_emmeans(f, "g", type = ty, weights = w)
    th <- as.data.frame(emmeans::emmeans(m, ~ g | y, weights = w,
                                         mode = if (ty == "response") "prob" else "latent"))
    th <- th[order(match(th$g, levels(d$g)), match(th$y, levels(d$y))), ]
    expect_equal(mine$estimate, th[[if (ty == "response") "prob" else "emmean"]],
                 tolerance = 1e-4, info = paste(w, ty))
    expect_equal(mine$se, th$SE, tolerance = 1e-4, info = paste(w, ty))
  }
  e <- ilm_emmeans(f, "g", type = "response")
  expect_identical(levels(e$category), c("lo", "mid", "hi"))
  expect_equal(as.numeric(tapply(e$estimate, e$g, sum)), c(1, 1, 1), tolerance = 1e-10)
  expect_true(all(e$lower > 0 & e$upper < 1))
  ## comparisons are made within a category, as differences in probability
  cc <- ilm_contrast(e, method = "trt.vs.ctrl", adjust = "none")
  ct <- as.data.frame(emmeans::contrast(emmeans::emmeans(m, ~ g | y, mode = "prob"),
                                        "trt.vs.ctrl", adjust = "none"))
  expect_identical(cc$contrast, paste0(ct$y, ": ", ct$contrast))
  expect_equal(cc$estimate, ct$estimate, tolerance = 1e-4)
  expect_equal(cc$se, ct$SE, tolerance = 1e-4)
  expect_output(print(e), "each group's sum")
})

test_that("an ordinal fit gets category probabilities, matching emmeans on polr", {
  skip_if_not_installed("emmeans")
  skip_if_not_installed("MASS")
  set.seed(8); n <- 500
  d <- data.frame(g = factor(sample(c("a", "b", "c"), n, TRUE, prob = c(.5, .3, .2))),
                  h = factor(sample(c("u", "v"), n, TRUE)), x = rnorm(n))
  d$yo <- cut(0.6 * d$x + 0.5 * (d$g == "b") - 0.4 * (d$g == "c") +
                0.3 * (d$h == "v") + stats::rlogis(n),
              c(-Inf, -0.7, 0.4, 1.3, Inf), labels = c("w", "x2", "y2", "z"),
              ordered_result = TRUE)
  f <- ilm_model(yo ~ g + h + x, data = d, family = "ordinal", verbose = FALSE)
  m <- MASS::polr(yo ~ g + h + x, data = d, Hess = TRUE, method = "logistic")
  for (w in c("equal", "cells")) {
    mine <- ilm_emmeans(f, "g", type = "response", weights = w)
    th <- as.data.frame(emmeans::emmeans(m, ~ g | yo, mode = "prob", weights = w))
    th <- th[order(match(th$g, levels(d$g)), match(th$yo, levels(d$yo))), ]
    ## the standard errors carry the thresholds' uncertainty, as polr's do
    expect_equal(mine$estimate, th$prob, tolerance = 1e-4, info = w)
    expect_equal(mine$se, th$SE, tolerance = 1e-4, info = w)
  }
  e <- ilm_emmeans(f, "g", type = "response")
  expect_equal(as.numeric(tapply(e$estimate, e$g, sum)), c(1, 1, 1), tolerance = 1e-10)
  cc <- ilm_contrast(e, method = "trt.vs.ctrl", adjust = "none")
  ct <- as.data.frame(emmeans::contrast(emmeans::emmeans(m, ~ g | yo, mode = "prob"),
                                        "trt.vs.ctrl", adjust = "none"))
  expect_equal(cc$estimate, ct$estimate, tolerance = 1e-4)
  expect_equal(cc$se, ct$SE, tolerance = 1e-4)
})

test_that("a multinomial fit gets a slope for every category, matching emtrends", {
  skip_if_not_installed("emmeans")
  skip_if_not_installed("nnet")
  set.seed(3); n <- 600
  d <- data.frame(g = factor(sample(c("a", "b"), n, TRUE)), x = rnorm(n))
  eta <- cbind(0, 0.3 + 0.6 * d$x + 0.5 * (d$g == "b") * d$x, -0.2 - 0.4 * d$x)
  P <- exp(eta) / rowSums(exp(eta))
  d$y <- factor(apply(P, 1, function(p) sample(c("lo", "mid", "hi"), 1, prob = p)),
                levels = c("lo", "mid", "hi"))
  f <- ilm_model(y ~ g * x, data = d, family = "multinomial", verbose = FALSE)
  m <- nnet::multinom(y ~ g * x, data = d, trace = FALSE, reltol = 1e-12,
                      maxit = 1000)
  tr <- ilm_trends(f, "g", "x")
  et <- as.data.frame(emmeans::emtrends(m, ~ g | y, var = "x", mode = "latent"))
  et <- et[order(match(et$g, levels(d$g)), match(et$y, levels(d$y))), ]
  expect_equal(tr$estimate, et$x.trend, tolerance = 1e-4)
  expect_equal(tr$se, et$SE, tolerance = 1e-4)
  ## compared within a category, never across
  cc <- ilm_contrast(tr)
  expect_identical(cc$contrast, c("lo: b - a", "mid: b - a", "hi: b - a"))
  expect_output(print(tr), "CENTRED log-odds")
})
