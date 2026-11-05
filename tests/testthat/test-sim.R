## ilm_sim(): whole simulated datasets from a fit or a scaffold.

quiet_fit <- function(...) suppressWarnings(suppressMessages(ilm_model(...)))

sim_grouped <- function(seed = 11) {
  set.seed(seed)
  d <- data.frame(id = factor(rep(1:12, each = 6)), x = stats::rnorm(72),
                  z = stats::rnorm(72))
  d$y <- 1 + 0.5 * d$x + stats::rnorm(12)[d$id] + stats::rnorm(72)
  d
}

test_that("the rows are the fitted ones, and the response is ilm_simulate()'s", {
  d <- sim_grouped()
  d$x[3] <- NA
  fit <- quiet_fit(y ~ x + (1 | id), data = d, family = "gaussian",
                   verbose = FALSE)
  s <- ilm_sim(fit, seed = 2)
  expect_s3_class(s, "data.frame")
  ## the model's columns only: z is not one of them
  expect_setequal(names(s), c("y", "x", "id"))
  expect_identical(nrow(s), 71L)
  expect_identical(s$x, d$x[-3])
  expect_identical(as.character(s$id), as.character(d$id[-3]))
  expect_equal(s$y, ilm_simulate(fit, 1, seed = 2)[, 1])
  l <- ilm_sim(fit, nsim = 3, seed = 2)
  expect_type(l, "list")
  expect_length(l, 3L)
  expect_equal(unname(sapply(l, `[[`, "y")),
               unname(ilm_simulate(fit, 3, seed = 2)))
  ## and the analysis runs on it as it runs on the data
  f2 <- quiet_fit(y ~ x + (1 | id), data = s, family = "gaussian",
                  verbose = FALSE)
  expect_s3_class(f2, "ilm_model")
})

test_that("the response keeps its type", {
  set.seed(3); n <- 150
  d <- data.frame(x = stats::rnorm(n))
  p <- stats::plogis(-0.2 + 0.8 * d$x)
  d$bf <- factor(ifelse(stats::runif(n) < p, "yes", "no"), levels = c("no", "yes"))
  d$bl <- stats::runif(n) < p
  d$k  <- as.integer(stats::rpois(n, exp(0.3 + 0.4 * d$x)))
  d$m  <- factor(sample(c("a", "b", "c"), n, TRUE))
  d$o  <- factor(sample(c("low", "mid", "high"), n, TRUE),
                 levels = c("low", "mid", "high"), ordered = TRUE)
  sb <- ilm_sim(quiet_fit(bf ~ x, data = d, family = "binomial", verbose = FALSE),
                seed = 1)$bf
  expect_s3_class(sb, "factor")
  expect_identical(levels(sb), c("no", "yes"))
  expect_type(ilm_sim(quiet_fit(bl ~ x, data = d, family = "binomial",
                                verbose = FALSE), seed = 1)$bl, "logical")
  expect_type(ilm_sim(quiet_fit(k ~ x, data = d, family = "poisson",
                                verbose = FALSE), seed = 1)$k, "integer")
  sm <- ilm_sim(quiet_fit(m ~ x, data = d, family = "multinomial",
                          verbose = FALSE), seed = 1)$m
  expect_identical(levels(sm), c("a", "b", "c"))
  so <- ilm_sim(quiet_fit(o ~ x, data = d, family = "ordinal", verbose = FALSE),
                seed = 1)$o
  expect_true(is.ordered(so))
  expect_identical(levels(so), c("low", "mid", "high"))
})

test_that("a transformed term's variable is there, not the computed column", {
  set.seed(4); n <- 100
  d <- data.frame(x = stats::rexp(n) + 0.1, e = stats::runif(n, 1, 3))
  d$k <- stats::rpois(n, d$e * exp(0.2 + 0.3 * log(d$x)))
  fit <- quiet_fit(k ~ log(x) + offset(log(e)), data = d, family = "poisson",
                   verbose = FALSE)
  s <- ilm_sim(fit, seed = 1)
  expect_setequal(names(s), c("k", "x", "e"))
  f2 <- quiet_fit(k ~ log(x) + offset(log(e)), data = s, family = "poisson",
                  verbose = FALSE)
  expect_s3_class(f2, "ilm_model")
})

test_that("a response that is a calculation stops and says why", {
  set.seed(5)
  d <- data.frame(x = stats::rnorm(40)); d$y <- exp(1 + d$x + stats::rnorm(40))
  fit <- quiet_fit(log(y) ~ x, data = d, family = "gaussian", verbose = FALSE)
  expect_error(ilm_sim(fit), "calculation on one")
})

test_that("a censored response comes with what refitting it needs", {
  set.seed(21); n <- 200
  d <- data.frame(x = stats::rnorm(n))
  tt <- exp(1.2 + 0.7 * d$x + 0.6 * log(stats::rexp(n)))
  ct <- stats::rexp(n, rate = 1 / (2 * stats::median(tt)))
  d$time <- pmin(tt, ct); d$event <- as.integer(tt <= ct)
  fit <- quiet_fit(time ~ x, data = d, family = "weibull",
                   censor = ilm_surv(d$time, d$event), verbose = FALSE)
  s <- ilm_sim(fit, seed = 3)
  expect_true(all(s$.event %in% c(0L, 1L)))
  ## a subject censored in the data is censored at the same time, or fails
  ## before it
  cz <- d$event == 0L
  expect_true(all(s$time[cz] <= d$time[cz] + 1e-12))
  expect_true(all(s$.event[cz & s$time >= d$time - 1e-12] == 0L))
  f2 <- quiet_fit(time ~ x, data = s, family = "weibull",
                  censor = ilm_surv(s$time, s$.event), verbose = FALSE)
  expect_s3_class(f2, "ilm_model")

  ## a floor needs no column: the drawn values are held at it
  d$yf <- pmax(0, 0.5 + d$x + stats::rnorm(n))
  ff <- quiet_fit(yf ~ x, data = d, family = "gaussian",
                  censor = ilm_censor(d$yf, lower = 0), verbose = FALSE)
  sf <- ilm_sim(ff, seed = 4)
  expect_false(any(c(".event", ".censored") %in% names(sf)))
  expect_gte(min(sf$yf), 0)
  expect_gt(sum(sf$yf == 0), 0)
})

test_that("a study of another size is the one ilm_power() draws", {
  d <- sim_grouped()
  fit <- quiet_fit(y ~ x + (1 | id), data = d, family = "gaussian",
                   verbose = FALSE)
  s <- ilm_sim(fit, n = 120, seed = 9)
  ## clusters of 6 drawn whole: 20 of them, renumbered by the copy drawn
  expect_identical(nrow(s), 120L)
  expect_s3_class(s$id, "factor")
  expect_identical(levels(s$id), as.character(1:20))
  expect_true(all(table(s$id) == 6L))
  ## the same draw, laid out as data and refitted through the formula, is
  ## ilm_power()'s replicate refitted from its stub
  set.seed(9)
  dr <- ilm_power_study(fit)$draw(120L, ilm_power_beta(fit, list(idx = integer(0)), 1))
  expect_equal(s$y, dr$y)
  r_stub <- ilm_refit_like(dr$stub, y = dr$y, restarts = 1L)
  r_data <- quiet_fit(y ~ x + (1 | id), data = s, family = "gaussian",
                      verbose = FALSE)
  expect_equal(coef(r_data), coef(r_stub), tolerance = 1e-6)
})

test_that("a nested grouping factor is relabelled with its copy, a crossed one is not", {
  set.seed(7)
  d <- expand.grid(item = factor(1:4), class = c("a", "b"), school = 1:8)
  d$y <- stats::rnorm(nrow(d))
  fit <- quiet_fit(y ~ 1 + (1 | school) + (1 | class) + (1 | item), data = d,
                   family = "gaussian", verbose = FALSE)
  s <- ilm_sim(fit, n = 64, seed = 2)
  ## class labels repeat across schools, so class is not nested in school here
  ## and keeps its labels, as an item does
  expect_type(s$school, "integer")
  expect_setequal(unique(s$item), factor(1:4))
  d$class2 <- paste(d$school, d$class)
  fit2 <- quiet_fit(y ~ 1 + (1 | school) + (1 | class2), data = d,
                    family = "gaussian", verbose = FALSE)
  s2 <- ilm_sim(fit2, n = 64, seed = 2)
  ## class2 is nested: each drawn copy of a school gets its own classes
  expect_true(all(startsWith(s2$class2, paste0(s2$school, "."))))
  expect_identical(length(unique(s2$class2)), 2L * length(unique(s2$school)))
})

test_that("a weighted row is drawn as that many rows, each of weight 1", {
  set.seed(8); n <- 60
  d <- data.frame(x = stats::rnorm(n), w = sample(1:3, n, TRUE))
  d$k <- stats::rpois(n, exp(0.2 + 0.3 * d$x))
  ## weights are found in `data` by name, so the call is made directly rather
  ## than through quiet_fit()'s dots
  fit <- suppressWarnings(suppressMessages(
    ilm_model(k ~ x, data = d, family = "poisson", weights = w, verbose = FALSE)))
  expect_identical(ilm_sim(fit, seed = 1)$w, d$w)
  s <- ilm_sim(fit, n = 50, seed = 1)
  expect_identical(nrow(s), 50L)
  expect_true(all(s$w == 1))
})

test_that("a scaffold gives its design, at its own size or another", {
  sc <- suppressWarnings(ilm_scaffold(
    y ~ arm * time + (1 | id),
    design = list(arm = c("control", "treatment"), time = c("pre", "post")),
    within = "time", n_unit = 40,
    cells = c(control.pre = 12, control.post = 12.2,
              treatment.pre = 12, treatment.post = 14.5),
    sd = 4, icc = 0.5, verbose = FALSE))
  s <- ilm_sim(sc, seed = 1)
  expect_identical(nrow(s), 80L)
  expect_setequal(names(s), c("y", "arm", "time", "id"))
  expect_equal(s$y, ilm_simulate(sc, 1, seed = 1)[, 1])
  s2 <- ilm_sim(sc, n = 160, seed = 1)
  expect_identical(nrow(s2), 160L)
  expect_identical(length(unique(s2$id)), 80L)
  expect_true(all(table(s2$arm) == 80L))
  expect_setequal(names(s2), names(s))
})

test_that("effect sets the coefficient the data are drawn at", {
  set.seed(9); n <- 80
  d <- data.frame(x = stats::rnorm(n), g = factor(sample(c("a", "b", "c"), n, TRUE)))
  d$y <- 1 + 0.5 * d$x + c(a = 0, b = 1, c = 2)[as.character(d$g)] + stats::rnorm(n)
  fit <- quiet_fit(y ~ x + g, data = d, family = "gaussian", verbose = FALSE)
  f1 <- fit; f1$beta[match("x", names(coef(fit)))] <- 2
  expect_equal(ilm_sim(fit, effect = 2, term = "x", seed = 3)$y,
               ilm_simulate(f1, 1, seed = 3)[, 1])
  ## a term of several coefficients is multiplied
  f2 <- fit; k <- grep("^g", names(coef(fit)))
  f2$beta[k] <- 0.5 * coef(fit)[k]
  expect_equal(ilm_sim(fit, effect = 0.5, term = "g", seed = 3)$y,
               ilm_simulate(f2, 1, seed = 3)[, 1])
})

test_that("a seed leaves the user's random stream as it was", {
  d <- sim_grouped()
  fit <- quiet_fit(y ~ x + (1 | id), data = d, family = "gaussian",
                   verbose = FALSE)
  set.seed(5); a <- stats::runif(1)
  set.seed(5); invisible(ilm_sim(fit, seed = 1)); b <- stats::runif(1)
  expect_identical(a, b)
  set.seed(5); invisible(ilm_sim(fit, n = 30, seed = 1)); b <- stats::runif(1)
  expect_identical(a, b)
})

test_that("arguments it cannot use stop with a reason", {
  d <- sim_grouped()
  fit <- quiet_fit(y ~ x + (1 | id), data = d, family = "gaussian",
                   verbose = FALSE)
  expect_error(ilm_sim(fit, nsim = 0), "whole number")
  expect_error(ilm_sim(fit, n = c(50, 100)), "one number")
  expect_error(ilm_sim(fit, effect = c(0.1, 0.2)), "one number")
  expect_error(ilm_sim(fit, term = "x"), "give `effect` too")
  expect_error(ilm_sim(fit, n = 5), "not a study")
  expect_error(ilm_sim(data.frame(y = 1)), "fitted ilm_model")
})
