bin_fit <- function(n = 600L, seed = 1L) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n), g = factor(sample(c("a", "b"), n, TRUE)),
                  z = rnorm(n))
  d$y <- rbinom(n, 1, plogis(0.8 * d$x + 0.6 * (d$g == "b") + 0.4 * d$z))
  ilm_model(y ~ x + g + z, data = d, family = "binomial", verbose = FALSE)
}

test_that("average marginal effects match marginaleffects", {
  skip_if_not_installed("marginaleffects")
  f <- bin_fit()
  a <- ilm_ame(f)
  expect_equal(nrow(a), 3L)
  expect_setequal(a$term, c("x", "g", "z"))
  expect_equal(a$kind[a$term == "g"], "contrast")
  expect_equal(a$kind[a$term == "x"], "slope")

  ilm_register_marginaleffects()
  m <- as.data.frame(marginaleffects::avg_slopes(f))
  for (v in c("x", "z")) {
    expect_equal(a$estimate[a$term == v], m$estimate[m$term == v],
                 tolerance = 1e-6)
    expect_equal(a$se[a$term == v], m$std.error[m$term == v], tolerance = 1e-4)
  }
  ## and the factor contrast, which is a different calculation on both sides
  expect_equal(a$estimate[a$term == "g"], m$estimate[m$term == "g"],
               tolerance = 1e-6)
})

test_that("a gaussian marginal effect is the coefficient itself", {
  set.seed(2); n <- 400
  d <- data.frame(x = rnorm(n)); d$y <- 0.7 * d$x + rnorm(n)
  f <- ilm_model(y ~ x, data = d, family = "gaussian", verbose = FALSE)
  a <- ilm_ame(f)
  expect_equal(a$estimate, unname(coef(f)["x"]), tolerance = 1e-5)
  expect_equal(a$se, unname(sqrt(diag(vcov(f)))["x"]), tolerance = 1e-4)
})

test_that("the interpretation names each effect on the response scale", {
  f <- bin_fit()
  r <- ilm_interpret(f)
  expect_s3_class(r, "ilm_interpretation")
  txt <- paste(unlist(r$sections), collapse = " ")
  expect_match(txt, "binomial model")
  ## one sentence per coefficient, each naming its own term
  expect_length(r$sections$effects, 3L)
  expect_match(txt, "odds ratio")
  expect_match(txt, "percentage points")
  ## a factor level is described against its reference, not as a bare name
  expect_match(txt, "being b rather than a")
  ## a numeric predictor reads as a change in it
  expect_match(txt, "a higher x")
  ## and the output is deterministic
  expect_identical(capture.output(print(r)), capture.output(print(ilm_interpret(f))))
})

test_that("causal language is withheld unless a design licenses it", {
  f <- bin_fit()
  plain <- paste(unlist(ilm_interpret(f)$sections), collapse = " ")
  expect_match(plain, "is associated with")
  expect_false(grepl("leads to", plain, fixed = TRUE))
  expect_match(plain, "These are associations")

  ## a DAG with a valid adjustment set licenses it
  set.seed(3); n <- 500
  z <- rnorm(n); x <- 0.5 * z + rnorm(n)
  d <- data.frame(x = x, z = z, y = 0.4 * x + 0.6 * z + rnorm(n))
  g <- ilm_dag("dag { x [exposure] ; y [outcome] ; z -> x -> y ; z -> y }")
  m <- ilm_dag_model(g, d, verbose = FALSE)
  causal <- paste(unlist(ilm_interpret(m, ame = FALSE)$sections), collapse = " ")
  expect_match(causal, "leads to")
  expect_match(causal, "minimal sufficient set")
  ## and it says the licence is an assumption, not a finding
  expect_match(causal, "assumption you supplied")

  ## forcing it the other way works too
  expect_match(paste(unlist(ilm_interpret(m, causal = FALSE,
                                          ame = FALSE)$sections), collapse = " "),
               "is associated with")
})

test_that("a non-identified DAG interpretation refuses to report an effect", {
  d <- data.frame(x = rnorm(200), y = rnorm(200))
  g <- ilm_dag("dag { x [exposure] ; y [outcome] ; u [latent]
                      u -> x ; u -> y ; x -> y }")
  m <- ilm_dag_model(g, d, verbose = FALSE)
  r <- ilm_interpret(m)
  txt <- paste(unlist(r$sections), collapse = " ")
  expect_match(txt, "NOT identified")
  expect_null(r$sections$effects)
  expect_false(r$causal)
})

test_that("the interpretation reports what the diagnostics found", {
  f <- bin_fit()
  r <- ilm_interpret(f)
  expect_match(paste(r$sections$diagnostics, collapse = " "),
               "checks passed")
  ## it never claims the model is right, only that the fit is sound
  expect_match(paste(r$sections$diagnostics, collapse = " "),
               "not that the model is right")
})

test_that("did and rdd interpretations name their identifying assumption", {
  set.seed(1)
  d <- expand.grid(unit = 1:40, time = 1:8)
  d$treated <- as.integer(d$unit <= 20); d$post <- as.integer(d$time >= 5)
  d$y <- 1 + 0.3 * d$time + rnorm(40)[d$unit] +
    0.8 * d$treated * d$post + rnorm(nrow(d))
  rd <- ilm_interpret(ilm_did(d, "y", "unit", "time", treated = "treated",
                              post = "post", verbose = FALSE))
  t1 <- paste(unlist(rd$sections), collapse = " ")
  expect_match(t1, "difference-in-differences")
  expect_match(t1, "[Pp]arallel trends")
  expect_match(t1, "would have moved together")

  set.seed(1); n <- 2000
  r <- runif(n, -1, 1)
  f2 <- ilm_rdd(data.frame(r = r, y = 0.5 * r + 0.8 * (r >= 0) + rnorm(n, 0, .5)),
                "y", "r", cutoff = 0, verbose = FALSE)
  t2 <- paste(unlist(ilm_interpret(f2)$sections), collapse = " ")
  expect_match(t2, "regression discontinuity")
  expect_match(t2, "local effect")
  expect_match(t2, "rdrobust")
})

test_that("a wide interval is not reported as no effect", {
  set.seed(9); n <- 200
  d <- data.frame(x = rnorm(n)); d$y <- 0.02 * d$x + rnorm(n)
  f <- ilm_model(y ~ x, data = d, family = "gaussian", verbose = FALSE)
  txt <- paste(unlist(ilm_interpret(f, ame = FALSE)$sections), collapse = " ")
  expect_match(txt, "consistent with no effect")
  expect_match(txt, "not that there is none")
  expect_false(grepl("no effect of", txt, fixed = TRUE))
})

test_that("insight methods let the easystats stack read an ilm_model", {
  skip_if_not_installed("insight")
  expect_true(ilm_register_insight())
  f <- bin_fit()
  expect_equal(insight::find_response(f), "y")
  expect_setequal(insight::find_predictors(f, flatten = TRUE), c("x", "g", "z"))
  expect_equal(insight::n_obs(f), stats::nobs(f))
  mi <- insight::model_info(f)
  expect_true(mi$is_binomial); expect_false(mi$is_linear)
  expect_equal(nrow(insight::get_data(f)), stats::nobs(f))
  ## the parameters insight reads are the ones illume reports
  gp <- insight::get_parameters(f)
  expect_equal(gp$Estimate, unname(coef(f)), tolerance = 1e-10)
  expect_equal(unname(diag(insight::get_varcov(f))),
               unname(diag(vcov(f))), tolerance = 1e-10)
})

test_that("an inferred family, an ordinal fit and the new objects are written up", {
  set.seed(8); n <- 400
  d <- data.frame(x = rnorm(n), g = factor(sample(c("a", "b", "c"), n, TRUE)))
  d$yc <- rpois(n, exp(0.3 + 0.4 * d$x))
  f1 <- suppressMessages(ilm_model(yc ~ x + g, data = d, verbose = FALSE))
  it1 <- ilm_interpret(f1, ame = FALSE)
  expect_match(it1$sections$model, "read off the response", fixed = TRUE)
  d$yo <- cut(0.8 * d$x + stats::rlogis(n), c(-Inf, -0.7, 0.4, 1.3, Inf),
              labels = c("w", "x2", "y2", "z"), ordered_result = TRUE)
  f2 <- ilm_model(yo ~ x, data = d, family = "ordinal", verbose = FALSE)
  it2 <- ilm_interpret(f2, ame = FALSE)
  expect_match(it2$sections$model, "^An ordinal model")
  expect_match(it2$sections$model, "cumulative logit", fixed = TRUE)
  expect_match(it2$sections$effects[1], "higher category", fixed = TRUE)
  expect_match(it2$sections$effects[1], "cumulative odds ratio", fixed = TRUE)
  ## comparisons, in percentage points when they are probabilities
  d$h <- factor(sample(c("u", "v"), n, TRUE))
  fm <- ilm_model(g ~ h + x, data = d, family = "multinomial", verbose = FALSE)
  itc <- ilm_interpret(ilm_contrast(ilm_emmeans(fm, "h", type = "response")))
  expect_length(itc$sections$effects, 3L)
  expect_match(itc$sections$effects[1], "percentage points", fixed = TRUE)
  expect_false(any(grepl(",.", itc$sections$effects, fixed = TRUE)))
  expect_output(print(itc), "INTERPRETATION (comparisons)", fixed = TRUE)
})

test_that("a power analysis is written up with its size and its uncertainty", {
  skip_on_cran()
  p <- ilm_power_design(y ~ arm, design = list(arm = c("control", "treatment")),
                        n_unit = c(40, 80, 160),
                        cells = c(control = 0, treatment = 0.5), sd = 1,
                        term = "arm", sims = 60, verbose = FALSE,
                        progress = FALSE)
  it <- ilm_interpret(p)
  expect_match(it$sections$model, "a t test", fixed = TRUE)
  expect_match(it$sections$model, "fresh draw of the planned design",
               fixed = TRUE)
  expect_match(it$sections$effects[1], "Reaching 80% power takes about",
               fixed = TRUE)
  expect_match(paste(it$sections$caveats, collapse = " "), "Monte Carlo",
               fixed = TRUE)
  expect_output(print(it), "INTERPRETATION (power analysis)", fixed = TRUE)
})
