## family = "auto": the family is read off the response, said, and written into
## the call. What the response cannot settle is asked about, not guessed.

sim_auto <- function(n = 240, seed = 11) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n))
  d$yg  <- 1 + 0.5 * d$x + rnorm(n)
  d$yb  <- rbinom(n, 1, plogis(0.4 * d$x))
  d$yl  <- d$yb == 1
  d$yf  <- factor(ifelse(d$yb == 1, "yes", "no"))
  d$yc  <- rpois(n, exp(0.3 + 0.4 * d$x))
  d$ym  <- factor(sample(c("a", "b", "c"), n, TRUE))
  d$yo  <- factor(sample(c("lo", "mid", "hi"), n, TRUE),
                  levels = c("lo", "mid", "hi"), ordered = TRUE)
  d$yp  <- plogis(0.3 * d$x + rnorm(n, 0, 0.5))
  d$ybp <- round(120 + 5 * d$x + rnorm(n, 0, 10))   # a blood pressure
  d
}

auto_msgs <- function(expr) {
  msgs <- character()
  val <- withCallingHandlers(expr, message = function(m) {
    msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage")
  })
  list(value = val, msgs = msgs)
}

test_that("each kind of response gets the family it needs, and says so", {
  d <- sim_auto()
  cases <- list(yg = "gaussian", yb = "binomial", yl = "binomial",
                yf = "binomial", yc = "poisson", ym = "multinomial",
                yo = "ordinal", yp = "beta", ybp = "gaussian")
  for (v in names(cases)) {
    r <- auto_msgs(ilm_model(stats::reformulate("x", v), data = d,
                             verbose = FALSE))
    expect_identical(r$value$family$name, cases[[v]], info = v)
    ## said once, naming the family and the variable
    hit <- grepl("inferred from", r$msgs, fixed = TRUE)
    expect_equal(sum(hit), 1L, info = v)
    expect_match(r$msgs[hit], sprintf("family = \"%s\"", cases[[v]]),
                 fixed = TRUE, info = v)
    expect_match(r$msgs[hit], sprintf("`%s`", v), fixed = TRUE, info = v)
    ## and recorded, so a refit through the call does not guess again
    expect_identical(r$value$call$family, cases[[v]], info = v)
    expect_type(r$value$family_inferred, "character")
  }
})

test_that("NULL means auto, and a family that is given is not second-guessed", {
  d <- sim_auto()
  r <- auto_msgs(ilm_model(yc ~ x, data = d, family = NULL, verbose = FALSE))
  expect_identical(r$value$family$name, "poisson")
  expect_true(any(grepl("inferred from", r$msgs, fixed = TRUE)))
  r2 <- auto_msgs(ilm_model(yc ~ x, data = d, family = "nbinom",
                            verbose = FALSE))
  expect_identical(r2$value$family$name, "nbinom")
  expect_false(any(grepl("inferred from", r2$msgs, fixed = TRUE)))
  expect_null(r2$value$family_inferred)
  ## the refit through the stored call is silent about the family
  f <- suppressMessages(ilm_model(yg ~ x, data = d, verbose = FALSE))
  r3 <- auto_msgs(eval(f$call))
  expect_false(any(grepl("inferred from", r3$msgs, fixed = TRUE)))
  expect_equal(coef(r3$value), coef(f))
})

test_that("the hint names the usual competitor", {
  d <- sim_auto()
  r <- auto_msgs(ilm_model(yc ~ x, data = d, verbose = FALSE))
  expect_match(paste(r$msgs, collapse = " "), "nbinom", fixed = TRUE)
  r <- auto_msgs(ilm_model(ym ~ x, data = d, verbose = FALSE))
  expect_match(paste(r$msgs, collapse = " "), "ordered factor", fixed = TRUE)
  r <- auto_msgs(ilm_model(ybp ~ x, data = d, verbose = FALSE))
  expect_match(paste(r$msgs, collapse = " "), "never come near zero",
               fixed = TRUE)
})

test_that("what the response cannot settle is asked about", {
  d <- sim_auto()
  d$y12 <- d$yb + 1
  expect_error(ilm_model(y12 ~ x, data = d, verbose = FALSE),
               "cannot tell which family .* two values, 1 and 2")
  d$yp01 <- pmin(1, pmax(0, d$yp + rnorm(nrow(d), 0, 0.3)))
  expect_error(ilm_model(yp01 ~ x, data = d, verbose = FALSE),
               "values at exactly 0 or 1")
  d$w <- 20
  expect_error(ilm_model(yp ~ x, data = d, weights = w, verbose = FALSE),
               "numbers of trials")
  d$when <- as.Date("2024-01-01") + seq_len(nrow(d))
  expect_error(ilm_model(when ~ x, data = d, verbose = FALSE),
               "a date, a time or a duration")
  d$k <- 3
  expect_error(ilm_model(k ~ x, data = d, verbose = FALSE),
               "single value")
  ## a survival time: which survival family is a modelling decision
  set.seed(2)
  s <- data.frame(x = rnorm(120))
  s$time <- rweibull(120, 1.5, exp(1 + 0.3 * s$x))
  s$event <- rbinom(120, 1, 0.8)
  expect_error(ilm_model(time ~ x, data = s,
                         censor = ilm_surv(s$time, s$event), verbose = FALSE),
               "family = \"rp\"")
})

test_that("a zero part and a floor change what is inferred", {
  d <- sim_auto()
  d$ypz <- ifelse(runif(nrow(d)) < 0.2, 0, d$yp)
  r <- auto_msgs(ilm_model(ypz ~ x, data = d, ziformula = ~ 1,
                           zi_type = "hurdle", verbose = FALSE))
  expect_identical(r$value$family$name, "beta")
  ## whole numbers above 1 are a count when the model has a zero part
  d$ycz <- ifelse(runif(nrow(d)) < 0.3, 0, d$yc + 2)
  r <- auto_msgs(ilm_model(ycz ~ x, data = d, ziformula = ~ 1,
                           verbose = FALSE))
  expect_identical(r$value$family$name, "poisson")
  ## a response censored at a floor is gaussian, in its Tobit form
  d$yfl <- pmax(d$yg, 0.5)
  r <- auto_msgs(ilm_model(yfl ~ x, data = d,
                           censor = ilm_censor(d$yfl, lower = 0.5),
                           verbose = FALSE))
  expect_identical(r$value$family$name, "gaussian")
  expect_match(paste(r$msgs, collapse = " "), "Tobit", fixed = TRUE)
})

test_that("the printed fit and its summary say the family was inferred", {
  d <- sim_auto()
  f <- suppressMessages(ilm_model(ym ~ x, data = d, verbose = FALSE))
  expect_output(print(f), "multinomial, inferred", fixed = TRUE)
  expect_output(print(summary(f)), "inferred from the response", fixed = TRUE)
  g <- ilm_model(ym ~ x, data = d, family = "multinomial", verbose = FALSE)
  expect_false(any(grepl("inferred", capture.output(print(summary(g))))))
})

test_that("multiple imputation infers the family once and holds it", {
  set.seed(4); n <- 150
  d <- data.frame(x = rnorm(n), z = rnorm(n))
  d$y <- rpois(n, exp(0.2 + 0.3 * d$x))
  d$x[sample(n, 25)] <- NA
  imp <- ilm_impute(d, m = 3, seed = 1, verbose = FALSE)
  r <- auto_msgs(ilm_mi_pool(imp, y ~ x + z))
  expect_equal(sum(grepl("inferred from", r$msgs, fixed = TRUE)), 1L)
  expect_true(all(vapply(attr(r$value, "fits"),
                         function(f) f$family$name, "") == "poisson"))
})

test_that("the narrated workflows use the same rules", {
  d <- sim_auto()
  g <- ilm_dag("dag { x [exposure] ; yo [outcome] ; x -> yo }")
  m <- ilm_dag_model(g, d[c("x", "yo")], verbose = FALSE)
  expect_identical(m$family, "ordinal")
  expect_true(any(grepl("family inferred as ordinal", m$steps)))
})
