## Two codings of one meaning (Craig's item 226): `censored` is the direction
## of censoring, `event` Surv()'s event indicator; and Surv() on the left of a
## formula is read as right-censored follow-up, fitted as ilm_surv() fits it.

cc_data <- function(n = 200, seed = 11) {
  set.seed(seed)
  d <- data.frame(x = stats::rnorm(n), g = factor(sample(letters[1:3], n, TRUE)))
  tt <- stats::rweibull(n, shape = 1.5, scale = exp(1 + 0.4 * d$x))
  cens <- stats::runif(n, 0.5, 6)
  d$time <- pmin(tt, cens)
  d$event <- as.integer(tt <= cens)
  d
}

test_that("censored and event code the same fit oppositely", {
  d <- cc_data()
  a <- suppressMessages(ilm_censor(d$time, censored = 1 - d$event))
  b <- ilm_censor(d$time, event = d$event)
  expect_identical(as.integer(a), as.integer(b))
  fa <- suppressMessages(ilm_model(time ~ x, data = d, family = "weibull",
                                   censor = a, verbose = FALSE))
  fb <- suppressMessages(ilm_model(time ~ x, data = d, family = "weibull",
                                   censor = b, verbose = FALSE))
  expect_equal(stats::coef(fa), stats::coef(fb), tolerance = 1e-12)
  expect_equal(as.numeric(logLik(fa)), as.numeric(logLik(fb)), tolerance = 1e-12)
  ## TRUE / FALSE is an event indicator too
  expect_identical(as.integer(ilm_censor(d$time, event = d$event == 1)), as.integer(b))
})

test_that("giving both codings is an error that names them", {
  expect_error(ilm_censor(1:3, censored = c(0, 1, 0), event = c(1, 0, 1)),
               "`censored` is the direction of censoring (-1 left, 0 observed, 1 right), `event` is Surv()'s event indicator",
               fixed = TRUE)
})

test_that("a censored vector of only 0s and 1s is said once a session", {
  rm(list = ls(ilm_censor_said), envir = ilm_censor_said)
  expect_message(ilm_censor(1:4, censored = c(0, 1, 0, 1)),
                 "censored codes 1 as right-censored (not observed); for an event indicator (1 = event observed, as in Surv) use event = or ilm_surv()",
                 fixed = TRUE)
  expect_silent(ilm_censor(1:4, censored = c(0, 1, 0, 1)))
  ## a -1 says it is a direction, so nothing is said of it
  rm(list = ls(ilm_censor_said), envir = ilm_censor_said)
  expect_silent(ilm_censor(1:4, censored = c(0, 1, -1, 1)))
  expect_silent(ilm_censor(1:4, event = c(0, 1, 0, 1)))
})

test_that("Surv() on the left is fitted exactly as ilm_surv()", {
  d <- cc_data()
  f1 <- suppressMessages(ilm_model(time ~ x + g, data = d, family = "weibull",
                                   censor = ilm_surv(d$time, d$event), verbose = FALSE))
  f2 <- suppressMessages(ilm_model(Surv(time, event) ~ x + g, data = d,
                                   family = "weibull", verbose = FALSE))
  f3 <- suppressMessages(ilm_model(survival::Surv(time, event == 1) ~ x + g, data = d,
                                   family = "weibull", verbose = FALSE))
  for (f in list(f2, f3)) {
    expect_equal(stats::coef(f), stats::coef(f1), tolerance = 1e-10)
    expect_equal(as.numeric(logLik(f)), as.numeric(logLik(f1)), tolerance = 1e-10)
    expect_equal(sqrt(diag(stats::vcov(f))), sqrt(diag(stats::vcov(f1))), tolerance = 1e-10)
    expect_identical(as.integer(f$censor), as.integer(f1$censor))
  }
  ## Surv()'s 1/2 coding, 2 = the event
  d$st <- d$event + 1L
  f4 <- suppressMessages(ilm_model(Surv(time, st) ~ x + g, data = d,
                                   family = "weibull", verbose = FALSE))
  expect_equal(stats::coef(f4), stats::coef(f1), tolerance = 1e-10)
  ## a row missing its event is dropped from the response with it
  d$event[3] <- NA
  f5 <- suppressMessages(ilm_model(Surv(time, event) ~ x, data = d,
                                   family = "weibull", verbose = FALSE))
  expect_identical(stats::nobs(f5), nrow(d) - 1L)
})

test_that("Surv() of any other kind stops and names the remedy", {
  d <- cc_data()
  d$t0 <- 0
  expect_error(ilm_model(Surv(time, event, type = "left") ~ x, data = d,
                         family = "weibull", verbose = FALSE),
               "ilm_censor(y, censored = )", fixed = TRUE)
  expect_error(ilm_model(Surv(t0, time, event) ~ x, data = d,
                         family = "weibull", verbose = FALSE),
               "survival::coxph()", fixed = TRUE)
  expect_error(ilm_model(Surv(time, time + 1, type = "interval2") ~ x, data = d,
                         family = "weibull", verbose = FALSE),
               "interval-censored", fixed = TRUE)
  expect_error(ilm_model(Surv(time, event) ~ x, data = d, family = "weibull",
                         censor = ilm_surv(d$time, d$event), verbose = FALSE),
               "not both", fixed = TRUE)
})

test_that("summary() says what was censored", {
  d <- cc_data()
  n_ev <- sum(d$event); n_c <- nrow(d) - n_ev
  f <- suppressMessages(ilm_model(Surv(time, event) ~ x, data = d,
                                  family = "weibull", verbose = FALSE))
  out <- capture.output(summary(f))
  expect_true(any(grepl(sprintf("Surv() read as right-censored: %d events, %d censored",
                                n_ev, n_c), out, fixed = TRUE)))
  expect_true(any(grepl(sprintf("Censoring: %d observed, %d right-censored", n_ev, n_c),
                        out, fixed = TRUE)))
  ## a floor and a ceiling: left-censored counted too
  set.seed(3); y <- stats::rnorm(150)
  dd <- data.frame(y = pmin(pmax(y, -1), 1.2), x = stats::rnorm(150))
  fc <- suppressMessages(ilm_model(y ~ x, data = dd, family = "gaussian",
                                   censor = ilm_censor(dd$y, lower = -1, upper = 1.2),
                                   verbose = FALSE))
  oc <- capture.output(summary(fc))
  expect_true(any(grepl(sprintf("Censoring: %d observed, %d right-censored, %d left-censored",
                                sum(dd$y > -1 & dd$y < 1.2), sum(dd$y >= 1.2),
                                sum(dd$y <= -1)), oc, fixed = TRUE)))
  expect_false(any(grepl("Surv() read", oc, fixed = TRUE)))
})
