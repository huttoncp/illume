## The optimizer check reads nlminb's stopping code beside the fit's own
## restarts from the solution: restarts that come back to the same optimum
## answer the code. And a restart that meets a non-finite gradient is set
## aside rather than stopping the fit.

q <- function(e) suppressMessages(suppressWarnings(e))

test_that("restarts that reach the same optimum answer a false-convergence code", {
  ## a negative binomial panel with a random walk and two random intercepts,
  ## where nlminb ends on "false convergence (8)" at an optimum its restarts
  ## reproduce (an external report showed the same in one row order of a
  ## data set and not the other). The code can differ by platform, so what is
  ## required holds either way: the line is OK, and says why when the code is
  ## not zero.
  set.seed(253)
  G <- 12; Tn <- 24
  d <- expand.grid(t = seq_len(Tn), g = factor(seq_len(G)))
  d$season <- factor((d$t - 1) %% 4 + 1)
  d$gs <- interaction(d$g, d$season, drop = TRUE)
  lat <- unlist(lapply(seq_len(G), function(i) cumsum(rnorm(Tn, 0, 0.05))))
  mu <- exp(4 + rnorm(G, 0, 0.8)[d$g] + rnorm(nlevels(d$gs), 0, 0.1)[d$gs] +
            c(0, 0.2, -0.1, 0.3)[d$season] + lat)
  d$y <- rnbinom(nrow(d), mu = mu, size = 30)
  f <- q(ilm_model(y ~ season + t + (1 | g) + (1 | gs), data = d, family = "nbinom",
                   ar = ilm_rw1(~ t | g), verbose = FALSE))
  r <- f$checks[f$checks$check == "optimizer", ]
  expect_identical(r$status, "OK")
  if (f$opt$convergence != 0L)
    expect_match(r$detail, "restarts from the solution reached the same optimum",
                 fixed = TRUE)
})

test_that("a restart that meets a non-finite gradient does not stop the fit", {
  ## twelve points of one AR(1) series whose correlation runs to -1, where a
  ## restart from the first solution evaluated to a NaN gradient and the whole
  ## fit stopped with "NA/NaN gradient evaluation"
  set.seed(279253)
  n <- 18
  a <- as.numeric(stats::arima.sim(list(ar = 0.7), n, sd = 0.8 * sqrt(1 - 0.49)))
  d <- data.frame(g = factor(1), t = seq_len(n))
  d$y <- 5 + a + stats::rnorm(n, 0, 0.5)
  d <- d[d$t <= 12, ]
  f <- NULL
  expect_no_error(f <- q(ilm_model(y ~ 1, data = d, family = "gaussian",
                                   ar = ilm_ar1(~ t | g), verbose = FALSE)))
  expect_s3_class(f, "ilm_model")
  ## and whatever it reached is graded by the checks, not hidden
  expect_true(all(f$checks$status[f$checks$check %in% c("optimizer", "gradient")] %in%
                    c("OK", "WARN", "FAIL")))
})
