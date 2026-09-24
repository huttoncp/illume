## Average marginal effects and the written interpretation of a multinomial
## fit -- the flagship model -- category by category. Both used to report ONE
## category without saying which: ilm_ame() took the last column of the
## predicted probabilities, and ilm_interpret() described only the first
## category's coefficients, quoting that last-category effect beside each.

mn_data <- function(n = 600L, seed = 3L) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n),
                  g = factor(sample(c("north", "central", "south"), n, TRUE)))
  eta <- cbind(0.8 * d$x, -0.5 * d$x + 0.7 * (d$g == "south"))
  P <- exp(cbind(eta, 0)); P <- P / rowSums(P)
  d$y <- factor(apply(P, 1, function(p) sample(c("lo", "mid", "hi"), 1, prob = p)),
                levels = c("lo", "mid", "hi"))
  d
}

test_that("a multinomial fit has an average marginal effect for every category", {
  d <- mn_data()
  f <- ilm_model(y ~ x + g, data = d, family = "multinomial", verbose = FALSE)
  a <- ilm_ame(f)
  expect_true("category" %in% names(a))
  ## every term and level, three categories each
  expect_equal(nrow(a), 3L * 3L)                 # x, g=north... one slope, two contrasts
  expect_setequal(unique(a$category), c("lo", "mid", "hi"))
  ## the probabilities sum to one, so the effects on them sum to zero
  key <- paste(a$term, a$level)
  sums <- tapply(a$estimate, key, sum)
  expect_true(all(abs(sums) < 1e-8))
  ## and each category's slope is the average derivative of ITS probability
  h <- 1e-4
  p1 <- predict(f, newdata = transform(d, x = x + h), type = "response")
  p0 <- predict(f, newdata = transform(d, x = x - h), type = "response")
  byhand <- colMeans((p1 - p0) / (2 * h))
  sl <- a[a$term == "x", ]
  expect_equal(sl$estimate[match(names(byhand), sl$category)], unname(byhand),
               tolerance = 1e-6)
})

test_that("the slopes agree with marginaleffects on the same model", {
  skip_if_not_installed("marginaleffects")
  skip_if_not_installed("nnet")
  d <- mn_data()
  f <- ilm_model(y ~ x + g, data = d, family = "multinomial", verbose = FALSE)
  m <- nnet::multinom(y ~ x + g, data = d, trace = FALSE)
  ref <- as.data.frame(marginaleffects::avg_slopes(m, variables = "x",
                                                  type = "probs"))
  a <- ilm_ame(f, "x")
  ## the two fits agree to about 1e-3 on the coefficients, so the effects do
  ## to about that as well
  expect_equal(a$estimate[match(ref$group, a$category)], ref$estimate,
               tolerance = 1e-3)
})

test_that("the interpretation describes every category, with its own effect", {
  d <- mn_data()
  f <- ilm_model(y ~ x + g, data = d, family = "multinomial", verbose = FALSE)
  txt <- paste(unlist(ilm_interpret(f)$sections), collapse = "\n")
  ## one sentence for x, with the predicted share of every category at both
  ## ends of its middle half -- not one category's figures quoted for all
  line <- grep("^x: ", strsplit(txt, "\n")[[1]], value = TRUE)
  expect_length(line, 1L)
  sh <- regmatches(line, regexec(paste0(
    "of 'lo' is ([0-9.]+)% against ([0-9.]+)%, of 'mid' ([0-9.]+)% against ",
    "([0-9.]+)% and of 'hi' ([0-9.]+)% against ([0-9.]+)%"), line))[[1]]
  expect_length(sh, 7L)
  ## one row per end, one column per category
  sh <- matrix(as.numeric(sh[-1]), 2L)
  ## as simulated, 'lo' rises with x and 'mid' falls
  expect_gt(sh[2, 1], sh[1, 1])
  expect_lt(sh[2, 2], sh[1, 2])
  ## and at either end the three shares make a whole, to rounding
  expect_true(all(abs(rowSums(sh) - 100) <= 2))
  ## the factor too: every category, level by level
  expect_match(txt, paste0("'lo' [0-9.]+%, [0-9.]+% and [0-9.]+%; 'mid' [0-9.]+%, ",
                           "[0-9.]+% and [0-9.]+%; and 'hi' [0-9.]+%, [0-9.]+% and [0-9.]+%"))
})
