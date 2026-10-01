## Separation: an outcome with no variation in a level has no finite
## coefficient there. The case behind the ruling: three all-zero rows in one
## level of a Poisson fit gave a coefficient near minus infinity, every check
## OK and a predicted count of 5e-8. Caught before the fit from the data and
## after it from a flat likelihood; the check FAILs, the interpretation says
## not estimable, and an empty level is announced rather than dropped quietly.

pois_sep <- function() {
  set.seed(3)
  d <- data.frame(g = factor(rep(c("a", "b", "c", "d"), c(40, 40, 40, 3))), x = stats::rnorm(123))
  d$y <- stats::rpois(123, exp(0.5 + 0.3 * d$x)); d$y[d$g == "d"] <- 0
  d
}
sep_row <- function(f) f$checks[f$checks$check == "separation", ]

test_that("a level whose counts are all zero is named, and the fit is not ok", {
  expect_message(f <- ilm_model(y ~ x + g, data = pois_sep(), family = "poisson", verbose = FALSE),
                 "does not vary within g 'd'")
  r <- sep_row(f)
  expect_identical(r$status, "FAIL")
  expect_match(r$detail, "g = d: every count zero (3 rows)", fixed = TRUE)
  expect_false(f$ok)
  expect_identical(f$separation$levels$level, "d")
  ## the other coefficients are still what they would be without the level
  f2 <- suppressMessages(ilm_model(y ~ x + g, data = droplevels(subset(pois_sep(), g != "d")),
                                   family = "poisson", verbose = FALSE))
  expect_equal(unname(coef(f)[c("x", "gb", "gc")]), unname(coef(f2)[c("x", "gb", "gc")]),
               tolerance = 1e-4)
  ## the remedies name the level, and the interpretation says not estimable
  rm <- ilm_remedies(f)
  expect_true(any(grepl("merge 'd'", rm$remedy, fixed = TRUE)))
  txt <- paste(capture.output(print(ilm_interpret(f))), collapse = " ")
  expect_match(gsub("\\s+", " ", txt), "g: not estimable for 'd'", fixed = TRUE)
  expect_match(gsub("\\s+", " ", txt), "estimated as usual", fixed = TRUE)
})

test_that("a separated reference level makes every comparison not estimable", {
  d <- pois_sep(); d$g <- stats::relevel(d$g, "d")
  f <- suppressMessages(ilm_model(y ~ x + g, data = d, family = "poisson", verbose = FALSE))
  expect_identical(sep_row(f)$status, "FAIL")
  txt <- gsub("\\s+", " ", paste(capture.output(print(ilm_interpret(f))), collapse = " "))
  expect_match(txt, "none of the comparisons is estimable", fixed = TRUE)
})

test_that("a binomial arm of all successes, and a cell of two factors, are named", {
  set.seed(4)
  b <- data.frame(arm = factor(rep(c("p", "q", "r"), each = 20)), x = stats::rnorm(60))
  b$y <- stats::rbinom(60, 1, 0.4); b$y[b$arm == "r"] <- 1
  f <- suppressMessages(ilm_model(y ~ x + arm, data = b, family = "binomial", verbose = FALSE))
  expect_match(sep_row(f)$detail, "arm = r: every trial a success (20 rows)", fixed = TRUE)
  b$s <- factor(rep(c("m", "f"), 30)); b$y <- stats::rbinom(60, 1, 0.5)
  b$y[b$arm == "q" & b$s == "f"] <- 0
  f2 <- suppressMessages(ilm_model(y ~ arm * s, data = b, family = "binomial", verbose = FALSE))
  expect_match(sep_row(f2)$detail, "arm:s = q:f: every trial a failure", fixed = TRUE)
})

test_that("separation by a numeric predictor is caught from the flat likelihood", {
  d <- data.frame(x = seq(-2, 2, length.out = 60)); d$y <- as.integer(d$x > 0)
  f <- suppressMessages(ilm_model(y ~ x, data = d, family = "binomial", verbose = FALSE))
  r <- sep_row(f)
  expect_identical(r$status, "FAIL")
  expect_match(r$detail, "coefficient x = ", fixed = TRUE)
  expect_true("x" %in% f$separation$flat$coef)
})

test_that("an ordinary fit has an OK separation check and no message", {
  set.seed(5)
  d <- data.frame(g = factor(rep(letters[1:3], each = 30)), x = stats::rnorm(90))
  d$y <- stats::rpois(90, exp(0.5 + 0.2 * d$x))
  expect_silent(f <- ilm_model(y ~ x + g, data = d, family = "poisson", verbose = FALSE))
  expect_identical(sep_row(f)$status, "OK")
  expect_true(f$ok)
  ## and a gaussian fit, where a level cannot be separated, has none
  d$z <- stats::rnorm(90)
  fg <- ilm_model(z ~ x + g, data = d, family = "gaussian", verbose = FALSE)
  expect_equal(nrow(sep_row(fg)), 0L)
})

test_that("a level no row uses is announced and recorded, not dropped quietly", {
  d <- pois_sep(); d <- d[d$g != "d", ]            # "d" stays a level with no rows
  expect_message(f <- ilm_model(y ~ x + g, data = d, family = "poisson", verbose = FALSE),
                 "`g` has a level with no rows (d)", fixed = TRUE)
  expect_identical(f$empty_levels, list(g = "d"))
  expect_identical(sep_row(f)$status, "OK")
})
