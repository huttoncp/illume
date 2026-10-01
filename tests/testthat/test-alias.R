## The whole fixed design, the smooths' unpenalised columns included, is
## checked for rank before the fit (item 4, as Craig ruled). Each of these
## used to fit and come back with a failed Hessian and NaN standard errors,
## or a FAIL after the fit, with nothing naming the overlap.

alias_data <- function(n = 300) {
  set.seed(1)
  d <- data.frame(x = stats::runif(n), z = stats::rnorm(n), w = stats::rnorm(n),
                  a = factor(sample(c("p", "q", "r"), n, TRUE)),
                  b = factor(sample(c("u", "v"), n, TRUE)))
  d$x1 <- stats::rnorm(n); d$x2 <- 2 * d$x1
  d$y <- sin(3 * d$x) + 0.5 * d$z + stats::rnorm(n, 0, 0.3)
  d
}
fit_or_msg <- function(f, d) tryCatch({
  suppressWarnings(suppressMessages(ilm_model(f, data = d, family = "gaussian", verbose = FALSE)))
  "fitted"
}, error = function(e) conditionMessage(e))

test_that("two smooths with the same numeric by stop, naming the overlap", {
  m <- fit_or_msg(y ~ s(x, by = z) + s(w, by = z), alias_data())
  expect_match(m, "not all separable", fixed = TRUE)
  expect_match(m, "the unpenalised part of s(w):z", fixed = TRUE)
  expect_match(m, "same numeric `by` variable", fixed = TRUE)
})

test_that("a parametric term inside a smooth's null space stops", {
  m <- fit_or_msg(y ~ x + s(x), alias_data())
  expect_match(m, "the unpenalised part of s(x) is a linear combination of the intercept and `x`",
               fixed = TRUE)
  m2 <- fit_or_msg(y ~ poly(z, 2) + s(x, by = z), alias_data())
  expect_match(m2, "not all separable", fixed = TRUE)
})

test_that("collinear ordinary columns and an empty interaction cell are named", {
  m <- fit_or_msg(y ~ x1 + x2, alias_data())
  expect_match(m, "`x2` is a linear combination of `x1`", fixed = TRUE)
  d <- alias_data(); d <- d[!(d$a == "r" & d$b == "v"), ]
  m2 <- fit_or_msg(y ~ a * b, d)
  expect_match(m2, "is zero in every row", fixed = TRUE)
  expect_match(m2, "combinations with no rows (r:v)", fixed = TRUE)
})

test_that("a design with full rank is untouched, and PR 21's case keeps its own words", {
  expect_identical(fit_or_msg(y ~ w + z:w + s(x, by = z), alias_data()), "fitted")
  expect_identical(fit_or_msg(y ~ x1 + z + s(x), alias_data()), "fitted")
  ## a column in dollars is not mistaken for a dependent one
  d <- alias_data(); d$big <- stats::rnorm(nrow(d)) * 1e6
  expect_identical(fit_or_msg(y ~ big + z, d), "fitted")
  m <- fit_or_msg(y ~ z + s(x, by = z), alias_data())
  expect_match(m, "is a term of the model and also the `by` variable", fixed = TRUE)
})

## aliased = "drop" (item 214): the dependent columns go, as lm() and lme4
## drop them, and the fit says so

test_that("the stop names aliased = \"drop\" where it applies, and only there", {
  m <- fit_or_msg(y ~ x1 + x2, alias_data())
  expect_match(m, "Or set aliased = \"drop\" to drop the dependent column, as lm() does.",
               fixed = TRUE)
  m2 <- fit_or_msg(y ~ x + s(x), alias_data())
  expect_no_match(m2, "aliased = \"drop\"", fixed = TRUE)
  ## asked for with a smooth involved, it stops all the same, and says why
  m3 <- tryCatch(ilm_model(y ~ x + s(x), data = alias_data(), family = "gaussian",
                           aliased = "drop", verbose = FALSE),
                 error = function(e) conditionMessage(e))
  expect_match(m3, "aliased = \"drop\" does not apply: a smooth's unpenalised part",
               fixed = TRUE)
})

test_that("collinear columns are dropped as lm() drops them, and said", {
  d <- alias_data()
  expect_message(f <- ilm_model(y ~ x1 + x2 + z, data = d, family = "gaussian",
                                aliased = "drop", verbose = FALSE),
                 "1 column was dropped as aliased, as aliased = \"drop\" asks: `x2` is a linear combination of `x1`",
                 fixed = TRUE)
  expect_identical(f$aliased$columns, "x2")
  expect_identical(names(stats::coef(f)), c("(Intercept)", "x1", "z"))
  ## the same estimates and standard errors as lm(), whose x2 is NA
  l <- stats::lm(y ~ x1 + x2 + z, data = d)
  expect_true(is.na(stats::coef(l)[["x2"]]))
  expect_equal(unname(stats::coef(f)), unname(stats::coef(l)[c(1, 2, 4)]), tolerance = 1e-6)
  expect_equal(unname(sqrt(diag(stats::vcov(f)))),
               unname(sqrt(diag(stats::vcov(l)))[c(1, 2, 4)]), tolerance = 1e-3)
  ## and the same model as one written without the column
  f0 <- suppressMessages(ilm_model(y ~ x1 + z, data = d, family = "gaussian", verbose = FALSE))
  expect_equal(as.numeric(logLik(f)), as.numeric(logLik(f0)), tolerance = 1e-8)
  ## new rows are built with the same columns
  nd <- d[1:5, ]
  expect_equal(as.vector(predict(f, nd)), as.vector(predict(f0, nd)), tolerance = 1e-8)
  ## said in summary()
  out <- paste(capture.output(summary(f)), collapse = " ")
  expect_match(out, "The fixed-effect columns were not all separable, so 1 column",
               fixed = TRUE)
  ## a term with no column left is not tested, and the table says so
  a <- ilm_anova(f)
  expect_identical(rownames(a), c("x1", "z"))
  expect_match(attr(a, "heading"), "Not tested, every column dropped as aliased: x2",
               fixed = TRUE, all = FALSE)
})

test_that("an interaction's empty cell is dropped, with a random effect beside it", {
  d <- alias_data(); d <- d[!(d$a == "r" & d$b == "v"), ]
  d$g <- factor(rep(1:20, length.out = nrow(d)))
  d$y <- d$y + stats::rnorm(20, 0, 0.5)[d$g]
  f <- suppressMessages(ilm_model(y ~ a * b + (1 | g), data = d, family = "gaussian",
                                  aliased = "drop", verbose = FALSE))
  expect_identical(f$aliased$columns, "ar:bv")
  expect_false("ar:bv" %in% names(stats::coef(f)))
  ## the same fit as the cells written as one factor
  d$ab <- interaction(d$a, d$b, drop = TRUE)
  f1 <- suppressMessages(ilm_model(y ~ ab + (1 | g), data = d, family = "gaussian",
                                   verbose = FALSE))
  expect_equal(as.numeric(logLik(f)), as.numeric(logLik(f1)), tolerance = 1e-6)
  ## predictions and marginal means at the observed cells
  nd <- unique(d[, c("a", "b")])
  nd$ab <- interaction(nd$a, nd$b, drop = TRUE)
  expect_equal(as.vector(predict(f, nd, groups = "typical")),
               as.vector(predict(f1, nd, groups = "typical")), tolerance = 1e-6)
  em <- ilm_emmeans(f, "a")
  expect_true(all(is.finite(em$emmean)))
})

test_that("an ordinal fit drops its intercept and an aliased column together", {
  d <- alias_data()
  d$o <- cut(d$y, c(-Inf, -0.3, 0.3, Inf), labels = c("lo", "mid", "hi"), ordered_result = TRUE)
  f <- suppressMessages(ilm_model(o ~ x1 + x2 + z, data = d, family = "ordinal",
                                  aliased = "drop", verbose = FALSE))
  f0 <- suppressMessages(ilm_model(o ~ x1 + z, data = d, family = "ordinal", verbose = FALSE))
  expect_identical(f$aliased$columns, "x2")
  expect_equal(stats::coef(f), stats::coef(f0), tolerance = 1e-6)
  nd <- d[1:4, ]
  expect_equal(predict(f, nd), predict(f0, nd), tolerance = 1e-6)
})
