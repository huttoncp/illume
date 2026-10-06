## Imputed values keep their column's type, and an imputation model that
## cannot be fitted is said (review findings 5.1 and 5.4). A binomial or
## multinomial draw used to be assigned into a numeric 0/1, logical or
## character column as the factor's integer codes; and a binary coded 1/2
## failed its model silently, leaving random starting draws in place.

imp_data <- function(seed = 1, n = 300) {
  set.seed(seed)
  d <- data.frame(x = stats::rnorm(n), z = stats::rnorm(n))
  d$b <- stats::rbinom(n, 1, stats::plogis(0.5 * d$x))
  d$y <- 0.5 * d$x + 0.3 * d$z + 0.4 * d$b + stats::rnorm(n)
  d$b[sample(n, 60)] <- NA
  d
}

test_that("a numeric 0/1, logical or factor binary is imputed in its own type and values", {
  d <- imp_data()
  imp <- ilm_impute(d, m = 2, seed = 1, verbose = FALSE)
  for (i in 1:2) {
    b <- imp$imputations[[i]]$b
    expect_true(is.numeric(b))
    expect_setequal(unique(b), c(0, 1))
  }
  d2 <- d; d2$b <- as.logical(d2$b)
  imp2 <- ilm_impute(d2, m = 2, seed = 1, verbose = FALSE)
  expect_true(is.logical(imp2$imputations[[1]]$b))
  expect_false(anyNA(imp2$imputations[[1]]$b))
  d3 <- d; d3$b <- factor(d3$b, labels = c("no", "yes"))
  imp3 <- ilm_impute(d3, m = 2, seed = 1, verbose = FALSE)
  expect_identical(levels(imp3$imputations[[1]]$b), c("no", "yes"))
  expect_false(anyNA(imp3$imputations[[1]]$b))
})

test_that("a character column is imputed with its own labels", {
  set.seed(1); n <- 300
  d <- data.frame(x = stats::rnorm(n), z = stats::rnorm(n))
  d$g <- sample(c("red", "green", "blue"), n, TRUE)
  d$s <- sample(c("m", "f"), n, TRUE)
  d$y <- 0.5 * d$x + stats::rnorm(n)
  d$g[sample(n, 50)] <- NA; d$s[sample(n, 50)] <- NA
  imp <- ilm_impute(d, m = 2, seed = 1, verbose = FALSE)
  g <- imp$imputations[[1]]$g; s <- imp$imputations[[1]]$s
  expect_true(is.character(g)); expect_true(is.character(s))
  expect_setequal(unique(g), c("red", "green", "blue"))
  expect_setequal(unique(s), c("m", "f"))
})

test_that("a binary coded 1/2 is imputed from its model, and the pooled fit recovers it", {
  set.seed(11); n <- 1000
  d <- data.frame(x = stats::rnorm(n), z = stats::rnorm(n))
  d$s <- 1 + stats::rbinom(n, 1, stats::plogis(2 * d$x))
  d$y <- 0.5 * d$x + 0.3 * d$z + 1.0 * d$s + stats::rnorm(n)
  full <- stats::coef(stats::lm(y ~ x + z + s, d))
  d$s[sample(n, 300)] <- NA
  mi <- which(is.na(d$s))
  imp <- expect_no_warning(ilm_impute(d, m = 10, seed = 1, verbose = FALSE))
  s1 <- imp$imputations[[1]]$s
  expect_setequal(unique(s1), c(1, 2))
  ## the imputed rows follow x as the observed ones do (it was -0.02 vs 0.58)
  expect_gt(stats::cor(d$x[mi], s1[mi]), 0.3)
  p <- ilm_mi_pool(imp, y ~ x + z + s, family = "gaussian")
  est <- if (is.data.frame(p)) p$estimate[p$term == "s"] else p$table$estimate[p$table$term == "s"]
  expect_lt(abs(est - full[["s"]]), 0.25)
})

test_that("an imputation model that cannot be fitted is said, not left silent", {
  d <- imp_data()
  local_mocked_bindings(ilm_model = function(...) stop("no fit"))
  expect_warning(ilm_impute(d, m = 2, seed = 1, verbose = FALSE),
                 "the imputation model for `b`", fixed = TRUE)
})
