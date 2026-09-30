## Which rows a model used: the counts ilm_rows_used() gives, stats'
## na.action() on a fit, and the line summary() and print() add when missing
## values took rows.

q <- function(e) suppressMessages(suppressWarnings(e))
mk <- function(n = 60, seed = 1) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n), z = rnorm(n), g = factor(rep(1:12, length.out = n)))
  d$y <- 1 + 0.5 * d$x + rnorm(12)[d$g] + rnorm(n)
  d
}

test_that("complete data: every row used, nothing to say", {
  d <- mk()
  f <- q(ilm_model(y ~ x + z, data = d, family = "gaussian", verbose = FALSE))
  r <- ilm_rows_used(f)
  expect_s3_class(r, "ilm_rows_used")
  expect_identical(c(r$n_input, r$n_used, r$n_dropped, r$n_zero_weight),
                   c(60L, 60L, 0L, 0L))
  expect_length(r$dropped_by, 0L)
  expect_null(stats::na.action(f))
  expect_no_match(capture.output(print(f)), "analysed")
  expect_no_match(capture.output(summary(f)), "^Rows:")
  expect_output(print(r), "Rows: 60 of 60 analysed", fixed = TRUE)
})

test_that("missing values: the counts, the columns, and stats' record agree", {
  d <- mk()
  d$y[c(2, 9)] <- NA          # the response
  d$x[c(9, 20, 31)] <- NA     # row 9 is missing in both
  d$z[40] <- NA
  f <- q(ilm_model(y ~ x + z, data = d, family = "gaussian", verbose = FALSE))
  r <- ilm_rows_used(f)
  expect_identical(c(r$n_input, r$n_used, r$n_dropped), c(60L, 55L, 5L))
  ## a row missing in two columns is dropped once and counted under each
  expect_identical(r$dropped_by, c(y = 2L, x = 3L, z = 1L))
  om <- stats::na.action(f)
  expect_s3_class(om, "omit")
  expect_identical(sort(as.integer(om)), c(2L, 9L, 20L, 31L, 40L))
  ## the rows used are the others, and they are the rows the model has
  expect_identical(nrow(f$X), length(setdiff(seq_len(r$n_input), as.integer(om))))
  ## the column counts add to 6 against 5 rows dropped, and the line says why
  out <- gsub("[[:space:]]+", " ", paste(capture.output(print(f)), collapse = " "))
  expect_match(out, paste0("Rows: 55 of 60 analysed; 5 dropped for missing values ",
                           "(y 2, x 3, z 1; a row missing in several columns ",
                           "counts in each)"), fixed = TRUE)
  expect_output(print(summary(f)), "Rows: 55 of 60 analysed; 5 dropped for missing values", fixed = TRUE)
  ## the line as ruled (item 164), exactly
  expect_identical(illume:::ilm_rows_line(r),
    "Rows: 55 of 60 analysed; 5 dropped for missing values (y 2, x 3, z 1; a row missing in several columns counts in each)")
})

test_that("a grouping factor's missing values count, in a mixed model", {
  d <- mk(); d$g[c(5, 6)] <- NA
  f <- q(ilm_model(y ~ x + (1 | g), data = d, family = "gaussian", verbose = FALSE))
  r <- ilm_rows_used(f)
  expect_identical(c(r$n_used, r$n_dropped), c(58L, 2L))
  expect_identical(r$dropped_by, c(g = 2L))
  expect_output(print(f), "Rows: 58 of 60 analysed; 2 dropped for missing values (g 2)",
                fixed = TRUE)
  expect_identical(sort(as.integer(stats::na.action(f))), c(5L, 6L))
})

test_that("a weight of zero is used, not dropped, and is counted", {
  d <- mk(); d$w <- rep(c(1, 2, 0), length.out = 60)
  f <- q(ilm_model(y ~ x, data = d, family = "gaussian", weights = w, verbose = FALSE))
  r <- ilm_rows_used(f)
  expect_identical(c(r$n_used, r$n_dropped, r$n_zero_weight), c(60L, 0L, 20L))
  expect_output(print(f), "20 with weight zero", fixed = TRUE)
})

test_that("a refit to other data counts its own rows", {
  d <- mk(); d$x[1:3] <- NA
  f <- q(ilm_model(y ~ x, data = d, family = "gaussian", verbose = FALSE))
  d2 <- mk(n = 40, seed = 2); d2$x[7] <- NA
  f2 <- q(ilm_refit(f, data = d2))
  expect_identical(c(ilm_rows_used(f)$n_input, ilm_rows_used(f)$n_dropped), c(60L, 3L))
  r2 <- ilm_rows_used(f2)
  expect_identical(c(r2$n_input, r2$n_used, r2$n_dropped), c(40L, 39L, 1L))
  expect_identical(as.integer(stats::na.action(f2)), 7L)
})

test_that("data that is not a data frame has no input count", {
  d <- mk()
  f <- q(ilm_model(y ~ x, data = as.list(d), family = "gaussian", verbose = FALSE))
  r <- ilm_rows_used(f)
  expect_true(is.na(r$n_input))
  expect_identical(r$n_used, 60L)
  expect_output(print(r), "Rows: 60 analysed", fixed = TRUE)
})

test_that("a DAG model's sets each count the rows they used", {
  set.seed(4); n <- 400
  z <- rnorm(n); w <- 0.8 * z + rnorm(n); x <- 0.5 * z + rnorm(n)
  d <- data.frame(x = x, z = z, w = w, y = 0.4 * x + 0.6 * w + rnorm(n))
  d$z[1:10] <- NA; d$w[11:15] <- NA
  g <- ilm_dag("dag { x [exposure] ; y [outcome] ; z -> x -> y ; z -> w -> y }")
  m <- q(ilm_dag_model(g, d, verbose = FALSE))
  expect_gte(length(m$fits), 2L)
  rr <- ilm_rows_used(m)
  expect_s3_class(rr, "data.frame")
  expect_identical(rr$set, seq_along(m$fits))
  ## each set loses the rows missing in what it adjusts for, and na.action()
  ## on each set's fit names them
  for (k in seq_along(m$fits)) {
    f <- m$fits[[k]]
    adj <- intersect(c("z", "w"), attr(stats::terms(f), "term.labels"))
    lost <- which(Reduce(`|`, lapply(adj, function(v) is.na(d[[v]]))))
    expect_identical(rr$n_dropped[k], length(lost))
    expect_identical(sort(as.integer(stats::na.action(f))), lost)
  }
})
