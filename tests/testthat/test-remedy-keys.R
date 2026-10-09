## Every remedy carries a key, `kind/target`, naming the change itself. The id
## numbers the rows of one table and changes when tables are combined or a
## standalone check is included; the key does not, so a script that applies
## "drop/g" applies the same remedy in another session or stops. The
## generics, c() and print() are illumex's; illume registers its methods on
## them and exports them again.

zero_fit <- function() {
  set.seed(1)
  d <- data.frame(g = factor(rep(1:20, each = 10)), x = rnorm(200), e = rnorm(200))
  d$x <- d$x - ave(d$x, d$g)
  d$y <- 1 + 0.5 * d$x + d$e - ave(d$e, d$g)
  ilm_model(y ~ x + (1 | g), data = d, verbose = FALSE)
}

over_fit <- function() {
  set.seed(2)
  d <- data.frame(x = stats::rnorm(300))
  d$y <- stats::rnbinom(300, mu = exp(1 + 0.5 * d$x), size = 1.5)
  ilm_model(y ~ x, data = d, family = "poisson", verbose = FALSE)
}

test_that("every remedy has a key, unique within its table", {
  f <- zero_fit()
  rem <- ilm_remedies(f)
  expect_true(nrow(rem) > 0L)
  expect_true(all(grepl("^[^/]+/.+", rem$key)))
  expect_false(anyDuplicated(gsub("\\s", "", rem$key)) > 0L)
  expect_true("drop/g" %in% rem$key)
  expect_true("boundary/avoid" %in% rem$key)
  g <- over_fit()
  dsp <- suppressMessages(ilm_check_dispersion(g, B = 100))
  rg <- ilm_remedies(g, dispersion = dsp)
  expect_false(anyDuplicated(rg$key) > 0L)
  expect_true("family/nbinom" %in% rg$key)
})

test_that("a key applies the same remedy as its id, and keeps doing so", {
  f <- zero_fit()
  rem <- ilm_remedies(f)
  i <- which(rem$key == "drop/g")
  a <- suppressMessages(ilm_apply_remedy(f, rem, rem$id[i]))
  b <- suppressMessages(ilm_apply_remedy(f, rem, "drop/g"))
  expect_identical(deparse(a$formula), deparse(b$formula))
  expect_identical(b$remedy_log$key, "drop/g")
  expect_identical(names(b$remedy_log)[1], "key")
  ## spaces in a key are ignored
  rk <- rem$key[grepl("^boundary/", rem$key)]
  expect_identical(illumex::ilm_remedy_find(rem, gsub("/", " / ", rk)),
                   which(rem$key == rk))
  ## a key the list does not have stops, naming the keys it has
  expect_error(ilm_apply_remedy(f, rem, "drop/h"), "drop/g")
})

test_that("a key survives what renumbers the ids", {
  g <- over_fit()
  dsp <- suppressMessages(ilm_check_dispersion(g, B = 100))
  zz <- suppressMessages(ilm_check_zeros(g, B = 100))
  r1 <- ilm_remedies(g, dispersion = dsp)
  r2 <- ilm_remedies(g, zeros = zz, dispersion = dsp)
  k <- "family/nbinom"
  expect_identical(r1$change[r1$key == k], r2$change[r2$key == k])
  ## another package's table, keyed by its arguments as illume keys them, so
  ## c() lists the change both name once
  ext <- ilm_remedy_table(g, "pit_shape", "FAIL", "structural",
                          "refit as a negative binomial",
                          args = list(list(family = "nbinom")))
  expect_identical(ext$key, k)
  all <- c(ext, r1)
  expect_identical(sum(all$key == k), 1L)
  expect_match(all$check[all$key == k], "pit_shape")
  expect_match(all$check[all$key == k], "ilm_check_dispersion")
  f2 <- suppressMessages(ilm_apply_remedy(g, all, k))
  expect_identical(f2$family$name, "nbinom")
})

test_that("another package's remedies made by hand get keys of their own", {
  g <- over_fit()
  ext <- ilm_remedy_table(g, c("pit_shape", "pit_shape", "coverage"),
                          c("FAIL", "FAIL", "WARN"),
                          c("structural", "estimand", "estimand"),
                          c("recalibrate", "use more folds", "widen"))
  expect_identical(sort(ext$key),
                   sort(c("by_hand/pit_shape:1", "by_hand/pit_shape:2",
                          "by_hand/coverage")))
  given <- ilm_remedy_table(g, "pit_shape", "FAIL", "structural", "a zero part",
                            args = list(list(ziformula = ~ 1)),
                            key = "zi/inflated")
  expect_identical(given$key, "zi/inflated")
  expect_error(ilm_remedy_table(g, "pit_shape", "FAIL", "structural", "x",
                                key = "nokind"), "kind/target")
})

test_that("illume's remedy generics are illumex's, so attaching masks nothing", {
  for (f in c("ilm_remedies", "ilm_apply_remedy", "ilm_remedy_table",
              "iml_remedies", "iml_apply_remedy", "iml_remedy_table"))
    expect_identical(getExportedValue("illume", f),
                     getExportedValue("illumex", f), label = f)
  ## and the methods are registered on illumex's generics
  expect_false(is.null(utils::getS3method("ilm_remedies", "ilm_model",
                                          envir = asNamespace("illumex"),
                                          optional = TRUE)))
  skip_on_cran()
  out <- system2(file.path(R.home("bin"), "Rscript"),
                 c("-e", shQuote(paste0(
                   ".libPaths(c(", paste(sprintf("'%s'", .libPaths()), collapse = ","),
                   ")); suppressPackageStartupMessages(library(illumex)); ",
                   "library(illume)"))),
                 stdout = TRUE, stderr = TRUE)
  expect_false(any(grepl("masked", out)),
               label = paste(out, collapse = "\n"))
})
