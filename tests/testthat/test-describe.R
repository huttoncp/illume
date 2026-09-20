# The describe layer. These check behaviour a user would notice, and in a few
# places they pin bugs that were found by testing rather than by inspection.

test_that("numeric output carries the documented columns", {
  d <- ilm_sim(n_id = 20)
  r <- ilm_describe(d, "score")
  expect_true(all(c("obs", "n", "na", "mean", "sd", "se", "p0", "p50", "p100",
                    "gauss", "gauss_note") %in% names(r)))
  expect_equal(nrow(r), 1L)
  # removed on purpose: p_na belongs to ilm_describe_na(), cases became obs
  expect_false(any(c("p_na", "cases") %in% names(r)))
  # off by default, since gauss and gauss_note say whether they are needed
  expect_false(any(c("skew", "kurt") %in% names(r)))
})

test_that("digits reaches the quantile columns", {
  # it did not: every other numeric column honoured digits and the percentiles
  # printed unrounded, which is only visible on a large-scale variable
  d <- ilm_sim(n_id = 20)
  r <- ilm_describe(d, "income", digits = 1)
  for (cn in c("p0", "p50", "p100"))
    expect_equal(r[[cn]], round(r[[cn]], 1), info = cn)
})

test_that("probs drives both the quantiles and the column names", {
  d <- ilm_sim(n_id = 20)
  r <- ilm_describe(d, "score", probs = c(0.025, 0.25, 0.5, 0.975))
  expect_true(all(c("p2.5", "p25", "p50", "p97.5") %in% names(r)))
  expect_false("p0" %in% names(r))
  expect_equal(unname(r$p50), unname(stats::median(d$score)), tolerance = 1e-6)
})

test_that("skew and kurt appear only when asked for", {
  d <- ilm_sim(n_id = 20)
  r <- ilm_describe(d, "score", skew = TRUE, kurt = TRUE)
  expect_true(all(c("skew", "kurt") %in% names(r)))
})

test_that("dispersion separates Poisson from negative binomial", {
  d <- ilm_sim()
  r <- ilm_describe_all(d, class = "numeric")
  disp <- setNames(r$dispersion, r$variable)
  expect_lt(abs(disp[["visits"]] - 1), 0.6)   # Poisson: variance ~ mean
  expect_gt(disp[["claims"]], 2)              # negative binomial: overdispersed
  expect_true(is.na(disp[["score"]]))         # continuous: not a count
})

test_that("a variable that is not a count gets no dispersion", {
  # an earlier rule keyed on support starting near zero gave 1-100 a ratio and
  # denied 101-200 one, though they are the same kind of variable
  d <- data.frame(a = as.integer(sample(1:100, 400, TRUE)),
                  b = as.integer(sample(101:200, 400, TRUE)),
                  c = rnorm(400))
  r <- ilm_describe_all(d, class = "numeric")
  disp <- setNames(r$dispersion, r$variable)
  expect_false(is.na(disp[["a"]]))
  expect_false(is.na(disp[["b"]]))            # treated the same as a
  expect_true(is.na(disp[["c"]]))
})

test_that("constant columns are split into their own section", {
  d <- ilm_sim(n_id = 20)
  r <- ilm_describe_all(d)
  expect_true("constant" %in% names(r))
  expect_true("cohort" %in% r$constant$variable)
  # and they do not appear in the ordinary sections
  expect_false("cohort" %in% r$categorical$variable)
  expect_true(all(c("variable", "class", "obs", "n", "na", "value") %in%
                    names(r$constant)))
})

test_that("categorical checks fire on the things that break models", {
  d <- ilm_sim()
  r <- ilm_describe_all(d, class = "categorical")
  g <- r[r$variable == "grp", ]
  s <- r[r$variable == "site", ]
  expect_equal(g$n_unused, 1L)          # a declared but unobserved level
  expect_gte(g$n_rare, 1L)              # a level too thin to estimate
  expect_gt(s$n_empty, 0L)              # "" is not NA
  expect_gte(s$case_variants, 1L)       # North / north / "North "
  expect_false("id_like" %in% names(r)) # moved to ilm_frame_issues()
})

test_that("a near-constant logical is flagged as a separation risk", {
  d <- ilm_sim()
  r <- ilm_describe_all(d, class = "logical")
  expect_match(r$note[r$variable == "consented"], "separation risk")
  expect_equal(r$note[r$variable == "flag"], "")
})

test_that("date summaries report spacing, gaps and duplicates", {
  d <- ilm_sim()
  r <- ilm_describe_all(d, class = "time")
  expect_equal(r$spacing, "monthly")
  expect_true(r$regular)
  expect_equal(r$n_gaps, 0L)
  # repeated dates are normal in long format, so the note describes rather
  # than warns
  expect_match(r$note, "long format")
})

test_that("calendar spacing is not called irregular", {
  # months run 28 to 31 days, so an exact test condemns every monthly series
  m <- data.frame(d = seq(as.Date("2020-01-01"), by = "month", length.out = 60))
  q <- data.frame(d = seq(as.Date("2020-01-01"), by = "quarter", length.out = 40))
  expect_true(ilm_describe(m, "d")$regular)
  expect_true(ilm_describe(q, "d")$regular)
})

test_that("real gaps are still caught", {
  g <- data.frame(d = as.Date("2020-01-01") + c(0:40, 60:98, 150:169))
  r <- ilm_describe(g, "d")
  expect_false(r$regular)
  expect_equal(r$n_gaps, 2L)
  expect_match(r$note, "AR\\(1\\)")
})

test_that("by works on describe and on describe_all", {
  d <- ilm_sim(n_id = 30)
  r1 <- ilm_describe(d, "score", by = "grp")
  expect_gt(nrow(r1), 1L)
  expect_true("grp" %in% names(r1))
  r2 <- ilm_describe_all(d, by = "grp", class = "numeric")
  expect_true("grp" %in% names(r2))
  expect_false("grp" %in% r2$variable)   # a grouping variable is not described
})

test_that("misspecified arguments name the valid options", {
  d <- ilm_sim(n_id = 20)
  expect_error(ilm_describe_all(d, class = "numerical"), "Options are")
  expect_error(ilm_describe_all(d, class = "numerical"), "'numeric'")
  expect_error(ilm_describe_all(d, by = "nope"), "not found")
  expect_error(ilm_describe(d, "score", probs = c(-1, 2)), "between 0 and 1")
})

test_that("ilm_frame_issues finds pairwise problems", {
  d <- data.frame(x = rnorm(100), k = 1L, id = seq_len(100),
                  code = paste0("a", 1:100), stringsAsFactors = FALSE)
  d$dup <- d$x
  r <- ilm_frame_issues(d)
  expect_true("constant" %in% r$issue)
  expect_true("duplicate_columns" %in% r$issue)
  expect_true("collinear" %in% r$issue)
  # a continuous column is unique per row by construction, not an identifier
  idl <- r$columns[r$issue == "id_like"]
  expect_true(all(c("id", "code") %in% idl))
  expect_false("x" %in% idl)
})
