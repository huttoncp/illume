# The display rules (R/shared-helpers.R, the same byte for byte in each of
# the family's packages), held to the hand-written cases in
# fixtures/format_cases.csv, the same file in every package.

cases_file <- function(f) utils::read.csv(test_path("fixtures", f),
                                          colClasses = "character", na.strings = NULL)

## Inputs beyond 10^22 or under 10^-22, and those of 17 figures, are written
## in the cases file as hex doubles, which parse exactly on every platform
## (on arm64 macOS a decimal literal is read in double precision only).
test_that("the cases' hex inputs read as the doubles they name", {
  expect_identical(as.numeric("0x1.fffffffffffffp+1023"), .Machine$double.xmax)
  expect_identical(as.numeric("0x1p-1022"), .Machine$double.xmin)
  expect_identical(as.numeric("0x0.0000000000001p-1022"), 2^-1074)
})

test_that("every display rule gives the hand-written text of every case", {
  cs <- cases_file("format_cases.csv")
  expect_identical(names(cs), c("input", "rule", "digits", "end", "trailing_zeros", "expected",
                                "note"))
  expect_true(all(cs$rule %in% c("signif", "fixed", "percent", "p", "clock", "ordinal",
                                 "decimals")))
  ## a group's values (rule "decimals", item 306) and their texts are
  ## separated by ";"
  got <- vapply(seq_len(nrow(cs)), function(i)
    paste(ilm_disp(as.numeric(strsplit(cs$input[i], ";", fixed = TRUE)[[1]]), cs$rule[i],
                   if (nzchar(cs$digits[i])) as.integer(cs$digits[i]),
                   end = identical(cs$end[i], "TRUE"),
                   ## kept unless a case says dropped (item 256)
                   trailing_zeros = !identical(cs$trailing_zeros[i], "FALSE"))$text,
          collapse = ";"), "")
  expect_identical(got, cs$expected)
})

test_that("a group's decimals follow its median and its smallest value, at most 6 (item 306)", {
  expect_identical(ilm_disp_group_decimals(c(12.3456, 0.5, 250.75)), 1L)
  expect_identical(ilm_disp_group_decimals(c(3, 120, -7, NA, Inf)), 0L)
  ## a value that will be scientific sets nothing (more than 6 decimals to
  ## show a figure, or 1e15 and up); one at exactly 6 stays in
  expect_identical(ilm_disp_group_decimals(c(1.25, 4e-8, 3.5)), 2L)
  expect_identical(ilm_disp_group_decimals(c(4e-8, 2e-8)), 0L)
  expect_identical(ilm_disp_group_decimals(c(1.25, 4e-6)), 6L)
  expect_identical(ilm_disp_group_decimals(c(0.5, 2e15)), 3L)
  expect_identical(ilm_disp_group_decimals(c(-0.034, 0.0125, -0.5)), 4L)
  expect_identical(ilm_disp_group_decimals(c(0, NA, NaN)), 0L)
  expect_identical(ilm_disp_group_decimals(c(1234.5, 2000.25)), 0L)  # big values need none
  r <- ilm_disp(c(0.5123, 0.4, 0.62), "decimals")
  expect_identical(r, list(text = c("0.512", "0.400", "0.620"), rule = "decimals", digits = 3L))
  expect_identical(ilm_disp(c(-0.0004, 0.5), "decimals")$text, c("-0.0004", "0.5000"))
  expect_identical(ilm_disp(c(1e-9, 2), "decimals")$text, c("1.00e-09", "2"))  # what is left is whole
  expect_identical(ilm_disp(c(0.5, 2e15), "decimals")$text, c("0.500", "2.00e+15"))
  expect_identical(ilm_disp(c(1.5, NA, -Inf), "decimals")$text, c("1.50", NA, "-Inf"))
  expect_error(ilm_disp(1.5, "decimals", c(1, 2)), "takes one `digits`")
})

test_that("a rule returns its text with the rule and its precision", {
  expect_identical(ilm_disp(0.2725, "percent", 1L), list(text = "27.3%", rule = "percent", digits = 1L))
  expect_identical(ilm_disp(0.0004, "p"), list(text = "< 0.001", rule = "p"))
  expect_identical(ilm_disp(2, "clock", end = TRUE), list(text = "02:59", rule = "clock", end = TRUE))
  expect_identical(ilm_disp(21, "clock"), list(text = "21:00", rule = "clock"))
  ## zeros kept is the rule and says nothing; dropped is the exception, said
  expect_identical(ilm_disp(2.5, "signif", 3L),
                   list(text = "2.50", rule = "signif", digits = 3L))
  expect_identical(ilm_disp(2.5, "signif", 3L, trailing_zeros = FALSE),
                   list(text = "2.5", rule = "signif", digits = 3L, trailing_zeros = FALSE))
  ## each value's text depends on it alone
  expect_identical(ilm_disp(c(1230, 12345.6, NA), "signif", 3L)$text, c("1,230", "12,300", NA))
  expect_error(ilm_disp(1, "signif"), "needs `digits`")
})
