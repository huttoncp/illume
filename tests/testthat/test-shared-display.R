# The display rules (R/shared-helpers.R, the same byte for byte in each of
# the family's packages), held to the hand-written cases in
# fixtures/format_cases.csv, the same file in every package.

cases_file <- function(f) utils::read.csv(test_path("fixtures", f),
                                          colClasses = "character", na.strings = NULL)

test_that("every display rule gives the hand-written text of every case", {
  cs <- cases_file("format_cases.csv")
  expect_identical(names(cs), c("input", "rule", "digits", "end", "trailing_zeros", "expected",
                                "note"))
  expect_true(all(cs$rule %in% c("signif", "fixed", "percent", "p", "clock", "ordinal")))
  got <- vapply(seq_len(nrow(cs)), function(i)
    ilm_disp(as.numeric(cs$input[i]), cs$rule[i],
             if (nzchar(cs$digits[i])) as.integer(cs$digits[i]),
             end = identical(cs$end[i], "TRUE"),
             ## kept unless a case says dropped (item 256)
             trailing_zeros = !identical(cs$trailing_zeros[i], "FALSE"))$text, "")
  expect_identical(got, cs$expected)
})
