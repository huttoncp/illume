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
