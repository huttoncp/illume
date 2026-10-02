## ilm_interpret()'s English, word for word, against the golden files in
## golden/interpret/ (written by dev/make_interpret_golden.R). A change here is
## a change users read: regenerate on purpose, name each changed sentence in
## the commit, and say so in NEWS.

test_that("every golden case says exactly what its file says", {
  cases <- golden_interpret_cases()
  dir <- test_path("golden", "interpret")
  expect_setequal(tools::file_path_sans_ext(list.files(dir)), names(cases))
  for (nm in names(cases)) {
    want <- readLines(file.path(dir, paste0(nm, ".txt")), encoding = "UTF-8", warn = FALSE)
    expect_identical(golden_interpret_text(cases[[nm]]), want, info = nm)
  }
})
