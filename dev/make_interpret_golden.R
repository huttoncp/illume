## The interpreter's golden English: one file per case in
## tests/testthat/helper-interpret-golden.R, every sentence as
## ilm_interpret() prints it, section by section. test-interpret-golden.R
## compares word for word. Rerun only for a deliberate change to the English,
## name each changed sentence in the commit, and say so in NEWS.
##
## Run from the package root:  Rscript dev/make_interpret_golden.R [outdir]
suppressMessages(pkgload::load_all(".", quiet = TRUE, helpers = FALSE))
source(file.path("tests", "testthat", "helper-interpret-golden.R"))
a <- commandArgs(TRUE)
out <- if (length(a)) a[1] else file.path("tests", "testthat", "golden", "interpret")
dir.create(out, showWarnings = FALSE, recursive = TRUE)
cases <- golden_interpret_cases()
for (nm in names(cases)) {
  writeLines(enc2utf8(golden_interpret_text(cases[[nm]])), file.path(out, paste0(nm, ".txt")),
             useBytes = TRUE)
  cat(nm, "\n")
}
