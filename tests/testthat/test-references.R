## inst/REFERENCES.bib holds the works the package's @references blocks cite,
## one entry each, under keys that documents cite. The two must not drift
## apart: a work cited with no entry cannot be cited by key, and an entry
## nothing cites is a work the package no longer relies on. They are matched
## by the first author's surname and the year, which every reference begins
## with.
##
## Most blocks become help pages. A few document internal functions (@noRd)
## and exist only in the source, so an entry cited only there is marked
## `keywords = {internal}`: under R CMD check, where the source is not
## installed, those are checked for their form only.

## the help pages: the source's man/ while developing, the installed
## package's Rd database under R CMD check
ref_rd <- function() {
  man <- test_path("..", "..", "man")
  if (dir.exists(man))
    lapply(list.files(man, "[.]Rd$", full.names = TRUE), tools::parse_Rd)
  else unname(as.list(tools::Rd_db("illume")))
}

ref_bib <- function() {
  src <- test_path("..", "..", "inst", "REFERENCES.bib")
  path <- if (file.exists(src)) src else
    system.file("REFERENCES.bib", package = "illume")
  readLines(path, encoding = "UTF-8", warn = FALSE)
}

## each reference paragraph: its first author's surname and year
ref_rows <- function(pars) {
  out <- lapply(pars, function(par) {
    par <- trimws(gsub("[[:space:]]+", " ", par))
    if (!nzchar(par)) return(NULL)
    yr <- regmatches(par, regexpr("\\(([0-9]{4})\\)", par))
    data.frame(surname = trimws(sub(",.*", "", par)),
               year = if (length(yr)) gsub("[()]", "", yr) else NA_character_,
               text = par, stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}

## the references on the help pages
rd_refs <- function(rds) {
  txt <- function(x) if (is.list(x)) paste(vapply(x, txt, ""), collapse = "")
                     else paste(as.character(x), collapse = "")
  pars <- character()
  for (rd in rds) {
    tags <- vapply(rd, function(s) attr(s, "Rd_tag") %||% "", "")
    for (sec in rd[tags == "\\references"])
      pars <- c(pars, strsplit(txt(sec), "\n[[:space:]]*\n")[[1L]])
  }
  ref_rows(pars)
}

## the @references blocks of the R source, internal ones included -- NULL
## where the source is not at hand (R CMD check)
src_refs <- function() {
  rdir <- test_path("..", "..", "R")
  if (!dir.exists(rdir)) return(NULL)
  pars <- character()
  for (f in list.files(rdir, "[.]R$", full.names = TRUE)) {
    l <- readLines(f, encoding = "UTF-8", warn = FALSE)
    for (s0 in grep("^#' @references", l)) {
      e <- s0
      while (e + 1L <= length(l) && grepl("^#'", l[e + 1L]) &&
             !grepl("^#' @", l[e + 1L])) e <- e + 1L
      blk <- sub("^#' ?", "", l[s0:e])
      blk[1L] <- sub("^@references ?", "", blk[1L])
      pars <- c(pars, strsplit(paste(blk, collapse = "\n"),
                               "\n[[:space:]]*\n")[[1L]])
    }
  }
  ref_rows(pars)
}

## each bib entry: its key, first author's surname, year, and whether only
## internal documentation cites it
bib_entries <- function(lines) {
  starts <- grep("^@[a-z]+\\{", lines)
  ends <- c(starts[-1L] - 1L, length(lines))
  field <- function(block, f) {
    l <- grep(paste0("^  ", f, " = \\{"), block, value = TRUE)
    if (!length(l)) NA_character_
    else sub(paste0("^  ", f, " = \\{(.*)\\},?$"), "\\1", l[1L])
  }
  do.call(rbind, lapply(seq_along(starts), function(i) {
    b <- lines[starts[i]:ends[i]]
    data.frame(key = sub("^@[a-z]+\\{([^,]+),.*$", "\\1", b[1L]),
               surname = gsub("[{}]", "", sub(",.*", "", field(b, "author"))),
               year = field(b, "year"),
               internal = identical(field(b, "keywords"), "internal"),
               stringsAsFactors = FALSE)
  }))
}

sy <- function(d) paste(d$surname, d$year)

## how many entries each reference matches: exactly one is right
n_entries <- function(refs, bib)
  vapply(sy(refs), function(k) sum(sy(bib) == k), 0L)

test_that("every work a help page cites has one entry in REFERENCES.bib", {
  refs <- rd_refs(ref_rd())
  bib <- bib_entries(ref_bib())
  expect_gt(nrow(refs), 40L)
  expect_false(anyNA(refs$year), label = "a reference with no (year)")
  expect_identical(refs$text[n_entries(refs, bib) != 1L], character(0),
                   label = "references with no single REFERENCES.bib entry")
})

test_that("every internal @references work has an entry too", {
  src <- src_refs()
  skip_if(is.null(src), "the R source is not at hand")
  bib <- bib_entries(ref_bib())
  expect_false(anyNA(src$year), label = "a reference with no (year)")
  expect_identical(src$text[n_entries(src, bib) != 1L], character(0),
                   label = "@references with no single REFERENCES.bib entry")
  ## an entry marked internal is cited in the source and on no help page
  onpage <- sy(bib) %in% sy(rd_refs(ref_rd()))
  expect_identical(bib$key[bib$internal & onpage], character(0),
                   label = "entries marked internal that a help page cites")
  expect_identical(bib$key[bib$internal & !sy(bib) %in% sy(src)],
                   character(0), label = "internal entries nothing cites")
})

test_that("every REFERENCES.bib entry is cited, under a stable key", {
  refs <- rd_refs(ref_rd())
  bib <- bib_entries(ref_bib())
  ## keys are unique, and are the surname's letters and the year
  expect_false(anyDuplicated(bib$key) > 0L)
  expect_true(all(grepl("^[A-Za-z]+[0-9]{4}[a-z]?$", bib$key)))
  expect_identical(sub("^[A-Za-z]+([0-9]{4})[a-z]?$", "\\1", bib$key), bib$year)
  ## an internal entry is cited in the source, which the test above reads
  cited <- sy(bib) %in% sy(refs) | bib$internal
  expect_identical(bib$key[!cited], character(0),
                   label = "REFERENCES.bib entries no help page cites")
})
