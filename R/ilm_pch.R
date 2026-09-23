## ---------------------------------------------------------------------------
## A copy of illumex's plotting-character helpers (illumex: R/ilm_plot.R).
## illume's own plots take a `pch` argument too, and one package should not
## reach into another's internals, so the helpers live in both.
## tests/testthat/test-shared-helpers.R fails if the two copies differ: change
## them together.
## ---------------------------------------------------------------------------

## ---- plotting characters by name -------------------------------------------
##
## Base R's `pch` is 26 integers nobody remembers. Nothing in the argument says
## that 16 is a filled circle, and the difference between 16, 19, 20 and 21 is
## not guessable from the numbers -- so a plot gets whichever code the author
## happened to recall, and a reader comparing two plots cannot tell whether a
## difference in the markers was meant.
##
## Names are also checkable, which numbers are not: a pch outside 0:25 is
## accepted by graphics and quietly draws nothing, while a name that is not in
## the table stops here and lists the ones that are.
##
## A SINGLE character is left alone, because base R draws it literally --
## pch = "x" means the letter x, and translating it would silently turn the
## plot into crosses.

#' Plotting characters, by name
#'
#' The lookup behind the `pch` argument of the plotting functions. Numbers pass
#' through after a range check, a single character passes through (base R draws
#' it literally), and a name becomes the number base R wants.
#'
#' @param x Numeric `pch` codes, single characters, or names such as
#'   `"filled circle"`. Case, spaces, underscores, hyphens and dots are all
#'   ignored, so `"filled_circle"` and `"Filled Circle"` are the same thing.
#' @return An integer `pch` vector, or the input unchanged where it was already
#'   numeric or a single character.
#' @keywords internal
#' @noRd
ilm_pch <- function(x) {
  if (is.null(x) || !length(x)) return(x)
  if (is.numeric(x)) {
    bad <- x[is.finite(x) & (x < 0 | x > 25 | x != as.integer(x))]
    if (length(bad))
      stop("`pch` codes run from 0 to 25; got ",
           paste(unique(bad), collapse = ", "),
           ". Names work too, such as \"filled circle\".", call. = FALSE)
    return(x)
  }
  if (!is.character(x)) return(x)
  tab <- ilm_pch_table()
  key <- tolower(gsub("[ _.-]+", "", trimws(x)))
  out <- vector("list", length(x))
  for (i in seq_along(x)) {
    if (is.na(x[i])) { out[[i]] <- NA_integer_; next }
    if (nchar(x[i]) == 1L) { out[[i]] <- x[i]; next }   # a literal glyph
    j <- match(key[i], names(tab))
    if (is.na(j))
      stop("'", x[i], "' is not a plotting character name. ",
           "The names are: ", paste(ilm_pch_names(), collapse = ", "), ".",
           call. = FALSE)
    out[[i]] <- tab[[j]]
  }
  ## a mix of glyphs and codes has to stay character, since that is the only
  ## vector that can carry both
  if (any(vapply(out, is.character, TRUE))) as.character(unlist(out))
  else as.integer(unlist(out))
}

## The names, in code order, primary name first. Aliases follow it, because a
## person reaching for "solid circle" should not have to find out that the package
## calls it something else.
#' @keywords internal
#' @noRd
ilm_pch_spec <- function() {
  list(
    c("0",  "open square", "square", "hollow square", "empty square"),
    c("1",  "open circle", "circle", "hollow circle", "empty circle"),
    c("2",  "open triangle", "triangle", "triangle up", "hollow triangle"),
    c("3",  "plus"),
    c("4",  "cross", "times"),
    c("5",  "open diamond", "diamond", "hollow diamond"),
    c("6",  "open triangle down", "triangle down", "down triangle"),
    c("7",  "square cross", "crossed square"),
    c("8",  "star", "asterisk"),
    c("9",  "diamond plus"),
    c("10", "circle plus"),
    c("11", "star of david", "double triangle"),
    c("12", "square plus"),
    c("13", "circle cross"),
    c("14", "square triangle"),
    c("15", "filled square", "solid square"),
    c("16", "filled circle", "solid circle", "point", "dot"),
    c("17", "filled triangle", "solid triangle"),
    c("18", "filled diamond", "solid diamond"),
    c("19", "large filled circle", "bold circle"),
    c("20", "small filled circle", "bullet", "small dot"),
    c("21", "circle fill", "filled circle outline", "bg circle"),
    c("22", "square fill", "filled square outline", "bg square"),
    c("23", "diamond fill", "filled diamond outline", "bg diamond"),
    c("24", "triangle fill", "filled triangle outline", "bg triangle"),
    c("25", "triangle down fill", "filled triangle down"))
}

#' @keywords internal
#' @noRd
ilm_pch_table <- function() {
  out <- list()
  for (s in ilm_pch_spec()) {
    code <- as.integer(s[1])
    for (nm in s[-1]) out[[gsub("[ _.-]+", "", nm)]] <- code
  }
  out
}

## The primary name of each code, spelled the way a person would write it --
## the lookup keys have had their spaces stripped and would read badly here.
#' @keywords internal
#' @noRd
ilm_pch_names <- function() {
  vapply(ilm_pch_spec(), function(s) sprintf("%s (%s)", s[2], s[1]), "")
}
