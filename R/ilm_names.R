## ---------------------------------------------------------------------------
## Column names that need backticks.
##
## A formula can name any column -- `my y` ~ `x 1` + (1 | `site id`) --
## because it is a language object, and a name inside it is a symbol however
## it is spelled. The trouble starts wherever code turns a name into TEXT and
## parses it back. all.vars() returns names bare, and deparse() of a lone
## symbol leaves its backticks off, so "x 1" goes back in as the symbol x
## followed by the number 1. That round trip was in the formula front end, in
## mgcv, in the prediction code and in half a dozen functions that take column
## names as strings, and every one of them failed on a name with a space in it.
##
## The rule these helpers keep: when a name becomes code, it is quoted; when
## it is used to find a column, it is not; and where a formula can be built
## from language objects instead of text, it is.
## ---------------------------------------------------------------------------

## Names as they must be written in R code: backticked when not syntactic,
## untouched when they are. deparse() of a symbol knows the rules, reserved
## words included, which a comparison against make.names() gets subtly wrong.
#' @keywords internal
#' @noRd
ilm_bq <- function(x) {
  if (!length(x)) return(character(0))
  vapply(as.character(x), function(v) deparse(as.name(v), backtick = TRUE), "",
         USE.NAMES = FALSE)
}

## An expression as TERM text, for reformulate(): backticks wherever needed,
## on one line however long it is. A bare deparse() splits at 60 characters,
## and pasting the pieces of a split deparse into a formula keeps only the
## first of them.
#' @keywords internal
#' @noRd
ilm_term_text <- function(e) deparse1(e, backtick = TRUE)

## An expression as the COLUMN NAME model.frame() gives it: a symbol bare, a
## call with its backticks. This is how stats names the columns of a model
## frame, so it is the key to look one up by.
#' @keywords internal
#' @noRd
ilm_mf_name <- function(e)
  deparse1(e, backtick = !is.symbol(e) && is.language(e))

## A term named the way a person writes it -- "x 1" -- or the way R writes it
## -- "`x 1`" -- resolved to the form the fit stores its labels in, which is
## R's. Anything already a label, or matching nothing, comes back unchanged,
## so the caller's own error about an unknown term still fires.
#' @keywords internal
#' @noRd
ilm_as_label <- function(v, labels) {
  v <- as.character(v)
  miss <- !(v %in% labels)
  if (any(miss)) {
    q <- ilm_bq(v[miss])
    hit <- q %in% labels
    v[miss][hit] <- q[hit]
  }
  v
}

## And back again: a label that is a single backticked name, as the plain
## name a data frame knows it by. A label that is a call -- log(`x 1`) --
## is code, and is left as it is.
#' @keywords internal
#' @noRd
ilm_unbq <- function(x) {
  x <- as.character(x)
  one <- grepl("^`[^`]+`$", x)
  x[one] <- substr(x[one], 2L, nchar(x[one]) - 1L)
  x
}

## `~ e` as a formula, built from the expression itself rather than from text,
## so nothing about how its names are spelled can be lost on the way.
#' @keywords internal
#' @noRd
ilm_one_sided <- function(e, env = parent.frame())
  stats::as.formula(call("~", e), env = env)

## Replace symbols in a formula or expression, by name. `map` is a named list
## of symbols; the attributes of a formula survive, so its environment does.
#' @keywords internal
#' @noRd
ilm_swap_names <- function(f, map) {
  if (!length(map)) return(f)
  if (is.symbol(f) || !is.language(f))
    return(if (is.symbol(f)) do.call(substitute, list(f, map)) else f)
  for (i in seq_along(f)[-1L])
    if (!is.null(f[[i]])) f[[i]] <- do.call(substitute, list(f[[i]], map))
  f
}

## mgcv::interpret.gam() rebuilds its formulas from TEXT, and fails on any
## name that needs backticks -- the response, a covariate or a smooth's
## argument alike, with "unexpected symbol" and nothing to say which. So the
## names that need them are swapped for syntactic stand-ins for that one call,
## and swapped back in everything it returns: the parametric formula, the
## fake formula, and the variable names and labels inside each smooth.
#' @keywords internal
#' @noRd
ilm_interpret_gam <- function(f) {
  v <- all.vars(f)
  bad <- v[ilm_bq(v) != v]
  if (!length(bad)) return(mgcv::interpret.gam(f))
  stand <- paste0(".ilm_nm", seq_along(bad), ".")
  while (any(stand %in% v)) stand <- paste0(stand, ".")
  fwd  <- stats::setNames(lapply(stand, as.name), bad)
  back <- stats::setNames(lapply(bad, as.name), stand)
  gp <- mgcv::interpret.gam(ilm_swap_names(f, fwd))
  gp$pf <- ilm_swap_names(gp$pf, back)
  if (!is.null(gp$fake.formula))
    gp$fake.formula <- ilm_swap_names(gp$fake.formula, back)
  gp$smooth.spec <- lapply(gp$smooth.spec, ilm_swap_spec, from = stand,
                           to = bad)
  for (k in intersect(c("response", "fake.names", "pred.names"), names(gp)))
    if (is.character(gp[[k]]))
      gp[[k]] <- ilm_swap_chr(gp[[k]], stand, bad)
  gp
}

## A smooth specification carries its variables as TEXT -- `term`, `by`, the
## `label` built from them, and the same again inside each margin of a
## tensor product. Names go back bare, because mgcv looks columns up by them;
## labels keep the name as written, which is how mgcv writes its own.
## `label = FALSE` leaves the label alone, for building on stand-ins while
## still reporting the smooth under the name the user wrote.
#' @keywords internal
#' @noRd
ilm_swap_spec <- function(sp, from, to, label = TRUE) {
  for (k in intersect(c("term", "by"), names(sp)))
    if (is.character(sp[[k]])) sp[[k]] <- ilm_swap_chr(sp[[k]], from, to)
  if (label && is.character(sp$label))
    for (i in seq_along(from))
      sp$label <- gsub(from[i], to[i], sp$label, fixed = TRUE)
  if (is.list(sp$margin))
    sp$margin <- lapply(sp$margin, ilm_swap_spec, from = from, to = to,
                        label = label)
  sp
}

## Every variable a smooth reads: its terms, its `by`, and its margins'.
#' @keywords internal
#' @noRd
ilm_spec_vars <- function(sp) {
  v <- c(sp$term, if (is.character(sp$by) && !identical(sp$by, "NA")) sp$by)
  if (is.list(sp$margin)) v <- c(v, unlist(lapply(sp$margin, ilm_spec_vars)))
  unique(v)
}

## mgcv cannot build or evaluate a smooth of a column whose name needs
## backticks: its get.var() parses the name as code whenever a lookup misses,
## outside the try() meant to catch that, and the lookup in the knots always
## misses. So such a smooth is built on syntactic stand-in columns, and the
## map (stand-in -> original) is kept with it; this puts the stand-ins onto
## whatever data the smooth is evaluated on.
#' @keywords internal
#' @noRd
ilm_add_standins <- function(data, map) {
  if (!length(map)) return(data)
  for (s in names(map)) data[[s]] <- data[[map[[s]]]]
  data
}

#' @keywords internal
#' @noRd
ilm_swap_chr <- function(x, from, to) {
  i <- match(x, from)
  x[!is.na(i)] <- to[i[!is.na(i)]]
  x
}
