## ---------------------------------------------------------------------------
## Which rows a model used.
##
## A table describing "the rows the model used" -- a Table 1 beside the model,
## an analysed-rows column -- needs the counts the fit made them from: the
## rows it was given, the rows it kept, the rows missing values took, and
## which columns took them. The fit keeps stats' own record of the dropped
## rows in `fit$na.action`, so stats::na.action(fit) finds them without
## illume; this adds the counts, which that record cannot give on its own
## (it does not know how many rows came in).
## ---------------------------------------------------------------------------

## The record, made where the model frame is built: `vars` are the data's
## columns the model reads, `w` the frequency weights of the rows kept (NULL
## for none). A row missing in two columns counts once in n_dropped and once
## under each column, so the per-column counts can sum to more.
#' @keywords internal
#' @noRd
ilm_rows_record <- function(data, mf, vars, w = NULL) {
  om <- attr(mf, "na.action")
  n_used <- nrow(mf)
  if (!is.data.frame(data))
    return(list(n_input = NA_integer_, n_used = n_used, n_dropped = NA_integer_,
                dropped_by = integer(0), n_zero_weight = ilm_rows_zero(w)))
  vars <- intersect(vars, names(data))
  idx <- as.integer(om)
  by <- vapply(vars, function(v) sum(is.na(data[[v]][idx])), 0L)
  list(n_input = nrow(data), n_used = n_used, n_dropped = nrow(data) - n_used,
       dropped_by = by[by > 0L], n_zero_weight = ilm_rows_zero(w))
}

#' @keywords internal
#' @noRd
ilm_rows_zero <- function(w) if (is.null(w)) 0L else sum(w == 0)

#' The rows a model used
#'
#' How many rows the model was given, how many it used, how many rows with
#' missing values it dropped, and which columns those were missing in.
#'
#' @details
#' `ilm_model()` drops a row missing in any column the model reads (the
#' default `na.action = na.omit`). Nothing else drops a row: there is no
#' `subset` argument, and a row with a frequency weight of zero stays in the
#' fit -- it counts as used, contributing nothing, and is counted as
#' `n_zero_weight`.
#'
#' Which rows went is stats' own record, `stats::na.action(fit)`: their
#' positions in `data`, with the row names. The rows used are the others.
#'
#' A row missing in two columns is dropped once and counted under each, so
#' `dropped_by` can sum to more than `n_dropped`.
#'
#' An [ilm_dag_model()] fits one model per adjustment set, and sets adjust
#' for different columns, so each can lose different rows: its method gives
#' one row per set.
#'
#' @param object A fitted `"ilm_model"`, or an `"ilm_dag_model"`.
#' @param ... Unused.
#' @return For a model, a list of class `"ilm_rows_used"`: `n_input` (the rows
#'   of `data`; `NA` when `data` was not a data frame), `n_used`, `n_dropped`,
#'   `dropped_by` (a named integer vector: for each column the model reads,
#'   the dropped rows missing in it) and `n_zero_weight`. For an
#'   `"ilm_dag_model"`, a data frame with one row per adjustment set, and the
#'   per-column counts as a list column.
#' @seealso [stats::na.action()], [ilm_model()].
#' @examples
#' set.seed(1)
#' d <- data.frame(x = rnorm(50), z = rnorm(50))
#' d$y <- 1 + d$x + rnorm(50)
#' d$x[1:4] <- NA; d$z[3:6] <- NA
#' fit <- ilm_model(y ~ x + z, data = d, family = "gaussian", verbose = FALSE)
#' ilm_rows_used(fit)
#' stats::na.action(fit)
#' @export
ilm_rows_used <- function(object, ...) UseMethod("ilm_rows_used")

#' @rdname ilm_rows_used
#' @export
ilm_rows_used.ilm_model <- function(object, ...) {
  r <- object$rows
  ## a fit made before the record was kept has the dropped rows and nothing
  ## to count them against
  if (is.null(r))
    r <- list(n_input = NA_integer_, n_used = nrow(object$X),
              n_dropped = if (is.null(object$n_dropped)) NA_integer_ else object$n_dropped,
              dropped_by = integer(0), n_zero_weight = ilm_rows_zero(object$weights))
  structure(r, class = "ilm_rows_used")
}

#' @rdname ilm_rows_used
#' @export
ilm_rows_used.ilm_dag_model <- function(object, ...) {
  if (!length(object$fits))
    return(data.frame(set = integer(0), n_input = integer(0), n_used = integer(0),
                      n_dropped = integer(0), n_zero_weight = integer(0)))
  rows <- lapply(seq_along(object$fits), function(k) {
    f <- object$fits[[k]]
    if (is.null(f)) return(NULL)
    r <- ilm_rows_used(f)
    out <- data.frame(set = k, n_input = r$n_input, n_used = r$n_used,
                      n_dropped = r$n_dropped, n_zero_weight = r$n_zero_weight)
    out$dropped_by <- list(r$dropped_by)
    out
  })
  do.call(rbind, rows)
}

#' @export
print.ilm_rows_used <- function(x, ...) {
  cat(ilm_rows_line(x), "\n", sep = "")
  invisible(x)
}

## The line's words, kept in one place: the wording is a ruling (item 164),
## and "used" may become "analysed" to match a Table 1 and illumex's print.
#' @keywords internal
#' @noRd
ilm_rows_words <- list(
  used = "used",
  dropped = "dropped for missing values",
  several = "a row missing in several columns counts in each",
  zero = "with weight zero, %s but contributing nothing")

## One line for print(), summary() and the record: "200 of 250 rows used; 50
## dropped for missing values (x 30, z 25)". The per-column counts can add up
## to more than the rows dropped, and the line says why when they do. NULL
## when nothing was dropped and no weight is zero, unless `always`.
#' @keywords internal
#' @noRd
ilm_rows_line <- function(r, always = TRUE) {
  w <- ilm_rows_words
  nd <- r$n_dropped; nz <- r$n_zero_weight
  quiet <- (is.na(nd) || nd == 0L) && (is.null(nz) || nz == 0L)
  if (quiet && !always) return(NULL)
  head <- if (is.na(r$n_input)) sprintf("%d rows %s", r$n_used, w$used)
          else sprintf("%d of %d rows %s", r$n_used, r$n_input, w$used)
  drop <- if (!is.na(nd) && nd > 0L) {
    by <- r$dropped_by
    cols <- if (!length(by)) "" else
      paste0(" (", paste(names(by), by, collapse = ", "),
             if (sum(by) > nd) paste0("; ", w$several) else "", ")")
    sprintf("; %d %s%s", nd, w$dropped, cols)
  } else ""
  zero <- if (!is.null(nz) && nz > 0L)
    sprintf(paste0("; %d ", w$zero), nz, w$used) else ""
  paste0(head, drop, zero)
}
