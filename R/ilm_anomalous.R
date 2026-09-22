## ---------------------------------------------------------------------------
## From "these rows are odd" to "and here is what they are".
##
## ilm_anomaly() ends where the interesting part starts. It says which rows do
## not fit the pattern; the next question is always whether those rows are
## alike -- one data-entry problem repeated, or a genuine subpopulation, or
## several unrelated oddities. Answering it means taking the flagged rows back
## to the exploration tools, and until now that meant the user writing the
## subset by hand from `row`, which is both a papercut and an invitation to
## line the wrong rows up.
##
## So the anomaly object carries the data it scored, and the exploration
## functions take the object directly. ilm_anomalous() is the extractor those
## use, exported because it composes with everything else too -- anything that
## takes a data frame takes its result.
## ---------------------------------------------------------------------------

#' The rows an anomaly scan flagged
#'
#' Pulls the flagged rows out of an [ilm_anomaly()] result, as a data frame,
#' so they can be explored with anything that takes one.
#'
#' The point of a scan is rarely the score. It is whether the odd rows are
#' alike: one mistake made repeatedly, a subpopulation the model does not
#' cover, or several unrelated things. That is a clustering or a description
#' question, and this is what hands the rows over.
#'
#' @param x An [ilm_anomaly()] result.
#' @param data The data frame that was scanned. Only needed if the result did
#'   not keep it -- see `keep_data` in [ilm_anomaly()].
#' @param flagged Return only the flagged rows (the default), or all of them
#'   with the score attached.
#' @param score Add the `score`, `p_adj` and `driver` columns.
#' @return A data frame of the flagged rows, in their original order, with
#'   `.row` giving their position in the input.
#' @seealso [ilm_anomaly()], [ilm_profile()], [ilm_cluster()],
#'   [ilm_describe_all()].
#' @examples
#' \donttest{
#' a <- ilm_anomaly(mtcars)
#' odd <- ilm_anomalous(a)
#' ilm_describe_all(odd)
#' }
#' @export
ilm_anomalous <- function(x, data = NULL, flagged = TRUE, score = TRUE) {
  if (!inherits(x, "ilm_anomaly"))
    stop("`x` must be an ilm_anomaly() result, not ", class(x)[1],
         call. = FALSE)
  d <- as.data.frame(x); class(d) <- "data.frame"
  dat <- if (!is.null(data)) data else attr(x, "data", exact = TRUE)
  if (is.null(dat))
    stop("this result did not keep the data it scored, so the rows cannot be ",
         "recovered. Pass `data =` -- the same frame ilm_anomaly() was given ",
         "-- or re-run it with keep_data = TRUE.", call. = FALSE)
  if (!is.data.frame(dat)) dat <- as.data.frame(dat)
  keep <- if (flagged) d$row[d$flag] else d$row
  if (!length(keep)) {
    ## an empty result is a legitimate answer, not a failure: it means nothing
    ## was flagged. It still has to carry every column the non-empty case
    ## would, or code that selects `.driver` breaks on exactly the data sets
    ## that behaved themselves.
    out <- dat[0, , drop = FALSE]
    out$.row <- integer(0)
    if (score) {
      out$.score <- numeric(0); out$.p_adj <- numeric(0)
      out$.driver <- character(0)
    }
    row.names(out) <- NULL
    return(out)
  }
  if (max(keep) > nrow(dat))
    stop("the result refers to row ", max(keep), " but `data` has ",
         nrow(dat), " rows. Pass the frame that was scanned.", call. = FALSE)
  ord <- order(keep)
  out <- dat[keep[ord], , drop = FALSE]
  out$.row <- keep[ord]
  if (score) {
    i <- match(keep[ord], d$row)
    out$.score <- d$score[i]
    out$.p_adj <- d$p_adj[i]
    out$.driver <- d$driver[i]
  }
  row.names(out) <- NULL
  out
}

## Used by ilm_profile(), ilm_cluster(), ilm_reduce() and the describe
## functions so that each of them takes an anomaly result directly.
##
## Whether the score columns come along depends on what happens next. For
## anything that measures distance -- reducing, clustering, profiling -- they
## are left out: they are properties of the scan rather than measurements on
## the rows, and a clustering built on them would recover the scan instead of
## the data. For describing they are kept, because `.driver` says which
## variable pushed each row over the line and a frequency table of it is
## often the whole answer.
#' @keywords internal
#' @noRd
ilm_from_anomaly <- function(x, fn, score = FALSE) {
  d <- ilm_anomalous(x, score = score)
  n <- nrow(d)
  if (!n)
    stop("no rows were flagged, so there is nothing for ", fn, "() to work ",
         "on. ilm_anomalous(x, flagged = FALSE) returns every row with its ",
         "score if you want to look at the whole scan.", call. = FALSE)
  if (n < 5L && !score)
    warning(fn, "(): only ", n, " row(s) were flagged. Describing them is ",
            "fine; clustering that few is not -- there is no structure to ",
            "find in a handful of points.", call. = FALSE)
  message(fn, "(): using the ", n, " row(s) ilm_anomaly() flagged, out of ",
          nrow(attr(x, "data", exact = TRUE) %||% d), ".")
  d$.row <- NULL
  d
}
