## ---------------------------------------------------------------------------
## Choosing columns without non-standard evaluation.
##
## Naming forty of two hundred columns as a character vector is the one place
## illume's string-based API genuinely hurts. The tidyselect answer to that is
## data masking, which buys the convenience at the cost of an rlang dependency,
## of making the package harder to drive FROM code -- illume's own study scripts
## loop over column names -- and of silent resolution when a bare name also
## exists in the calling environment.
##
## None of that is necessary. A `cols` argument can take a character vector, a
## regular expression, or a predicate function, and all three are ordinary
## values: assignable to a variable, passable through `...`, and printable. The
## expressive part of tidyselect is the SELECTION, not the masking.
## ---------------------------------------------------------------------------

#' How columns can be chosen
#'
#' Wherever illume takes a `cols` or `by` argument, it accepts any of four
#' things. All are ordinary values, so any of them can be held in a variable and
#' passed along -- which a bare-name interface cannot.
#'
#' * `NULL`, meaning every eligible column.
#' * A **character vector** of names: `c("age", "income")`. Every name must
#'   exist, and a typo is an error naming the columns that do exist.
#' * A **regular expression**, as a single string that is not itself a column
#'   name: `"^score_"` takes every column whose name starts with `score_`.
#'   Matching nothing is an error rather than a silent empty selection.
#' * A **predicate function** applied to each column: `is.numeric`, or
#'   `function(v) is.numeric(v) && !anyNA(v)`.
#'
#' The one ambiguity is a single string, which could be a name or a pattern.
#' A name wins: `cols = "score"` selects the column called `score` even if
#' `scores_2024` also exists. Anything that is not a column name is read as a
#' pattern.
#'
#' @name ilm_selection
#' @examples
#' d <- data.frame(id = 1:5, score_a = rnorm(5), score_b = rnorm(5),
#'                 label = letters[1:5])
#' ilm_outliers_all(d, cols = "^score_")
#' ilm_outliers_all(d, cols = is.numeric)
NULL

#' @keywords internal
#' @noRd
ilm_resolve_cols <- function(data, cols, exclude = character(),
                             arg = "cols", eligible = NULL) {
  nms <- setdiff(names(data), exclude)
  if (!is.null(eligible)) nms <- intersect(nms, eligible)
  if (is.null(cols)) return(nms)

  if (is.function(cols)) {
    keep <- vapply(nms, function(v) {
      ok <- tryCatch(cols(data[[v]]), error = function(e) FALSE)
      isTRUE(ok)
    }, TRUE)
    out <- nms[keep]
    if (!length(out))
      stop("`", arg, "` is a function that matched no column", call. = FALSE)
    return(out)
  }

  if (!is.character(cols))
    stop("`", arg, "` must be column names, a regular expression, or a ",
         "function; it is ", class(cols)[1], call. = FALSE)

  ## every element a real column name: take them as names
  known <- cols %in% names(data)
  if (all(known)) {
    out <- intersect(cols, nms)
    if (!length(out))
      stop("`", arg, "` named only columns that are excluded here: ",
           paste(cols, collapse = ", "), call. = FALSE)
    return(out)
  }
  ## a single string that is not a column name: read it as a pattern
  if (length(cols) == 1L) {
    out <- tryCatch(grep(cols, nms, value = TRUE),
                    error = function(e)
                      stop("`", arg, "` is neither a column name nor a valid ",
                           "regular expression: ", cols, call. = FALSE))
    ## A single unknown string is ambiguous: a typo'd name and a pattern that
    ## matches nothing look identical from here, so the message covers both
    ## readings rather than guessing which the user meant.
    if (!length(out))
      stop("column not found in the data, and matched nothing as a pattern ",
           "either: ", cols, ". Available: ",
           paste(utils::head(nms, 12), collapse = ", "), call. = FALSE)
    return(out)
  }
  stop("column(s) not found in the data: ",
       paste(cols[!known], collapse = ", "), ". Available: ",
       paste(utils::head(names(data), 12), collapse = ", "), call. = FALSE)
}
