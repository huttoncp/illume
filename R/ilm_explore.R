## ---------------------------------------------------------------------------
## Counts, duplicate detection, cleaning and recoding.
##
## Built on collapse primitives (fcount, group) and base R.
## ---------------------------------------------------------------------------

## ---- counts ----------------------------------------------------------------

#' Frequency counts of a vector's unique values
#'
#' @param y A vector.
#' @param n Number of rows to return, or `"all"`.
#' @param order `"d"` by count descending, `"a"` by count ascending, or `"i"`
#'   by the value itself (factor level order, not label order).
#' @param na.rm Drop missing values before counting.
#' @return A data frame with `value` and `n`.
#' @seealso [ilm_counts_tb()] for both ends at once, [ilm_counts_all()] for a
#'   whole data frame.
#' @examples
#' d <- ilm_sim()
#' ilm_counts(d$grp)
#' ilm_counts(d$site, n = 3)
#' @export
ilm_counts <- function(y, n = "all", order = c("d", "a", "i"), na.rm = TRUE) {
  order <- match.arg(order)
  if (is.data.frame(y))
    stop("`ilm_counts()` summarises a vector; use `ilm_counts_all()` for a data frame",
         call. = FALSE)
  keep <- if (na.rm) !is.na(y) else rep(TRUE, length(y))
  yv <- y[keep]
  tb <- fcount(yv)                       # columns: x, N
  names(tb) <- c("value", "n")
  ## "i" sorts by the value itself, which for a factor means level order, not
  ## alphabetical order of the labels
  idx <- switch(order,
                d = order(-tb$n, tb$value),
                a = order(tb$n, tb$value),
                i = order(tb$value))
  tb <- tb[idx, , drop = FALSE]
  ## `n` was documented as the number of rows to return but had no effect in
  ## the previous implementation; it is honoured here
  if (!identical(n, "all")) {
    n <- as.integer(n)
    if (is.na(n) || n < 1L)
      stop("`n` must be \"all\" or a positive whole number", call. = FALSE)
    tb <- utils::head(tb, n)
  }
  rownames(tb) <- NULL
  tb
}

## top and bottom of the frequency distribution side by side: the two ends are
## where the problems live (a dominant level, and levels too thin to model)
#' The most and least frequent values, side by side
#'
#' The two ends of a frequency distribution are where the problems live: a
#' level that dominates, and levels too thin to estimate.
#'
#' @inheritParams ilm_counts
#' @param n How many values from each end.
#' @return A data frame with `top_value`, `top_n`, `bot_value` and `bot_n`.
#' @examples
#' ilm_counts_tb(ilm_sim()$grp, n = 2)
#' @export
ilm_counts_tb <- function(y, n = 10L, na.rm = TRUE) {
  full <- ilm_counts(y, n = "all", order = "d", na.rm = na.rm)
  k <- min(as.integer(n), nrow(full))
  if (k < 1L) return(data.frame(top_value = character(0), top_n = integer(0),
                                bot_value = character(0), bot_n = integer(0)))
  top <- full[seq_len(k), , drop = FALSE]
  bot <- full[rev(seq(nrow(full) - k + 1L, nrow(full))), , drop = FALSE]
  data.frame(top_value = top$value, top_n = top$n,
             bot_value = bot$value, bot_n = bot$n,
             stringsAsFactors = FALSE, row.names = NULL)
}

## values are stacked across variables of different classes, so they are
## rendered as character; the counts stay integer
#' Frequency counts for every column
#'
#' Values are stacked across columns of different classes, so they are rendered
#' as character; the counts stay integer.
#'
#' @inheritParams ilm_counts
#' @param data A data frame.
#' @param by Optional grouping columns.
#' @return A data frame with `variable`, `value` and `n`.
#' @examples
#' ilm_counts_all(ilm_sim()[, c("grp", "site")], n = 2)
#' @export
ilm_counts_all <- function(data, by = NULL, n = "all", order = c("d", "a", "i"),
                           na.rm = TRUE) {
  order <- match.arg(order)
  miss <- setdiff(by, names(data))
  if (length(miss))
    stop("`by` variable(s) not found in the data: ", paste(miss, collapse = ", "),
         call. = FALSE)
  cols <- setdiff(names(data), by)
  one <- function(d, cols) do.call(rbind, lapply(cols, function(cn) {
    tb <- ilm_counts(d[[cn]], n = n, order = order, na.rm = na.rm)
    data.frame(variable = cn, value = as.character(tb$value), n = tb$n,
               stringsAsFactors = FALSE)
  }))
  if (is.null(by)) return(one(data, cols))
  g <- interaction(data[by], drop = TRUE)
  parts <- lapply(split(seq_len(nrow(data)), g),
                  function(i) one(data[i, , drop = FALSE], cols))
  res <- do.call(rbind, parts)
  cbind(setNames(data.frame(rep(names(parts), vapply(parts, nrow, 1L)),
                            stringsAsFactors = FALSE), paste(by, collapse = ".")),
        res, row.names = NULL)
}

#' Most and least frequent values for every column
#'
#' @inheritParams ilm_counts_tb
#' @param data A data frame.
#' @return A data frame with `variable` and the top/bottom columns.
#' @examples
#' ilm_counts_tb_all(ilm_sim()[, c("grp", "site")], n = 2)
#' @export
ilm_counts_tb_all <- function(data, n = 10L, na.rm = TRUE) {
  do.call(rbind, lapply(names(data), function(cn) {
    tb <- ilm_counts_tb(data[[cn]], n = n, na.rm = na.rm)
    if (!nrow(tb)) return(NULL)
    ## assigned by name rather than with transform(), whose non-standard
    ## evaluation reads as an undefined global to R CMD check
    tb$top_value <- as.character(tb$top_value)
    tb$bot_value <- as.character(tb$bot_value)
    cbind(variable = cn, tb, stringsAsFactors = FALSE, row.names = NULL)
  }))
}

## ---- copies / dupes --------------------------------------------------------

#' Find copied or duplicated rows
#'
#' @param data A data frame.
#' @param ... Columns defining a duplicate. All columns if none are given.
#' @param filter `"all"` appends `copy_number` and `n_copies`; `"dupes"` keeps
#'   only repeated rows; `"first"`, `"last"` and `"unique"` filter rows.
#' @param na_last Sort missing values last.
#' @param sort_by Sort the result by the key columns.
#' @return A data frame; the columns depend on `filter`.
#' @seealso [ilm_dupes()] for the common case.
#' @examples
#' d <- data.frame(a = c(1, 1, 2), b = c("x", "x", "y"))
#' ilm_copies(d)
#' ilm_copies(d, filter = "unique")
#' @export
ilm_copies <- function(data, ..., filter = c("all", "dupes", "first", "last", "unique"),
                       na_last = TRUE, sort_by = TRUE) {
  filter <- match.arg(filter)
  vars <- as.character(unlist(list(...)))
  if (!length(vars)) vars <- names(data)
  miss <- setdiff(vars, names(data))
  if (length(miss))
    stop("column(s) not found in the data: ", paste(miss, collapse = ", "),
         call. = FALSE)

  key <- data[vars]
  g <- group(key)                       # integer group id per row
  sizes <- tabulate(g, attr(g, "N.groups"))
  n_copies <- sizes[g]
  ## position within the group, in row order
  copy_number <- ilm_ave_seq(g)

  out <- switch(filter,
    all    = cbind(data, copy_number = copy_number, n_copies = n_copies),
    dupes  = cbind(data, n_copies = n_copies)[n_copies > 1L, , drop = FALSE],
    first  = data[copy_number == 1L, , drop = FALSE],
    last   = data[copy_number == n_copies, , drop = FALSE],
    unique = data[n_copies == 1L, , drop = FALSE])

  if (sort_by && filter %in% c("all", "dupes")) {
    ord <- do.call(order, c(unname(as.list(out[vars])),
                            list(na.last = na_last)))
    out <- out[ord, , drop = FALSE]
  }
  rownames(out) <- NULL
  out
}

## sequence position within each group, without a split-apply round trip
#' @keywords internal
#' @noRd
ilm_ave_seq <- function(g) {
  o <- order(g)
  r <- integer(length(g))
  gs <- g[o]
  r[o] <- sequence(tabulate(gs, attr(g, "N.groups")))
  r
}

#' Duplicated rows only
#'
#' A shorthand for `ilm_copies(data, ..., filter = "dupes")`.
#'
#' @inheritParams ilm_copies
#' @return A data frame of repeated rows with `n_copies` appended.
#' @examples
#' ilm_dupes(data.frame(a = c(1, 1, 2), b = c("x", "x", "y")))
#' @export
ilm_dupes <- function(data, ..., na_last = TRUE) {
  ilm_copies(data, ..., filter = "dupes", na_last = na_last)
}

## ---- wash_df ---------------------------------------------------------------

## snake_case without depending on janitor: strip non-alphanumerics to
## underscores, split camelCase, collapse repeats, and de-duplicate
#' @keywords internal
#' @noRd
ilm_snake <- function(x) {
  x <- gsub("([a-z0-9])([A-Z])", "\\1_\\2", as.character(x))
  x <- tolower(x)
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  x[!nzchar(x)] <- "v"
  make.unique(x, sep = "_")
}

## A column read as text may still be numeric or logical. Convert only when
## EVERY non-missing value converts cleanly, so a single stray "n/a" cannot
## silently turn a column into NAs.
#' @keywords internal
#' @noRd
ilm_retype <- function(v) {
  if (!is.character(v)) return(v)
  s <- trimws(v)
  s[!nzchar(s)] <- NA_character_
  obs <- s[!is.na(s)]
  if (!length(obs)) return(s)
  lu <- tolower(obs)
  if (all(lu %in% c("true", "false", "t", "f"))) {
    r <- rep(NA, length(s))
    r[!is.na(s)] <- tolower(trimws(s[!is.na(s)])) %in% c("true", "t")
    return(r)
  }
  num <- suppressWarnings(as.numeric(obs))
  if (!anyNA(num)) {
    r <- rep(NA_real_, length(s))
    r[!is.na(s)] <- num
    if (all(abs(num - round(num)) < 1e-8) && all(abs(num) < .Machine$integer.max))
      return(as.integer(r))
    return(r)
  }
  s
}

#' Clean up a messy data frame
#'
#' Removes empty rows and columns, standardises names to snake_case, and
#' converts text columns that are really numeric or logical.
#'
#' Conversion happens only when every non-missing value converts cleanly, so a
#' single stray `"n/a"` cannot silently turn a column into `NA`s.
#'
#' @param data A data frame.
#' @param clean_names Standardise column names to snake_case.
#' @param retype Convert character columns that are really numeric or logical.
#' @param drop_empty Drop rows and columns that are entirely missing or blank.
#' @param column_to_rownames Use a column's values as row names.
#' @param names_col The column to use when `column_to_rownames = TRUE`.
#' @return A cleaned data frame.
#' @examples
#' m <- data.frame("Col One" = c("1", "2", ""), someFlag = c("TRUE", "FALSE", ""),
#'                 check.names = FALSE)
#' ilm_wash_df(m)
#' @export
ilm_wash_df <- function(data, clean_names = TRUE, retype = TRUE,
                        drop_empty = TRUE, column_to_rownames = FALSE,
                        names_col = NULL) {
  if (!is.data.frame(data)) stop("`data` must be a data frame", call. = FALSE)
  d <- as.data.frame(data, stringsAsFactors = FALSE)

  if (column_to_rownames) {
    if (is.null(names_col) || !names_col %in% names(d))
      stop("`names_col` must name a column of `data` when `column_to_rownames = TRUE`",
           call. = FALSE)
    rn <- as.character(d[[names_col]])
    d[[names_col]] <- NULL
  } else rn <- NULL

  blank <- function(v) is.na(v) | (is.character(v) & !nzchar(trimws(as.character(v))))
  if (drop_empty && nrow(d) && ncol(d)) {
    ## a row or column is empty when every cell is missing or whitespace
    keep_col <- !vapply(d, function(v) all(blank(v)), TRUE)
    d <- d[, keep_col, drop = FALSE]
    if (!is.null(rn)) rn <- rn
    if (ncol(d)) {
      keep_row <- !Reduce(`&`, lapply(d, blank))
      d <- d[keep_row, , drop = FALSE]
      if (!is.null(rn)) rn <- rn[keep_row]
    }
  }
  if (retype && ncol(d)) d[] <- lapply(d, ilm_retype)
  if (clean_names && ncol(d)) names(d) <- ilm_snake(names(d))
  if (!is.null(rn)) rownames(d) <- make.unique(rn) else rownames(d) <- NULL
  d
}

## ---- translate -------------------------------------------------------------

## Recode against a dictionary held as two parallel vectors -- the shape a
## lookup table takes when it comes out of a spreadsheet. Values with no entry
## in `old` become NA, which is deliberate: a silent pass-through would hide
## codes the dictionary does not cover.
#' Recode a variable against a dictionary held as two vectors
#'
#' The shape a lookup table takes when it comes out of a spreadsheet: one
#' vector of old codes, one of new, matched by position.
#'
#' Values with no entry in `old` become `NA` by default, which is deliberate: a
#' silent pass-through would hide codes the dictionary does not cover. Use
#' `default` to change that.
#'
#' @param y A vector to recode.
#' @param old Values of the old coding scheme.
#' @param new Replacements, the same length as `old`.
#' @param default Value for entries not found in `old`. `NA` by default.
#' @return A vector of the same length as `y`.
#' @seealso [match()], [ilm_recode_errors()].
#' @examples
#' ilm_translate(c(1, 2, 3, 99), old = 1:3, new = c("low", "mid", "high"))
#' ilm_translate(c(1, 2, 99), old = 1:2, new = c(10, 20), default = -1)
#' @export
ilm_translate <- function(y, old, new, default = NA) {
  if (length(old) != length(new))
    stop("`old` and `new` must be the same length (they are ",
         length(old), " and ", length(new), ")", call. = FALSE)
  if (anyDuplicated(old))
    warning("`old` contains duplicate entries; the first match is used",
            call. = FALSE)
  i <- fmatch(as.character(y), as.character(old))
  out <- new[i]
  if (!identical(default, NA)) out[is.na(i) & !is.na(y)] <- default
  out
}

## ---- recode_errors ---------------------------------------------------------

## Replace known-bad codes (999, -1, "unknown", ...) with NA or another value.
## Restricting to rows/cols matters: a sentinel that means "missing" in one
## column can be a legitimate measurement in another.

#' Replace known-bad values in a vector
#'
#' The vector path of [ilm_recode_errors()], exposed for use inside other
#' pipelines. Factor levels that are recoded away are dropped.
#'
#' @param x A vector.
#' @param errors Values to recode.
#' @param replacement What to put in their place.
#' @return A vector of the same length as `x`.
#' @examples
#' ilm_recode_errors_vec(factor(c("a", "b", "unknown")), errors = "unknown")
#' @export
ilm_recode_errors_vec <- function(x, errors, replacement = NA) {
  if (missing(errors) || !length(errors))
    stop("`errors` must contain at least one value to recode", call. = FALSE)
  hit <- if (is.factor(x)) as.character(x) %in% as.character(errors)
         else x %in% errors
  if (is.factor(x)) {
    lv <- setdiff(levels(x), as.character(errors))
    x <- factor(ifelse(hit, NA_character_, as.character(x)), levels = lv)
    if (!is.na(replacement)) {
      x <- factor(ifelse(hit, as.character(replacement), as.character(x)),
                  levels = unique(c(lv, as.character(replacement))))
    }
    return(x)
  }
  x[hit] <- replacement
  x
}

#' Replace known-bad values with NA or another value
#'
#' Sentinels such as 999, -1 or `"unknown"` mean missing in one column and can
#' be a legitimate measurement in another, which is why the replacement can be
#' restricted to particular rows and columns.
#'
#' @param data A vector, data frame or matrix.
#' @param errors Values to recode.
#' @param replacement What to put in their place. `NA` by default.
#' @param rows,cols Restrict the replacement (data frame or matrix input).
#' @param ind Restrict the replacement (vector input).
#' @return An object of the same shape as `data`.
#' @examples
#' ilm_recode_errors(c(1, 2, 999, -1), errors = c(999, -1))
#' d <- data.frame(x = c(1, 999, 3), y = c(999, 2, 3))
#' ilm_recode_errors(d, errors = 999, cols = "x")
#' @export
ilm_recode_errors <- function(data, errors, replacement = NA,
                              rows = NULL, cols = NULL, ind = NULL) {
  if (missing(errors) || !length(errors))
    stop("`errors` must contain at least one value to recode", call. = FALSE)

  if (is.data.frame(data) || is.matrix(data)) {
    if (!is.null(ind))
      stop("`ind` applies to vector input; use `rows` and `cols` for a ",
           if (is.matrix(data)) "matrix" else "data frame", call. = FALSE)
    nms <- if (is.data.frame(data)) names(data) else colnames(data)
    ci <- if (is.null(cols)) seq_along(nms) else {
      if (is.character(cols)) {
        miss <- setdiff(cols, nms)
        if (length(miss))
          stop("column(s) not found: ", paste(miss, collapse = ", "),
               ". Available: ", paste(nms, collapse = ", "), call. = FALSE)
        match(cols, nms)
      } else {
        bad <- cols[cols < 1 | cols > length(nms)]
        if (length(bad))
          stop("column index out of range: ", paste(bad, collapse = ", "),
               " (data has ", length(nms), " columns)", call. = FALSE)
        as.integer(cols)
      }
    }
    ri <- if (is.null(rows)) seq_len(nrow(data)) else {
      bad <- rows[rows < 1 | rows > nrow(data)]
      if (length(bad))
        stop("row index out of range: ", paste(bad, collapse = ", "),
             " (data has ", nrow(data), " rows)", call. = FALSE)
      as.integer(rows)
    }
    if (is.matrix(data)) {
      sub <- data[ri, ci, drop = FALSE]
      sub[sub %in% errors] <- replacement
      data[ri, ci] <- sub
      return(data)
    }
    for (j in ci) {
      v <- data[[j]]
      v[ri] <- ilm_recode_errors_vec(v[ri], errors, replacement)
      data[[j]] <- v
    }
    return(data)
  }

  if (!is.null(rows) || !is.null(cols))
    stop("`rows` and `cols` apply to data frame or matrix input; use `ind` ",
         "for a vector", call. = FALSE)
  if (is.null(ind)) return(ilm_recode_errors_vec(data, errors, replacement))
  bad <- ind[ind < 1 | ind > length(data)]
  if (length(bad))
    stop("index out of range: ", paste(bad, collapse = ", "),
         " (vector has ", length(data), " elements)", call. = FALSE)
  data[ind] <- ilm_recode_errors_vec(data[ind], errors, replacement)
  data
}
