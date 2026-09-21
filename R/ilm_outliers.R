## ---------------------------------------------------------------------------
## Flagging unusual values.
##
## A boxplot shows outliers through its whiskers and leaves you squinting at
## them. What is wanted alongside that is the programmatic question: which
## values are flagged, in which rows, and by how far. All three rules are base R.
##
## None of them decides whether a value is WRONG. An outlier is a value far from
## the others by some rule, which may mean a typo, a different population, or a
## genuinely heavy tail that the model should be handling rather than the
## cleaning step. Deleting flagged rows is almost never the right response, and
## the documentation says so rather than leaving it implied.
## ---------------------------------------------------------------------------

#' Flag unusual values in a numeric vector
#'
#' Scores every value by how far it sits from the centre, in the units of
#' whichever rule is chosen, and flags those past a threshold.
#'
#' @section The three rules:
#'
#' * `"iqr"` -- distance beyond the nearer quartile in interquartile ranges.
#'   The classic Tukey fence, and the same rule a boxplot's whiskers draw.
#' * `"mad"` -- distance from the median in median absolute deviations. More
#'   robust than `"zscore"` on skewed or heavy-tailed data, because neither the
#'   median nor the MAD is itself dragged around by the values being detected,
#'   where the mean and standard deviation both are.
#' * `"zscore"` -- distance from the mean in standard deviations. Included
#'   because it is expected, but it is the weakest of the three for exactly
#'   that reason: a large outlier inflates the standard deviation and so hides
#'   itself.
#'
#' @section What a flag is not:
#'
#' It is not a verdict that the value is wrong. A flagged value may be a
#' recording error, a member of a different population, or an ordinary draw
#' from a heavy tail -- and the third is common. Dropping flagged rows because
#' they are flagged changes the estimand and biases whatever is fitted next; if
#' the tail is real, the remedy is a model that expects it, which for illume
#' means a different family or [ilm_model(dispformula = )][ilm_model].
#'
#' @param y A numeric vector.
#' @param method `"iqr"`, `"mad"` or `"zscore"`.
#' @param threshold Flagging threshold. Defaults per rule: 1.5 for `"iqr"`
#'   (the boxplot convention), 3.5 for `"mad"` (Iglewicz and Hoaglin), 3 for
#'   `"zscore"`.
#' @param na.rm Compute the reference statistics with missing values removed.
#'   The result always has one row per element of `y`, with `NA` for `score`
#'   and `is_outlier` wherever `y` is `NA`.
#' @return A data frame with `value`, `method`, `threshold`, `score` (distance
#'   from the centre in the rule's own units) and `is_outlier`.
#' @references
#' Iglewicz, B. and Hoaglin, D. C. (1993). How to Detect and Handle Outliers.
#' ASQC Quality Press.
#' @seealso [ilm_outliers_all()] for a whole data frame, [ilm_describe()] for
#'   the distribution a flag should be read against.
#' @examples
#' table(ilm_outliers(mtcars$hp)$is_outlier)
#' head(ilm_outliers(mtcars$hp, method = "mad"))
#' @export
ilm_outliers <- function(y, method = c("iqr", "mad", "zscore"),
                         threshold = NULL, na.rm = TRUE) {
  method <- match.arg(method)
  if (!is.numeric(y))
    stop("`y` must be numeric; it is ", class(y)[1], call. = FALSE)
  if (is.null(threshold))
    threshold <- switch(method, iqr = 1.5, mad = 3.5, zscore = 3)

  ref <- if (na.rm) y[!is.na(y)] else y
  if (length(ref) < 2L)
    stop("at least 2 non-missing values are needed to say what is unusual",
         call. = FALSE)

  score <- if (method == "iqr") {
    ## Tukey's HINGES, via fivenum(), not quantile(type = 7). The two differ on
    ## small samples -- on mtcars$wt the fences land at 5.311 and 5.153 -- and
    ## the boxplot whisker, which is what this rule claims to reproduce, is
    ## drawn from the hinges. With type-7 quantiles a value of 5.25 was flagged
    ## here while sitting inside the whisker, so the documented equivalence to
    ## ilm_plot_box() was false. Checked against grDevices::boxplot.stats().
    q <- stats::fivenum(ref)[c(2L, 4L)]
    iqr <- q[2] - q[1]
    ## Distance beyond the nearer quartile in IQR units, with NO threshold
    ## baked in: `threshold` is applied exactly once, below. An earlier version
    ## subtracted threshold * iqr into the fence AND THEN compared the score
    ## against the threshold, silently doubling the fence -- caught by a test
    ## that expected a lower threshold to flag more values and saw it flag
    ## fewer.
    if (iqr == 0) rep(0, length(y))
    else pmax((q[1] - y) / iqr, (y - q[2]) / iqr, 0)
  } else if (method == "mad") {
    md <- stats::mad(ref)
    if (md == 0) rep(0, length(y)) else abs(y - stats::median(ref)) / md
  } else {
    s <- stats::sd(ref)
    if (s == 0) rep(0, length(y)) else abs(y - mean(ref)) / s
  }

  is_outlier <- score > threshold
  is_outlier[is.na(y)] <- NA
  data.frame(value = y, method = method, threshold = threshold,
             score = round(score, 3), is_outlier = is_outlier,
             stringsAsFactors = FALSE, row.names = NULL)
}

#' Flag unusual values across a data frame
#'
#' Runs [ilm_outliers()] over every numeric column, optionally within groups,
#' and returns one row per flagged value with the row it came from.
#'
#' Grouping matters more than it looks. A value can be perfectly ordinary for
#' its own group and extreme against the pooled distribution, so flagging
#' without `by` on data that has groups mostly rediscovers the groups.
#'
#' @param data A data frame.
#' @param by Grouping column(s), as a character vector. Reference statistics,
#'   and so the flags, are computed separately within each group.
#' @param cols Columns to use. A character vector of names, a
#'   regular expression, a predicate function such as `is.numeric`, or
#'   `NULL` for all of them -- see [ilm_selection].
#' @param flagged_only Return only the flagged values. `FALSE` returns every
#'   value's score.
#' @inheritParams ilm_outliers
#' @return A data frame with the `by` columns, `row_id` (the row's position in
#'   `data`), `variable`, `value`, `score` and `is_outlier`.
#' @seealso [ilm_outliers()], [ilm_plot_box()] to see them.
#' @examples
#' head(ilm_outliers_all(mtcars))
#' head(ilm_outliers_all(mtcars, by = "cyl", method = "mad"))
#' @export
ilm_outliers_all <- function(data, by = NULL, cols = NULL,
                             method = c("iqr", "mad", "zscore"),
                             threshold = NULL, flagged_only = TRUE,
                             na.rm = TRUE) {
  method <- match.arg(method)
  if (!is.data.frame(data))
    stop("`data` must be a data frame; it is ", class(data)[1], call. = FALSE)
  g <- ilm_resolve_cols(data, by, arg = "by")
  if (is.null(by)) g <- character()
  num <- ilm_resolve_cols(data, cols, exclude = g, arg = "cols",
                          eligible = names(data)[vapply(data, is.numeric, TRUE)])
  if (!length(num))
    stop("no numeric columns to check for unusual values", call. = FALSE)

  ## one grouping key, or a single group when none was asked for
  key <- if (!length(g)) rep("", nrow(data))
         else interaction(data[g], drop = TRUE, sep = "\r")
  ## Accumulate into plain vectors and build ONE frame at the end. Binding a
  ## frame per (column, group) means 28 columns x 8 groups = 224 rbind calls,
  ## each copying everything accumulated so far; that alone was most of the
  ## 1.58s a grouped run took on 20,000 rows.
  ri <- vi <- va <- sc <- io <- vector("list", length(num) * length(unique(key)))
  k <- 0L
  for (v in num) {
    for (lv in unique(key)) {
      idx <- which(key == lv)
      if (length(idx) < 2L) next
      o <- tryCatch(ilm_outliers(data[[v]][idx], method = method,
                                 threshold = threshold, na.rm = na.rm),
                    error = function(e) NULL)
      if (is.null(o)) next
      k <- k + 1L
      ri[[k]] <- idx; vi[[k]] <- rep(v, length(idx))
      va[[k]] <- o$value; sc[[k]] <- o$score; io[[k]] <- o$is_outlier
    }
  }
  out <- if (!k) NULL else {
    rid <- unlist(ri[seq_len(k)], use.names = FALSE)
    o2 <- data.frame(row_id = rid,
                     variable = unlist(vi[seq_len(k)], use.names = FALSE),
                     value = unlist(va[seq_len(k)], use.names = FALSE),
                     score = unlist(sc[seq_len(k)], use.names = FALSE),
                     is_outlier = unlist(io[seq_len(k)], use.names = FALSE),
                     stringsAsFactors = FALSE)
    if (length(g)) cbind(data[rid, g, drop = FALSE], o2) else o2
  }
  if (is.null(out)) {
    out <- data.frame(row_id = integer(), variable = character(),
                      value = numeric(), score = numeric(),
                      is_outlier = logical(), stringsAsFactors = FALSE)
    if (length(g)) out <- cbind(data[0, g, drop = FALSE], out)
    return(out)
  }
  if (flagged_only) out <- out[!is.na(out$is_outlier) & out$is_outlier, ,
                               drop = FALSE]
  out$variable <- factor(out$variable, levels = num)
  ord <- do.call(order, c(lapply(g, function(k) out[[k]]),
                          list(out$variable, out$row_id)))
  out <- out[ord, , drop = FALSE]
  out$variable <- as.character(out$variable)
  rownames(out) <- NULL
  out[, c(g, "row_id", "variable", "value", "score", "is_outlier"),
      drop = FALSE]
}
