## ---------------------------------------------------------------------------
## Named EDA plots, one per geometry.
##
## `ilm_plot()` picks a geometry for you and takes `geom =` to override it.
## These are the same drawing done the other way round: you name the plot you
## want. Both sit on tinyplot, which is base graphics with no dependency cost
## beyond what illume already carries, and columns are named as strings
## throughout, the way every other illume function takes them.
##
## No interactivity and no pie chart, both dropped deliberately.
## ---------------------------------------------------------------------------

## Resolve a column argument that may be NULL, returning the vector.
#' @keywords internal
#' @noRd
ilm_col_vec <- function(data, nm, arg, required = FALSE, discrete = FALSE) {
  if (is.null(nm)) {
    if (required)
      stop("`", arg, "` must name a column of `data`", call. = FALSE)
    return(NULL)
  }
  if (!is.character(nm) || length(nm) != 1L)
    stop("`", arg, "` must be a single column name, as a string", call. = FALSE)
  if (!nm %in% names(data))
    stop("column not found in the data: ", nm, ". Available: ",
         paste(utils::head(names(data), 12), collapse = ", "), call. = FALSE)
  v <- data[[nm]]
  ## A grouping column for a discrete geometry is a factor whatever it is
  ## stored as. Handing tinyplot a numeric one makes it try a continuous
  ## legend, fail, and warn on every call.
  if (discrete && !is.null(v) && !is.factor(v)) v <- factor(v)
  v
}

#' @keywords internal
#' @noRd
ilm_plot_frame_check <- function(data) {
  if (!is.data.frame(data))
    stop("`data` must be a data frame; it is ", class(data)[1], call. = FALSE)
}

#' Histogram
#'
#' @param data A data frame.
#' @param x Name of the numeric column to plot.
#' @param by Optional grouping column, overlaying one histogram per level.
#' @param breaks Passed to [tinyplot::type_histogram()].
#' @param ... Passed to [tinyplot::tinyplot()].
#' @return `NULL`, invisibly.
#' @seealso [ilm_plot()], which picks a geometry for you.
#' @examples
#' ilm_plot_histogram(mtcars, "mpg")
#' ilm_plot_histogram(mtcars, "mpg", by = "cyl")
#' @export
ilm_plot_histogram <- function(data, x, by = NULL, breaks = "Sturges", ...) {
  ilm_plot_frame_check(data)
  tinyplot::tinyplot(x = ilm_col_vec(data, x, "x", TRUE),
                     by = ilm_col_vec(data, by, "by", discrete = TRUE),
                     type = tinyplot::type_histogram(breaks = breaks),
                     xlab = x, ...)
  invisible(NULL)
}

#' Density plot
#'
#' @inheritParams ilm_plot_histogram
#' @return `NULL`, invisibly.
#' @seealso [ilm_plot_histogram()].
#' @examples
#' ilm_plot_density(mtcars, "mpg", by = "cyl")
#' @export
ilm_plot_density <- function(data, x, by = NULL, ...) {
  ilm_plot_frame_check(data)
  tinyplot::tinyplot(x = ilm_col_vec(data, x, "x", TRUE),
                     by = ilm_col_vec(data, by, "by", discrete = TRUE), type = "density",
                     xlab = x, ...)
  invisible(NULL)
}

#' Boxplot
#'
#' The whiskers use the same 1.5 interquartile-range fence that
#' [ilm_outliers()] scores against, so the points beyond them are the values
#' that function flags under `method = "iqr"`.
#'
#' @param data A data frame.
#' @param y Name of the numeric column.
#' @param x Optional categorical column to split by.
#' @param by Optional grouping column.
#' @param ... Passed to [tinyplot::tinyplot()].
#' @param pch Plotting character. Takes a NAME as well as a number:
#'   `"filled circle"` is 16, and every code from 0 to 25 has one. Case,
#'   spaces, underscores and hyphens are ignored. A single character is
#'   drawn literally, so `pch = "x"` is still the letter x.
#' @return `NULL`, invisibly.
#' @seealso [ilm_outliers()], [ilm_plot_violin()].
#' @examples
#' ilm_plot_box(mtcars, "mpg", x = "cyl")
#' @export
ilm_plot_box <- function(data, y, x = NULL, by = NULL, ..., pch = NULL) {
  ilm_plot_frame_check(data)
  xv <- ilm_col_vec(data, x, "x")
  tinyplot::tinyplot(x = if (is.null(xv)) "" else xv,
                     y = ilm_col_vec(data, y, "y", TRUE),
                     by = ilm_col_vec(data, by, "by", discrete = TRUE),
                     type = "boxplot", ylab = y, pch = ilm_pch(pch), ...)
  invisible(NULL)
}

#' Violin plot
#'
#' @inheritParams ilm_plot_box
#' @return `NULL`, invisibly.
#' @seealso [ilm_plot_box()], which shows the quartiles rather than the shape.
#' @examples
#' ilm_plot_violin(mtcars, "mpg", x = "cyl")
#' @export
ilm_plot_violin <- function(data, y, x = NULL, by = NULL, ..., pch = NULL) {
  ilm_plot_frame_check(data)
  xv <- ilm_col_vec(data, x, "x")
  tinyplot::tinyplot(x = if (is.null(xv)) "" else xv,
                     y = ilm_col_vec(data, y, "y", TRUE),
                     by = ilm_col_vec(data, by, "by", discrete = TRUE),
                     type = "violin", ylab = y, pch = ilm_pch(pch), ...)
  invisible(NULL)
}

#' Scatter plot
#'
#' @param data A data frame.
#' @param y,x Names of the numeric columns.
#' @param by Optional grouping column.
#' @param trend `"none"`, `"lm"` or `"loess"`. A trend line here is a
#'   description of the two columns shown and nothing more -- it holds nothing
#'   else fixed, so it is not the effect [ilm_model()] would estimate.
#' @param ... Passed to [tinyplot::tinyplot()].
#' @param pch Plotting character. Takes a NAME as well as a number:
#'   `"filled circle"` is 16, and every code from 0 to 25 has one. Case,
#'   spaces, underscores and hyphens are ignored. A single character is
#'   drawn literally, so `pch = "x"` is still the letter x.
#' @return `NULL`, invisibly.
#' @seealso [ilm_plot_var_pairs()] for every pair at once.
#' @examples
#' ilm_plot_scatter(mtcars, "mpg", "wt", by = "cyl", trend = "lm")
#' @export
ilm_plot_scatter <- function(data, y, x, by = NULL,
                             trend = c("none", "lm", "loess"), ...,
                             pch = NULL) {
  trend <- match.arg(trend)
  ilm_plot_frame_check(data)
  tinyplot::tinyplot(x = ilm_col_vec(data, x, "x", TRUE),
                     y = ilm_col_vec(data, y, "y", TRUE),
                     by = ilm_col_vec(data, by, "by"),
                     type = if (trend == "none") "points" else trend,
                     xlab = x, ylab = y, pch = ilm_pch(pch), ...)
  invisible(NULL)
}

#' Bar plot
#'
#' @param data A data frame.
#' @param x Name of the categorical column.
#' @param by Optional grouping column.
#' @param ... Passed to [tinyplot::tinyplot()].
#' @return `NULL`, invisibly.
#' @seealso [ilm_counts()] for the same information as a table.
#' @examples
#' ilm_plot_bar(mtcars, "cyl", by = "am")
#' @export
ilm_plot_bar <- function(data, x, by = NULL, ...) {
  ilm_plot_frame_check(data)
  tinyplot::tinyplot(x = ilm_col_vec(data, x, "x", TRUE),
                     by = ilm_col_vec(data, by, "by", discrete = TRUE), type = "barplot",
                     xlab = x, ...)
  invisible(NULL)
}

#' Line plot
#'
#' @param data A data frame.
#' @param y,x Column names; `x` is usually a date or a sequence.
#' @param by Optional grouping column.
#' @param ... Passed to [tinyplot::tinyplot()].
#' @param pch Plotting character. Takes a NAME as well as a number:
#'   `"filled circle"` is 16, and every code from 0 to 25 has one. Case,
#'   spaces, underscores and hyphens are ignored. A single character is
#'   drawn literally, so `pch = "x"` is still the letter x.
#' @return `NULL`, invisibly.
#' @seealso [ilm_plot_acf()] for what a line plot of residuals cannot show.
#' @examples
#' d <- data.frame(t = 1:40, v = cumsum(rnorm(40)))
#' ilm_plot_line(d, "v", "t")
#' @export
ilm_plot_line <- function(data, y, x, by = NULL, ..., pch = NULL) {
  ilm_plot_frame_check(data)
  tinyplot::tinyplot(x = ilm_col_vec(data, x, "x", TRUE),
                     y = ilm_col_vec(data, y, "y", TRUE),
                     by = ilm_col_vec(data, by, "by"), type = "lines",
                     xlab = x, ylab = y, pch = ilm_pch(pch), ...)
  invisible(NULL)
}

#' Group means or medians with an error bar
#'
#' A summary per group with a measure of spread around it.
#'
#' The two `stat` options show different things and are not interchangeable.
#' `"mean"` draws the standard error, which is about where the *mean* is;
#' `"median"` draws the quartiles, which is about where the *data* are. The
#' first shrinks as the sample grows and the second does not.
#'
#' Overlapping error bars are a poor test of a difference -- they are
#' conservative and lossy. [ilm_boot_diff()] gives the difference itself with
#' its own interval.
#'
#' @param data A data frame.
#' @param y Name of the numeric column to summarise.
#' @param x Name of the column defining the groups.
#' @param by Optional secondary grouping column.
#' @param stat `"mean"` with its standard error, or `"median"` with the
#'   quartiles.
#' @param ... Passed to [tinyplot::tinyplot()].
#' @param pch Plotting character. Takes a NAME as well as a number:
#'   `"filled circle"` is 16, and every code from 0 to 25 has one. Case,
#'   spaces, underscores and hyphens are ignored. A single character is
#'   drawn literally, so `pch = "x"` is still the letter x.
#' @return `NULL`, invisibly.
#' @seealso [ilm_boot_diff()] for the comparison this plot invites.
#' @examples
#' ilm_plot_stat_error(mtcars, "mpg", "cyl")
#' @export
ilm_plot_stat_error <- function(data, y, x, by = NULL,
                                stat = c("mean", "median"), ...,
                                pch = NULL) {
  stat <- match.arg(stat)
  ilm_plot_frame_check(data)
  yv <- ilm_col_vec(data, y, "y", TRUE)
  xv <- ilm_col_vec(data, x, "x", TRUE)
  bv <- ilm_col_vec(data, by, "by", discrete = TRUE)
  key <- if (is.null(bv)) data.frame(x = xv, stringsAsFactors = FALSE)
         else data.frame(x = xv, by = bv, stringsAsFactors = FALSE)
  sp <- split(seq_along(yv), interaction(key, drop = TRUE, sep = "\r"))
  rows <- lapply(sp, function(i) {
    v <- yv[i]; v <- v[!is.na(v)]
    if (!length(v)) return(NULL)
    ct <- if (stat == "mean") {
      se <- stats::sd(v) / sqrt(length(v))
      c(mean(v), mean(v) - se, mean(v) + se)
    } else unname(stats::quantile(v, c(0.5, 0.25, 0.75)))
    data.frame(x = key$x[i[1]],
               by = if (is.null(bv)) NA else key$by[i[1]],
               center = ct[1], lo = ct[2], hi = ct[3],
               stringsAsFactors = FALSE)
  })
  s <- do.call(rbind, rows)
  if (is.null(s) || !nrow(s))
    stop("no non-missing values to summarise", call. = FALSE)
  tinyplot::tinyplot(x = s$x, y = s$center, ymin = s$lo, ymax = s$hi,
                     by = if (is.null(bv)) NULL else s$by, type = "pointrange",
                     xlab = x, ylab = paste(stat, y), pch = ilm_pch(pch), ...)
  invisible(NULL)
}
