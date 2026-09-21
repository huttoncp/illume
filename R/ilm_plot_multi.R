## ---------------------------------------------------------------------------
## Missingness plots, auto-dispatching plots, pairs, and multi-panel figures.
##
## ilm_plot_c() is a different shape from a patchwork-style combinator, and the
## reason is base graphics rather than a gap here: base draws immediately and
## keeps no retained plot object, so there is nothing to combine after the fact.
## It therefore takes plot-producing EXPRESSIONS and evaluates each into a panel
## of a par(mfrow) grid.
## ---------------------------------------------------------------------------

## Counts of missing and present per column, optionally within groups.
#' @keywords internal
#' @noRd
ilm_na_summary <- function(data, cols, g = character()) {
  key <- if (!length(g)) factor(rep("", nrow(data)))
         else interaction(data[g], drop = TRUE, sep = " / ")
  rows <- list()
  for (lv in levels(key)) {
    idx <- which(key == lv)
    for (v in cols) {
      z <- data[[v]][idx]
      rows[[length(rows) + 1L]] <- data.frame(
        group = lv, variable = v, n = sum(!is.na(z)), na = sum(is.na(z)),
        p_na = if (length(z)) mean(is.na(z)) else NA_real_,
        stringsAsFactors = FALSE)
    }
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' Missing values by column
#'
#' A bar per column. The first question to ask of an unfamiliar data set, and
#' the one [ilm_check_missing()] then turns into advice.
#'
#' @param data A data frame.
#' @param by Optional grouping column(s), as a character vector -- bars are
#'   drawn per group, which is how you see whether missingness is concentrated
#'   somewhere.
#' @param stat `"p_na"` (proportion missing), `"na"` (count missing) or `"n"`
#'   (count present).
#' @param ... Passed to [tinyplot::tinyplot()].
#' @return `NULL`, invisibly.
#' @seealso [ilm_check_missing()] for whether it matters, [ilm_profile_na()]
#'   for which columns go missing together.
#' @examples
#' ilm_plot_na_all(airquality)
#' ilm_plot_na_all(airquality, by = "Month")
#' @export
ilm_plot_na_all <- function(data, by = NULL, stat = c("p_na", "na", "n"), ...) {
  stat <- match.arg(stat)
  ilm_plot_frame_check(data)
  g <- if (is.null(by)) character() else as.character(by)
  miss <- setdiff(g, names(data))
  if (length(miss))
    stop("column(s) not found in the data: ", paste(miss, collapse = ", "),
         call. = FALSE)
  cols <- setdiff(names(data), g)
  if (!length(cols))
    stop("no columns left to plot once `by` is set aside", call. = FALSE)
  s <- ilm_na_summary(data, cols, g)
  tinyplot::tinyplot(x = s$variable, y = s[[stat]],
                     by = if (length(g)) s$group else NULL,
                     type = if (length(g)) tinyplot::type_barplot(beside = TRUE)
                            else "barplot",
                     xlab = "variable", ylab = stat, ...)
  invisible(NULL)
}

#' Missing values in one column, across groups
#'
#' Where [ilm_plot_na_all()] compares columns, this compares groups within one
#' column -- which is how a pattern in *who* is missing shows itself.
#'
#' @param data A data frame.
#' @param x Name of the column whose missingness to show.
#' @param by Grouping column(s) to split by, as a character vector. Required:
#'   without it there is only one bar, which is [ilm_plot_na_all()]'s job.
#' @param stat `"p_na"`, `"na"` or `"n"`.
#' @param ... Passed to [tinyplot::tinyplot()].
#' @return `NULL`, invisibly.
#' @seealso [ilm_plot_na_all()], [ilm_check_missing()].
#' @examples
#' ilm_plot_na(airquality, "Ozone", by = "Month")
#' @export
ilm_plot_na <- function(data, x, by = NULL, stat = c("p_na", "na", "n"), ...) {
  stat <- match.arg(stat)
  ilm_plot_frame_check(data)
  if (is.null(by))
    stop("`by` is required: without a grouping variable there is one bar. ",
         "Use ilm_plot_na_all() to compare columns instead.", call. = FALSE)
  g <- as.character(by)
  miss <- setdiff(c(x, g), names(data))
  if (length(miss))
    stop("column(s) not found in the data: ", paste(miss, collapse = ", "),
         call. = FALSE)
  s <- ilm_na_summary(data, x, g)
  tinyplot::tinyplot(x = s$group, y = s[[stat]], type = "barplot",
                     xlab = paste(g, collapse = " / "),
                     ylab = paste0(stat, " (", x, ")"), ...)
  invisible(NULL)
}

#' @keywords internal
#' @noRd
ilm_classify_var <- function(v)
  if (is.numeric(v) || inherits(v, "Date") || inherits(v, "POSIXt"))
    "continuous" else "categorical"

#' Plot one or two variables, choosing the geometry
#'
#' Picks from the columns' types: a continuous variable alone gets a density, a
#' categorical one a bar chart, two continuous ones a scatter, one of each a
#' boxplot, two categorical ones a grouped bar chart.
#'
#' [ilm_plot()] does the same job through a `geom = "auto"` argument. This is
#' the same choice made by naming variables rather than a geometry, and it says
#' which it picked when asked.
#'
#' @param data A data frame.
#' @param var1 Name of the primary column.
#' @param var2 Optional name of a second column.
#' @param by Optional grouping column.
#' @param verbose Say which plot was chosen and why.
#' @param ... Passed to the underlying plot function.
#' @return `NULL`, invisibly.
#' @seealso [ilm_plot()], [ilm_plot_var_all()].
#' @examples
#' ilm_plot_var(mtcars, "mpg")
#' ilm_plot_var(mtcars, "mpg", "cyl", verbose = TRUE)
#' @export
ilm_plot_var <- function(data, var1, var2 = NULL, by = NULL, verbose = FALSE,
                         ...) {
  ilm_plot_frame_check(data)
  miss <- setdiff(c(var1, var2, by), names(data))
  if (length(miss))
    stop("column(s) not found in the data: ", paste(miss, collapse = ", "),
         call. = FALSE)
  say <- function(...) if (verbose) message(...)
  c1 <- ilm_classify_var(data[[var1]])
  if (is.null(var2)) {
    if (c1 == "continuous") {
      say("`", var1, "` is continuous -- density")
      ilm_plot_density(data, var1, by = by, ...)
    } else {
      say("`", var1, "` is categorical -- bar")
      ilm_plot_bar(data, var1, by = by, ...)
    }
    return(invisible(NULL))
  }
  c2 <- ilm_classify_var(data[[var2]])
  if (c1 == "continuous" && c2 == "continuous") {
    say("both continuous -- scatter")
    ilm_plot_scatter(data, var1, var2, by = by, ...)
  } else if (c1 != c2) {
    num <- if (c1 == "continuous") var1 else var2
    cat_ <- if (c1 == "categorical") var1 else var2
    say("`", num, "` continuous by `", cat_, "` categorical -- boxplot")
    ilm_plot_box(data, num, x = cat_, by = by, ...)
  } else {
    say("both categorical -- grouped bar")
    ilm_plot_bar(data, var1, by = var2, ...)
  }
  invisible(NULL)
}

## Work out a panel grid that is as square as it can be.
#' @keywords internal
#' @noRd
ilm_panel_grid <- function(n, nrow, ncol) {
  if (is.null(nrow) && is.null(ncol)) {
    ncol <- ceiling(sqrt(n)); nrow <- ceiling(n / ncol)
  } else if (is.null(nrow)) nrow <- ceiling(n / ncol)
  else if (is.null(ncol)) ncol <- ceiling(n / nrow)
  c(nrow, ncol)
}

#' Plot every column of a data frame
#'
#' [ilm_plot_var()] once per column, arranged in a grid.
#'
#' @inheritParams ilm_plot_var
#' @param cols Columns to use. A character vector of names, a
#'   regular expression, a predicate function such as `is.numeric`, or
#'   `NULL` for all of them -- see [ilm_selection].
#' @param nrow,ncol Panel grid. Default is as square as it goes.
#' @return `NULL`, invisibly.
#' @seealso [ilm_plot_all()], [ilm_plot_var_pairs()].
#' @examples
#' ilm_plot_var_all(mtcars, cols = c("mpg", "cyl", "wt", "gear"))
#' @export
ilm_plot_var_all <- function(data, var2 = NULL, by = NULL, cols = NULL,
                             nrow = NULL, ncol = NULL, verbose = FALSE, ...) {
  ilm_plot_frame_check(data)
  g <- if (is.null(by)) character() else as.character(by)
  miss <- setdiff(c(var2, g, cols), names(data))
  if (length(miss))
    stop("column(s) not found in the data: ", paste(miss, collapse = ", "),
         call. = FALSE)
  target <- ilm_resolve_cols(data, cols, exclude = c(var2, g))
  if (!length(target)) stop("no columns to plot", call. = FALSE)
  ## several grouping columns become one real column on a local copy, since the
  ## single-variable plotters each take one `by` column by name
  by_col <- if (length(g) == 1L) g else if (length(g) > 1L) {
    data[[".ilm_by"]] <- interaction(data[g], drop = TRUE, sep = " / ")
    ".ilm_by"
  } else NULL
  gr <- ilm_panel_grid(length(target), nrow, ncol)
  op <- graphics::par(mfrow = gr); on.exit(graphics::par(op))
  for (cn in target)
    ilm_plot_var(data, cn, var2 = var2, by = by_col, verbose = verbose, ...)
  invisible(NULL)
}

#' Pairwise plots
#'
#' Every pair of columns at once, handling mixed numeric and categorical
#' columns rather than only numeric ones as a classic scatterplot matrix does.
#'
#' @param data A data frame.
#' @param cols Columns to use. A character vector of names, a
#'   regular expression, a predicate function such as `is.numeric`, or
#'   `NULL` for all of them -- see [ilm_selection].
#' @param by Optional grouping column.
#' @param ... Passed to [tinyplot::tinypairs()].
#' @return `NULL`, invisibly.
#' @seealso [ilm_plot_scatter()] for one pair,
#'   [ilm_check_collinearity()] for what a pairs plot cannot show about a fit.
#' @examples
#' ilm_plot_var_pairs(mtcars, cols = c("mpg", "wt", "hp"), by = "cyl")
#' @export
ilm_plot_var_pairs <- function(data, cols = NULL, by = NULL, ...) {
  ilm_plot_frame_check(data)
  g <- if (is.null(by)) character() else as.character(by)
  miss <- setdiff(c(g, cols), names(data))
  if (length(miss))
    stop("column(s) not found in the data: ", paste(miss, collapse = ", "),
         call. = FALSE)
  target <- ilm_resolve_cols(data, cols, exclude = g)
  if (length(target) < 2L)
    stop("at least 2 columns are needed to plot pairs of them; ",
         length(target), " given", call. = FALSE)
  bv <- if (length(g)) interaction(data[g], drop = TRUE, sep = " / ") else NULL
  ## tinypairs() arrived in tinyplot 0.7.0 and is the only part of this package
  ## that needs it -- everything else works from 0.6.1. Say so, and name what
  ## to do instead, rather than failing on a missing object.
  pairs_fn <- tryCatch(getExportedValue("tinyplot", "tinypairs"),
                       error = function(e) NULL)
  if (is.null(pairs_fn))
    stop("a pairs plot needs tinyplot >= 0.7.0, which supplies tinypairs(); ",
         "this is tinyplot ", utils::packageVersion("tinyplot"),
         ". Either update tinyplot, or use ilm_plot_scatter() one pair at a ",
         "time, or ilm_plot_var_all() for each column on its own.",
         call. = FALSE)
  pairs_fn(data[target], by = bv, ...)
  invisible(NULL)
}

#' The bootstrap distribution behind a group difference
#'
#' Draws the replicate differences from one row of an [ilm_boot_diff()] result,
#' with zero and the interval marked. The summary says where the difference is;
#' this says what the resampling actually produced -- whether it is symmetric,
#' skewed, or piled against a boundary, which the interval alone cannot show.
#'
#' @param x An [ilm_boot_diff()] result.
#' @param row Which comparison to draw, when the result has several. The
#'   default draws the first, and having more than one is normal now that
#'   [ilm_boot_diff()] compares every pair.
#' @param type `"density"` or `"histogram"`.
#' @param ref_line Where to draw the reference line.
#' @param ... Passed to [tinyplot::tinyplot()].
#' @return `NULL`, invisibly.
#' @seealso [ilm_boot_diff()].
#' @examples
#' d <- ilm_sim()
#' b <- ilm_boot_diff(d, "score", "grp", R = 300, seed = 1)
#' ilm_plot_boot_diff(b, row = 1)
#' @export
ilm_plot_boot_diff <- function(x, row = 1L, type = c("density", "histogram"),
                               ref_line = 0, ...) {
  type <- match.arg(type)
  dr <- attr(x, "draws")
  if (!is.data.frame(x) || is.null(dr))
    stop("`x` must be an ilm_boot_diff() result, which carries the replicate ",
         "draws; it is ", class(x)[1], call. = FALSE)
  row <- as.integer(row)
  if (length(row) != 1L || is.na(row) || row < 1L || row > nrow(x))
    stop("`row` must be one of the ", nrow(x), " comparison(s) in `x`",
         call. = FALSE)
  d <- dr[, row]
  tinyplot::tinyplot(d,
    type = if (type == "histogram") tinyplot::type_histogram() else "density",
    xlab = paste0(x$to[row], " - ", x$from[row], " (", x$stat[row], ")"),
    main = sprintf("%d replicates; %.0f%% above %g", length(d),
                   100 * mean(d > ref_line), ref_line), ...)
  tinyplot::tinyplot_add(x = ref_line, type = "vline", lty = 2)
  tinyplot::tinyplot_add(x = c(x$lower[row], x$upper[row]), type = "vline",
                         lty = 3, col = "gray40")
  invisible(NULL)
}

#' Combine several plots into one figure
#'
#' Sets up a panel grid and evaluates each plot-producing expression into it.
#'
#' It takes **expressions**, not plot objects, and that is base graphics rather
#' than a shortcut: a base plot draws immediately and leaves nothing behind to
#' combine, so the drawing has to happen inside the grid.
#'
#' @param ... Plot-producing expressions, for instance
#'   `ilm_plot_histogram(mtcars, "mpg")`.
#' @param nrow,ncol Panel grid. Default is as square as it goes.
#' @return `NULL`, invisibly.
#' @seealso [ilm_plot_var_all()].
#' @examples
#' ilm_plot_c(
#'   ilm_plot_histogram(mtcars, "mpg"),
#'   ilm_plot_box(mtcars, "mpg", x = "cyl"),
#'   nrow = 1
#' )
#' @export
ilm_plot_c <- function(..., nrow = NULL, ncol = NULL) {
  exprs <- as.list(substitute(list(...)))[-1]
  n <- length(exprs)
  if (!n) stop("give at least one plot expression", call. = FALSE)
  gr <- ilm_panel_grid(n, nrow, ncol)
  op <- graphics::par(mfrow = gr); on.exit(graphics::par(op))
  env <- parent.frame()
  for (e in exprs) eval(e, envir = env)
  invisible(NULL)
}
