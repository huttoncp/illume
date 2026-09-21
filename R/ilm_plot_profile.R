## ---------------------------------------------------------------------------
## Plots for the profiling set, on the same tinyplot backend the rest of the
## exploration layer uses.
##
## The `_na` variants are type-checked wrappers and nothing more: an
## ilm_reduce_na() result already inherits from ilm_reduce(), so the base
## functions work on one unchanged. They exist so that the missingness surface
## is discoverable under a consistent name, not because anything new was needed.
## ---------------------------------------------------------------------------

#' Map the observations from a reduction
#'
#' Each row placed on two of the dimensions -- the usual individuals map of a
#' PCA, MCA or mixed analysis.
#'
#' @param x An [ilm_reduce()] result.
#' @param dims Which two dimensions, by number.
#' @param by Optional vector to colour by, the same length as the data the
#'   reduction was built from -- a cluster assignment, or a column of the
#'   original data.
#' @param ... Passed to [tinyplot::tinyplot()].
#' @return `NULL`, invisibly.
#' @seealso [ilm_plot_cluster()], [ilm_plot_reduce_scree()].
#' @examples
#' r <- ilm_reduce(mtcars)
#' ilm_plot_reduce(r)
#' ilm_plot_reduce(r, by = factor(mtcars$cyl))
#' @export
ilm_plot_reduce <- function(x, dims = c(1, 2), by = NULL, ...) {
  if (!inherits(x, "ilm_reduce"))
    stop("`x` must be an ilm_reduce() result; it is ", class(x)[1], call. = FALSE)
  if (length(dims) != 2L)
    stop("`dims` must name exactly 2 dimensions", call. = FALSE)
  d1 <- paste0("dim", dims[1]); d2 <- paste0("dim", dims[2])
  if (!all(c(d1, d2) %in% names(x$ind_coord)))
    stop("those dimensions were not retained; the reduction kept ", x$ndim,
         call. = FALSE)
  pct <- x$eig$pct_var
  tinyplot::tinyplot(x = x$ind_coord[[d1]], y = x$ind_coord[[d2]], by = by,
                     type = "points",
                     xlab = sprintf("Dim %d (%.1f%%)", dims[1], pct[dims[1]]),
                     ylab = sprintf("Dim %d (%.1f%%)", dims[2], pct[dims[2]]),
                     ...)
  tinyplot::tinyplot_add(x = 0, type = "vline", lty = 3, col = "gray70")
  tinyplot::tinyplot_add(y = 0, type = "hline", lty = 3, col = "gray70")
  invisible(NULL)
}

#' Scree plot for a reduction
#'
#' How much variance each dimension accounts for, which is what decides how
#' many are worth keeping.
#'
#' @param x An [ilm_reduce()] result.
#' @param ... Passed to [tinyplot::tinyplot()].
#' @return `NULL`, invisibly.
#' @seealso [ilm_reduce()].
#' @examples
#' ilm_plot_reduce_scree(ilm_reduce(mtcars))
#' @export
ilm_plot_reduce_scree <- function(x, ...) {
  if (!inherits(x, "ilm_reduce"))
    stop("`x` must be an ilm_reduce() result; it is ", class(x)[1], call. = FALSE)
  tinyplot::tinyplot(x = factor(x$eig$dim), y = x$eig$pct_var, type = "barplot",
                     xlab = "dimension", ylab = "% of variance explained", ...)
  invisible(NULL)
}

#' Which variables a dimension is made of
#'
#' A sorted dot plot of squared loadings. Deliberately not the classic
#' correlation circle: that plot asks more visual literacy of a reader than a
#' ranked list does, and gives no more here.
#'
#' @param x An [ilm_reduce()] result.
#' @param dim Which dimension, by number.
#' @param top_n How many variables to show.
#' @param ... Passed to [tinyplot::tinyplot()].
#' @return `NULL`, invisibly.
#' @seealso [ilm_reduce()].
#' @examples
#' ilm_plot_reduce_contrib(ilm_reduce(mtcars), dim = 1)
#' @export
ilm_plot_reduce_contrib <- function(x, dim = 1, top_n = 10, ...) {
  if (!inherits(x, "ilm_reduce"))
    stop("`x` must be an ilm_reduce() result; it is ", class(x)[1], call. = FALSE)
  ## `dim_want`, not a bare `dim`: var_contrib has a column of that name, so a
  ## bare `dim` in the filter below would resolve to the column and turn the
  ## comparison into dim == dim, which is always TRUE. Caught by a test that
  ## expected an out-of-range dimension to error and got a plot instead.
  dim_want <- dim
  d <- x$var_contrib[x$var_contrib$dim == dim_want, , drop = FALSE]
  if (!nrow(d))
    stop("dimension ", dim, " was not retained; the reduction kept ", x$ndim,
         call. = FALSE)
  d <- d[order(-d$sqload), , drop = FALSE]
  d <- d[seq_len(min(top_n, nrow(d))), , drop = FALSE]
  d$variable <- factor(d$variable, levels = rev(d$variable))
  tinyplot::tinyplot(x = d$sqload, y = d$variable, type = "points",
                     xlab = "squared loading", ylab = "",
                     main = sprintf("Dimension %d", dim), ...)
  invisible(NULL)
}

#' Map the clusters
#'
#' The observations on two dimensions, coloured by which cluster they landed in.
#'
#' @param x An [ilm_cluster()] or [ilm_profile()] result.
#' @param dims Which two coordinates, by number.
#' @param ... Passed to [tinyplot::tinyplot()].
#' @return `NULL`, invisibly.
#' @seealso [ilm_cluster()], [ilm_plot_cluster_gap()].
#' @examples
#' ilm_plot_cluster(ilm_cluster(ilm_reduce(mtcars), k_max = 5, B = 25, seed = 1))
#' @export
ilm_plot_cluster <- function(x, dims = c(1, 2), ...) {
  if (inherits(x, "ilm_profile")) x <- x$cluster
  if (!inherits(x, "ilm_cluster"))
    stop("`x` must be an ilm_cluster() or ilm_profile() result; it is ",
         class(x)[1], call. = FALSE)
  if (ncol(x$coords) < max(dims))
    stop("`dims` asks for more coordinates than there are (", ncol(x$coords),
         ")", call. = FALSE)
  tinyplot::tinyplot(x = x$coords[, dims[1]], y = x$coords[, dims[2]],
                     by = factor(x$ind_cluster$cluster), type = "points",
                     xlab = paste0("Dim ", dims[1]),
                     ylab = paste0("Dim ", dims[2]),
                     legend = list(title = "cluster"), ...)
  invisible(NULL)
}

#' The gap statistic across every k considered
#'
#' With the chosen `k` marked. Only available when `k` was searched for rather
#' than given.
#'
#' @param x An [ilm_cluster()] or [ilm_profile()] result.
#' @param ... Passed to [tinyplot::tinyplot()].
#' @return `NULL`, invisibly.
#' @seealso [ilm_cluster()].
#' @examples
#' ilm_plot_cluster_gap(ilm_cluster(ilm_reduce(mtcars), k_max = 5, B = 25,
#'                                  seed = 1))
#' @export
ilm_plot_cluster_gap <- function(x, ...) {
  if (inherits(x, "ilm_profile")) x <- x$cluster
  if (!inherits(x, "ilm_cluster"))
    stop("`x` must be an ilm_cluster() or ilm_profile() result; it is ",
         class(x)[1], call. = FALSE)
  if (is.null(x$gap))
    stop("there is no gap statistic to plot: `k` was given rather than ",
         "searched for", call. = FALSE)
  tab <- x$gap$Tab
  tinyplot::tinyplot(x = seq_len(nrow(tab)), y = tab[, "gap"],
                     ymin = tab[, "gap"] - tab[, "SE.sim"],
                     ymax = tab[, "gap"] + tab[, "SE.sim"],
                     type = "pointrange", xlab = "number of clusters (k)",
                     ylab = "gap statistic", ...)
  tinyplot::tinyplot_add(x = x$k, type = "vline", lty = 2, col = "firebrick")
  invisible(NULL)
}

#' Map the clusters from a profile
#'
#' Convenience wrapper for [ilm_plot_cluster()] when what you have is a whole
#' [ilm_profile()] rather than the clustering on its own.
#'
#' @inheritParams ilm_plot_cluster
#' @param x An [ilm_profile()] result.
#' @return `NULL`, invisibly.
#' @seealso [ilm_profile()].
#' @examples
#' ilm_plot_profile(ilm_profile(mtcars, k_max = 5, B = 25, seed = 1))
#' @export
ilm_plot_profile <- function(x, dims = c(1, 2), ...) {
  if (!inherits(x, "ilm_profile"))
    stop("`x` must be an ilm_profile() result; it is ", class(x)[1],
         call. = FALSE)
  ilm_plot_cluster(x$cluster, dims = dims, ...)
}

## ---- the missingness surface ----------------------------------------------

#' @rdname ilm_plot_reduce
#' @export
ilm_plot_reduce_na <- function(x, dims = c(1, 2), by = NULL, ...) {
  if (!inherits(x, "ilm_reduce_na"))
    stop("`x` must be an ilm_reduce_na() result; it is ", class(x)[1],
         call. = FALSE)
  ilm_plot_reduce(x, dims = dims, by = by, ...)
}

#' @rdname ilm_plot_reduce_scree
#' @export
ilm_plot_reduce_scree_na <- function(x, ...) {
  if (!inherits(x, "ilm_reduce_na"))
    stop("`x` must be an ilm_reduce_na() result; it is ", class(x)[1],
         call. = FALSE)
  ilm_plot_reduce_scree(x, ...)
}

#' @rdname ilm_plot_reduce_contrib
#' @export
ilm_plot_reduce_contrib_na <- function(x, dim = 1, top_n = 10, ...) {
  if (!inherits(x, "ilm_reduce_na"))
    stop("`x` must be an ilm_reduce_na() result; it is ", class(x)[1],
         call. = FALSE)
  ilm_plot_reduce_contrib(x, dim = dim, top_n = top_n, ...)
}

#' @rdname ilm_plot_cluster
#' @export
ilm_plot_cluster_na <- function(x, dims = c(1, 2), ...) {
  if (!inherits(x, "ilm_cluster_na") && !inherits(x, "ilm_profile_na"))
    stop("`x` must be an ilm_cluster_na() or ilm_profile_na() result; it is ",
         class(x)[1], call. = FALSE)
  ilm_plot_cluster(x, dims = dims, ...)
}

#' @rdname ilm_plot_cluster_gap
#' @export
ilm_plot_cluster_gap_na <- function(x, ...) {
  if (!inherits(x, "ilm_cluster_na") && !inherits(x, "ilm_profile_na"))
    stop("`x` must be an ilm_cluster_na() or ilm_profile_na() result; it is ",
         class(x)[1], call. = FALSE)
  ilm_plot_cluster_gap(x, ...)
}

#' @rdname ilm_plot_profile
#' @export
ilm_plot_profile_na <- function(x, dims = c(1, 2), ...) {
  if (!inherits(x, "ilm_profile_na"))
    stop("`x` must be an ilm_profile_na() result; it is ", class(x)[1],
         call. = FALSE)
  ilm_plot_profile(x, dims = dims, ...)
}
