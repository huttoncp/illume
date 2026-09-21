## ---------------------------------------------------------------------------
## Dimension reduction over mixed column types.
##
## PCA when every column is numeric, MCA when every column is categorical, and
## a mixed method when both are present -- chosen from the data rather than
## asked for, because the choice is forced by the column types and making the
## user name it only invites getting it wrong.
##
## The mixed case goes through PCAmixdata::PCAmix(), which handles all three
## natively. That was checked against FactoMineR::FAMD() on the same data before
## being relied on: eigenvalues agreed to two decimals and the individual
## coordinates correlated at |r| = 1.0000 on every dimension, while scaling
## equal or better at every size tried, up to 20,000 rows and 120 columns.
## PCAmixdata brings nothing beyond base R's `graphics`, where FactoMineR pulls
## in 119 transitive packages -- which is why it stays in Suggests behind a
## guard rather than becoming a hard dependency.
## ---------------------------------------------------------------------------

#' @keywords internal
#' @noRd
ilm_require_pcamixdata <- function() {
  if (!requireNamespace("PCAmixdata", quietly = TRUE))
    stop("the PCAmixdata package is needed for ilm_reduce(). Install it with ",
         'install.packages("PCAmixdata").', call. = FALSE)
}

#' Reduce a data frame's variables to a few dimensions
#'
#' Takes the columns of a data frame and summarises them as a small number of
#' dimensions. Which method that means is decided by the column types: PCA when
#' they are all numeric, multiple correspondence analysis when they are all
#' categorical, and a mixed method when both are present. A date or date-time
#' column cannot be used by any of the three and is dropped with a message;
#' convert it to something numeric first if it should count.
#'
#' The mixed method is Chavent et al.'s, which belongs to the same
#' generalised-PCA family as the FAMD of Pages without being a
#' reimplementation of it. The two were checked against each other directly and
#' agree for this purpose, matching on eigenvalues and on individual
#' coordinates.
#'
#' @param data A data frame.
#' @param cols Columns to include, as a character vector. Default is all of
#'   them.
#' @param ndim Number of dimensions to keep.
#' @return An object of class `"ilm_reduce"`: `method` (`"pca"`, `"mca"` or
#'   `"famd"`), `eig` (dimension, eigenvalue, percent of variance and its
#'   cumulative total), `ind_coord` (`row_id` and one column per retained
#'   dimension -- this is what [ilm_cluster()] takes), `var_contrib`
#'   (`variable`, `dim`, `sqload`: how strongly each original variable relates
#'   to each dimension, on a 0 to 1 scale, for numeric and categorical
#'   variables alike), `n`, and `fit`, the underlying `PCAmixdata::PCAmix()`
#'   object for anyone who wants to go past this wrapper.
#' @seealso [ilm_cluster()] to group the rows, [ilm_profile()] for the whole
#'   pipeline, [ilm_reduce_na()] for the same thing applied to missingness.
#' @references
#' Chavent, M., Kuentz-Simonet, V., Labenne, A. and Saracco, J. (2014).
#' Multivariate analysis of mixed data: the PCAmixdata R package. arXiv
#' 1411.4911.
#' @examples
#' r <- ilm_reduce(mtcars)
#' r
#' head(r$var_contrib[order(-r$var_contrib$sqload), ])
#' @export
ilm_reduce <- function(data, cols = NULL, ndim = 5) {
  ilm_require_pcamixdata()
  if (!is.data.frame(data))
    stop("`data` must be a data frame; it is ", class(data)[1], call. = FALSE)
  keep <- if (is.null(cols)) names(data) else as.character(cols)
  miss <- setdiff(keep, names(data))
  if (length(miss))
    stop("column(s) not found in the data: ", paste(miss, collapse = ", "),
         ". Available: ", paste(utils::head(names(data), 12), collapse = ", "),
         call. = FALSE)
  sub <- data[keep]

  is_num <- vapply(sub, is.numeric, TRUE)
  is_cat <- vapply(sub, function(x)
    is.factor(x) || is.character(x) || is.logical(x), TRUE)
  dropped <- names(sub)[!is_num & !is_cat]
  if (length(dropped))
    message("ilm_reduce(): dropping column(s) that are neither numeric nor ",
            "categorical, so no method here can use them: ",
            paste(dropped, collapse = ", "))
  if (!any(is_num) && !any(is_cat))
    stop("no numeric or categorical columns to reduce", call. = FALSE)

  quanti <- if (any(is_num)) sub[is_num] else NULL
  quali <- if (any(is_cat)) as.data.frame(lapply(sub[is_cat], as.factor)) else NULL
  method <- if (!is.null(quanti) && !is.null(quali)) "famd"
            else if (!is.null(quanti)) "pca" else "mca"

  ## MCA's dimensionality is the number of levels less the number of variables,
  ## not the number of variables: asking for more than that is degenerate
  max_dim <- if (method == "mca")
               sum(vapply(quali, nlevels, 1L)) - ncol(quali) else ncol(sub) - 1L
  ndim <- max(1L, min(ndim, max_dim, nrow(sub) - 1L))

  fit <- PCAmixdata::PCAmix(X.quanti = quanti, X.quali = quali, ndim = ndim,
                            rename.level = TRUE, graph = FALSE)

  eig <- data.frame(dim = seq_len(nrow(fit$eig)),
                    eigenvalue = round(fit$eig[, 1], 4),
                    pct_var = round(fit$eig[, 2], 3),
                    cum_pct_var = round(fit$eig[, 3], 3),
                    stringsAsFactors = FALSE)
  rownames(eig) <- NULL

  ic <- as.data.frame(fit$ind$coord)
  names(ic) <- paste0("dim", seq_len(ncol(ic)))
  ind_coord <- cbind(row_id = seq_len(nrow(ic)), ic)
  rownames(ind_coord) <- NULL

  sq <- fit$sqload
  var_contrib <- data.frame(
    variable = rep(rownames(sq), ncol(sq)),
    dim = rep(seq_len(ncol(sq)), each = nrow(sq)),
    sqload = round(as.vector(sq), 4), stringsAsFactors = FALSE)
  var_contrib <- var_contrib[order(var_contrib$dim, -var_contrib$sqload), ,
                             drop = FALSE]
  rownames(var_contrib) <- NULL

  structure(list(method = method, eig = eig, ind_coord = ind_coord,
                 var_contrib = var_contrib, n = nrow(sub), ndim = ndim,
                 cols = keep, fit = fit), class = "ilm_reduce")
}

#' @export
print.ilm_reduce <- function(x, ...) {
  tag <- if (inherits(x, "ilm_reduce_na")) "ilm_reduce_na" else "ilm_reduce"
  cat(sprintf("<%s> method = %s, n = %d, %d dimension(s) retained\n",
              tag, x$method, x$n, x$ndim))
  cat(sprintf("  first %d dimension(s) explain %.1f%% of the variance\n",
              min(3L, nrow(x$eig)), x$eig$cum_pct_var[min(3L, nrow(x$eig))]))
  cat("\n  strongest variable per dimension (squared loading)\n")
  for (d in sort(unique(x$var_contrib$dim))) {
    sl <- x$var_contrib[x$var_contrib$dim == d, , drop = FALSE]
    k <- which.max(sl$sqload)
    cat(sprintf("    dim %-3d %-24s %.3f\n", d, sl$variable[k], sl$sqload[k]))
  }
  invisible(x)
}

## ---- the same thing, applied to what is missing ----------------------------
##
## A present/missing marker is a two-level categorical variable, so building a
## frame of those markers and handing it to ilm_reduce() -- which then picks MCA,
## every column being categorical -- reuses the whole pipeline rather than
## duplicating it. Checked on synthetic data with two indicators built to go
## missing together: the resulting dimension recovered them at squared loadings
## of about 0.85 each, while a third, independently missing indicator loaded on
## a separate dimension.

#' @keywords internal
#' @noRd
ilm_build_na_indicator <- function(data, cols) {
  keep0 <- if (is.null(cols) || !length(cols)) names(data) else as.character(cols)
  p_na <- vapply(data[keep0], function(v) mean(is.na(v)), 1)
  degenerate <- p_na == 0 | p_na == 1
  if (any(degenerate))
    message("ilm_reduce_na(): dropping column(s) whose missingness never ",
            "varies (always or never missing): ",
            paste(keep0[degenerate], collapse = ", "))
  keep <- keep0[!degenerate]
  if (length(keep) < 2L)
    stop("at least 2 columns need some -- but not all -- values missing ",
         "before there is a missingness pattern to profile; ", length(keep),
         " qualif", if (length(keep) == 1L) "ies" else "y", call. = FALSE)
  out <- as.data.frame(lapply(keep, function(cn)
    factor(ifelse(is.na(data[[cn]]), "missing", "present"),
           levels = c("present", "missing"))))
  names(out) <- keep
  out
}

#' Reduce a data frame's missingness pattern to a few dimensions
#'
#' The missingness counterpart to [ilm_reduce()]. Builds a present/missing
#' marker for every column, drops any whose missingness never varies, and
#' reduces the markers -- which, being two-level categorical variables, always
#' takes the MCA route. The dimensions that come back describe which columns
#' tend to go missing *together*, which is what separates a block of variables
#' lost to one skipped section from values that went missing independently.
#'
#' @param data A data frame.
#' @param cols Columns to consider, as a character vector. Default is all.
#' @param ndim Number of dimensions to keep.
#' @return An object of class `"ilm_reduce_na"`, which is also an
#'   `"ilm_reduce"` -- see [ilm_reduce()] for the shared structure. In
#'   `var_contrib`, `sqload` means how strongly a column's *missingness*
#'   relates to a dimension, not its values.
#' @seealso [ilm_profile_na()] for the whole pipeline,
#'   [ilm_check_missing()] for whether any of it matters to your model.
#' @examples
#' r <- ilm_reduce_na(airquality)
#' r
#' @export
ilm_reduce_na <- function(data, cols = NULL, ndim = 5) {
  if (!is.data.frame(data))
    stop("`data` must be a data frame; it is ", class(data)[1], call. = FALSE)
  if (!is.null(cols)) {
    miss <- setdiff(as.character(cols), names(data))
    if (length(miss))
      stop("column(s) not found in the data: ", paste(miss, collapse = ", "),
           call. = FALSE)
  }
  out <- ilm_reduce(ilm_build_na_indicator(data, cols), ndim = ndim)
  class(out) <- c("ilm_reduce_na", class(out))
  out
}
