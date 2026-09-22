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
#' @param cols Columns to use. A character vector of names, a
#'   regular expression, a predicate function such as `is.numeric`, or
#'   `NULL` for all of them -- see [ilm_selection].
#' @param method `"pcamix"` (the default) for PCA, MCA or FAMD depending on
#'   the column types, in closed form. `"glrm"` fits a generalized low rank
#'   model instead, which uses a loss appropriate to each column's type rather
#'   than squared error on one-hot indicators, and reconstructs a category as a
#'   category. It costs an iterative fit, and on all-numeric data the two are
#'   the same model -- see [ilm_glrm()] for when it is worth that.
#' @param ... Passed to [ilm_glrm()] when `method = "glrm"`.
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
ilm_reduce <- function(data, cols = NULL, ndim = 5,
                       method = c("pcamix", "glrm"), ...) {
  method <- match.arg(method)
  ## An ilm_anomaly() result is accepted directly: the flagged rows are what
  ## the user wants to look at next, and rebuilding that subset by hand from
  ## `row` is both a papercut and a chance to line the wrong rows up.
  if (inherits(data, "ilm_anomaly"))
    data <- ilm_from_anomaly(data, "ilm_reduce")
  if (method == "glrm") return(ilm_reduce_glrm(data, cols, ndim, ...))
  ilm_require_pcamixdata()
  if (!is.data.frame(data))
    stop("`data` must be a data frame; it is ", class(data)[1], call. = FALSE)
  keep <- ilm_resolve_cols(data, cols)
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

  ## as.data.frame(), not just the subset. `sub[is_num]` on a TIBBLE returns a
  ## tibble, and PCAmix checks its columns with `is.numeric(X.quanti[, j])` --
  ## which for a tibble is a one-column tibble rather than a vector, so every
  ## column looks non-numeric and the whole call dies on "All variables in
  ## X.quanti must be numeric". Nothing in that message mentions tibbles, and
  ## a tibble is what anyone gets from readr, dplyr or gapminder, so this hit
  ## the most ordinary way to arrive at the function. `quali` was already
  ## coerced on the next line, which is why only the numeric half broke.
  quanti <- if (any(is_num)) as.data.frame(sub[is_num]) else NULL
  quali <- if (any(is_cat)) as.data.frame(lapply(sub[is_cat], as.factor)) else NULL
  method <- if (!is.null(quanti) && !is.null(quali)) "famd"
            else if (!is.null(quanti)) "pca" else "mca"

  ## How many dimensions exist at all. For MCA that is the number of levels
  ## less the number of variables. For PCA and the mixed case it is
  ## min(n - 1, p) -- NOT p - 1, which is a different quantity and one short:
  ## two columns have two components, not one, and asking PCAmix for one is an
  ## error rather than a smaller answer, so ilm_reduce() used to fail outright
  ## on any two-column selection.
  max_dim <- if (method == "mca")
               sum(vapply(quali, nlevels, 1L)) - ncol(quali)
             else min(nrow(sub) - 1L, ncol(sub))
  ## PCAmix needs at least two, and with fewer than two available there is
  ## nothing to reduce
  ndim <- min(ndim, max_dim)
  if (ndim < 2L) {
    if (max_dim < 2L)
      stop("there are only ", max_dim, " dimension(s) available from these ",
           "columns, and a reduction needs at least 2", call. = FALSE)
    ndim <- 2L
  }

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
#' @param cols Columns to use. A character vector of names, a
#'   regular expression, a predicate function such as `is.numeric`, or
#'   `NULL` for all of them -- see [ilm_selection].
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
  ## An ilm_anomaly() result is accepted directly: the flagged rows are what
  ## the user wants to look at next, and rebuilding that subset by hand from
  ## `row` is both a papercut and a chance to line the wrong rows up.
  if (inherits(data, "ilm_anomaly"))
    data <- ilm_from_anomaly(data, "ilm_reduce_na")
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

## ---- the generalized low rank route ----------------------------------------
##
## The same SHAPE of answer as the PCAmix one -- scores per row, a contribution
## per variable per dimension -- so ilm_cluster() and ilm_profile() need to
## know nothing about which produced it. What differs is underneath: a loss per
## column type rather than squared error on scaled indicators.
##
## "Variance explained" here is the share of the fitted linear predictor each
## dimension carries. That coincides with the usual quantity when every loss is
## quadratic and generalises it when they are not. It is not an eigenvalue, and
## nothing downstream treats it as one.

#' @keywords internal
#' @noRd
ilm_reduce_glrm <- function(data, cols, ndim, ...) {
  g <- ilm_glrm(data, cols = cols, rank = ndim, ...)
  k <- g$rank
  n <- nrow(g$scores)
  ss <- vapply(seq_len(k), function(j)
    sum((as.matrix(g$scores)[, j, drop = FALSE] %*%
           g$archetypes[j, , drop = FALSE])^2), 0)
  ## Against the TOTAL variation in the encoded data, not against the part
  ## these dimensions already account for -- normalising by the latter makes
  ## the retained dimensions explain 100% of it by construction, whatever the
  ## rank, which tells the reader nothing.
  A <- do.call(cbind, lapply(g$encoding$blocks, function(b) {
    if (b$loss %in% c("quadratic", "poisson")) {
      v <- b$target; v[is.na(v)] <- 0; matrix(v, ncol = 1L)
    } else {
      yi <- b$target; M <- matrix(0, n, length(b$cols))
      okr <- !is.na(yi); M[cbind(which(okr), yi[okr])] <- 1
      sweep(M, 2L, colMeans(M), "-")
    }
  }))
  tot <- max(sum(A^2), .Machine$double.eps)
  pct <- 100 * pmin(ss / tot, 1)
  ord <- order(-pct)
  eig <- data.frame(dim = seq_len(k), eigenvalue = ss[ord],
                    pct_var = pct[ord], cum_pct_var = cumsum(pct[ord]),
                    row.names = NULL)
  ## a variable contributes to a dimension through every column its block
  ## occupies, so the squared loadings are summed over the block
  vc <- do.call(rbind, lapply(seq_len(k), function(j) {
    sq <- vapply(g$encoding$blocks, function(b)
      sum(g$archetypes[ord[j], b$cols]^2), 0)
    data.frame(dim = j, variable = names(g$encoding$blocks),
               sqload = sq / max(sum(sq), .Machine$double.eps),
               row.names = NULL)
  }))
  co <- as.data.frame(as.matrix(g$scores)[, ord, drop = FALSE])
  names(co) <- paste0("dim", seq_len(k))
  structure(list(method = "glrm", eig = eig, ind_coord = co,
                 var_contrib = vc, n = n, ndim = k,
                 cols = g$columns, fit = g), class = "ilm_reduce")
}
