## ---------------------------------------------------------------------------
## Reduce, cluster, then say what each cluster IS.
##
## A cluster number tells you nothing. What is wanted is "cluster 4 is the
## heavy, thirsty cars", and the route to that runs through the dimensions: find
## which dimensions a cluster sits unusually far along, then name the variables
## that load most heavily on those dimensions.
##
## The unusualness measure is the v-test, the same device FactoMineR's catdes()
## and condes() use: for cluster c on dimension d,
## (mean_c - mean_overall) / (sd_overall / sqrt(n_c)). It is a standardised
## statement of how surprising the cluster's average position would be if
## membership had nothing to do with that dimension. It is NOT a p-value from a
## model fitted to the whole design -- the clusters were found from the same
## coordinates being tested -- so the 1.96 default is a useful cutoff rather
## than a formal test, and is documented as one.
##
## The characterisation logic is shared with the missingness variant. Only the
## labelling differs: high/low for values, "missing: x" for indicators, since a
## present/missing marker has no meaningful high or low.
## ---------------------------------------------------------------------------

## For values: the direction combines the v-test's sign with the variable's own
## correlation sign on that dimension, because a dimension can run either way.
#' @keywords internal
#' @noRd
ilm_label_var_direction <- function(v, vtest, d_idx, reduce_res) {
  qc <- reduce_res$fit$quanti.cor          # NULL when the method is MCA
  if (!is.null(qc) && v %in% rownames(qc))
    paste0(if (sign(vtest) * sign(qc[v, d_idx]) > 0) "high " else "low ", v)
  else v
}

## For missingness: every variable is a present/missing marker, so direction
## does not apply and the fact itself is the label.
#' @keywords internal
#' @noRd
ilm_label_var_missing <- function(v, vtest, d_idx, reduce_res)
  paste0("missing: ", v)

#' @keywords internal
#' @noRd
ilm_characterize_clusters <- function(reduce_res, cluster_res, vtest_threshold,
                                      top_n_vars, label_fn) {
  coords <- reduce_res$ind_coord
  dims <- setdiff(names(coords), "row_id")
  cl <- cluster_res$ind_cluster$cluster
  overall_mean <- vapply(coords[dims], mean, 1)
  overall_sd <- vapply(coords[dims], stats::sd, 1)

  rows <- list()
  for (cc in cluster_res$clusters$cluster) {
    members <- which(cl == cc)
    n_c <- length(members)
    if (!n_c) next
    for (d in dims) {
      d_idx <- as.integer(sub("dim", "", d))
      if (!is.finite(overall_sd[[d]]) || overall_sd[[d]] == 0) next
      vtest <- (mean(coords[[d]][members]) - overall_mean[[d]]) /
        (overall_sd[[d]] / sqrt(n_c))
      if (!is.finite(vtest) || abs(vtest) < vtest_threshold) next
      tv <- reduce_res$var_contrib[reduce_res$var_contrib$dim == d_idx, ,
                                   drop = FALSE]
      tv <- tv[order(-tv$sqload), , drop = FALSE]
      tv <- tv[seq_len(min(top_n_vars, nrow(tv))), , drop = FALSE]
      labs <- vapply(tv$variable, label_fn, "", vtest = vtest, d_idx = d_idx,
                     reduce_res = reduce_res)
      rows[[length(rows) + 1L]] <- data.frame(
        cluster = cc, dim = d_idx, vtest = round(vtest, 2),
        direction = if (vtest > 0) "high" else "low",
        top_variables = paste(labs, collapse = ", "),
        stringsAsFactors = FALSE)
    }
  }
  out <- if (length(rows)) do.call(rbind, rows) else data.frame(
    cluster = integer(0), dim = integer(0), vtest = numeric(0),
    direction = character(0), top_variables = character(0),
    stringsAsFactors = FALSE)
  if (nrow(out)) out <- out[order(out$cluster, -abs(out$vtest)), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' @keywords internal
#' @noRd
ilm_profile_summary_lines <- function(characterization, cluster_res) {
  ct <- cluster_res$clusters
  ind <- cluster_res$ind_cluster
  vapply(ct$cluster, function(cc) {
    rc <- characterization[characterization$cluster == cc, , drop = FALSE]
    si <- ct[ct$cluster == cc, , drop = FALSE]
    n_amb <- sum(ind$cluster == cc & ind$is_ambiguous)
    head_txt <- sprintf("Cluster %d (n = %d, %.1f%% of the data, %s)", cc,
                        si$size, si$pct, as.character(si$stability))
    ## Two dimensions can have the same top-loading variables -- which happens
    ## whenever there are few variables to go round, as in a missingness
    ## profile of a frame with two incomplete columns. Naming the same pair
    ## twice in one sentence says nothing the first mention did not, so only
    ## the first dimension to use a given set is reported.
    rc <- rc[!duplicated(rc$top_variables), , drop = FALSE]
    body <- if (!nrow(rc))
      paste0(head_txt, ": no dimension stood out past the threshold.")
    else sprintf("%s is characterised by %s.", head_txt,
                 paste0("dim ", rc$dim, " (", rc$top_variables, ")",
                        collapse = "; "))
    notes <- character()
    if (isTRUE(si$anomalous))
      notes <- c(notes, sprintf(
        "It is a small cluster, %.1f%% of observations: possibly a real minority pattern, possibly a data problem, but worth looking at either way.",
        si$pct))
    if (n_amb > 0)
      ## "N of its members" is plural whatever N is; only the verb agrees
      notes <- c(notes, sprintf(
        "%d of its members sit%s close enough to another cluster to be uncertain.",
        n_amb, if (n_amb == 1L) "s" else ""))
    paste(c(body, notes), collapse = " ")
  }, "")
}

#' Profile a data set: reduce, cluster, and describe the clusters
#'
#' Runs [ilm_reduce()], then [ilm_cluster()] on the dimensions it produces,
#' then says what distinguishes each cluster -- which dimensions it sits
#' unusually far along, and which original variables those dimensions are made
#' of. The result is a sentence per cluster rather than a column of numbers.
#'
#' @section How a cluster gets characterised:
#'
#' By a v-test: for cluster `c` and dimension `d`,
#' `(mean_c - mean_overall) / (sd_overall / sqrt(n_c))`. It measures how
#' surprising the cluster's average position on that dimension would be if
#' membership had nothing to do with it.
#'
#' It is **not a p-value**. The clusters were found from the very coordinates
#' being tested, so the usual sampling argument does not apply and the 1.96
#' default is a threshold that behaves sensibly rather than a 5% test. Treat a
#' characterisation as a description of the partition you have, not as evidence
#' that the partition is real -- the stability and silhouette figures in
#' [ilm_cluster()] are what speak to that.
#'
#' For numeric variables the reported direction combines the v-test's sign with
#' the variable's own correlation on that dimension, since a dimension can run
#' either way. Categorical variables are named without a direction.
#'
#' @inheritParams ilm_reduce
#' @param ... Passed to [ilm_cluster()], for instance `k`, `k_max`, `method`,
#'   `B` or `seed`.
#' @param vtest_threshold Smallest `|vtest|` for a dimension to count toward a
#'   cluster's description.
#' @param top_n_vars How many top-loading variables to name per dimension.
#' @return An object of class `"ilm_profile"`: `reduce`, `cluster`,
#'   `characterization` (one row per dimension that characterises a cluster)
#'   and `summary`, one sentence per cluster.
#' @seealso [ilm_reduce()], [ilm_cluster()], [ilm_plot_profile()],
#'   [ilm_profile_na()].
#' @examples
#' p <- ilm_profile(mtcars, k_max = 5, B = 25, seed = 1)
#' cat(p$summary, sep = "\n")
#' @export
ilm_profile <- function(data, cols = NULL, ndim = 5, ...,
                        vtest_threshold = 1.96, top_n_vars = 2) {
  rr <- ilm_reduce(data, cols = cols, ndim = ndim)
  cr <- ilm_cluster(rr, ...)
  ch <- ilm_characterize_clusters(rr, cr, vtest_threshold, top_n_vars,
                                  ilm_label_var_direction)
  structure(list(reduce = rr, cluster = cr, characterization = ch,
                 summary = ilm_profile_summary_lines(ch, cr)),
            class = "ilm_profile")
}

#' Profile which values are missing, and for whom
#'
#' The missingness counterpart to [ilm_profile()]: runs [ilm_reduce_na()] then
#' [ilm_cluster_na()], and describes each cluster of rows by which columns'
#' *missingness* sets it apart -- "cluster 2 is characterised by dim 1 (missing:
#' income, missing: age)" rather than by those columns' values.
#'
#' This is how a structured gap shows itself: a block of variables that go
#' missing together points at a shared cause, such as a section of a form
#' everyone in one group skipped, which is a different problem from values
#' going missing one at a time. [ilm_check_missing()] then says whether any of
#' it threatens the model you intend to fit.
#'
#' @inheritParams ilm_reduce_na
#' @param ... Passed to [ilm_cluster_na()].
#' @param vtest_threshold Smallest `|vtest|` for a dimension to count.
#' @param top_n_vars How many top-loading indicators to name per dimension.
#' @return An object of class `"ilm_profile_na"`, which is also an
#'   `"ilm_profile"`.
#' @seealso [ilm_check_missing()], [ilm_impute()], [ilm_profile()].
#' @examples
#' p <- ilm_profile_na(airquality, k_max = 4, B = 25, seed = 1)
#' cat(p$summary, sep = "\n")
#' @export
ilm_profile_na <- function(data, cols = NULL, ndim = 5, ...,
                           vtest_threshold = 1.96, top_n_vars = 2) {
  rr <- ilm_reduce_na(data, cols = cols, ndim = ndim)
  cr <- ilm_cluster_na(rr, ...)
  ch <- ilm_characterize_clusters(rr, cr, vtest_threshold, top_n_vars,
                                  ilm_label_var_missing)
  structure(list(reduce = rr, cluster = cr, characterization = ch,
                 summary = ilm_profile_summary_lines(ch, cr)),
            class = c("ilm_profile_na", "ilm_profile"))
}

#' @export
print.ilm_profile <- function(x, ...) {
  cat(if (inherits(x, "ilm_profile_na"))
        "<ilm_profile_na>  (profiling the pattern of missing values)"
      else "<ilm_profile>", "\n\n")
  print(x$reduce)
  cat("\n")
  print(x$cluster)
  cat("\n  what each cluster is\n")
  for (s in x$summary) cat(ilm_wrap(s, 76L, "    "), "\n\n")
  invisible(x)
}
