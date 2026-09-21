## ---------------------------------------------------------------------------
## Grouping observations, with the number of groups chosen rather than assumed.
##
## k comes from the gap statistic unless it is supplied. Both the search and the
## final clustering take k-means or hierarchical clustering; the hierarchical
## default is Euclidean distance with Ward's minimum-variance linkage as
## "ward.D2", which implements Ward's criterion correctly on ordinary
## (unsquared) distances. The older "ward.D" needs pre-squared distances to be
## right and survives in base R only for backwards compatibility.
##
## Two things are reported that a bare cluster assignment does not give you.
## Cluster STABILITY, by bootstrap Jaccard: resample the rows, recluster, and
## for each original cluster find its best-matching resampled one, averaged over
## B resamples. That is Hennig's method, written out here rather than taken from
## fpc::clusterboot(), which would cost 21 transitive packages for the one
## function; the computation itself -- intersection over union of two label sets
## -- is unambiguous, so the risk in writing it directly is low. And per
## observation AMBIGUITY, from the silhouette width, which is independent of how
## big a cluster is: a point can sit on a boundary inside a large, perfectly
## stable cluster, and that is worth knowing separately.
## ---------------------------------------------------------------------------

#' @keywords internal
#' @noRd
ilm_require_cluster <- function() {
  if (!requireNamespace("cluster", quietly = TRUE))
    stop("the cluster package is needed to choose k by the gap statistic. It ",
         "ships with R, so this should already be available; otherwise pass ",
         "`k` directly and no search is needed.", call. = FALSE)
}

#' @keywords internal
#' @noRd
ilm_cluster_funcluster <- function(method, dist_method, hclust_method, nstart) {
  if (method == "kmeans") {
    function(x, k) {
      ## kmeans()'s default Hartigan-Wong algorithm emits "Quick-TRANSfer stage
      ## steps exceeded maximum" on data with many tied points, which
      ## missingness-indicator coordinates inherently are: p binary columns
      ## admit at most 2^p distinct positions, so far fewer than there are
      ## rows. Checked directly -- the warning fired more than 70 times across
      ## one routine gap search while the resulting cluster sizes matched the
      ## data's true missingness combinations exactly. It reports the
      ## algorithm's reassignment bookkeeping hitting its cap on ties, not a
      ## failure to reach a valid partition, so it is silenced rather than
      ## raised on every call.
      list(cluster = withCallingHandlers(
        stats::kmeans(x, centers = k, nstart = nstart)$cluster,
        warning = function(w) {
          if (grepl("Quick-TRANSfer", conditionMessage(w)))
            invokeRestart("muffleWarning")
        }))
    }
  } else {
    function(x, k) {
      d <- stats::dist(x, method = dist_method)
      list(cluster = stats::cutree(stats::hclust(d, method = hclust_method),
                                   k = k))
    }
  }
}

#' @keywords internal
#' @noRd
ilm_cluster_stability <- function(coords, FUNcluster, k, orig_cluster, B,
                                  progress = NULL) {
  n <- nrow(coords)
  jac_sum <- numeric(k); jac_n <- integer(k)
  pb <- ilm_progress(B, progress); on.exit(pb$done(), add = TRUE)
  for (b in seq_len(B)) {
    pb$tick(b)
    idx <- sample.int(n, n, replace = TRUE)
    drawn <- unique(idx)
    boot_cluster <- tryCatch(FUNcluster(coords[idx, , drop = FALSE], k)$cluster,
                             error = function(e) NULL)
    if (is.null(boot_cluster)) next
    for (kk in seq_len(k)) {
      ## Restrict the ORIGINAL cluster to the rows this resample actually drew,
      ## so both sides of the comparison share one universe. Without it, the
      ## original cluster is compared against a set that structurally cannot
      ## contain the ~36.8% of rows a bootstrap leaves out, which caps Jaccard
      ## near 0.632 even for a perfect clustering -- caught by a test that
      ## expected stability near 1.0 on trivially separated synthetic clusters
      ## and saw 0.63.
      orig_members <- intersect(which(orig_cluster == kk), drawn)
      if (!length(orig_members)) next
      best <- 0
      for (bb in seq_len(k)) {
        bm <- unique(idx[boot_cluster == bb])
        uni <- length(union(orig_members, bm))
        j <- if (uni == 0) 0 else length(intersect(orig_members, bm)) / uni
        if (j > best) best <- j
      }
      jac_sum[kk] <- jac_sum[kk] + best
      jac_n[kk] <- jac_n[kk] + 1L
    }
  }
  ifelse(jac_n > 0, jac_sum / jac_n, NA_real_)
}

#' Cluster observations, choosing the number of clusters
#'
#' Groups the rows of an [ilm_reduce()] result, or any numeric coordinates, by
#' k-means or hierarchical clustering. When `k` is not given it is searched for
#' over `1..k_max` by the gap statistic. Every cluster gets a stability score
#' from bootstrap resampling, every observation gets a silhouette width, and
#' clusters holding less than `small_cluster_frac` of the data are flagged.
#'
#' @section Two different things called anomalous:
#'
#' A **small cluster** is a property of the cluster: a handful of rows that sit
#' apart from everything else. It may be a real minority pattern or it may be a
#' data problem, and either way it is worth looking at.
#'
#' An **ambiguous observation** is a property of the row, taken from its
#' silhouette width -- its average distance to its own cluster against its
#' average distance to the nearest other one. Near 1 means clearly placed, near
#' 0 means sitting on a boundary, negative means it is closer to another cluster
#' than to its own. This is independent of cluster size: a point can be
#' ambiguous inside a large and otherwise stable cluster, which a size-based
#' flag alone would never show. Kaufman and Rousseeuw treat widths below 0.25 as
#' no substantial structure; the default threshold here is stricter at 0.1.
#'
#' @param x An [ilm_reduce()] result, or a numeric matrix or data frame of
#'   coordinates with one row per observation.
#' @param k Number of clusters. Omit to choose it by the gap statistic.
#' @param k_max Largest `k` to consider in that search. Ignored when `k` is
#'   given.
#' @param method `"kmeans"` or `"hclust"`.
#' @param dist_method Distance for `method = "hclust"`. Euclidean by default,
#'   which is right here because the clustering runs on continuous
#'   dimension-reduction coordinates rather than raw mixed-type columns.
#' @param hclust_method Linkage for `method = "hclust"`; `"ward.D2"` by
#'   default.
#' @param nstart Random starts for k-means, guarding against a poor local
#'   optimum.
#' @param B Bootstrap resamples, used both for the gap statistic and for the
#'   stability score.
#' @param gap_method Rule for reading `k` off the gap curve; see
#'   [cluster::maxSE()].
#' @param small_cluster_frac Clusters holding less than this fraction of the
#'   data are flagged.
#' @param ambiguous_threshold Silhouette width below which an observation is
#'   flagged individually.
#' @param seed Random seed.
#' @param progress Show a progress bar. Defaults to [interactive()], so a
#'   bar appears when someone is watching and nothing is written in a
#'   script or a knitted document. See [ilm_progress_arg].
#' @return An object of class `"ilm_cluster"`: `method`, `k`, `gap` (the
#'   search, or `NULL` when `k` was given), `clusters` (one row per cluster,
#'   with `size`, `pct`, `jaccard`, `stability` labelled `"stable"` above 0.75,
#'   `"moderate"` from 0.5 to 0.75 and `"unstable"` below, following Hennig,
#'   plus `mean_silhouette` and `anomalous`), `ind_cluster` (one row per
#'   observation, with `silhouette`, `is_small_cluster`, `is_ambiguous` and
#'   `is_anomalous`), and `coords`.
#' @references
#' Hennig, C. (2007). Cluster-wise assessment of cluster stability.
#' Computational Statistics and Data Analysis 52(1).
#'
#' Tibshirani, R., Walther, G. and Hastie, T. (2001). Estimating the number of
#' clusters in a data set via the gap statistic. JRSS B 63(2).
#' @seealso [ilm_reduce()], [ilm_profile()], [ilm_plot_cluster()].
#' @examples
#' cl <- ilm_cluster(ilm_reduce(mtcars), k_max = 5, B = 25, seed = 1)
#' cl
#' @export
ilm_cluster <- function(x, k = NULL, k_max = 10, method = c("kmeans", "hclust"),
                        dist_method = "euclidean", hclust_method = "ward.D2",
                        nstart = 25, B = 100, gap_method = "firstSEmax",
                        small_cluster_frac = 0.05, ambiguous_threshold = 0.1,
                        seed = NULL, progress = NULL) {
  method <- match.arg(method)
  coords <- if (inherits(x, "ilm_reduce"))
    as.matrix(x$ind_coord[, setdiff(names(x$ind_coord), "row_id"), drop = FALSE])
  else if (is.matrix(x) || is.data.frame(x)) as.matrix(x)
  else stop("`x` must be an ilm_reduce() result or a numeric matrix or data ",
            "frame of coordinates; it is ", class(x)[1], call. = FALSE)
  if (!is.numeric(coords))
    stop("the coordinates to cluster must be numeric", call. = FALSE)
  n <- nrow(coords)
  if (!is.null(seed)) set.seed(seed)
  FUNcluster <- ilm_cluster_funcluster(method, dist_method, hclust_method, nstart)

  gap <- NULL
  if (is.null(k)) {
    ilm_require_cluster()
    ## cap at the number of DISTINCT points, not just n - 1: asking for more
    ## clusters than there are distinct positions is always degenerate, and
    ## low-cardinality input such as missingness indicators can have far fewer
    ## distinct rows than observations
    k_max <- max(1L, min(k_max, n - 1L, nrow(unique(coords))))
    gap <- cluster::clusGap(coords, FUNcluster = FUNcluster, K.max = k_max, B = B)
    k <- cluster::maxSE(gap$Tab[, "gap"], gap$Tab[, "SE.sim"], method = gap_method)
  }

  cluster_assign <- FUNcluster(coords, k)$cluster
  stability <- ilm_cluster_stability(coords, FUNcluster, k, cluster_assign, B,
                                     progress)

  ## the silhouette's distance matches whatever the clustering used: k-means is
  ## implicitly Euclidean, hclust uses what it was told
  sil_dist <- if (method == "hclust") dist_method else "euclidean"
  sil_width <- if (k > 1L) {
    ilm_require_cluster()
    cluster::silhouette(cluster_assign,
                        stats::dist(coords, method = sil_dist))[, "sil_width"]
  } else rep(NA_real_, n)     # undefined for a single cluster

  sizes <- as.integer(table(factor(cluster_assign, levels = seq_len(k))))
  is_small <- sizes < max(1L, ceiling(small_cluster_frac * n))
  mean_sil <- if (k > 1L)
    vapply(seq_len(k), function(kk) mean(sil_width[cluster_assign == kk]), 1)
  else NA_real_

  clusters <- data.frame(
    cluster = seq_len(k), size = sizes, pct = round(100 * sizes / n, 1),
    jaccard = round(stability, 3),
    stability = cut(stability, breaks = c(-Inf, 0.5, 0.75, Inf),
                    labels = c("unstable", "moderate", "stable")),
    mean_silhouette = round(mean_sil, 3), anomalous = is_small,
    stringsAsFactors = FALSE)
  rownames(clusters) <- NULL

  is_ambig <- !is.na(sil_width) & sil_width < ambiguous_threshold
  ind_cluster <- data.frame(
    row_id = seq_len(n), cluster = cluster_assign,
    silhouette = round(sil_width, 3),
    is_small_cluster = is_small[cluster_assign], is_ambiguous = is_ambig,
    is_anomalous = is_small[cluster_assign] | is_ambig,
    stringsAsFactors = FALSE)
  rownames(ind_cluster) <- NULL

  structure(list(method = method, k = k,
                 k_max = if (is.null(gap)) NA_integer_ else k_max, gap = gap,
                 clusters = clusters, ind_cluster = ind_cluster, coords = coords,
                 dist_method = if (method == "hclust") dist_method else NA_character_,
                 hclust_method = if (method == "hclust") hclust_method else NA_character_,
                 small_cluster_frac = small_cluster_frac,
                 ambiguous_threshold = ambiguous_threshold),
            class = "ilm_cluster")
}

#' Cluster observations by which values they are missing
#'
#' The missingness counterpart to [ilm_cluster()]: groups the dimensions that
#' [ilm_reduce_na()] produced, so that rows missing the same columns end up
#' together. A type-checked wrapper -- the clustering is identical, only the
#' input differs.
#'
#' @param x An [ilm_reduce_na()] result.
#' @param ... Passed to [ilm_cluster()].
#' @return An object of class `"ilm_cluster_na"`, which is also an
#'   `"ilm_cluster"`.
#' @seealso [ilm_profile_na()].
#' @examples
#' cl <- ilm_cluster_na(ilm_reduce_na(airquality), k_max = 4, B = 25, seed = 1)
#' cl
#' @export
ilm_cluster_na <- function(x, ...) {
  if (!inherits(x, "ilm_reduce_na"))
    stop("`x` must be an ilm_reduce_na() result; it is ", class(x)[1],
         call. = FALSE)
  out <- ilm_cluster(x, ...)
  class(out) <- c("ilm_cluster_na", class(out))
  out
}

#' @export
print.ilm_cluster <- function(x, ...) {
  tag <- if (inherits(x, "ilm_cluster_na")) "ilm_cluster_na" else "ilm_cluster"
  cat(sprintf("<%s> method = %s, k = %d %s\n\n", tag, x$method, x$k,
              if (!is.null(x$gap)) "(chosen by gap statistic)" else "(as given)"))
  cl <- x$clusters
  cat(sprintf("  %7s %6s %6s %8s %10s %7s\n", "cluster", "size", "pct",
              "jaccard", "stability", "sil"))
  for (i in seq_len(nrow(cl)))
    cat(sprintf("  %7d %6d %6.1f %8.3f %10s %7s%s\n", cl$cluster[i],
                cl$size[i], cl$pct[i], cl$jaccard[i],
                as.character(cl$stability[i]),
                ifelse(is.na(cl$mean_silhouette[i]), "-",
                       sprintf("%.3f", cl$mean_silhouette[i])),
                if (isTRUE(cl$anomalous[i])) "  small" else ""))
  ns <- sum(x$ind_cluster$is_small_cluster)
  if (ns)
    cat(sprintf("\n  %d observation(s) in a cluster holding less than %.0f%% of the data\n",
                ns, 100 * x$small_cluster_frac))
  na_ <- sum(x$ind_cluster$is_ambiguous)
  if (na_)
    cat(sprintf("  %d observation(s) sit between clusters (silhouette below %.2f) and may be misassigned\n",
                na_, x$ambiguous_threshold))
  invisible(x)
}
