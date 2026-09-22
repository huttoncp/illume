## ---------------------------------------------------------------------------
## Which variables are carrying the clustering, and which are only adding
## distance.
##
## WHY THIS EXISTS: neither FAMD nor any distance in this package does variable
## selection, and an irrelevant variable is not neutral -- it contributes its
## full share of the geometry and none of the structure. Measured on data with
## three known clusters and two informative numeric columns, adding two PURE
## NOISE factors took recovery from 0.301 to 0.006. That is not degradation, it
## is annihilation, and it is the largest single effect measured anywhere in
## this comparison: larger than the choice between FAMD and Gower, larger than
## k-means against PAM, larger than the choice of k.
##
## A mixture model that selects variables (VarSelLCM) survives that case --
## 0.397 where FAMD scored 0.006 -- but it assumes the variables are
## independent within a cluster, and when that fails it collapses the other
## way: at a within-cluster correlation of 0.75 it scored 0.047 against FAMD's
## 0.265. So swapping the estimator trades one failure for another.
##
## Measuring the variables instead does not. This reports what each one
## contributes to the separation that was found, so a variable that is only
## adding distance can be seen and dropped -- which is the package's usual
## bargain: name the problem and name a remedy that exists here.
## ---------------------------------------------------------------------------

#' What each variable contributes to a clustering
#'
#' Reports, for every variable, how strongly it separates the clusters that
#' were found. A variable near the bottom is contributing distance without
#' contributing structure, and dropping it usually improves the result rather
#' than losing information.
#'
#' @details
#' The statistic is the **between-cluster share of variance** -- eta squared
#' for a numeric variable, Cramer's V for a categorical one -- so both are on
#' a 0 to 1 scale and directly comparable. A permutation reference is computed
#' by shuffling the cluster labels, which preserves each variable's own
#' distribution and destroys only its relationship with the clustering, and the
#' variable is called uninformative when it does not beat that.
#'
#' Neither FAMD nor any distance used here selects variables, and an
#' uninformative one is not harmless: on three known clusters with two
#' informative columns, adding two pure-noise factors took recovery from 0.301
#' to 0.006. Checking is cheap and the remedy is `cols`.
#'
#' @section What this cannot tell you:
#'
#' The clustering was fitted to these variables, so a variable the clustering
#' **used** will separate the clusters it helped make, whether or not it means
#' anything. On three real clusters plus one pure-noise numeric column, the
#' noise column scored 0.441 against the informative columns' 0.576 and 0.550
#' -- visibly weaker, and still "better than chance", because k-means had split
#' on it.
#'
#' So this ranks contribution; it does not prove relevance. Two readings are
#' safe and are what it is for: a variable the clustering **ignored** (no
#' better than a shuffled label) is adding distance and nothing else, and a
#' variable that separates the clusters almost perfectly while nothing else
#' does has not explained the clustering -- it *is* the clustering, and the
#' partition is its levels under another name. That second case is reported
#' separately because from inside a clustering it looks like the best possible
#' result.
#'
#' Settling relevance properly means refitting without each variable and
#' seeing whether the partition survives, which costs one clustering per
#' variable and is not done here.
#'
#' @param x An [ilm_profile()] result, or an [ilm_cluster()] result that
#'   carries its reduction.
#' @param data The data frame the clustering came from.
#' @param B Permutations for the reference.
#' @param seed Random seed.
#' @return A data frame, one row per variable, with the separation statistic,
#'   the permutation mean, a p-value and a verdict. Sorted strongest first.
#' @seealso [ilm_profile()], [ilm_cluster()], [ilm_reduce()].
#' @examples
#' \donttest{
#' p <- ilm_profile(mtcars, k = 3)
#' ilm_var_contrib(p, mtcars)
#' }
#' @export
ilm_var_contrib <- function(x, data, B = 199L, seed = 1L) {
  cl <- ilm_var_contrib_labels(x)
  if (!is.data.frame(data))
    stop("`data` must be the data frame the clustering came from; it is ",
         class(data)[1], call. = FALSE)
  if (length(cl) != nrow(data)) {
    ## the clustering may have dropped incomplete rows, and lining the labels
    ## up with the wrong rows would score every variable against noise
    stop("the clustering has ", length(cl), " labels and `data` has ",
         nrow(data), " rows. Pass the same data frame the clustering was ",
         "built from; if rows were dropped for missing values, ",
         "ilm_impute() keeps them and the two then match.", call. = FALSE)
  }
  use <- names(data)[vapply(data, function(v)
    is.numeric(v) || is.factor(v) || is.character(v) || is.logical(v), TRUE)]
  if (!length(use)) stop("no usable columns in `data`", call. = FALSE)
  set.seed(seed)
  res <- do.call(rbind, lapply(use, function(v) {
    obs <- ilm_var_sep(data[[v]], cl)
    nul <- vapply(seq_len(B), function(b)
      ilm_var_sep(data[[v]], sample(cl)), 0)
    p <- (1 + sum(nul >= obs)) / (B + 1)
    data.frame(variable = v,
               type = if (is.numeric(data[[v]])) "numeric" else "categorical",
               separation = obs, permuted = mean(nul), p = p,
               stringsAsFactors = FALSE)
  }))
  res$p_adj <- stats::p.adjust(res$p, "BH")
  res$verdict <- ifelse(res$separation >= 0.9 & res$p_adj < 0.05,
                        "defines the clusters",
                 ifelse(res$p_adj >= 0.05, "no better than chance",
                 ifelse(res$separation >= 0.15, "carries the clustering",
                        "weak")))
  res <- res[order(-res$separation), , drop = FALSE]
  row.names(res) <- NULL
  ## A variable that separates the clusters PERFECTLY while nothing else does
  ## has not explained the clustering, it IS the clustering: the partition is
  ## that variable's levels wearing a different name. This is the dangerous
  ## case and it is not obvious from the table, because from inside a
  ## clustering a variable that defines it looks like the best variable.
  dom <- nrow(res) > 1L && !is.na(res$separation[1]) &&
    res$separation[1] >= 0.9 &&
    (is.na(res$separation[2]) || res$separation[2] < 0.25)
  structure(res, class = c("ilm_var_contrib", "data.frame"), B = B,
            k = length(unique(cl)), dominated_by = if (dom) res$variable[1])
}

## Between-cluster share of variance, on one 0-1 scale for both types so the
## two can be ranked against each other: eta squared for a numeric variable,
## Cramer's V for a categorical one.
#' @keywords internal
#' @noRd
ilm_var_sep <- function(v, cl) {
  ok <- !is.na(v) & !is.na(cl)
  if (sum(ok) < 3L || length(unique(cl[ok])) < 2L) return(NA_real_)
  v <- v[ok]; cl <- factor(cl[ok])
  if (is.numeric(v)) {
    gm <- mean(v)
    ssb <- sum(vapply(split(v, cl), function(z) length(z) * (mean(z) - gm)^2, 0))
    sst <- sum((v - gm)^2)
    if (!is.finite(sst) || sst <= 0) return(NA_real_)
    return(ssb / sst)
  }
  tb <- table(factor(v), cl)
  if (any(dim(tb) < 2L)) return(NA_real_)
  cs <- suppressWarnings(stats::chisq.test(tb, correct = FALSE)$statistic)
  sqrt(as.numeric(cs) / (sum(tb) * (min(dim(tb)) - 1L)))
}

#' @keywords internal
#' @noRd
ilm_var_contrib_labels <- function(x) {
  if (inherits(x, "ilm_profile")) return(x$cluster$ind_cluster$cluster)
  if (inherits(x, "ilm_cluster")) return(x$ind_cluster$cluster)
  stop("`x` must be an ilm_profile() or ilm_cluster() result; it is ",
       class(x)[1], call. = FALSE)
}

#' @export
print.ilm_var_contrib <- function(x, ...) {
  d <- as.data.frame(x); class(d) <- "data.frame"
  cat(sprintf("Variable contribution to %d clusters (%d permutations)\n",
              attr(x, "k"), attr(x, "B")))
  d$separation <- round(d$separation, 3)
  d$permuted <- round(d$permuted, 3)
  d$p <- signif(d$p, 3); d$p_adj <- signif(d$p_adj, 3)
  print(d, row.names = FALSE)
  dom <- attr(x, "dominated_by")
  dead <- d$variable[d$verdict == "no better than chance"]
  ## Order matters here. When one variable defines the partition, every OTHER
  ## variable scores low against it however informative it is -- so "drop the
  ## low ones" is exactly the wrong advice in precisely the case where the
  ## table looks most decisive. It has to be pre-empted, not printed beside.
  if (!is.null(dom)) {
    cat("\n  This clustering is a re-labelling of `", dom, "`.\n",
        "  It separates the clusters almost perfectly while nothing else\n",
        "  does, which means the partition found IS its levels. Every other\n",
        "  variable then scores low against it however informative it is, so\n",
        "  do NOT read the bottom of this table as a list to drop.\n",
        "  If `", dom, "` is not the grouping you were looking for, exclude\n",
        "  it with cols = and cluster again.\n", sep = "")
  } else if (length(dead)) {
    cat("\n  ", paste(dead, collapse = ", "), " separate",
        if (length(dead) == 1L) "s" else "",
        " the clusters no better than a shuffled label does.\n",
        "  Nothing here selects variables, so ",
        if (length(dead) == 1L) "it is" else "they are",
        " still contributing distance.\n",
        "  Refit with cols = to drop ", if (length(dead) == 1L) "it" else "them",
        ": on known clusters, two irrelevant factors\n",
        "  took recovery from 0.301 to 0.006.\n", sep = "")
  } else {
    cat("\n  Every variable separates the clusters better than chance.\n",
        "  Read the ranking, not the verdict: the clustering was fitted to\n",
        "  these variables, so one it USED separates the clusters it helped\n",
        "  make whether or not it means anything. A pure-noise column scored\n",
        "  0.441 here against 0.576 for a real one.\n", sep = "")
  }
  invisible(x)
}
