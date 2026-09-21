## ---------------------------------------------------------------------------
## Multivariate anomaly detection from low-rank reconstruction error.
##
## ilm_outliers() asks whether a value is extreme for its own column. This asks
## a different question: whether a ROW is a plausible combination. Someone 150
## cm tall is unremarkable, someone weighing 110 kg is unremarkable, and
## someone who is both is not. Nothing in either column's distribution says so.
##
## The method: most of the variation in a set of correlated columns lies in a
## few directions. Fit those directions, project each row onto them, and
## measure what is left over. A row that respects the correlations is
## reconstructed from k numbers; one that does not is not.
##
## Three things had to be measured rather than assumed, and two of them
## overturned the first design.
##
## THE RANK IS NOT THE IMPUTATION RANK. ilm_impute() chooses it by
## cross-validating held-out cells, which is the right question for filling a
## value in and the wrong one here. On a rank-2 structure in eight columns that
## chose 6 or 7 -- on clean data as well as contaminated -- and the extra
## components spanned the very directions the anomalies departed along.
## Detection fell from 0.975 to 0.560. Parallel analysis chose 2 every time.
##
## THE FIT IS TRIMMED, NOT HELD OUT. The first design scored every row against
## a fit built from other folds of rows, on the reasoning that an anomaly
## inside the fit bends the directions towards itself. The reasoning is sound
## and that remedy did nothing: measured at four contamination levels the two
## agreed to three decimals, because four fifths of the anomalies are still in
## each training fold. Fitting on the rows with the smallest scores and
## iterating does work, because it removes them from the fit rather than a
## fifth of them.
##
## THE REFERENCE IS SIMULATED AND CALIBRATED. The score is not chi-squared: the
## noise scale is estimated, the rank was chosen from the same data, and the
## residual is taken against an estimated subspace. Simulating from the fitted
## structure needs the noise PER COLUMN -- standardising gives every column
## unit variance, so one that loads weakly on the shared directions keeps
## proportionally more of its variance off them -- and needs its level matched
## to the observed median rather than estimated from shrunken in-sample
## residuals. Getting either wrong flagged 38.9% of the rows of data containing
## no anomalies at all.
##
## References:
##   Hawkins, D. M. (1974). The detection of errors in multivariate data using
##     principal components. JASA 69, 340-344.
##   Horn, J. L. (1965). A rationale and test for the number of factors in
##     factor analysis. Psychometrika 30, 179-185.
##   Hubert, M., Rousseeuw, P. J. and Vanden Branden, K. (2005). ROBPCA: a new
##     approach to robust principal component analysis. Technometrics 47, 64-79.
## ---------------------------------------------------------------------------

#' How many directions are real structure?
#'
#' Horn's parallel analysis. Each column is permuted independently, which
#' destroys everything the columns share while leaving each one's own
#' distribution alone, and a component is kept when it beats what the permuted
#' data produces.
#'
#' @keywords internal
#' @noRd
ilm_anom_rank <- function(Z, B = 30L, q = 0.95) {
  dobs <- svd(Z, nu = 0L, nv = 0L)$d
  perm <- vapply(seq_len(B), function(b)
    svd(apply(Z, 2L, sample), nu = 0L, nv = 0L)$d, numeric(length(dobs)))
  thr <- apply(perm, 1L, stats::quantile, q)
  max(1L, min(sum(dobs > thr), ncol(Z) - 1L))
}

#' Score rows against a trimmed low-rank fit
#'
#' The directions are estimated from the rows that fit them best, so a row far
#' off the structure is scored against a structure it did not help define.
#' Every row is then scored, including the trimmed ones.
#'
#' @param Z Standardised numeric matrix with no missing values.
#' @param k Rank.
#' @param trim Share of rows excluded from the fit at each pass; `0` fits all.
#' @param iter Refitting passes.
#' @return A list with the per-row score, the per-cell residuals and the basis.
#' @keywords internal
#' @noRd
ilm_anom_score <- function(Z, k, trim = 0.25, iter = 3L) {
  n <- nrow(Z)
  keep <- seq_len(n)
  nkeep <- max(k + 2L, floor(n * (1 - trim)))
  reps <- if (trim > 0 && nkeep < n) iter else 1L
  sc <- numeric(n); R <- Z; V <- NULL
  for (it in seq_len(reps)) {
    sv <- svd(Z[keep, , drop = FALSE], nu = 0L,
              nv = min(k, length(keep) - 1L))
    V <- sv$v[, seq_len(min(k, ncol(sv$v))), drop = FALSE]
    R <- Z - Z %*% V %*% t(V)
    sc <- rowSums(R^2)
    keep <- order(sc)[seq_len(nkeep)]
  }
  list(score = sc, residual = R, V = V)
}

#' Rows that do not fit the pattern the other rows make
#'
#' Finds observations that are implausible as a **combination** of values, even
#' when no single value is extreme. [ilm_outliers()] asks whether a number is
#' far out in its own column; this asks whether a row is far from the structure
#' the columns share.
#'
#' A set of correlated columns puts most of its variation in a few directions.
#' Those directions are estimated, each row is projected onto them, and the
#' score is what is left over.
#'
#' @section How many directions:
#'
#' By parallel analysis, and deliberately not by the cross-validation
#' [ilm_impute()] uses on the same decomposition, because the two answer
#' different questions. Imputation wants the rank that best predicts a missing
#' cell. This wants the number of directions that are real shared structure, so
#' that what is left over is residual rather than signal it failed to fit.
#'
#' The difference is large. On a rank-2 structure in eight columns,
#' cross-validation chose 6 or 7 -- on clean data as well as contaminated --
#' and those extra components span the directions the anomalies depart along.
#' Detection fell from 0.975 to 0.560. Parallel analysis chose 2 every time.
#'
#' @section Why the fit is trimmed:
#'
#' The anomalies sit in the same data the directions are estimated from, so
#' they pull the directions towards themselves and are then reconstructed well.
#' `trim` excludes the worst-fitting rows from the fit and refits, so a row is
#' scored against a structure it did not help define. Every row is still
#' scored, the trimmed ones included.
#'
#' It earns its place only where the anomalies share a direction, which is
#' where they can form a component between them. Detection with `trim = 0`
#' against the default 0.25, over 12 datasets each:
#'
#' ```
#'   anomalies along random directions      along ONE shared direction
#'    2%   0.925  vs  0.917                  2%   0.950  vs  0.942
#'    5%   0.967  vs  0.977                  5%   0.797  vs  0.887
#'   10%   0.970  vs  0.978                 10%   0.587  vs  0.670
#' ```
#'
#' Holding rows out in folds instead, which was the first design here, does not
#' work and is not offered: measured at four contamination levels it matched
#' in-sample scoring to three decimals, because most of the anomalies remain in
#' every training fold.
#'
#' @section What it cannot do:
#'
#' When a large share of rows depart along the **same** direction they are not
#' anomalies, they are a subpopulation, and a rank-k fit of the whole data
#' legitimately includes their direction. Detection degrades accordingly: with
#' anomalies sharing one direction, 0.89 of them were found at 5%
#' contamination, 0.67 at 10% and 0.41 at 20%. That is the method reaching its
#' limit rather than failing quietly, and [ilm_cluster()] is the tool for a
#' second group, since finding one is what it is for.
#'
#' @section What the score is compared with:
#'
#' Not a chi-squared distribution. Datasets with the same structure and no
#' anomalies are simulated and scored the same way. `B` can be modest because
#' each simulated dataset contributes `n` null scores rather than one.
#'
#' Measured on 100 matrices of 400 by 8 with a rank-2 structure and **no
#' anomalies at all**: 11 rows of 40,000 were flagged, a rate of 0.00028, with
#' 95 of the 100 datasets producing none. With anomalies present, the rows that
#' were not planted were flagged at 0.0015 to 0.0030 across every design tried
#' above.
#'
#' `flag` uses `p_adj` and is what those numbers describe. The raw `p` runs a
#' little hot -- about 0.064 of clean rows fall below 0.05 rather than 0.05 of
#' them -- so read it as a ranking rather than as a test.
#'
#' @param data A data frame or matrix.
#' @param cols Columns to use; see [ilm_selection]. Numeric columns only --
#'   the reconstruction is a projection, and a category has no residual along a
#'   direction. Anything else is dropped with a note.
#' @param rank Number of directions; `NULL` uses parallel analysis.
#' @param trim Share of the worst-fitting rows held out of the fit. `0` fits
#'   every row, which lets the anomalies define the structure they are scored
#'   against.
#' @param B Simulated null datasets for the reference. A Benjamini-Hochberg
#'   adjusted p-value cannot fall much below `1 / B`, so a single anomaly among
#'   many rows needs `B` comfortably above `1 / alpha`.
#' @param alpha Flagging level, after a Benjamini-Hochberg adjustment across
#'   rows -- one row in twenty at 0.05 would be 50 rows in a thousand, which is
#'   a list nobody reads.
#' @param seed Random seed.
#' @param progress Show a progress bar; see [ilm_progress_arg].
#' @return A data frame with one row per observation: `row`, `score`, `p`,
#'   `p_adj`, `flag`, and `driver`, the column contributing most to the score.
#'   The rank and the residual matrix are attributes.
#' @seealso [ilm_outliers()] for the one-column-at-a-time question,
#'   [ilm_cluster()] when the unusual rows turn out to be a group,
#'   [ilm_impute()] which fits the same decomposition to fill values in.
#' @references Hawkins, D. M. (1974). The detection of errors in multivariate
#'   data using principal components. *Journal of the American Statistical
#'   Association* 69, 340-344.
#'
#'   Horn, J. L. (1965). A rationale and test for the number of factors in
#'   factor analysis. *Psychometrika* 30, 179-185.
#' @examples
#' set.seed(1)
#' n <- 300
#' f <- rnorm(n)
#' d <- data.frame(a = f + rnorm(n, 0, .3), b = 2 * f + rnorm(n, 0, .3),
#'                 c = -f + rnorm(n, 0, .3))
#' ## a row that is ordinary in every column but breaks the pattern
#' d[1, ] <- c(1.2, -2.4, 1.2)
#' head(ilm_anomaly(d), 3)
#' @export
ilm_anomaly <- function(data, cols = NULL, rank = NULL, trim = 0.25,
                        B = 39L, alpha = 0.05, seed = 1L, progress = NULL) {
  if (is.matrix(data)) data <- as.data.frame(data)
  if (!is.data.frame(data))
    stop("`data` must be a data frame or matrix, not ", class(data)[1],
         call. = FALSE)
  if (!is.numeric(trim) || length(trim) != 1L || trim < 0 || trim >= 0.5)
    stop("`trim` must be a single number in [0, 0.5): it is the share of rows ",
         "held out of the FIT, and trimming half of them leaves the structure ",
         "defined by whichever half happened to fit first.", call. = FALSE)
  sel <- ilm_resolve_cols(data, cols)
  num <- sel[vapply(sel, function(v) is.numeric(data[[v]]), TRUE)]
  drop <- setdiff(sel, num)
  if (length(drop))
    message("ilm_anomaly(): using the numeric columns only; a reconstruction ",
            "is a projection and a category has no residual along a ",
            "direction. Dropped: ", paste(drop, collapse = ", "),
            ". ilm_reduce() handles mixed types, and a cluster far from every ",
            "centre there is the categorical analogue of this.")
  if (length(num) < 3L)
    stop("at least 3 numeric columns are needed: with two there is one ",
         "direction to fit and one left over, which is a scatterplot rather ",
         "than a multivariate question.", call. = FALSE)

  X <- as.matrix(data[num])
  keep <- stats::complete.cases(X)
  if (sum(keep) < 20L)
    stop("only ", sum(keep), " complete rows across these columns. Fill them ",
         "in with ilm_impute() first, or choose columns with fewer gaps -- ",
         "scoring a row against a structure fitted to a fifth of the data is ",
         "not worth doing.", call. = FALSE)
  Xc <- X[keep, , drop = FALSE]
  n <- nrow(Xc); p <- ncol(Xc)
  if (B * alpha <= 1)
    warning("B = ", B, " puts the smallest achievable adjusted p-value at ",
            "about ", signif(1 / B, 2), ", which is not below alpha = ", alpha,
            ", so a lone anomaly among many rows cannot be flagged however ",
            "extreme it is. Use B >= ", ceiling(2 / alpha), ".", call. = FALSE)

  ctr <- colMeans(Xc)
  scl <- apply(Xc, 2L, stats::sd)
  scl[!is.finite(scl) | scl <= 0] <- 1
  Z <- sweep(sweep(Xc, 2L, ctr, "-"), 2L, scl, "/")

  set.seed(seed)
  k <- if (is.null(rank)) ilm_anom_rank(Z) else as.integer(rank)
  k <- max(1L, min(k, p - 1L))
  obs <- ilm_anom_score(Z, k, trim)

  ## The null: the same structure, the same noise, no anomalies.
  sv <- svd(Z, nu = k, nv = k)
  V <- sv$v[, seq_len(k), drop = FALSE]
  d <- sv$d[seq_len(k)]
  Rin <- Z - Z %*% V %*% t(V)
  dfc <- (n * p) / max(n * p - (n + p) * k, 1)
  rvar <- apply(Rin, 2L, function(z) stats::mad(z)^2) * dfc
  rvar[!is.finite(rvar) | rvar <= 0] <- .Machine$double.eps

  sim <- function(rv, reps) {
    unlist(lapply(seq_len(reps), function(b) {
      U <- matrix(stats::rnorm(n * k), n, k)
      E <- vapply(seq_len(p),
                  function(j) stats::rnorm(n, 0, sqrt(rv[j])), numeric(n))
      ## standardised and scored exactly as the observed matrix was, so the
      ## two are the same pipeline and not merely similar ones
      ilm_anom_score(scale(U %*% diag(d / sqrt(n), k, k) %*% t(V) + E),
                     k, trim)$score
    }))
  }
  pb <- ilm_progress(B + 3L, progress, "simulating the reference")
  cal <- sim(rvar, 3L); pb$tick(3L)
  mo <- stats::median(obs$score); mc <- stats::median(cal)
  if (is.finite(mc) && mc > 0) rvar <- rvar * (mo / mc)
  null <- numeric(0)
  for (b in seq_len(B)) { null <- c(null, sim(rvar, 1L)); pb$tick(3L + b) }
  pb$done()

  pv <- (1 + vapply(obs$score, function(s) sum(null >= s), 0L)) /
    (length(null) + 1)
  padj <- stats::p.adjust(pv, "BH")
  drv <- colnames(Z)[max.col(abs(obs$residual), ties.method = "first")]

  out <- data.frame(row = which(keep), score = obs$score, p = pv,
                    p_adj = padj, flag = padj < alpha, driver = drv,
                    stringsAsFactors = FALSE, row.names = NULL)
  out <- out[order(-out$score), , drop = FALSE]
  structure(out, class = c("ilm_anomaly", "data.frame"), rank = k,
            residual = obs$residual, columns = num, trim = trim,
            n_null = length(null), alpha = alpha, dropped = drop)
}

#' @export
print.ilm_anomaly <- function(x, n = 10L, ...) {
  d <- as.data.frame(x)
  class(d) <- "data.frame"
  nf <- sum(d$flag)
  cat(sprintf("Multivariate anomalies: %d of %d rows flagged at %.3g\n",
              nf, nrow(d), attr(x, "alpha")))
  cat(sprintf("  rank %d over %d columns, %s, %d null scores\n",
              attr(x, "rank"), length(attr(x, "columns")),
              if (attr(x, "trim") > 0)
                sprintf("fitted on the best %.0f%% of rows",
                        100 * (1 - attr(x, "trim")))
              else "fitted on every row",
              attr(x, "n_null")))
  if (attr(x, "trim") <= 0)
    cat("  With trim = 0 the anomalous rows help define the structure they\n",
        "  are then scored against, and hide themselves in it.\n", sep = "")
  show <- utils::head(d, max(n, nf))
  show$score <- round(show$score, 3)
  show$p <- signif(show$p, 3); show$p_adj <- signif(show$p_adj, 3)
  print(show, row.names = FALSE)
  if (nrow(d) > nrow(show))
    cat(sprintf("  ... %d more rows\n", nrow(d) - nrow(show)))
  if (nf)
    cat("\n  `driver` is the column contributing most to each row's score.\n",
        "  A flagged row is a combination the other rows do not make. It is\n",
        "  not necessarily an error, and deleting it because a method said so\n",
        "  is how real effects get removed.\n", sep = "")
  invisible(x)
}
