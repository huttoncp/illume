## illume: parametric-bootstrap likelihood-ratio tests.
##
## WHY.  There is no Kenward-Roger or Satterthwaite correction for a multinomial
## GLMM, so the Wald and asymptotic-LR p-values are purely asymptotic.  Two
## things make that uncomfortable here: the Laplace approximation, and the
## Hauck-Donner effect (measured earlier: Wald 586 vs LR 1114 on the same term).
## The parametric bootstrap replaces the chi-square reference with the LR
## statistic's ACTUAL null distribution under this model and this design.
##
## Procedure (Halekoh & Hojsgaard):
##   1. fit the full model M1 and the reduced model M0 (term dropped)
##   2. observed statistic LR = 2 * (logLik(M1) - logLik(M0))
##   3. simulate B datasets FROM M0 -- the null is what generates the reference
##   4. refit BOTH models on each, recording LR_b
##   5. p = (1 + #{LR_b >= LR}) / (1 + #usable)
## The +1 in numerator and denominator keeps p strictly positive: with B
## replicates you cannot resolve a p-value below 1/(B+1), and reporting 0 would
## overstate the evidence.

#' Parametric-bootstrap likelihood-ratio test
#'
#' Tests one fixed-effect term by comparing the observed likelihood-ratio
#' statistic to its **simulated** null distribution, instead of to a chi-square
#' approximation.
#'
#' @section Why bother:
#' The usual p-value for a likelihood-ratio test assumes the statistic follows a
#' chi-square distribution, which is a large-sample result. Linear mixed models
#' have corrections for this (Kenward-Roger, Satterthwaite) but **no such
#' correction exists for a multinomial GLMM**, so the chi-square reference is all
#' the asymptotic route offers.
#'
#' How wrong can it be? On a design with 22 groups, the test that should reject a
#' true null 5% of the time rejected it 8.5% of the time. The bootstrap replaces
#' the assumed reference with the statistic's actual behaviour under this model
#' and this design.
#'
#' @section How it works:
#' Fit the full model and the model without the term. Simulate many datasets
#' **from the reduced model**, because that is what the null hypothesis says
#' generated the data. Refit both models to each simulated dataset to build the
#' null distribution, then see where the observed statistic falls.
#'
#' @section Reading the output:
#' With `B` replicates no p-value can fall below `1/(B+1)`, so the smallest
#' reportable value at `B = 200` is about 0.005. The function says when a result
#' has hit that floor, rather than implying more precision than the simulation
#' can deliver. It also reports the implied type-I error of the asymptotic test,
#' which tells you directly whether the chi-square shortcut would have been
#' adequate for your data.
#'
#' @param object A fitted `"ilm_model"` object, fitted through the formula interface.
#' @param term Character name or integer index of the term to test.
#' @param B Integer. Bootstrap replicates. 200 or more is advisable; each costs
#'   two model fits.
#' @param ncores Integer. Worker processes.
#' @param seed Integer. Random seed.
#' @param restarts Integer. Optimiser restarts within replicates.
#' @param verbose Logical. Print the result.
#' @return Invisibly, a list with `LR`, `df`, `p_boot`, `p_chisq`, the simulated
#'   null distribution, and calibration diagnostics.
#' @references
#' Halekoh, U., & Hojsgaard, S. (2014). A Kenward-Roger approximation and
#' parametric bootstrap methods for tests in linear mixed models: the R package
#' pbkrtest. *Journal of Statistical Software*, 59(9), 1--30.
#' @seealso [ilm_anova()].
#' @export
ilm_pb_lrt <- function(object, term, B = 200L, ncores = 1L, seed = 1L,
                        restarts = 1L, verbose = TRUE) {
  if (is.null(object$assign))
    stop("ilm_pb_lrt() needs a model fitted through the formula interface")
  ilm_stop_reml_lrt(object, "ilm_pb_lrt()")
  j <- if (is.numeric(term)) as.integer(term) else which(object$term_labels == term)
  if (!length(j) || is.na(j)) stop("term '", term, "' not found")
  lab <- object$term_labels[j]
  df <- sum(!is.na(object$assign) & object$assign == j) * object$C
  if (!df) stop("term '", lab, "' contributes no columns")

  t0 <- proc.time()[3]
  f0 <- ilm_refit_drop(object, j, restarts = 2L)
  lr_obs <- 2 * ((-object$opt$objective) - (-f0$opt$objective))
  if (lr_obs < 0) warning("observed LR is negative: a refit did not converge", call. = FALSE)

  ## simulate under the NULL
  ys <- ilm_simulate(f0, B, seed)

  keep <- is.na(object$assign) | object$assign != j
  X1 <- object$X; X0 <- object$X[, keep, drop = FALSE]
  J <- object$J; rl <- ilm_re_list_of(object); rs <- object$re_struct
  arr <- object$ar; yl <- object$ylevels
  fm <- object$family; wt <- object$weights
  cl <- ilm_pool(ncores)
  on.exit(if (!is.null(cl)) try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
  one <- function(b) {
    yb <- ys[, b]
    a <- try(ilm_fit(X1, yb, J, rl, rs, arr, ylevels = yl, weights = wt,
                      family = fm, verbose = FALSE,
                      restarts = restarts), silent = TRUE)
    if (inherits(a, "try-error") || a$opt$convergence != 0) return(NA_real_)
    z <- try(ilm_fit(X0, yb, J, rl, rs, arr, ylevels = yl, weights = wt,
                      family = fm, verbose = FALSE,
                      restarts = restarts), silent = TRUE)
    if (inherits(z, "try-error") || z$opt$convergence != 0) return(NA_real_)
    2 * ((-a$opt$objective) - (-z$opt$objective))
  }
  lr_null <- unlist(ilm_lapply(cl, seq_len(B), one))
  ok <- is.finite(lr_null) & lr_null >= 0
  nok <- sum(ok)
  p_pb <- (1 + sum(lr_null[ok] >= lr_obs)) / (1 + nok)
  p_chisq <- stats::pchisq(lr_obs, df, lower.tail = FALSE)
  el <- proc.time()[3] - t0

  ## Calibration of the asymptotic reference.  Judge it in the UPPER TAIL, not at
  ## the mean: p-values are governed by the tail, and the two can disagree.
  ## Measured on a 22-group design, the null LR mean ratio was 0.94-1.13 (looks
  ## fine) while the 95th-percentile ratio was 1.21 and the actual type-I error
  ## at nominal 0.05 was 0.085.  A mean-based check would have missed that.
  mnull <- mean(lr_null[ok]); ratio <- mnull / df
  q95 <- unname(stats::quantile(lr_null[ok], 0.95))
  tail_ratio <- q95 / stats::qchisq(0.95, df)
  size05 <- mean(lr_null[ok] > stats::qchisq(0.95, df))
  out <- list(term = lab, df = df, LR = lr_obs, p_boot = p_pb, p_chisq = p_chisq,
              n_ok = nok, B = B, null = lr_null[ok], mean_null = mnull,
              ratio = ratio, tail_ratio = tail_ratio, size_05 = size05,
              seconds = el)
  if (verbose) {
    cat(sprintf("\nparametric-bootstrap LRT: %s\n", lab))
    cat(sprintf("  LR = %.3f on %d df | %d of %d replicates usable | %.0f s\n",
                lr_obs, df, nok, B, el))
    cat(sprintf("  p (bootstrap) = %.4f%s | p (asymptotic chi-square) = %.4g\n",
                p_pb, if (p_pb <= 1 / (nok + 1)) sprintf(" (floor: < %.4f)", 1 / (nok + 1)) else "",
                p_chisq))
    cat(sprintf("  null LR: mean %.3f (df %d) | 95th pct %.3f vs chi-square %.3f (ratio %.2f)\n",
                mnull, df, q95, stats::qchisq(0.95, df), tail_ratio))
    cat(sprintf("  implied type-I error of the asymptotic test at nominal 0.05: %.3f%s\n",
                size05,
                if (size05 > 0.075 || size05 < 0.03)
                  "  <- asymptotic chi-square is MISCALIBRATED here; prefer p_boot" else ""))
    if (nok < 0.8 * B)
      cat(sprintf("  >> WARNING: %d of %d replicates failed to refit; the bootstrap\n     reference is conditioned on convergence and may be optimistic.\n",
                  B - nok, B))
  }
  invisible(out)
}
