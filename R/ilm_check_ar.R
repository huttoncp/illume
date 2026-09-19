## illume: temporal-autocorrelation diagnostics.
##
## The check the original AR experiment lacked.  Residual autocorrelation is
## computed WITHIN group, by lag, on the log-score quantile residuals, and
## compared against an envelope built by simulating from the fit and REFITTING.
## The refit-based reference matters for the same reason as everywhere else: the
## in-sample residuals are not exactly uniform, so a theoretical +/- 2/sqrt(n)
## band (what acf() draws) is the wrong reference and will over-reject.
##
## Note what this can and cannot see.  It detects autocorrelation the model does
## not ACCOUNT for -- so it fires on a model fitted without an AR term to data
## that have one.  It cannot, by construction, tell you an AR term you DID fit is
## mis-estimated: simulating from that fit reproduces its own autocorrelation, so
## the misspecification cancels out of the comparison (the same limitation
## established for omitted covariates).
##
## WHICH RESIDUAL, and why not the log-score one.  Measured separation of lag-1
## residual ACF between data with rho = 0.8 and data with none, model omitting
## AR in both cases:
##   one-vs-rest Pearson  -0.025 -> +0.095   separation  0.119   <- works
##   Pearson (rms)         0.087 ->  0.078              -0.009
##   mid-PIT               0.033 ->  0.024              -0.009
##   randomised quantile   0.021 ->  0.014              -0.007   <- useless here
##   log score             0.083 ->  0.078              -0.005
## The log score is a function of p_hat at the OBSERVED category only, so it is
## direction-blind, and autocorrelation is directional -- both magnitude-style
## statistics fail identically.  Signed one-vs-rest Pearson residuals keep the
## direction and are what this check uses.  (The log-score quantile residual is
## still the right tool for the omnibus QQ and calibration panels.)

#' Signed one-against-the-rest Pearson residuals
#'
#' For each category, the difference between whether that category occurred and
#' its predicted probability, scaled by its standard deviation. Unlike the
#' log-score residual these keep the **direction** of the discrepancy, which is
#' what makes them able to detect autocorrelation.
#'
#' @param object A fitted `"ilm_model"` object.
#' @return A matrix with one column per category.
#' @keywords internal
#' @noRd
ilm_pearson_ovr <- function(object) {
  P <- ilm_fitted(object, TRUE); N <- nrow(P)
  Y <- matrix(0, N, ncol(P)); Y[cbind(seq_len(N), object$y)] <- 1
  (Y - P) / sqrt(pmax(P * (1 - P), 1e-8))
}

#' Within-group autocorrelation of residuals
#'
#' Correlates each residual with the residual a given number of time steps
#' earlier in the same group.
#'
#' @param z Numeric residual vector.
#' @param grp Integer group codes.
#' @param tim Integer time indices.
#' @param maxlag Integer. Largest lag.
#' @param min_pairs Integer. Minimum pairs needed to report a lag.
#' @return A numeric vector of correlations, one per lag.
#' @keywords internal
#' @noRd
ilm_resid_acf <- function(z, grp, tim, maxlag = 8L, min_pairs = 20L) {
  ord <- order(grp, tim)
  g <- grp[ord]; t <- tim[ord]; zz <- z[ord]; n <- length(zz)
  out <- rep(NA_real_, maxlag)
  for (L in seq_len(maxlag)) {
    if (n <= L) break
    i1 <- seq_len(n - L); i2 <- i1 + L
    ok <- g[i1] == g[i2] & (t[i2] - t[i1]) == L
    if (sum(ok) >= min_pairs) {
      a <- zz[i1][ok]; b <- zz[i2][ok]
      if (stats::sd(a) > 0 && stats::sd(b) > 0) out[L] <- stats::cor(a, b)
    }
  }
  out
}

#' Check for leftover temporal autocorrelation
#'
#' Tests whether residuals from observations close together in time within the
#' same group are more similar than the model predicts. If they are, the model
#' is missing temporal structure.
#'
#' @section Which residual, and why it matters:
#' This uses signed one-against-the-rest Pearson residuals rather than the
#' log-score quantile residuals used elsewhere. The reason is worth
#' understanding: the log score depends only on the probability assigned to the
#' category that actually occurred, so it registers the *size* of a discrepancy
#' but not its **direction**. Autocorrelation is directional -- consecutive
#' observations being pushed the same way -- so a direction-blind residual cannot
#' see it. Measured side by side, the log-score residual had essentially no
#' power here while the signed Pearson residual separated the cases cleanly.
#'
#' @section What it cannot tell you:
#' It detects autocorrelation the model does not **account for**, so it fires on
#' a model fitted without an AR term to data that have one. It cannot tell you
#' that an AR term you *did* fit is poorly estimated, because simulating from
#' that fit reproduces its own autocorrelation and the discrepancy cancels out.
#'
#' @section The comparison band:
#' Built by refitting simulated data rather than using the usual
#' `plus or minus 2/sqrt(n)` lines, which assume independent observations and
#' would over-flag here. Because the statistic takes a maximum across categories,
#' its reference is calibrated by leave-one-out to avoid an artificial bias
#' toward significance.
#'
#' @param object A fitted `"ilm_model"` object.
#' @param time Integer time index, one per observation. Must be evenly spaced.
#' @param group Grouping variable, one per observation.
#' @param maxlag Integer. Largest lag to examine.
#' @param B Integer. Simulated datasets. No p-value can fall below `1/(B+1)`, so
#'   `B` of about 100 is needed before the strongest verdict becomes reachable.
#' @param ncores Integer. Worker processes.
#' @param seed Integer. Random seed.
#' @param verbose Logical. Print the table.
#' @param plot Logical. Draw the autocorrelation and its band.
#' @return Invisibly, a list with the per-lag table, per-category z scores and
#'   the simulated null.
#' @references
#' Rue, H., & Held, L. (2005). *Gaussian Markov Random Fields: Theory and
#' Applications*. Chapman & Hall/CRC. (Background on latent autoregressive
#' structures.)
#' @seealso [ilm_appraise()], [ilm_check_omitted()].
#' @export
ilm_check_ar <- function(object, time, group, maxlag = 8L, B = 30L,
                          ncores = 1L, seed = 1L, verbose = TRUE, plot = FALSE) {
  tim <- as.integer(time); grp <- as.integer(factor(group))
  if (length(tim) != nrow(object$X) || length(grp) != nrow(object$X))
    stop("time and group must have one entry per observation")
  acf_all <- function(fit) {
    R <- ilm_pearson_ovr(fit)
    vapply(seq_len(ncol(R)), function(j) ilm_resid_acf(R[, j], grp, tim, maxlag),
           numeric(maxlag))            # maxlag x J
  }
  obs <- acf_all(object)

  ys <- ilm_sim_cond(object, B, seed + 1L)
  X <- object$X; J <- object$J; rl <- ilm_re_list_of(object); rs <- object$re_struct
  arr <- object$ar; yl <- object$ylevels; asg <- object$assign; tl <- object$term_labels
  fm <- object$family; wt <- object$weights
  cl <- ilm_pool(ncores)
  on.exit(if (!is.null(cl)) try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
  if (!is.null(cl))
    try(parallel::clusterExport(cl, c("grp", "tim", "maxlag"), envir = environment()),
        silent = TRUE)
  one <- function(b) {
    f <- try(ilm_fit(X, ys[, b], J, rl, rs, arr, ylevels = yl,
                      weights = wt, family = fm,
                      verbose = FALSE, restarts = 1L), silent = TRUE)
    if (inherits(f, "try-error") || f$opt$convergence != 0) return(rep(NA_real_, maxlag))
    f$assign <- asg; f$term_labels <- tl; f$y <- ys[, b]
    acf_all(f)
  }
  reps <- ilm_lapply(cl, seq_len(B), one)
  arr <- array(unlist(reps), c(maxlag, object$J, length(reps)))
  nok <- sum(!is.na(arr[1, 1, ]))
  mu <- apply(arr, 1:2, mean, na.rm = TRUE); sdv <- apply(arr, 1:2, stats::sd, na.rm = TRUE)
  zmat <- (obs - mu) / sdv
  ## one-vs-rest residuals across categories are linearly dependent, so the J
  ## per-lag tests are not independent; summarise each lag by its worst category
  jbest <- apply(abs(zmat), 1, function(v) if (all(is.na(v))) NA_integer_ else which.max(v))
  z <- vapply(seq_len(maxlag), function(L)
    if (is.na(jbest[L])) NA_real_ else zmat[L, jbest[L]], 0)
  obsv <- vapply(seq_len(maxlag), function(L)
    if (is.na(jbest[L])) NA_real_ else obs[L, jbest[L]], 0)
  muv <- vapply(seq_len(maxlag), function(L)
    if (is.na(jbest[L])) NA_real_ else mu[L, jbest[L]], 0)
  lo <- apply(arr, 1:2, stats::quantile, 0.025, na.rm = TRUE)
  hi <- apply(arr, 1:2, stats::quantile, 0.975, na.rm = TRUE)
  lov <- vapply(seq_len(maxlag), function(L) if (is.na(jbest[L])) NA_real_ else lo[L, jbest[L]], 0)
  hiv <- vapply(seq_len(maxlag), function(L) if (is.na(jbest[L])) NA_real_ else hi[L, jbest[L]], 0)
  obs <- obsv; mu <- muv; lo <- lov; hi <- hiv

  ## z above is a MAX over J categories, so judging it against a normal threshold
  ## over-rejects (an extreme-value statistic against the wrong reference).
  ## Calibrate against the null distribution of the SAME max statistic.
  ##
  ## LEAVE-ONE-OUT IS NOT OPTIONAL HERE.  Standardising a null replicate by a
  ## mean and sd computed from a set that INCLUDES it shrinks its own deviation,
  ## while the observed value -- which is not in that set -- is not shrunk.  The
  ## comparison is then rigged: measured, every lag came back at the p floor,
  ## true positives and false positives alike.
  pmax_lag <- rep(NA_real_, maxlag)
  for (L in seq_len(maxlag)) {
    M <- arr[L, , ]                                   # J x B
    if (is.null(dim(M))) M <- matrix(M, nrow = object$J)
    nb <- ncol(M)
    zz <- matrix(NA_real_, nrow(M), nb)
    for (bb in seq_len(nb)) {
      Mo <- M[, -bb, drop = FALSE]
      mm <- apply(Mo, 1, mean, na.rm = TRUE)
      ss <- pmax(apply(Mo, 1, stats::sd, na.rm = TRUE), 1e-12)
      zz[, bb] <- (M[, bb] - mm) / ss
    }
    nullmax <- apply(abs(zz), 2, max, na.rm = TRUE)
    nullmax <- nullmax[is.finite(nullmax)]
    if (!is.na(z[L]) && length(nullmax) >= 10L)
      pmax_lag[L] <- (1 + sum(nullmax >= abs(z[L]))) / (1 + length(nullmax))
  }
  st <- if (nok < 10L) rep("INCONCLUSIVE", maxlag)
        else ifelse(is.na(pmax_lag), "INCONCLUSIVE",
             ifelse(pmax_lag <= 0.01, "FAIL", ifelse(pmax_lag <= 0.05, "WARN", "OK")))
  tab <- data.frame(lag = seq_len(maxlag),
                    category = object$ylevels[jbest], acf = round(obs, 4),
                    null_mean = round(mu, 4), lo = round(lo, 4), hi = round(hi, 4),
                    maxz = round(z, 2), p = round(pmax_lag, 4), status = st)
  if (verbose) {
    cat(sprintf("\nresidual autocorrelation vs simulated envelope (%d refits)\n", nok))
    print(tab, row.names = FALSE)
    nb <- sum(st %in% c("WARN", "FAIL"))
    ## Lags are correlated with one another, so this is a rough guide, not a
    ## multiple-testing correction.
    if (nb) {
      w <- which(st %in% c("WARN", "FAIL"))
      cat(sprintf("\n>> %d of %d lags outside the envelope (lag %s).\n",
                  nb, maxlag, paste(w, collapse = ", ")))
      cat("   why: residual autocorrelation the model does not account for\n")
      cat(sprintf("   try: %s\n", if (is.null(object$ar))
        "add an AR(1) term over time within group -- but check obs-per-AR-latent first, since one latent per observation is not estimable"
        else "raise the AR order, or coarsen the AR time grid"))
      cat("   note: lags are not independent, so treat the count as indicative\n")
      ## With B replicates no p-value can fall below 1/(B+1), so a FAIL verdict
      ## (p <= 0.01) is simply unreachable for small B -- an enormous z can then
      ## only ever reach WARN.  Say so rather than let it look like weak evidence.
      if (1 / (nok + 1) > 0.01)
        cat(sprintf("   note: p cannot fall below %.3f with %d replicates, so FAIL is unreachable;\n         raise B to about 100 if you need FAIL-level verdicts\n",
                    1 / (nok + 1), nok))
    } else if (nok >= 10L)
      cat("\n>> residual autocorrelation is consistent with the fitted model.\n")
  }
  if (plot) {
    plot(seq_len(maxlag), obs, type = "n", ylim = range(c(obs, lo, hi), na.rm = TRUE),
         xlab = "lag", ylab = "residual ACF", main = "Autocorrelation vs envelope")
    graphics::polygon(c(seq_len(maxlag), rev(seq_len(maxlag))), c(lo, rev(hi)),
                      col = "grey85", border = NA)
    graphics::abline(h = 0, col = "grey50", lty = 2)
    graphics::lines(seq_len(maxlag), mu, col = "grey40")
    graphics::points(seq_len(maxlag), obs, pch = 16,
                     col = ifelse(tab$status == "OK", "black", "red"))
  }
  invisible(list(table = tab, z_by_category = zmat, null = arr, n_ok = nok))
}
