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
  fam <- if (!is.null(object$family)) object$family$name else "gaussian"
  P <- ilm_fitted(object, TRUE); N <- nrow(P)
  ## For a univariate family the Pearson residual is the ordinary one: the
  ## deviation scaled by the standard deviation the family implies at that
  ## mean. Only the multinomial needs the one-vs-rest indicator construction.
  if (!identical(fam, "multinomial")) {
    mu <- as.numeric(P[, 1]); y <- as.numeric(object$y)
    d <- ilm_disp_vec(object)
    if (is.null(d)) d <- numeric(0)
    ## For an accelerated failure time model the residual is the normal score
    ## of the quantile residual. The obvious alternative -- the error on the
    ## log-time scale, (log t - eta) / scale -- is centred and unit-scaled only
    ## when nothing is censored: a censored row is evaluated at the censoring
    ## time rather than at the event, so it sits systematically low. Measured
    ## at 40% censoring that pulled the mean to -0.64, -0.84 and -0.56 for the
    ## three families. The normal score fills the censoring interval in and is
    ## centred and unit-scaled whatever the censoring.
    if (!is.null(object$family$cdf))
      return(matrix(stats::qnorm(pmin(pmax(ilm_rqr(object, TRUE, 1L), 1e-8),
                                      1 - 1e-8)), ncol = 1L))
    w <- if (is.null(object$weights)) rep(1, N) else pmax(1, round(object$weights))
    v <- switch(fam,
      gaussian = if (length(d)) unname(d)^2 else rep(1, N),
      poisson  = mu,
      nbinom   = mu + mu^2 / unname(d),
      binomial = mu * (1 - mu) / w,
      stop("no Pearson residual is defined for family ", fam, call. = FALSE))
    return(matrix((y - mu) / sqrt(pmax(v, 1e-8)), ncol = 1L))
  }
  Y <- matrix(0, N, ncol(P)); Y[cbind(seq_len(N), as.integer(object$y))] <- 1
  (Y - P) / sqrt(pmax(P * (1 - P), 1e-8))
}

#' Within-group autocorrelation of residuals
#'
#' Correlates each residual with the residual a given number of time steps
#' earlier in the same group. Pairs that straddle a group boundary, or that are
#' not exactly `L` time units apart, are dropped rather than quietly treated as
#' adjacent.
#'
#' @param z Numeric residual vector.
#' @param grp Integer group codes.
#' @param tim Integer time indices.
#' @param maxlag Integer. Largest lag.
#' @param min_pairs Integer. Minimum pairs needed to report a lag.
#' @return A numeric vector of correlations, one per lag, carrying an
#'   `"n_pairs"` attribute giving how many pairs each rests on.
#' @keywords internal
#' @noRd
ilm_resid_acf <- function(z, grp, tim, maxlag = 8L, min_pairs = 20L) {
  ord <- order(grp, tim)
  g <- grp[ord]; t <- tim[ord]; zz <- z[ord]; n <- length(zz)
  out <- rep(NA_real_, maxlag); np <- rep(0L, maxlag)
  for (L in seq_len(maxlag)) {
    if (n <= L) break
    i1 <- seq_len(n - L); i2 <- i1 + L
    ok <- g[i1] == g[i2] & (t[i2] - t[i1]) == L
    np[L] <- sum(ok)
    if (np[L] >= min_pairs) {
      a <- zz[i1][ok]; b <- zz[i2][ok]
      if (stats::sd(a) > 0 && stats::sd(b) > 0) out[L] <- stats::cor(a, b)
    }
  }
  attr(out, "n_pairs") <- np
  out
}

#' Partial autocorrelation from an autocorrelation sequence
#'
#' Durbin-Levinson recursion. The partial autocorrelation at lag k is what is
#' left at that lag once lags 1 to k-1 have been accounted for, which is why the
#' two panels answer different questions: for an AR(p) process the PACF cuts off
#' after lag p while the ACF decays, and for a moving-average process it is the
#' other way round.
#'
#' A sample autocorrelation sequence need not be a valid (positive definite)
#' one, and the recursion diverges once it is not. Rather than report the
#' divergence as a statistic, the recursion stops there and the remaining lags
#' come back `NA`.
#'
#' @param r Numeric vector of autocorrelations at lags 1, 2, and so on.
#' @return A numeric vector of partial autocorrelations of the same length.
#' @keywords internal
#' @noRd
ilm_pacf_dl <- function(r) {
  m <- length(r); phi <- rep(NA_real_, m)
  if (!m || !is.finite(r[1])) return(phi)
  phi[1] <- r[1]
  if (abs(r[1]) >= 1) return(phi)
  prev <- r[1]
  for (k in seq_len(m)[-1]) {
    if (!is.finite(r[k])) break
    back <- rev(seq_len(k - 1L))
    den <- 1 - sum(prev * r[seq_len(k - 1L)])
    if (!is.finite(den) || abs(den) < 1e-10) break
    pk <- (r[k] - sum(prev * r[back])) / den
    if (!is.finite(pk)) break
    phi[k] <- pk
    if (abs(pk) >= 1) break
    prev <- c(prev - pk * prev[back], pk)
  }
  phi
}

#' Per-lag verdict for a statistic computed per residual component
#'
#' The value shown at each lag is the worst of the `C` residual components,
#' which is a maximum over correlated tests. Judging a maximum against a
#' per-component percentile over-rejects, so both the band and the p-value come
#' from the null distribution of the same max statistic. The band is therefore
#' the line the verdict is actually made on: a point outside it is a point the
#' verdict flags, which is not true of a conventional 2/sqrt(n) band.
#'
#' Leave-one-out is not optional here. Standardising a null replicate by a mean
#' and standard deviation computed from a set that includes it shrinks its own
#' deviation, while the observed value, which is not in that set, is not shrunk.
#' Measured, that put every lag at the p floor, true positives and false
#' positives alike.
#'
#' @param obs A `maxlag` by `C` matrix of observed values.
#' @param arr A `maxlag` by `C` by `B` array of null values.
#' @param labels Component names.
#' @param n_pairs Pairs behind each lag, for the table.
#' @return A list with the per-lag table, the full z matrix, the chosen
#'   component per lag and the critical value the band is drawn at.
#' @keywords internal
#' @noRd
ilm_lag_verdict <- function(obs, arr, labels, n_pairs = NULL) {
  maxlag <- nrow(obs); C <- ncol(obs)
  mu <- apply(arr, 1:2, mean, na.rm = TRUE)
  sdv <- pmax(apply(arr, 1:2, stats::sd, na.rm = TRUE), 1e-12)
  zmat <- (obs - mu) / sdv
  pick <- apply(abs(zmat), 1, function(v)
    if (all(is.na(v))) NA_integer_ else which.max(v))
  take <- function(M) vapply(seq_len(maxlag), function(L)
    if (is.na(pick[L])) NA_real_ else M[L, pick[L]], 0)
  z <- take(zmat); o <- take(obs); m <- take(mu); s <- take(sdv)

  p <- crit <- rep(NA_real_, maxlag)
  for (L in seq_len(maxlag)) {
    M <- matrix(arr[L, , ], nrow = C)
    zz <- matrix(NA_real_, C, ncol(M))
    for (b in seq_len(ncol(M))) {
      Mo <- M[, -b, drop = FALSE]
      mm <- apply(Mo, 1, mean, na.rm = TRUE)
      ss <- pmax(apply(Mo, 1, stats::sd, na.rm = TRUE), 1e-12)
      zz[, b] <- (M[, b] - mm) / ss
    }
    nullmax <- suppressWarnings(apply(abs(zz), 2, max, na.rm = TRUE))
    nullmax <- nullmax[is.finite(nullmax)]
    if (length(nullmax) >= 10L) {
      crit[L] <- unname(stats::quantile(nullmax, 0.95))
      if (!is.na(z[L]))
        p[L] <- (1 + sum(nullmax >= abs(z[L]))) / (1 + length(nullmax))
    }
  }
  st <- ifelse(is.na(p), "INCONCLUSIVE",
        ifelse(p <= 0.01, "FAIL", ifelse(p <= 0.05, "WARN", "OK")))
  tab <- data.frame(lag = seq_len(maxlag),
                    component = labels[pick],
                    estimate = round(o, 4),
                    null_mean = round(m, 4),
                    lo = round(m - crit * s, 4),
                    hi = round(m + crit * s, 4),
                    maxz = round(z, 2),
                    p = round(p, 4),
                    status = st,
                    stringsAsFactors = FALSE)
  if (!is.null(n_pairs)) tab$n_pairs <- as.integer(n_pairs)
  list(table = tab, z = zmat, pick = pick, crit = crit)
}

## ---- the shared refit loop -------------------------------------------------

## Simulate from the fit, refit each replicate, and recompute a statistic on it.
## Every envelope diagnostic in the package needs exactly this, and it is the
## expensive part; keeping one copy means the autocorrelation check and the
## variogram cannot drift apart in how they build their reference.
##
## `statfun` takes a fitted model and returns a matrix of the shape `dims`.
#' @keywords internal
#' @noRd
ilm_refit_stat <- function(object, statfun, dims, B, ncores, seed,
                           exports = NULL, where = parent.frame()) {
  ys <- ilm_sim_cond(object, B, seed + 1L)
  X <- object$X; J <- object$J; rl <- ilm_re_list_of(object)
  rs <- object$re_struct; arstr <- object$ar; yl <- object$ylevels
  asg <- object$assign; tl <- object$term_labels
  fm <- object$family; wt <- object$weights; cnsr <- object$censor
  n <- prod(dims)
  cl <- ilm_pool(ncores)
  on.exit(if (!is.null(cl)) try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
  if (!is.null(cl) && length(exports))
    try(parallel::clusterExport(cl, exports, envir = where), silent = TRUE)
  one <- function(b) {
    ## A replicate that fails to fit is dropped and reduces n_ok, so the
    ## optimiser's complaints about it are noise the user cannot act on.
    f <- suppressWarnings(try(ilm_fit(X, ys[, b], J, rl, rs, arstr, ylevels = yl,
                                      weights = wt, censor = cnsr, family = fm,
                                      verbose = FALSE, restarts = 1L),
                              silent = TRUE))
    if (inherits(f, "try-error") || f$opt$convergence != 0)
      return(rep(NA_real_, n))
    f$assign <- asg; f$term_labels <- tl; f$y <- ys[, b]
    as.numeric(statfun(f))
  }
  reps <- ilm_lapply(cl, seq_len(B), one)
  arr <- array(unlist(reps), c(dims, length(reps)))
  list(null = arr, n_ok = sum(!is.na(arr[1, 1, ])))
}

## ---- the shared engine -----------------------------------------------------

## Simulating from the fit and refitting is far and away the expensive part of
## this diagnostic, and both ilm_check_ar() and ilm_plot_acf() need exactly the
## same thing from it, so it happens once here. The result can be handed to
## either, which is why ilm_plot_acf() accepts one in place of a model.

#' @keywords internal
#' @noRd
ilm_ar_envelope <- function(object, time, group, maxlag = 8L, B = 30L,
                            ncores = 1L, seed = 1L) {
  if (!inherits(object, "ilm_model"))
    stop("`object` must be a fitted ilm_model object, not ", class(object)[1],
         call. = FALSE)
  N <- nrow(object$X)
  if (missing(time) || missing(group))
    stop("`time` and `group` are both required: autocorrelation is measured ",
         "between observations of the same unit at known time points. Supply ",
         "one value of each per observation (the model has ", N, " rows).",
         call. = FALSE)
  if (length(time) != N)
    stop("`time` has ", length(time), " values but the model has ", N,
         " rows. If the model dropped incomplete cases, subset `time` the ",
         "same way.", call. = FALSE)
  if (length(group) != N)
    stop("`group` has ", length(group), " values but the model has ", N,
         " rows. If the model dropped incomplete cases, subset `group` the ",
         "same way.", call. = FALSE)
  maxlag <- suppressWarnings(as.integer(maxlag)[1])
  if (is.na(maxlag) || maxlag < 1L)
    stop("`maxlag` must be a positive whole number", call. = FALSE)
  ## as.numeric() on a factor returns level codes, not labels, which would
  ## silently redefine what a lag means. Make that conversion the caller's.
  if (is.factor(time))
    stop("`time` is a factor. Converting it here would use the level order, ",
         "not the labels, and silently change what one lag means. If periods ",
         "are equally spaced, pass as.integer(time); otherwise pass the real ",
         "times, e.g. as.integer(as.character(time)).", call. = FALSE)
  tnum <- suppressWarnings(as.numeric(time))
  if (anyNA(tnum))
    stop("`time` must be numeric, integer or Date, not ", class(time)[1],
         ". Convert it first; lag k means k units of whatever you supply.",
         call. = FALSE)
  if (any(abs(tnum - round(tnum)) > 1e-8))
    warning("`time` is not whole-numbered, so it has been rounded. Lag k means ",
            "k units of `time`, so rescale it first if one unit is not one ",
            "sampling interval.", call. = FALSE)
  tim <- as.integer(round(tnum)); grp <- as.integer(factor(group))
  if (anyDuplicated(cbind(grp, tim)))
    warning("some group and time combinations are repeated. Lagged pairs are ",
            "matched on exact time differences, so repeats are dropped from ",
            "the pairing rather than averaged.", call. = FALSE)
  span <- max(tapply(tim, grp, function(v) diff(range(v))))
  if (maxlag > span)
    stop("`maxlag` is ", maxlag, " but the longest series spans only ", span,
         " time unit", if (span == 1L) "" else "s",
         ", so no pair of observations can be that far apart. Use maxlag <= ",
         span, ".", call. = FALSE)
  B <- suppressWarnings(as.integer(B)[1])
  if (is.na(B) || B < 2L)
    stop("`B` must be at least 2; the envelope is the spread across ",
         "simulated replicates", call. = FALSE)

  ## one ACF series per residual component: J one-against-the-rest columns for
  ## a multinomial, a single plain Pearson residual for every other family
  acf_all <- function(fit) {
    R <- ilm_pearson_ovr(fit)
    m <- matrix(NA_real_, maxlag, ncol(R)); np <- NULL
    for (j in seq_len(ncol(R))) {
      a <- ilm_resid_acf(R[, j], grp, tim, maxlag)
      if (is.null(np)) np <- attr(a, "n_pairs")
      m[, j] <- as.numeric(a)
    }
    attr(m, "n_pairs") <- np
    m
  }
  obs <- acf_all(object)
  npair <- attr(obs, "n_pairs")
  C <- ncol(obs)
  labels <- if (C == 1L) "residual" else object$ylevels

  env <- ilm_refit_stat(object, acf_all, c(maxlag, C), B, ncores, seed,
                        exports = c("grp", "tim", "maxlag"),
                        where = environment())
  nullarr <- env$null; nok <- env$n_ok

  ## the same lags in partial form, derived from the autocorrelations rather
  ## than refitted again: the Durbin-Levinson map is deterministic, so the null
  ## replicates transform straight through and cost nothing
  pobs <- matrix(apply(obs, 2, ilm_pacf_dl), maxlag, C)
  pnull <- array(apply(nullarr, c(2, 3), ilm_pacf_dl),
                 c(maxlag, C, dim(nullarr)[3]))

  vacf <- ilm_lag_verdict(obs, nullarr, labels, npair)
  vpacf <- ilm_lag_verdict(pobs, pnull, labels, npair)
  if (nok < 10L) {
    vacf$table$status <- "INCONCLUSIVE"
    vpacf$table$status <- "INCONCLUSIVE"
  }
  structure(list(acf = vacf$table, pacf = vpacf$table,
                 obs_acf = obs, obs_pacf = pobs,
                 z_acf = vacf$z, z_pacf = vpacf$z,
                 null_acf = nullarr, null_pacf = pnull,
                 labels = labels, n_pairs = npair, n_ok = nok, B = B,
                 maxlag = maxlag, has_ar = !is.null(object$ar)),
            class = "ilm_ar_envelope")
}

## Telling a cycle apart from autoregression.
##
## Both fill the ACF with flagged lags, but they need OPPOSITE remedies, and
## AR(1) cannot represent a cycle at all. The distinguishing feature is not the
## sign change -- an AR(1) fitted without its AR term goes negative at the
## longer lags too, once the group mean is absorbed -- it is the RETURN. A
## seasonal autocorrelation comes back up to a positive peak at the period; an
## autoregressive one decays through zero and stays there.
##
## Measured on a 50 x 36 monthly panel with an annual cycle and no
## autoregression, the ACF ran +0.45 +0.28 -0.02 -0.27 -0.46 -0.55 -0.43 -0.26
## +0.01 +0.25 +0.45 +0.50 -- a clean peak at lag 12 -- while the same rule
## found nothing on AR(1) data.
#' @keywords internal
#' @noRd
ilm_period_hint <- function(est, bad, maxlag) {
  ok <- is.finite(est)
  if (length(bad) < 2L || sum(ok) < 4L) return(NULL)
  neg <- which(ok & est < 0)
  if (!length(neg) || min(neg) >= maxlag) return(NULL)
  after <- seq.int(min(neg) + 1L, maxlag)
  after <- after[ok[after] & est[after] > 0 & after %in% bad]
  if (!length(after)) return(NULL)
  s <- after[which.max(est[after])]
  ## a period of 2 is alternation, which is a different animal and usually a
  ## coding artefact rather than a season
  if (s < 3L) return(NULL)
  s
}

## When the cycle is longer than the lags examined, the return never appears and
## the pattern looks autoregressive. Rising back toward zero at the last lag is
## the tell that there is more to see.
#' @keywords internal
#' @noRd
ilm_period_truncated <- function(est, maxlag) {
  if (maxlag < 3L) return(FALSE)
  e <- est[c(maxlag - 1L, maxlag)]
  isTRUE(all(is.finite(e)) && e[1] < 0 && e[2] > e[1])
}

## Shared verdict prose, so the printed report and the plot annotation cannot
## drift apart.
#' @keywords internal
#' @noRd
ilm_ar_advice <- function(spec, which = "acf") {
  tab <- spec[[which]]
  bad <- which(tab$status %in% c("WARN", "FAIL"))
  if (!length(bad)) return(NULL)
  why <- if (identical(which, "acf"))
    "residual autocorrelation the model does not account for"
  else "correlation at these lags beyond what earlier lags explain"
  per <- ilm_period_hint(tab$estimate, bad, spec$maxlag)
  if (!is.null(per)) {
    why <- paste0("a cycle of period ", per,
                  ", not autoregression: the correlation returns to a peak ",
                  "there rather than decaying away")
    fix <- paste0("add a seasonal term of period ", per,
                  " to the fixed effects -- ilm_fourier(time, ", per,
                  ") for a plain rise and fall, ilm_cyclic(time, ", per,
                  ") for a cycle with a sharp local peak, or a factor for the ",
                  "phase if the shape jumps. An AR term cannot represent a ",
                  "cycle, so adding one here will not help")
    fix_short <- paste0("periodic at lag ", per, ": try ilm_fourier or ",
                        "ilm_cyclic(time, ", per, "), not AR")
  } else {
    fix <- if (!spec$has_ar)
      paste0("add an AR(1) term over time within group, but check observations ",
             "per AR latent first, since one latent per observation is not ",
             "estimable")
    else "raise the AR order, or coarsen the AR time grid"
    ## the plot gets a shorter form: a caption that has to be truncated to fit
    ## is worse than a caption that was written to fit
    fix_short <- if (!spec$has_ar) "try an AR(1) term over time within group"
                 else "try raising the AR order or coarsening its time grid"
    if (ilm_period_truncated(tab$estimate, spec$maxlag)) {
      fix <- paste0(fix, ". The correlation is also turning back up at lag ",
                    spec$maxlag, ", which can be the start of a cycle: raise ",
                    "maxlag past one suspected period before settling on AR")
      fix_short <- paste0("still turning at lag ", spec$maxlag,
                          ": raise maxlag to check for a cycle")
    }
  }
  ## The order reading, and the one caveat that keeps it honest. Residual
  ## autocorrelation computed from a model that OMITS the structure is
  ## distorted at the higher lags by the omission itself -- an AR(1) process
  ## fitted without an AR term leaves a flagged lag 2 as an artefact -- so a
  ## contiguous run is reported as what it looks like, not asserted as the
  ## order.
  ord <- NULL
  if (identical(which, "pacf") && is.null(per) && 1L %in% bad) {
    k <- 1L
    while (k < spec$maxlag && (k + 1L) %in% bad) k <- k + 1L
    ord <- if (k == 1L)
      "lag 1 alone is outside the envelope, which is the AR(1) signature"
    else if (!spec$has_ar)
      paste0("lags 1 to ", k, " are outside the envelope. That looks like AR(",
             k, "), but a model with no AR term distorts its own residual ",
             "autocorrelation at the higher lags, so fit AR(1) and look again ",
             "before reading the order as above 1")
    else
      paste0("lags 1 to ", k, " are still outside the envelope with AR(1) ",
             "fitted, which points to an order above the AR(1) illume offers")
  }
  ord_short <- if (is.null(ord)) NULL
    else if (grepl("lag 1 alone", ord, fixed = TRUE))
      "lag 1 alone: the AR(1) signature"
    else paste0(strsplit(ord, ". ", fixed = TRUE)[[1]][1],
                ", so fit AR(1) and look again")
  list(lags = bad, why = why, fix = fix, order = ord,
       fix_short = fix_short, order_short = ord_short)
}

#' @keywords internal
#' @noRd
ilm_ar_report <- function(spec, which = "acf") {
  lab <- if (identical(which, "acf")) "autocorrelation" else
    "partial autocorrelation"
  cat(sprintf("\nresidual %s vs simulated envelope (%d refits)\n", lab,
              spec$n_ok))
  print(spec[[which]], row.names = FALSE)
  ad <- ilm_ar_advice(spec, which)
  if (is.null(ad)) {
    if (spec$n_ok >= 10L)
      cat(sprintf("\n>> residual %s is consistent with the fitted model.\n", lab))
    else
      cat(sprintf("\n>> only %d of %d replicates refitted, so no verdict.\n",
                  spec$n_ok, spec$B))
    return(invisible(NULL))
  }
  cat(sprintf("\n>> %d of %d lags outside the envelope (lag %s).\n",
              length(ad$lags), spec$maxlag, paste(ad$lags, collapse = ", ")))
  cat("   why: ", ad$why, "\n", sep = "")
  if (!is.null(ad$order)) cat("   read: ", ad$order, "\n", sep = "")
  cat("   try: ", ad$fix, "\n", sep = "")
  ## Lags are correlated with one another, so the count is a guide rather than
  ## a multiple-testing correction.
  cat("   note: lags are not independent, so treat the count as indicative\n")
  ## With B replicates no p-value can fall below 1/(B+1), so a FAIL verdict
  ## (p <= 0.01) is simply unreachable for small B, and an enormous z can then
  ## only ever reach WARN. Say so rather than let it look like weak evidence.
  if (1 / (spec$n_ok + 1) > 0.01)
    cat(sprintf("   note: p cannot fall below %.3f with %d replicates, so FAIL is\n         unreachable; raise B to about 100 if you need FAIL-level verdicts\n",
                1 / (spec$n_ok + 1), spec$n_ok))
  invisible(NULL)
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
#' would over-flag here. Because the statistic takes a maximum across residual
#' components, its reference is calibrated by leave-one-out to avoid an
#' artificial bias toward significance.
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
#' @return Invisibly, a list with the per-lag `table`, the matching `pacf`
#'   table, the observed values and per-component z scores, the simulated null
#'   and the number of replicates that refitted. The whole thing can be handed
#'   straight to [ilm_plot_acf()] to avoid refitting.
#' @references
#' Rue, H., & Held, L. (2005). *Gaussian Markov Random Fields: Theory and
#' Applications*. Chapman & Hall/CRC. (Background on latent autoregressive
#' structures.)
#' @seealso [ilm_plot_acf()] for the picture, [ilm_appraise()],
#'   [ilm_check_omitted()].
#' @examples
#' set.seed(1)
#' d <- ilm_sim(n_id = 20, n_period = 8)
#' f <- ilm_model(score ~ income + (1 | id), data = d, family = "gaussian",
#'                verbose = FALSE)
#' ilm_check_ar(f, time = as.integer(factor(d$date)), group = d$id,
#'              maxlag = 3, B = 12)
#' @export
ilm_check_ar <- function(object, time, group, maxlag = 8L, B = 30L,
                         ncores = 1L, seed = 1L, verbose = TRUE, plot = FALSE) {
  spec <- ilm_ar_envelope(object, time, group, maxlag, B, ncores, seed)
  if (verbose) ilm_ar_report(spec, "acf")
  if (plot) ilm_acf_panel(spec, "acf")
  invisible(list(table = spec$acf, pacf = spec$pacf,
                 acf_by_component = spec$obs_acf,
                 z_by_category = spec$z_acf, null = spec$null_acf,
                 n_ok = spec$n_ok, envelope = spec))
}
