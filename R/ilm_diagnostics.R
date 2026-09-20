## illume: residual diagnostics for a NOMINAL multinomial mixed model.
##
## Quantile residuals need a CDF, and a CDF needs an ordering.  Nominal
## categories have none, which is why DHARMa cannot be applied directly.  The
## way out: you do not need an ordering of categories, only of some scalar, and
## the model itself supplies one.  For observation i the log score
##      s_i = -log p_hat[i, y_i]
## is a real-valued discrete random variable taking the value -log p_hat[i, j]
## with probability p_hat[i, j].  Its null CDF is therefore fully known, and the
## standard randomised PIT applies:
##      u_i = F_i(s^-) + v_i * (F_i(s) - F_i(s^-)),   v_i ~ U(0, 1)
## giving u_i ~ U(0, 1) under a correct model, with no arbitrary ordering
## anywhere -- the fitted probabilities induce it.

#' Randomised quantile residuals for a nominal outcome
#'
#' Turns categorical outcomes into residuals that should be uniformly
#' distributed when the model is right.
#'
#' @section The obstacle, and the way around it:
#' Quantile residuals need a cumulative distribution function, and a CDF needs an
#' ordering. Nominal categories have none -- "red", "green", "blue" cannot be put
#' in order -- which is why `DHARMa` cannot be applied directly to this model.
#'
#' The resolution is that you do not need an ordering of *categories*, only of
#' some number, and the model supplies one. For each observation the **log
#' score** `-log(p)` at the observed category is a real-valued quantity whose
#' distribution under the model is completely known: it takes the value
#' `-log(p_j)` with probability `p_j`. The standard randomised transform then
#' applies, and the fitted probabilities themselves provide the ordering.
#'
#' @section Do not test these against a uniform distribution:
#' Verified by simulation: when the probabilities are known exactly the residuals
#' are uniform, as intended. But **in sample** they are not, because the random
#' effects were fitted to the same data, which pulls the observed log score down.
#' A Kolmogorov-Smirnov test against uniform therefore rejects a correctly
#' specified model. Use [ilm_rqr_test()], which builds its reference by
#' refitting simulated data and so accounts for this automatically.
#'
#' @param object A fitted `"ilm_model"` object.
#' @param conditional Logical. Evaluate at the fitted random effects.
#' @param seed Integer. Random seed; the transform uses randomisation.
#' @return A numeric vector on the unit interval.
#' @references
#' Dunn, P. K., & Smyth, G. K. (1996). Randomized quantile residuals. *Journal of
#' Computational and Graphical Statistics*, 5(3), 236--244.
#'
#' Hartig, F. (2022). DHARMa: Residual diagnostics for hierarchical regression
#' models. R package. (The simulation-based approach this borrows from.)
#' @seealso [ilm_rqr_test()], [ilm_appraise()].
#' @export
ilm_rqr <- function(object, conditional = TRUE, seed = 1L) {
  set.seed(seed)
  fam <- if (!is.null(object$family)) object$family$name else "gaussian"
  y <- as.numeric(object$y)
  N <- length(y)

  ## Randomised quantile residuals (Dunn & Smyth 1996). For a continuous
  ## response the CDF at the observation is already uniform. For a discrete one
  ## the CDF jumps, so the residual is drawn uniformly across the jump, which
  ## restores uniformity. Only the multinomial needs the log-score construction
  ## below, because its outcomes have no natural order to take a CDF along.
  if (!identical(fam, "multinomial")) {
    mu <- as.numeric(ilm_fitted(object, conditional)[, 1])
    w <- if (is.null(object$weights)) rep(1, N) else object$weights
    disp <- object$dispersion
    ## An accelerated failure time family supplies its own distribution
    ## function, so this does not need a third switch over families to keep in
    ## step with the two the likelihood already has.
    u <- if (isTRUE(object$family$aft))
      object$family$cdf(y, as.numeric(ilm_eta_hat(object, conditional)[, 1]),
                        unname(disp[1]))
    else switch(fam,
      gaussian = stats::pnorm(y, mu, if (length(disp)) unname(disp[1]) else 1),
      poisson  = {
        lo <- stats::ppois(y - 1, mu); hi <- stats::ppois(y, mu)
        lo + stats::runif(N) * (hi - lo)
      },
      nbinom = {
        k <- unname(disp[1])
        lo <- stats::pnbinom(y - 1, size = k, mu = mu)
        hi <- stats::pnbinom(y, size = k, mu = mu)
        lo + stats::runif(N) * (hi - lo)
      },
      binomial = {
        ## y is a proportion when weights give the number of trials, and 0/1
        ## otherwise; both are counts of successes out of size
        size <- pmax(1, round(w)); k <- round(y * size)
        lo <- stats::pbinom(k - 1, size, mu); hi <- stats::pbinom(k, size, mu)
        lo + stats::runif(N) * (hi - lo)
      },
      stop("no quantile residual is defined for family ", fam, call. = FALSE))

    ## A censored observation is known only as an interval, so its residual is
    ## drawn uniformly across the probability of that interval -- the same
    ## randomisation Dunn and Smyth use for a discrete response, applied to the
    ## censoring interval instead of to a jump. Without it every censored row
    ## lands on the same quantile and the QQ plot shows a spike at the limit
    ## that has nothing to do with model fit.
    cz <- ilm_censor_for(object$censor, y)
    if (!is.null(cz) && any(cz != 0L)) {
      l <- cz < 0L; r <- cz > 0L
      if (any(l)) u[l] <- stats::runif(sum(l), 0, u[l])
      if (any(r)) u[r] <- stats::runif(sum(r), u[r], 1)
    }
    return(pmin(pmax(u, 0), 1))
  }

  P <- ilm_fitted(object, conditional)
  S <- -log(pmax(P, .Machine$double.eps))       # atom values, N x J
  obs <- S[cbind(seq_len(N), as.integer(y))]    # observed log score
  tol <- 1e-10
  ## a length-N vector recycles down the columns of an N x J matrix, so each
  ## element is compared against its own row's observed score
  lo <- rowSums(P * (S <  obs - tol))
  hi <- rowSums(P * (S <= obs + tol))           # ties merged into one atom
  lo + stats::runif(N) * (hi - lo)
}

## DO NOT TEST THESE AGAINST THE THEORETICAL UNIFORM.  Established by simulation:
##   oracle (probabilities KNOWN, nothing estimated): mean u = 0.5046, and a KS
##     test against U(0,1) rejected 0 of 40 replicates at alpha = 0.05.  The
##     construction is exactly right.
##   in-sample (conditional modes fitted to the same data): mean u = 0.4869
##     across refits, sd 0.0039.  Fitting pulls the observed log score down.
## A naive KS test against U(0,1) therefore gives p = 0.002 on a CORRECTLY
## specified model.  The reference distribution must come from refitting data
## simulated from the fit -- which is what ilm_rqr_test() does, and why the
## verdict framework is simulation-calibrated throughout.

#' Draw replicate outcomes from fitted probabilities
#' @param P Matrix of fitted probabilities.
#' @param B Integer. Number of replicates.
#' @param seed Integer. Random seed.
#' @return An integer matrix of category codes, one column per replicate.
#' @keywords internal
#' @noRd
ilm_sim_cond <- function(object, B, seed = 1L) {
  set.seed(seed)
  fam <- object$family
  eta <- ilm_eta_hat(object, TRUE)
  N <- nrow(eta)
  if (is.null(fam) || identical(fam$name, "multinomial")) {
    Tct <- t(contr.sum(object$J))
    return(vapply(seq_len(B),
                  function(b) fam$sim(eta, NULL, numeric(0), Tct = Tct),
                  integer(N)))
  }
  w <- if (is.null(object$weights)) rep(1, N) else object$weights
  ## the family simulators take the dispersion on the log scale, which is how
  ## it is parameterised internally
  disp <- if (!is.null(object$dispersion)) log(unname(object$dispersion))
          else numeric(0)
  ## A replicate has to be censored the way the data were, or every envelope
  ## built on it is calibrated against a process that is not the fitted one.
  ## Whether a given draw lands beyond the limit is itself random, which is
  ## part of the variability the envelope is meant to carry.
  vapply(seq_len(B), function(b)
    ilm_censor_apply(object$censor, as.numeric(fam$sim(eta, w, disp))),
    numeric(N))
}

#' Simulate category draws from a fitted probability matrix
#'
#' Multinomial-only. [ilm_sim_cond()] is the family-aware entry point.
#'
#' @param P An N x J matrix of fitted probabilities.
#' @param B Integer. Number of simulated data sets.
#' @param seed Integer. Random seed.
#' @return An N x B integer matrix of category indices.
#' @keywords internal
#' @noRd
ilm_sim_from_P <- function(P, B, seed = 1L) {
  set.seed(seed); N <- nrow(P); cp <- t(apply(P, 1, cumsum))
  vapply(seq_len(B), function(b) as.integer(rowSums(runif(N) > cp)) + 1L, integer(N))
}

#' Calibration of predicted probabilities
#'
#' Asks a direct question: among observations the model gave roughly a 30%
#' chance of category A, did about 30% turn out to be category A? Predictions
#' that pass this are said to be well calibrated, and for a categorical outcome
#' this is often more informative than any residual plot.
#'
#' The comparison band is simulated from the model itself rather than taken from
#' a formula, so it reflects how much scatter is expected at this sample size.
#'
#' @param object A fitted `"ilm_model"` object.
#' @param nbins Integer. Number of bins of predicted probability.
#' @param B Integer. Simulated datasets for the band.
#' @param seed Integer. Random seed.
#' @return A list with one element per category, each holding mean predicted and
#'   observed proportions per bin plus the simulated interval.
#' @export
ilm_calibration <- function(object, nbins = 10L, B = 200L, seed = 1L) {
  fam <- if (!is.null(object$family)) object$family$name else "gaussian"
  ## Calibration compares predicted PROBABILITIES against observed frequencies,
  ## so it needs a model that predicts a probability. For gaussian, Poisson and
  ## negative binomial the honest analogue is not a calibration curve but the
  ## randomised quantile residuals, which ilm_rqr() already provides.
  if (!fam %in% c("multinomial", "binomial"))
    stop("ilm_calibration() applies to the binomial and multinomial families, ",
         "which predict probabilities. For ", fam,
         " use ilm_rqr() or the residual panels in ilm_appraise().",
         call. = FALSE)

  P <- ilm_fitted(object, TRUE)
  ## A binomial fit predicts one probability, of the modelled outcome. Putting
  ## it in the two-column form the multinomial path already uses lets one
  ## implementation serve both, with the second column the complement.
  if (identical(fam, "binomial")) {
    p1 <- as.numeric(P[, 1])
    P <- cbind(1 - p1, p1)
    lab <- if (!is.null(object$ylevels) && length(object$ylevels) == 2L)
      object$ylevels else c("0", "1")
    colnames(P) <- lab
    y <- as.integer(object$y) + 1L        # 0/1 becomes column index 1/2
    J <- 2L
  } else {
    y <- object$y; J <- object$J
  }

  ## The envelope's 2.5% and 97.5% points are order statistics of B draws. At
  ## B = 60 they are so poorly estimated that the envelope comes out too narrow
  ## and roughly a fifth of correct models pick up a warning, against the 8.6%
  ## the binomial null predicts. Above ~150 the two agree.
  if (B < 150L)
    warning("B = ", B, " gives an unstable envelope: the 2.5% and 97.5% points ",
            "are estimated from too few draws, which inflates the apparent ",
            "number of bins outside it. Use B >= 200.", call. = FALSE)
  ys <- ilm_sim_cond(object, B, seed)
  if (identical(fam, "binomial")) ys <- ys + 1L

  out <- lapply(seq_len(J), function(j) {
    br <- unique(stats::quantile(P[, j], seq(0, 1, length.out = nbins + 1L)))
    ## a prediction that barely varies cannot be binned into a curve
    if (length(br) < 3L) return(NULL)
    bin <- cut(P[, j], br, include.lowest = TRUE)
    mp  <- tapply(P[, j], bin, mean)
    ob  <- tapply(y == j, bin, mean)
    sim <- vapply(seq_len(B),
                  function(b) tapply(ys[, b] == j, bin, mean),
                  numeric(nlevels(bin)))
    ## The envelope is what makes this worth plotting. An in-sample curve on
    ## its own is close to self-fulfilling: fitting with an intercept forces
    ## the mean residual to zero, so calibration-in-the-large is guaranteed.
    ## What the envelope tests is the SHAPE -- systematic bending points at a
    ## misspecified link, and misfit at the extremes at missing nonlinearity.
    lo <- apply(sim, 1, stats::quantile, 0.025, na.rm = TRUE)
    hi <- apply(sim, 1, stats::quantile, 0.975, na.rm = TRUE)
    outside <- sum(ob < lo | ob > hi, na.rm = TRUE)
    ## With 95% envelopes, the number of bins falling outside is Binomial(bins,
    ## 0.05) when the model is correct -- so at 10 bins there is a 40% chance
    ## of at least one, and flagging on that would mean warning about most
    ## correctly specified models. The status is graded against that null.
    nb <- nlevels(bin)
    ptail <- stats::pbinom(outside - 1L, nb, 0.05, lower.tail = FALSE)
    list(mean_p = mp, obs = ob, lo = lo, hi = hi,
         n_bins = nb, n_outside = outside, p_outside = round(ptail, 4),
         status = if (ptail < 0.05) "FAIL" else if (ptail < 0.20) "WARN" else "OK")
  })
  names(out) <- if (!is.null(colnames(P))) colnames(P) else
    if (!is.null(object$ylevels)) object$ylevels else as.character(seq_len(J))
  out
}

#' Distances of the fitted random effects from zero
#'
#' The model assumes random effects are multivariate normal. Mahalanobis
#' distance reduces each group's vector of effects to one number measuring how
#' far it lies from the centre, allowing for their covariance, so the set can be
#' compared against a chi-square reference.
#'
#' @section An important caveat:
#' Fitted random effects are **shrunk** toward zero -- that is what makes them
#' useful -- so their spread is smaller than the true random effects. Comparing
#' them to the theoretical chi-square distribution is therefore
#' anti-conservative: a misspecified covariance can still look acceptable. Treat
#' this panel as indicative rather than as a test. The plot produced by
#' [ilm_appraise()] is annotated accordingly.
#'
#' @param object A fitted `"ilm_model"` object.
#' @param term Optional name of the grouping term; defaults to the first.
#' @return A list with the distances, degrees of freedom and term name, or
#'   `NULL` if the model has no grouping terms.
#' @export
ilm_re_mahalanobis <- function(object, term = NULL) {
  gi <- which(vapply(object$re, function(e) e$kind != "basis", TRUE))
  if (!length(gi)) return(NULL)
  k <- if (is.null(term)) gi[1] else which(names(object$re) == term)
  B <- ilm_Bhat_term(object, k); nl <- object$nlk[k]; w <- object$wk[k]
  Bi <- B[seq_len(nl), , drop = FALSE]                 # first RE dimension
  isrr <- identical(object$re_struct[[k]]$type, "rr")
  d <- if (isrr) rowSums(Bi^2) else {                  # rr prior is identity
    S <- object$Sigma[[k]]
    Si <- tryCatch(solve(S + diag(1e-10, nrow(S))), error = function(e) NULL)
    if (is.null(Si)) return(NULL)
    rowSums((Bi %*% Si) * Bi)
  }
  list(d = d, df = w, term = names(object$re)[k])
}

#' Test residuals against a simulated reference
#'
#' Takes any single-number summary of the residuals and compares it to the
#' distribution that summary takes when the model is true, obtained by simulating
#' from the fit and refitting each replicate.
#'
#' The refitting matters: it builds in the fact that in-sample residuals are not
#' exactly uniform (see [ilm_rqr()]), so a correctly specified model is not
#' flagged for a discrepancy that fitting itself created.
#'
#' @section What this can and cannot detect:
#' Broad summaries of the residuals have essentially **no power** against an
#' omitted covariate. This is structural rather than a matter of tuning: the
#' reference distribution is "this model refitting data it generated", and a
#' misspecified model reproduces its own behaviour faithfully, so the
#' misspecification cancels from both sides of the comparison.
#'
#' A statistic aimed at specific structure does have power. Use
#' [ilm_check_covariate()] or [ilm_check_omitted()] to point the test at
#' particular variables, including ones the model does not contain.
#'
#' @param object A fitted `"ilm_model"` object.
#' @param B Integer. Simulated datasets.
#' @param ncores Integer. Worker processes.
#' @param seed Integer. Random seed.
#' @param stat Function taking the residual vector and returning one number.
#' @param verbose Logical. Print the result.
#' @return Invisibly, a list with `status`, the observed statistic, the simulated
#'   null, and a z score and p-value.
#' @export
ilm_rqr_test <- function(object, B = 30L, ncores = 1L, seed = 1L,
                          stat = function(u) mean(u), verbose = TRUE) {
  obs <- stat(ilm_rqr(object, TRUE, seed))
  ys <- ilm_sim_cond(object, B, seed + 1L)
  X <- object$X; J <- object$J; rl <- ilm_re_list_of(object); rs <- object$re_struct
  arr <- object$ar; yl <- object$ylevels; asg <- object$assign; tl <- object$term_labels
  fm <- object$family; wt <- object$weights; cnsr <- object$censor
  cl <- ilm_pool(ncores)
  on.exit(if (!is.null(cl)) try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
  ## `stat` may close over data the workers do not have (e.g. a covariate the
  ## model omits, which is exactly the case with power -- see below), so ship
  ## its enclosing environment too.
  if (!is.null(cl)) {
    src <- environment(stat); if (is.null(src)) src <- globalenv()
    free <- tryCatch(codetools::findGlobals(stat, merge = FALSE)$variables,
                     error = function(e) character(0))
    have <- free[vapply(free, function(n) exists(n, envir = src), TRUE)]
    if (length(have)) try(parallel::clusterExport(cl, have, envir = src), silent = TRUE)
  }
  one <- function(b) {
    f <- try(ilm_fit(X, ys[, b], J, rl, rs, arr, ylevels = yl,
                      weights = wt, censor = cnsr, family = fm,
                      verbose = FALSE, restarts = 1L), silent = TRUE)
    if (inherits(f, "try-error") || f$opt$convergence != 0) return(NA_real_)
    f$assign <- asg; f$term_labels <- tl
    stat(ilm_rqr(f, TRUE, seed = seed + 1000L + b))
  }
  null <- unlist(ilm_lapply(cl, seq_len(B), one))
  nok <- sum(is.finite(null))
  if (nok < 10L) {
    if (verbose) cat("INCONCLUSIVE: only", nok, "of", B, "replicates refitted.\n")
    return(invisible(list(status = "INCONCLUSIVE", observed = obs, null = null)))
  }
  z <- (obs - mean(null, na.rm = TRUE)) / sd(null, na.rm = TRUE)
  p <- mean(abs(null - mean(null, na.rm = TRUE)) >=
              abs(obs - mean(null, na.rm = TRUE)), na.rm = TRUE)
  st <- if (abs(z) >= 4) "FAIL" else if (abs(z) >= 2) "WARN" else "OK"
  if (verbose)
    cat(sprintf("residual statistic %.4f | simulated null %.4f (sd %.4f, %d refits) | z = %.2f | p = %.3f -> %s\n",
                obs, mean(null, na.rm = TRUE), sd(null, na.rm = TRUE), nok, z, p, st))
  invisible(list(status = st, observed = obs, null = null, z = z, p = p))
}

#' Diagnostic plots for a fitted model
#'
#' Six panels: a quantile-residual normal plot, residuals against fitted
#' probability, per-category calibration with a simulated band, observed against
#' simulated category frequencies, a random-effects distance plot, and the model
#' check verdicts.
#'
#' These are built for a nominal categorical outcome rather than adapted from
#' tools designed for continuous responses, because the usual residual plots have
#' no clear meaning here. `performance::check_model()` routes to this function.
#'
#' @param object A fitted `"ilm_model"` object.
#' @param nbins Integer. Bins for the calibration panel.
#' @param B Integer. Simulated datasets for the bands.
#' @param seed Integer. Random seed.
#' @param ... Unused.
#' @return Invisibly, the residuals, calibration data and random-effect
#'   distances.
#' @seealso [ilm_check_ar()] for temporal correlation,
#'   [ilm_check_omitted()] for omitted variables.
#' @export
ilm_appraise <- function(object, nbins = 10L, B = 200L, seed = 1L, ...) {
  op <- par(mfrow = c(2, 3), mar = c(4, 4, 3, 1), cex = 0.8); on.exit(par(op))
  fam <- if (!is.null(object$family)) object$family$name else "gaussian"
  mn <- identical(fam, "multinomial")
  u <- ilm_rqr(object, TRUE, seed); z <- qnorm(pmin(pmax(u, 1e-6), 1 - 1e-6))
  P <- ilm_fitted(object, TRUE)

  ## 1 QQ of the log-score randomised quantile residuals
  qqnorm(z, main = "Quantile residuals", pch = 16, cex = 0.4,
         col = "#00000055"); qqline(z, col = "red", lwd = 2)

  ## 2 residuals against the fitted value. For a multinomial that is the
  ## predicted probability of the category actually observed; for every other
  ## family it is the fitted mean.
  py <- if (mn) P[cbind(seq_len(nrow(P)), as.integer(object$y))]
        else as.numeric(P[, 1])
  plot(py, z, pch = 16, cex = 0.4, col = "#00000055",
       xlab = if (mn) "fitted p(observed)" else "fitted",
       ylab = "quantile residual", main = "Residual vs fitted")
  abline(h = 0, col = "red", lwd = 2)
  lines(lowess(py, z), col = "blue", lwd = 2)

  ## 3 calibration where a probability is predicted, scale-location otherwise.
  ## A calibration curve needs a predicted probability to bin by, which a
  ## gaussian or count model does not have; what those models can get wrong
  ## instead is the spread, so that is what the panel shows.
  cal <- NULL
  if (mn || identical(fam, "binomial")) {
    cal <- ilm_calibration(object, nbins, B, seed)
    plot(0:1, 0:1, type = "n", xlab = "predicted", ylab = "observed",
         main = "Calibration (envelope = model)")
    abline(0, 1, col = "grey60", lty = 2)
    cols <- seq_along(cal) + 1L
    for (j in seq_along(cal)) {
      cj <- cal[[j]]; if (is.null(cj)) next
      arrows(cj$mean_p, cj$lo, cj$mean_p, cj$hi, angle = 90, code = 3,
             length = 0.02, col = adjustcolor(cols[j], 0.4))
      points(cj$mean_p, cj$obs, col = cols[j], pch = 16)
      lines(cj$mean_p, cj$obs, col = cols[j])
    }
    legend("topleft", legend = names(cal), col = cols, pch = 16, bty = "n",
           cex = 0.7)
  } else {
    plot(py, sqrt(abs(z)), pch = 16, cex = 0.4, col = "#00000055",
         xlab = "fitted", ylab = expression(sqrt(abs(residual))),
         main = "Scale-location")
    if (length(py) > 20L) lines(lowess(py, sqrt(abs(z))), col = "red", lwd = 2)
    mtext("flat is homoscedastic; see ilm_check_variance()", side = 3,
          line = -1, cex = 0.6, col = "grey30")
  }

  ## 4 observed against simulated: category frequencies for a multinomial,
  ## the response distribution otherwise
  ys <- ilm_sim_cond(object, B, seed)
  if (mn) {
    obs <- tabulate(as.integer(object$y), object$J) / length(object$y)
    sim <- vapply(seq_len(B), function(b)
      tabulate(ys[, b], object$J) / nrow(ys), numeric(object$J))
    rng <- range(c(obs, sim))
    plot(seq_len(object$J), obs, ylim = rng, pch = 16, xaxt = "n", xlab = "",
         ylab = "proportion", main = "Category frequencies")
    axis(1, seq_len(object$J), object$ylevels, las = 2, cex.axis = 0.7)
    for (j in seq_len(object$J))
      arrows(j, quantile(sim[j, ], .025), j, quantile(sim[j, ], .975),
             angle = 90, code = 3, length = 0.03, col = "grey50")
    points(seq_len(object$J), obs, pch = 16, col = "red")
  } else {
    yv <- as.numeric(object$y)
    dobs <- density(yv[is.finite(yv)])
    plot(dobs, main = "Observed against simulated", xlab = "response",
         lwd = 2.5, ylim = c(0, max(dobs$y) * 1.3))
    for (b in seq_len(min(B, 60L))) {
      db <- try(density(ys[is.finite(ys[, b]), b]), silent = TRUE)
      if (!inherits(db, "try-error"))
        lines(db, col = adjustcolor("grey50", 0.3))
    }
    lines(dobs, lwd = 2.5)
  }

  ## 5 random-effect Mahalanobis QQ
  rm_ <- ilm_re_mahalanobis(object)
  if (is.null(rm_)) { plot.new(); title("No grouping random effects") } else {
    qs <- qchisq(ppoints(length(rm_$d)), rm_$df)
    plot(qs, sort(rm_$d), pch = 16, cex = 0.5, col = "#00000077",
         xlab = sprintf("chi-square(%d) quantiles", rm_$df), ylab = "observed",
         main = sprintf("RE Mahalanobis: %s", rm_$term))
    abline(0, 1, col = "red", lwd = 2)
    mtext("modes are shrunk: anti-conservative", side = 3, line = -1, cex = 0.6, col = "grey30")
  }

  ## 6 verdicts from the fit's own check layer
  plot.new(); ck <- object$checks
  bad <- ck[ck$status != "OK", , drop = FALSE]
  title("Model checks")
  if (!nrow(bad)) text(0.5, 0.5, "all checks passed", cex = 1.1, col = "darkgreen")
  else {
    txt <- paste0("[", bad$status, "] ", bad$check)
    text(0, seq(0.95, by = -0.11, length.out = min(8L, nrow(bad))),
         head(txt, 8L), adj = 0, cex = 0.75,
         col = ifelse(head(bad$status, 8L) == "FAIL", "red", "darkorange"))
  }
  invisible(list(rqr = u, calibration = cal, re = rm_))
}

#' performance::check_model method
#'
#' Routes `performance::check_model()` to [ilm_appraise()], so the panels are
#' the ones built for this model rather than generic ones that would not apply.
#'
#' @param x A fitted `"ilm_model"` object.
#' @param ... Passed to [ilm_appraise()].
#' @return Invisibly, the result of [ilm_appraise()].
#' @export
check_model.ilm_model <- function(x, ...) ilm_appraise(x, ...)

## ---------------------------------------------------- targeted covariate checks
## WHY TARGETED.  Global summaries of the residuals have essentially NO power
## against an omitted covariate.  Measured on a model missing a real effect:
##   mean(u)  z = 0.24 -> OK      ks(u)  z = -0.06 -> OK      sd(u)  z = 0.05 -> OK
## The reason is structural, not a tuning failure: the reference distribution is
## "this model refitting data it generated", and a misspecified model reproduces
## its own behaviour faithfully, so the misspecification cancels out of the
## comparison.  A statistic aimed at the candidate structure does have power --
## residual spread across levels of the omitted variable gave z = 4.04 -> FAIL on
## the same pair of models.  So residual checks must be POINTED AT SPECIFIC
## COVARIATES, including ones the model does not contain.

#' Test whether residuals shift across a variable
#'
#' Asks whether residuals are systematically higher or lower at some values of a
#' variable than others. If they are, the model is missing something about that
#' variable -- it may need to be added, or added in a more flexible form.
#'
#' The variable does **not** have to be in the model; pointing the test at
#' variables you left out is exactly the useful case, and is the one situation
#' where residual testing has real power here. Broad summaries of residuals
#' detect almost nothing; a test aimed at specific structure detects a great
#' deal. See [ilm_rqr_test()] for why.
#'
#' @param object A fitted `"ilm_model"` object.
#' @param x A variable, of the same length as the data. Numeric variables are
#'   binned; factors are used as-is.
#' @param name Character label for printing.
#' @param nbin Integer. Bins for a numeric variable.
#' @param B Integer. Simulated datasets.
#' @param ncores Integer. Worker processes.
#' @param seed Integer. Random seed.
#' @param verbose Logical. Print the result.
#' @return Invisibly, a list with `status`, the observed statistic and the test.
#' @seealso [ilm_check_omitted()].
#' @export
ilm_check_covariate <- function(object, x, name = NULL, nbin = 5L, B = 30L,
                                 ncores = 1L, seed = 1L, verbose = TRUE) {
  if (is.null(name)) name <- deparse(substitute(x))
  xv <- if (is.numeric(x) && length(unique(x)) > nbin)
    cut(x, stats::quantile(x, seq(0, 1, length.out = nbin + 1L)), include.lowest = TRUE)
  else factor(x)
  gv <- as.integer(xv)
  if (length(unique(gv)) < 2L) {
    if (verbose) cat(sprintf("%-18s INCONCLUSIVE (only one level)\n", name))
    return(invisible(list(status = "INCONCLUSIVE", term = name)))
  }
  stat <- function(u) { m <- tapply(u, gv, mean); max(m) - min(m) }
  r <- ilm_rqr_test(object, B = B, ncores = ncores, seed = seed,
                     stat = stat, verbose = FALSE)
  r$term <- name
  if (verbose) {
    cat(sprintf("[%s] %-16s residual spread %.4f vs null %.4f (sd %.4f) | z = %5.2f | p = %.3f\n",
                switch(r$status, OK = "  ok  ", WARN = " WARN ", FAIL = " FAIL ", "  ??  "),
                name, r$observed,
                if (is.null(r$null)) NA else mean(r$null, na.rm = TRUE),
                if (is.null(r$null)) NA else sd(r$null, na.rm = TRUE),
                if (is.null(r$z)) NA else r$z, if (is.null(r$p)) NA else r$p))
    if (!is.null(r$status) && r$status %in% c("WARN", "FAIL")) {
      ## suggestion has to match the variable's type: s() on a factor is nonsense,
      ## and a factor with many levels usually belongs in the random structure
      nlev <- length(unique(x))
      tr <- if (is.numeric(x))
        sprintf("add %s to the fixed effects, or s(%s) if the shift is nonlinear", name, name)
      else if (nlev > 10L)
        sprintf("add %s to the fixed effects, or (1 | %s) if its %d levels are exchangeable",
                name, name, nlev)
      else sprintf("add %s to the fixed effects", name)
      cat(sprintf("        why: residuals shift systematically across %s, which the model does not account for\n        try: %s\n",
                  name, tr))
    }
  }
  invisible(r)
}

#' Screen variables the model does not use
#'
#' Runs [ilm_check_covariate()] over every variable in a data frame that the
#' model's formula does not mention, and reports any whose residual pattern
#' suggests it should have been included.
#'
#' Because one test is run per variable, some flags are expected by chance; the
#' output states how many. Treat a single borderline flag among many variables
#' with appropriate caution.
#'
#' @param object A fitted `"ilm_model"` object.
#' @param data A data frame containing candidate variables.
#' @param vars Optional character vector naming variables to test; defaults to
#'   everything not already in the model.
#' @param B Integer. Simulated datasets.
#' @param ncores Integer. Worker processes.
#' @param seed Integer. Random seed.
#' @param verbose Logical. Print results.
#' @return Invisibly, a named list of per-variable results.
#' @export
ilm_check_omitted <- function(object, data, vars = NULL, B = 30L, ncores = 1L,
                               seed = 1L, verbose = TRUE) {
  used <- all.vars(object$formula)
  cand <- if (is.null(vars)) setdiff(names(data), used) else vars
  cand <- cand[vapply(cand, function(v) {
    z <- data[[v]]; !is.null(z) && (is.numeric(z) || is.factor(z) || is.character(z))
  }, TRUE)]
  if (!length(cand)) {
    if (verbose) cat("no candidate omitted variables found\n"); return(invisible(NULL))
  }
  if (verbose)
    cat("checking", length(cand), "variable(s) not in the model:",
        paste(cand, collapse = ", "), "\n")
  res <- lapply(cand, function(v)
    ilm_check_covariate(object, data[[v]], name = v, B = B, ncores = ncores,
                         seed = seed, verbose = verbose))
  st <- vapply(res, function(r) if (is.null(r$status)) "INCONCLUSIVE" else r$status, "")
  if (verbose) {
    nb <- sum(st %in% c("WARN", "FAIL"))
    ## multiplicity: one test per candidate variable
    if (nb) cat(sprintf("\n>> %d of %d variables flagged. With %d tested at |z| >= 2,\n   about %.1f flags are expected by chance.\n",
                        nb, length(cand), length(cand), 0.0455 * length(cand)))
    else cat("\n>> no omitted-variable structure detected in the residuals.\n")
  }
  invisible(setNames(res, cand))
}
