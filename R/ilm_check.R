## Model checks that close specific gaps against performance::check_model()
## and DHARMa.
##
## Each follows the pattern the rest of illume uses: a status, the reason, and
## what to do about it. A check that reports a number without saying whether it
## is a problem has left the hard part to the reader.

## ---- collinearity ----------------------------------------------------------

#' Variance inflation among fixed effects
#'
#' [ilm_frame_issues()] catches pairs of columns that are nearly identical, but
#' the case that actually breaks a model is a predictor collinear with a
#' *combination* of the others, which no pairwise correlation reveals. That is
#' what variance inflation measures.
#'
#' For a term with more than one degree of freedom -- a factor, a spline basis
#' -- the plain VIF is not comparable across terms, so the generalised form
#' (GVIF) is reported alongside `gvif^(1/(2*df))`, which is on the scale of a
#' standard error inflation and can be compared directly.
#'
#' @param object A fitted `"ilm_model"` object.
#' @param warn,fail Thresholds on the standard-error inflation factor. The
#'   defaults correspond to the conventional VIF cut-offs of 5 and 10.
#' @return A data frame with one row per fixed-effect term: `gvif`, `df`,
#'   `se_inflation`, and a `status` of `"OK"`, `"WARN"` or `"FAIL"`.
#' @references
#' Fox, J. and Monette, G. (1992). Generalized collinearity diagnostics.
#' Journal of the American Statistical Association, 87(417), 178-183.
#' @seealso [ilm_frame_issues()] for collinearity in the raw data.
#' @examples
#' set.seed(1)
#' d <- ilm_sim()
#' f <- ilm_model(score ~ income + grp + (1 | id), data = d,
#'                family = "gaussian", verbose = FALSE)
#' ilm_check_collinearity(f)
#' @export
ilm_check_collinearity <- function(object, warn = sqrt(5), fail = sqrt(10)) {
  if (!inherits(object, "ilm_model"))
    stop("`object` must be a fitted ilm_model object, not ", class(object)[1],
         call. = FALSE)
  V <- suppressWarnings(stats::vcov(object))
  if (is.null(V) || !all(is.finite(V)))
    stop("this model has no usable coefficient covariance, so variance ",
         "inflation cannot be computed; see object$checks", call. = FALSE)
  tl <- object$term_labels
  if (is.null(tl) || !length(tl))
    stop("variance inflation needs a model fitted through the formula ",
         "interface", call. = FALSE)

  nm <- names(stats::coef(object))
  keep <- !grepl("(Intercept)", nm, fixed = TRUE)
  if (sum(keep) < 2L)
    stop("variance inflation needs at least two non-intercept coefficients; ",
         "this model has ", sum(keep), call. = FALSE)
  R <- stats::cov2cor(V[keep, keep, drop = FALSE])
  detR <- det(R)
  if (!is.finite(detR) || detR <= 0)
    stop("the coefficient correlation matrix is singular, which means exact ",
         "collinearity among the fixed effects; drop a redundant term",
         call. = FALSE)

  ## positions of each term among the retained coefficients
  pos <- lapply(seq_along(tl), function(j) {
    ix <- ilm_term_idx(object, j)
    match(ix, which(keep))
  })
  ok <- vapply(pos, function(z) length(z) && !anyNA(z), TRUE)

  out <- do.call(rbind, lapply(which(ok), function(j) {
    i <- pos[[j]]
    if (length(i) >= ncol(R)) return(NULL)   # nothing left to compare against
    g <- det(R[i, i, drop = FALSE]) * det(R[-i, -i, drop = FALSE]) / detR
    dfj <- length(i)
    sei <- g^(1 / (2 * dfj))
    data.frame(term = tl[j], df = dfj, gvif = round(g, 3),
               se_inflation = round(sei, 3),
               status = if (!is.finite(sei)) "INCONCLUSIVE"
                        else if (sei >= fail) "FAIL"
                        else if (sei >= warn) "WARN" else "OK",
               stringsAsFactors = FALSE)
  }))
  if (is.null(out))
    return(data.frame(term = character(0), df = integer(0), gvif = numeric(0),
                      se_inflation = numeric(0), status = character(0)))
  rownames(out) <- NULL
  attr(out, "note") <- paste(
    "se_inflation is how much wider a coefficient's standard error is than it",
    "would be with uncorrelated predictors.")
  out
}

## ---- predictive check ------------------------------------------------------

#' Compare the observed response with data the model would generate
#'
#' The frequentist counterpart of a posterior predictive check: simulate
#' datasets from the fitted model and see whether the real response looks like
#' one of them. It is a blunt instrument, and that is the point -- it catches
#' gross misspecification (the wrong family, unmodelled zero inflation, a
#' bounded response fitted as gaussian) at a glance, where a residual plot
#' makes you work for it.
#'
#' @param object A fitted `"ilm_model"` object.
#' @param B Simulated datasets.
#' @param seed Random seed.
#' @param plot Draw the overlay.
#' @return Invisibly, a list with the observed and simulated summaries and a
#'   `status`.
#' @seealso [ilm_check_dispersion()], [ilm_check_zeros()] for the specific
#'   failures this can only hint at.
#' @examples
#' set.seed(1)
#' d <- ilm_sim()
#' f <- ilm_model(visits ~ income + (1 | id), data = d, family = "poisson",
#'                verbose = FALSE)
#' ilm_check_predictive(f, B = 30)
#' @export
ilm_check_predictive <- function(object, B = 50L, seed = 1L, plot = TRUE) {
  if (!inherits(object, "ilm_model"))
    stop("`object` must be a fitted ilm_model object, not ", class(object)[1],
         call. = FALSE)
  fam <- if (!is.null(object$family)) object$family$name else "gaussian"
  if (identical(fam, "multinomial"))
    stop("a density overlay has no meaning for an unordered categorical ",
         "response; use ilm_calibration() instead", call. = FALSE)
  y <- as.numeric(object$y)
  ys <- ilm_sim_cond(object, B, seed)

  ## Where does the observed mean sit among the simulated means? A value near
  ## 0 or 1 says the model generates data centred somewhere else entirely.
  sm <- colMeans(ys, na.rm = TRUE)
  q_mean <- mean(sm <= mean(y, na.rm = TRUE))
  sv <- apply(ys, 2L, stats::var, na.rm = TRUE)
  q_var <- mean(sv <= stats::var(y, na.rm = TRUE))
  extreme <- function(q) min(q, 1 - q) < 0.025

  if (plot) {
    disc <- all(abs(y - round(y)) < 1e-8, na.rm = TRUE)
    if (disc) {
      rng <- 0:max(y, na.rm = TRUE)
      obs <- tabulate(match(y, rng), length(rng)) / length(y)
      sim <- vapply(seq_len(B), function(b)
        tabulate(match(ys[, b], rng), length(rng)) / nrow(ys), numeric(length(rng)))
      graphics::matplot(rng, sim, type = "l", lty = 1,
                        col = grDevices::adjustcolor("grey50", 0.35),
                        xlab = "value", ylab = "proportion",
                        main = "Observed against simulated")
      graphics::lines(rng, obs, col = "black", lwd = 2.5)
    } else {
      do <- stats::density(y[is.finite(y)])
      graphics::plot(do, main = "Observed against simulated", xlab = "value",
                     lwd = 2.5, ylim = c(0, max(do$y) * 1.35))
      for (b in seq_len(B)) {
        db <- try(stats::density(ys[is.finite(ys[, b]), b]), silent = TRUE)
        if (!inherits(db, "try-error"))
          graphics::lines(db, col = grDevices::adjustcolor("grey50", 0.35))
      }
      graphics::lines(do, lwd = 2.5)
    }
    graphics::legend("topright", c("observed", "simulated"), lwd = c(2.5, 1),
                     col = c("black", "grey50"), bty = "n", cex = 0.8)
  }

  status <- if (extreme(q_mean) || extreme(q_var)) "FAIL"
            else if (min(q_mean, 1 - q_mean, q_var, 1 - q_var) < 0.10) "WARN"
            else "OK"
  res <- list(q_mean = q_mean, q_var = q_var, B = B, status = status,
              note = if (status == "OK") ""
                     else paste("the observed", if (extreme(q_mean)) "mean" else "variance",
                                "is outside what this model generates; check the",
                                "family and the linear predictor"))
  if (status != "OK") message("ilm_check_predictive: ", res$note)
  invisible(res)
}

## ---- dispersion ------------------------------------------------------------

#' Is the response more variable than the model allows?
#'
#' The Pearson statistic divided by its degrees of freedom, compared against
#' the same quantity computed on data simulated from the fitted model. The
#' reference distribution is simulated rather than assumed, so the test does
#' not rely on the chi-square approximation that fails for sparse counts.
#'
#' @param object A fitted `"ilm_model"` object.
#' @param B Simulated datasets. A simulated p-value cannot fall below
#'   `1 / (B + 1)`, so B must exceed 100 for the strongest verdict to be
#'   reachable at all.
#' @param seed Random seed.
#' @return Invisibly, a list with the observed ratio, the simulated
#'   distribution, a p-value and a `status`.
#' @seealso [ilm_check_zeros()] when the excess variance is concentrated at zero.
#' @examples
#' set.seed(1)
#' d <- ilm_sim()
#' f <- ilm_model(claims ~ income + (1 | id), data = d, family = "poisson",
#'                verbose = FALSE)
#' ilm_check_dispersion(f, B = 200)
#' @export
ilm_check_dispersion <- function(object, B = 200L, seed = 1L) {
  if (!inherits(object, "ilm_model"))
    stop("`object` must be a fitted ilm_model object, not ", class(object)[1],
         call. = FALSE)
  fam <- if (!is.null(object$family)) object$family$name else "gaussian"
  if (!fam %in% c("poisson", "nbinom", "binomial"))
    stop("a dispersion check applies to poisson, nbinom and binomial models; ",
         "for ", fam, " the residual variance is estimated rather than fixed, ",
         "so there is nothing to test", call. = FALSE)

  ## A simulated p-value cannot fall below 1 / (B + 1). At B = 60 the floor is
  ## 0.016, so the FAIL threshold of 0.01 is unreachable and even gross
  ## misspecification can only ever report WARN.
  if (1 / (B + 1) > 0.01)
    warning("B = ", B, " puts the smallest achievable p-value at ",
            signif(1 / (B + 1), 2), ", so a FAIL verdict is unreachable. ",
            "Use B >= 100.", call. = FALSE)
  mu <- as.numeric(ilm_fitted(object, TRUE)[, 1])
  vf <- switch(fam,
    poisson  = mu,
    binomial = mu * (1 - mu),
    nbinom   = { k <- unname(object$dispersion[1]); mu + mu^2 / k })
  pear <- function(yy) sum((yy - mu)^2 / vf, na.rm = TRUE)
  dfres <- length(mu) - length(stats::coef(object))
  obs <- pear(as.numeric(object$y)) / dfres

  ys <- ilm_sim_cond(object, B, seed)
  sim <- apply(ys, 2L, function(z) pear(as.numeric(z)) / dfres)
  p <- (1 + sum(sim >= obs)) / (B + 1)      # one-sided: overdispersion

  status <- if (p < 0.01) "FAIL" else if (p < 0.05) "WARN" else "OK"
  res <- list(ratio = round(obs, 3), sim_median = round(stats::median(sim), 3),
              p = round(p, 4), B = B, status = status,
              note = if (status == "OK") "" else
                paste0("more variable than ", fam, " allows",
                       if (fam == "poisson") "; try family = \"nbinom\"" else "",
                       ", or check for an omitted predictor"))
  if (status != "OK") message("ilm_check_dispersion: ", res$note)
  invisible(res)
}

## ---- zero inflation --------------------------------------------------------

#' Are there more zeros than the model expects?
#'
#' More zeros than a count model expects is the one misspecification that
#' cannot be read off a residual plot, because a zero is a perfectly ordinary
#' value of a count. This simulates from the fit and asks how often it produces
#' as many zeros as were seen.
#'
#' The remedy is [ilm_model(ziformula = )][ilm_model], which adds a second
#' linear predictor for the probability of an excess zero. Whether that should
#' be a mixture (`zi_type = "inflated"`) or two processes (`"hurdle"`) is a
#' question about what the zeros mean rather than one this check can answer.
#'
#' @param object A fitted `"ilm_model"` object.
#' @param B Simulated datasets. A simulated p-value cannot fall below
#'   `1 / (B + 1)`, so B must exceed 100 for the strongest verdict to be
#'   reachable at all.
#' @param seed Random seed.
#' @return Invisibly, a list with the observed and expected zero counts, a
#'   p-value and a `status`.
#' @examples
#' set.seed(1)
#' d <- ilm_sim()
#' f <- ilm_model(downtime ~ income + (1 | id), data = d, family = "poisson",
#'                verbose = FALSE)
#' ilm_check_zeros(f, B = 200)
#' @export
ilm_check_zeros <- function(object, B = 200L, seed = 1L) {
  if (!inherits(object, "ilm_model"))
    stop("`object` must be a fitted ilm_model object, not ", class(object)[1],
         call. = FALSE)
  fam <- if (!is.null(object$family)) object$family$name else "gaussian"
  ## A zero means something when the response distribution gives it its own
  ## probability -- a count -- or when the model has a part devoted to it. A
  ## plain beta fit is neither: it could not have been fitted to data
  ## containing a zero in the first place.
  if (!fam %in% c("poisson", "nbinom") && is.null(object$Zzi))
    stop("a zero-inflation check applies to count models (poisson, nbinom), ",
         "or to any fit with a zero part. For a ", fam, " model without one, ",
         "a zero carries no special meaning -- and a beta fit cannot contain ",
         "one at all, since its density has no mass there.", call. = FALSE)
  ## A simulated p-value cannot fall below 1 / (B + 1). At B = 60 the floor is
  ## 0.016, so the FAIL threshold of 0.01 is unreachable and even gross
  ## misspecification can only ever report WARN.
  if (1 / (B + 1) > 0.01)
    warning("B = ", B, " puts the smallest achievable p-value at ",
            signif(1 / (B + 1), 2), ", so a FAIL verdict is unreachable. ",
            "Use B >= 100.", call. = FALSE)
  y <- as.numeric(object$y)
  obs <- sum(y == 0, na.rm = TRUE)
  ys <- ilm_sim_cond(object, B, seed)
  sim <- colSums(ys == 0, na.rm = TRUE)
  p <- (1 + sum(sim >= obs)) / (B + 1)      # one-sided: excess zeros

  status <- if (p < 0.01) "FAIL" else if (p < 0.05) "WARN" else "OK"
  res <- list(observed = obs, expected = round(mean(sim), 1),
              p = round(p, 4), B = B, status = status,
              note = if (status == "OK") "" else paste0(
                "observed ", obs, " zeros against about ", round(mean(sim)),
                " expected; refit with ilm_model(ziformula = ~ 1), or ",
                "ziformula = ~ x if the excess depends on a predictor. Use ",
                "zi_type = \"hurdle\" when every zero comes from one process ",
                "and \"inflated\" when some units were never at risk"))
  if (status != "OK") message("ilm_check_zeros: ", res$note)
  invisible(res)
}

## ---- binned residuals ------------------------------------------------------

#' Binned residuals for a binary response
#'
#' A raw residual from a 0/1 outcome takes one of two values and tells you
#' almost nothing. Averaging residuals within bins of fitted probability makes
#' the pattern visible: points outside the band mark regions of the predictor
#' space where the model is systematically wrong.
#'
#' @param object A fitted `"ilm_model"` object with a binomial family.
#' @param nbins Number of bins. The default follows the usual `sqrt(n)` rule.
#' @param plot Draw the plot.
#' @return Invisibly, a data frame of bin midpoints, mean residuals, bounds and
#'   a `status`.
#' @references
#' Gelman, A. and Hill, J. (2007). Data Analysis Using Regression and
#' Multilevel/Hierarchical Models. Cambridge University Press, ch. 5.
#' @examples
#' set.seed(1)
#' d <- ilm_sim()
#' f <- ilm_model(flag ~ income + (1 | id), data = d, family = "binomial",
#'                verbose = FALSE)
#' ilm_binned_residuals(f)
#' @export
ilm_binned_residuals <- function(object, nbins = NULL, plot = TRUE) {
  if (!inherits(object, "ilm_model"))
    stop("`object` must be a fitted ilm_model object, not ", class(object)[1],
         call. = FALSE)
  fam <- if (!is.null(object$family)) object$family$name else "gaussian"
  if (!identical(fam, "binomial"))
    stop("binned residuals are for a binary response; for ", fam,
         " use ilm_rqr() or the panels in ilm_appraise()", call. = FALSE)
  p <- as.numeric(ilm_fitted(object, TRUE)[, 1])
  y <- as.numeric(object$y)
  r <- y - p
  n <- length(p)
  if (is.null(nbins)) nbins <- max(5L, min(100L, floor(sqrt(n))))

  br <- unique(stats::quantile(p, seq(0, 1, length.out = nbins + 1L),
                               names = FALSE))
  if (length(br) < 3L)
    stop("the fitted probabilities barely vary, so there is nothing to bin",
         call. = FALSE)
  bin <- cut(p, br, include.lowest = TRUE)
  mp <- tapply(p, bin, mean)
  mr <- tapply(r, bin, mean)
  nb <- tapply(r, bin, length)
  ## two standard errors of a bin mean, the conventional band
  se <- 2 * sqrt(tapply(r, bin, stats::var) / nb)
  outside <- sum(abs(mr) > se, na.rm = TRUE)
  k <- sum(!is.na(mr))
  ptail <- stats::pbinom(outside - 1L, k, 0.05, lower.tail = FALSE)

  out <- data.frame(mean_p = as.numeric(mp), mean_resid = as.numeric(mr),
                    n = as.integer(nb), lower = -as.numeric(se),
                    upper = as.numeric(se), row.names = NULL)
  status <- if (ptail < 0.05) "FAIL" else if (ptail < 0.20) "WARN" else "OK"
  attr(out, "status") <- status
  attr(out, "n_outside") <- outside

  if (plot) {
    graphics::plot(out$mean_p, out$mean_resid, pch = 19,
                   ylim = range(c(out$lower, out$upper, out$mean_resid),
                                na.rm = TRUE),
                   xlab = "fitted probability", ylab = "mean residual",
                   main = sprintf("Binned residuals (%d of %d bins outside)",
                                  outside, k))
    graphics::abline(h = 0, lty = 2, col = "grey50")
    graphics::lines(out$mean_p, out$lower, col = "grey40")
    graphics::lines(out$mean_p, out$upper, col = "grey40")
  }
  invisible(out)
}

## ---- homoscedasticity ------------------------------------------------------

## Two ways the variance can be wrong, each with a different remedy:
## a spread that trends with the fitted value, and a spread that differs
## between groups. Both are measured on randomised quantile residuals, which
## are standard normal when the model -- INCLUDING its variance -- is right, so
## any remaining pattern in their spread is variance misspecification.
ilm_var_stats <- function(fit, by_vec, seed) {
  u <- ilm_rqr(fit, TRUE, seed)
  z <- stats::qnorm(pmin(pmax(u, 1e-6), 1 - 1e-6))
  mu <- as.numeric(ilm_fitted(fit, TRUE)[, 1])
  ok <- is.finite(z) & is.finite(mu)
  ## Spearman, so a curved mean-variance relationship is still caught and no
  ## scale assumption is smuggled in
  trend <- if (sum(ok) > 10L)
    suppressWarnings(stats::cor(mu[ok], abs(z[ok]), method = "spearman")) else NA_real_
  ratio <- NA_real_
  if (!is.null(by_vec)) {
    v <- tapply(z[ok], by_vec[ok], stats::var)
    v <- v[is.finite(v) & v > 0]
    if (length(v) > 1L) ratio <- max(v) / min(v)
  }
  c(trend = trend, ratio = ratio)
}

#' Is the residual spread constant?
#'
#' The assumption illume had no check for. Heteroscedasticity does not usually
#' bias the coefficients much, but it does bias their standard errors, which is
#' what an inference package is for.
#'
#' Two patterns are tested, because they call for different fixes: spread that
#' trends with the fitted value, and spread that differs between the levels of
#' a grouping variable. Both are measured on randomised quantile residuals,
#' which are standard normal when the model and its variance are both correct,
#' and both are compared against a null simulated from the fitted model rather
#' than against an assumed distribution.
#'
#' @param object A fitted `"ilm_model"` object.
#' @param by Optional name of a column in the model frame, or a vector, whose
#'   levels might have different variances.
#' @param B Simulated datasets. A simulated p-value cannot fall below
#'   `1 / (B + 1)`, so B must exceed 100 for the strongest verdict to be
#'   reachable.
#' @param seed Random seed.
#' @param ncores Worker processes for the refits.
#' @param plot Draw the scale-location panel.
#' @param verbose Print the verdict.
#' @param progress Show a progress bar. Defaults to [interactive()], so a
#'   bar appears when someone is watching and nothing is written in a
#'   script or a knitted document. See [ilm_progress_arg].
#' @return Invisibly, a list with the observed statistics, their simulated
#'   nulls, p-values, a `status` and a suggested remedy.
#' @seealso [ilm_rqr_test()] for other residual statistics,
#'   [ilm_check_dispersion()] when the whole response is over-dispersed rather
#'   than unevenly dispersed.
#' @examples
#' set.seed(1)
#' d <- ilm_sim(n_id = 25)
#' f <- ilm_model(score ~ income + (1 | id), data = d, family = "gaussian",
#'                verbose = FALSE)
#' ilm_check_variance(f, B = 20, plot = FALSE)
#' @export
ilm_check_variance <- function(object, by = NULL, B = 100L, seed = 1L,
                               ncores = 1L, plot = TRUE, verbose = TRUE,
                               progress = NULL) {
  if (!inherits(object, "ilm_model"))
    stop("`object` must be a fitted ilm_model object, not ", class(object)[1],
         call. = FALSE)
  if (1 / (B + 1) > 0.01)
    warning("B = ", B, " puts the smallest achievable p-value at ",
            signif(1 / (B + 1), 2), ", so a FAIL verdict is unreachable. ",
            "Use B >= 100.", call. = FALSE)

  by_vec <- NULL; by_lab <- NULL
  if (!is.null(by)) {
    if (length(by) == 1L && is.character(by)) {
      mf <- object$model
      if (is.null(mf) || !by %in% names(mf))
        stop("`by` (", by, ") is not a column of the model frame. Available: ",
             paste(names(mf), collapse = ", "), call. = FALSE)
      by_vec <- factor(mf[[by]]); by_lab <- by
    } else {
      if (length(by) != nrow(object$X))
        stop("`by` has ", length(by), " values but the model has ",
             nrow(object$X), " rows", call. = FALSE)
      by_vec <- factor(by); by_lab <- "by"
    }
  }

  obs <- ilm_var_stats(object, by_vec, seed)
  ys <- ilm_sim_cond(object, B, seed + 1L)
  ## the same model, every part of it -- see ilm_refit_stub()
  stub <- ilm_refit_stub(object)
  asg <- object$assign; tl <- object$term_labels
  cl <- ilm_pool(ncores)
  on.exit(if (!is.null(cl)) try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
  if (!is.null(cl)) try(parallel::clusterExport(cl, "by_vec", envir = environment()),
                        silent = TRUE)
  one <- function(b) {
    f <- try(ilm_refit_like(stub, y = ys[, b], restarts = 1L), silent = TRUE)
    if (inherits(f, "try-error") || f$opt$convergence != 0)
      return(c(trend = NA_real_, ratio = NA_real_))
    f$assign <- asg; f$term_labels <- tl
    ilm_var_stats(f, by_vec, seed + 1000L + b)
  }
  nullm <- do.call(rbind,
                   ilm_lapply_progress(cl, seq_len(B), one,
                                       progress = progress))
  nok <- sum(is.finite(nullm[, "trend"]))
  if (nok < 10L) {
    if (verbose) cat("INCONCLUSIVE: only", nok, "of", B, "replicates refitted.\n")
    return(invisible(list(status = "INCONCLUSIVE", observed = obs, null = nullm)))
  }

  ## two-sided for the trend, one-sided for the ratio, which cannot be small
  pv <- function(o, nl, sided) {
    nl <- nl[is.finite(nl)]
    if (!is.finite(o) || !length(nl)) return(NA_real_)
    if (sided == 2) {
      m <- stats::median(nl)
      (1 + sum(abs(nl - m) >= abs(o - m))) / (length(nl) + 1)
    } else (1 + sum(nl >= o)) / (length(nl) + 1)
  }
  p_trend <- pv(obs[["trend"]], nullm[, "trend"], 2)
  p_ratio <- pv(obs[["ratio"]], nullm[, "ratio"], 1)
  worst <- suppressWarnings(min(c(p_trend, p_ratio), na.rm = TRUE))
  if (!is.finite(worst)) worst <- NA_real_
  status <- if (is.na(worst)) "INCONCLUSIVE" else
            if (worst < 0.01) "FAIL" else if (worst < 0.05) "WARN" else "OK"

  ## the remedy depends on which pattern fired
  fix <- character(0)
  ## Both remedies are a dispersion formula, which is why the wording names it.
  if (isTRUE(p_trend < 0.05))
    fix <- c(fix, "spread changes with the fitted value: dispformula = ~ mu, or a family whose variance grows with the mean")
  if (isTRUE(p_ratio < 0.05))
    fix <- c(fix, sprintf("spread differs across %s: dispformula = ~ %s",
                          if (is.null(by_lab)) "groups" else by_lab,
                          if (is.null(by_lab)) "group" else by_lab))

  if (plot) {
    u <- ilm_rqr(object, TRUE, seed)
    z <- stats::qnorm(pmin(pmax(u, 1e-6), 1 - 1e-6))
    mu <- as.numeric(ilm_fitted(object, TRUE)[, 1])
    ok <- is.finite(z) & is.finite(mu)
    graphics::plot(mu[ok], sqrt(abs(z[ok])), pch = 19,
                   col = grDevices::adjustcolor("black", 0.35),
                   xlab = "fitted", ylab = expression(sqrt(abs(residual))),
                   main = sprintf("Scale-location (%s)", status))
    if (sum(ok) > 20L) {
      lw <- try(stats::lowess(mu[ok], sqrt(abs(z[ok]))), silent = TRUE)
      if (!inherits(lw, "try-error")) graphics::lines(lw, col = "red", lwd = 2)
    }
  }

  res <- list(trend = unname(obs[["trend"]]), p_trend = round(p_trend, 4),
              ratio = unname(obs[["ratio"]]), p_ratio = round(p_ratio, 4),
              n_refits = nok, B = B, status = status,
              suggestion = paste(fix, collapse = "; "))
  if (verbose) {
    cat(sprintf("spread vs fitted: rho = %.3f, p = %s\n", res$trend,
                format(res$p_trend)))
    if (is.finite(res$ratio))
      cat(sprintf("variance ratio across %s: %.2f, p = %s\n",
                  by_lab, res$ratio, format(res$p_ratio)))
    cat(status, if (nzchar(res$suggestion)) paste0(" -- ", res$suggestion) else "", "\n", sep = "")
  }
  invisible(res)
}
