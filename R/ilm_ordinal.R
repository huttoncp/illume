## ---------------------------------------------------------------------------
## Ordered outcomes.
##
## An ordered response treated as unordered throws away the ordering and pays
## for it in parameters: a multinomial spends J - 1 coefficients per predictor
## where a cumulative link model spends one. Treated as numeric it invents a
## spacing the categories do not have -- "agree" is not twice "neutral".
##
## The cumulative link model sits between the two. See ilm_ordinal_family().
## ---------------------------------------------------------------------------

#' Names for the thresholds of a cumulative link fit
#'
#' @keywords internal
#' @noRd
ilm_ord_cut_names <- function(ylevels, J) {
  lv <- if (!is.null(ylevels) && length(ylevels) == J) ylevels else
    as.character(seq_len(J))
  paste(lv[-J], lv[-1L], sep = "|")
}

#' Thresholds of a cumulative link fit
#'
#' The cut points on the latent scale, with their standard errors. There are
#' `J - 1` of them for `J` categories and they are increasing by construction:
#' the fit estimates the first and the logarithms of the gaps, so an ordering
#' violation is not something the optimiser can reach.
#'
#' Read them against the linear predictor, which is **subtracted**:
#' `P(Y <= j) = F(threshold_j - eta)`. A positive coefficient therefore pushes
#' probability towards the higher categories.
#'
#' @param object A fitted [ilm_model()] with an ordinal family.
#' @param level Confidence level.
#' @return A data frame of `cut`, `estimate`, `se`, `lower` and `upper`.
#' @seealso [ilm_check_proportional()] for the assumption that makes one
#'   coefficient per predictor enough.
#' @examples
#' set.seed(1); n <- 400
#' d <- data.frame(x = rnorm(n))
#' z <- 0.8 * d$x + rlogis(n)
#' d$y <- factor(cut(z, c(-Inf, -1, 1, Inf), labels = c("lo", "mid", "hi")),
#'               ordered = TRUE)
#' f <- ilm_model(y ~ x, data = d, family = "ordinal", verbose = FALSE)
#' ilm_thresholds(f)
#' @export
ilm_thresholds <- function(object, level = 0.95) {
  if (!inherits(object, "ilm_model"))
    stop("`object` must be a fitted ilm_model, not ", class(object)[1],
         call. = FALSE)
  if (is.null(object$zeta))
    stop("this model has no thresholds: they belong to an ordered response. ",
         "Fit one with family = \"ordinal\".", call. = FALSE)
  est <- object$zeta
  se <- object$zeta_se
  if (is.null(se)) se <- rep(NA_real_, length(est))
  crit <- stats::qnorm(1 - (1 - level) / 2)
  data.frame(cut = names(est), estimate = unname(est), se = unname(se),
             lower = unname(est) - crit * se, upper = unname(est) + crit * se,
             row.names = NULL)
}

#' Category probabilities from a cumulative link fit
#'
#' @param object A fitted [ilm_model()] with an ordinal family.
#' @param newdata Optional data to predict for.
#' @return An `N x J` matrix of probabilities whose rows sum to one.
#' @keywords internal
#' @noRd
ilm_ord_pmat <- function(object, eta) {
  ilm_ord_probs(eta, object$zeta, object$family$pfun)
}

#' Drop the intercept column from a rebuilt design, for an ordinal fit
#'
#' A cumulative link fit has thresholds where every other family has an
#' intercept, so its design carries no column of ones. Anything that rebuilds
#' a design from the model's terms -- predictions for new data, a reference
#' grid, an average marginal effect -- gets one back from `model.matrix()` and
#' has to drop it, or the columns no longer line up with the coefficients.
#'
#' @keywords internal
#' @noRd
ilm_drop_intercept <- function(X, object) {
  if (!isTRUE(object$ordinal)) return(X)
  ic <- match("(Intercept)", colnames(X), nomatch = 0L)
  if (ic > 0L) X <- X[, -ic, drop = FALSE]
  X
}

#' Does one coefficient per predictor describe every cut?
#'
#' A cumulative link model says a predictor shifts all `J - 1` cut points by
#' the same amount. Under a logit link that is the proportional-odds
#' assumption, and it is what makes the model economical -- and wrong, when it
#' does not hold, in a way that is invisible in the coefficient table.
#'
#' This is the Brant idea, done by refitting rather than by a score test. Each
#' cut becomes a binary outcome (`y > j`), a separate model is fitted to each,
#' and the statistic is how far the resulting coefficients spread apart
#' relative to their standard errors.
#'
#' @section Why the reference is simulated:
#'
#' Those separate fits share their data, so their estimates are correlated,
#' and the spread statistic is not chi-squared however tempting the shape.
#' Treating it as one was measured here at a **31.5% false-alarm rate** for a
#' binary predictor whose effect was perfectly proportional -- with two cuts
#' the statistic is `(b1 - b2)^2` over an average variance, where the honest
#' denominator is `v1 + v2 - 2c`, and the covariance `c` that is missing
#' depends on the predictor.
#'
#' The reference is therefore simulated from the fitted proportional-odds
#' model, which carries that correlation because it arises from the same
#' refitting. Over 400 replicates per design, with the assumption holding:
#'
#' ```
#'                    flagged at 0.05        p-values under the null
#'   n = 400   x 0.032, g 0.040          p<.05  0.042 / 0.056
#'   n = 1500  x 0.052, g 0.072          p<.25  0.255 / 0.244
#'                                       p<.50  0.510 / 0.484
#' ```
#'
#' and with `x` acting on the first cut only, the offending term is flagged in
#' 400 of 400 samples at both sizes, while the term that does act
#' proportionally is flagged 1.5% and 2.8% of the time -- it localises the
#' violation rather than condemning the model. The 7.2% at the largest size is
#' a little above nominal and is reported here rather than rounded away.
#'
#' @section What to do when it fails:
#'
#' The remedy is `family = "multinomial"`, which spends `J - 1` coefficients
#' per predictor and assumes nothing about the ordering. It nests this model,
#' so [AIC()] will say whether the extra parameters earn their place. When
#' only one or two terms offend, that is the trade to weigh: the ordering is
#' real information, and giving all of it up to accommodate one predictor may
#' cost more than it buys.
#'
#' @param object A fitted [ilm_model()] with an ordinal family.
#' @param B Simulated datasets for the reference. A simulated p-value cannot
#'   fall below `1 / (B + 1)`.
#' @param alpha Level for the per-term verdicts.
#' @param seed Random seed.
#' @param progress Show a progress bar; see [illumex::ilm_progress_arg].
#' @return Invisibly, a data frame with one row per fixed-effect term giving
#'   the spread statistic, a simulated p-value, the largest gap between cuts
#'   and a `status`, with the per-cut coefficients in a `coefs` attribute.
#' @references Brant, R. (1990). Assessing proportionality in the proportional
#'   odds model for ordinal logistic regression. *Biometrics* 46, 1171-1178.
#' @seealso [ilm_thresholds()], [ilm_model()].
#' @examples
#' set.seed(1); n <- 500
#' d <- data.frame(x = rnorm(n))
#' z <- 0.8 * d$x + rlogis(n)
#' d$y <- factor(cut(z, c(-Inf, -1, 1, Inf), labels = c("lo", "mid", "hi")),
#'               ordered = TRUE)
#' f <- ilm_model(y ~ x, data = d, family = "ordinal", verbose = FALSE)
#' ilm_check_proportional(f, B = 99)
#' @export
ilm_check_proportional <- function(object, B = 199L, alpha = 0.05, seed = 1L,
                                   progress = NULL) {
  if (!inherits(object, "ilm_model"))
    stop("`object` must be a fitted ilm_model, not ", class(object)[1],
         call. = FALSE)
  if (is.null(object$zeta))
    stop("the proportional-odds assumption belongs to an ordered response; ",
         "this model has no thresholds. See family = \"ordinal\".",
         call. = FALSE)
  if (length(object$re))
    stop("this check refits the model one cut at a time as a binary outcome, ",
         "and does not carry the random effects across. Test the assumption ",
         "on a fixed-effects fit of the same mean structure, or compare with ",
         "family = \"multinomial\" using AIC().", call. = FALSE)
  ## A simulated p-value cannot fall below 1 / (B + 1), and the strongest
  ## verdict needs it below alpha / 5. At B = 99 that floor is exactly 0.01
  ## against a threshold of 0.01, so a FAIL is unreachable and a flagrant
  ## violation can only ever read WARN.
  if (1 / (B + 1) >= alpha / 5)
    warning("B = ", B, " puts the smallest achievable p-value at ",
            signif(1 / (B + 1), 2), ", so a FAIL verdict is unreachable ",
            "(it needs p < ", signif(alpha / 5, 2),
            ") and even a flagrant violation can only read WARN. Use B >= ",
            ceiling(5 / alpha), ".", call. = FALSE)

  J <- object$J; X <- object$X
  w <- if (is.null(object$weights)) rep(1, nrow(X)) else object$weights
  asg <- object$assign; tl <- object$term_labels
  rows <- lapply(seq_along(tl), function(k) which(!is.na(asg) & asg == k))
  keep <- which(vapply(rows, length, 0L) > 0L)
  if (!length(keep))
    stop("no parametric term to test: every column here belongs to a smooth, ",
         "and a smooth has no single coefficient to compare across cuts.",
         call. = FALSE)

  obs <- ilm_prop_stat(as.integer(object$y), X, w, J, rows[keep])
  if (all(is.na(obs$stat)))
    stop("no cut could be fitted as a binary outcome, so there is nothing to ",
         "compare. A category with very few observations is the usual cause.",
         call. = FALSE)

  ## The null is the fitted proportional-odds model itself, so draws from it
  ## carry both the correlation between cuts and whatever the design does to
  ## it -- neither of which a chi-squared reference knows about.
  ys <- ilm_sim_cond(object, B, seed)
  pb <- ilm_progress(B, progress, "simulating the proportional-odds reference")
  ref <- matrix(NA_real_, B, length(keep))
  for (b in seq_len(B)) {
    ref[b, ] <- ilm_prop_stat(ys[, b], X, w, J, rows[keep])$stat
    pb$tick(b)
  }
  pb$done()
  p <- vapply(seq_along(keep), function(i) {
    r <- ref[, i]; r <- r[is.finite(r)]
    if (!length(r) || !is.finite(obs$stat[i])) return(NA_real_)
    (1 + sum(r >= obs$stat[i])) / (length(r) + 1)
  }, 0)

  res <- data.frame(term = tl[keep], stat = obs$stat, p = p,
                    max_gap = obs$gap,
                    status = ifelse(is.na(p), "INCONCLUSIVE",
                             ifelse(p < alpha / 5, "FAIL",
                             ifelse(p < alpha, "WARN", "OK"))),
                    stringsAsFactors = FALSE, row.names = NULL)
  bad <- res[res$status %in% c("FAIL", "WARN"), , drop = FALSE]
  if (nrow(bad))
    message("ilm_check_proportional: ", paste(bad$term, collapse = ", "),
            if (nrow(bad) > 1L) " do not act " else " does not act ",
            "the same way on every cut (smallest p = ",
            signif(min(bad$p), 3), "). One coefficient per predictor is then ",
            "an average over effects that differ by category. The remedy is ",
            "family = \"multinomial\", which nests this model -- compare them ",
            "with AIC(), since dropping the ordering costs ", J - 2L,
            " extra coefficient(s) per predictor.")
  structure(res, class = c("ilm_prop_check", "data.frame"),
            coefs = obs$B, B = B)
}

#' Spread of the per-cut binary coefficients, one number per term
#'
#' Splitting an ordered response at each cut and fitting a binary model to
#' each gives `J - 1` estimates of every coefficient. Under proportional odds
#' they estimate the same thing.
#'
#' @keywords internal
#' @noRd
ilm_prop_stat <- function(y, X, w, J, rows) {
  p <- ncol(X)
  Bm <- matrix(NA_real_, J - 1L, p)
  Vm <- matrix(NA_real_, J - 1L, p)
  Xi <- cbind(`(Intercept)` = 1, X)
  ## Each cut is an ordinary logistic regression with no random effects, so
  ## it goes through IRLS rather than through ilm_fit(): the estimates are the
  ## same maximum likelihood ones to eight figures, and building a TMB tape
  ## for each of 2 * B refits would make the reference below cost minutes
  ## instead of seconds. That difference decides whether the check gets run.
  for (j in seq_len(J - 1L)) {
    yj <- as.numeric(y > j)
    ## a cut with no variation carries no information about any coefficient
    if (length(unique(yj)) < 2L) next
    fj <- try(suppressWarnings(
      stats::glm.fit(Xi, yj, weights = w,
                     family = stats::binomial())), silent = TRUE)
    if (inherits(fj, "try-error") || !isTRUE(fj$converged)) next
    Bm[j, ] <- unname(fj$coefficients)[-1L]
    ## the inverse of the weighted cross-product is the usual covariance, and
    ## the QR factor glm.fit already has is the cheap way to it
    vv <- try(chol2inv(fj$qr$qr[seq_len(fj$rank), seq_len(fj$rank),
                                drop = FALSE]), silent = TRUE)
    if (!inherits(vv, "try-error") && fj$rank == ncol(Xi)) {
      dv <- diag(vv)[order(fj$qr$pivot[seq_len(fj$rank)])]
      Vm[j, ] <- dv[-1L]
    }
  }
  ok <- which(!is.na(Bm[, 1]))
  stat <- rep(NA_real_, length(rows)); gap <- rep(NA_real_, length(rows))
  if (length(ok) >= 2L) for (i in seq_along(rows)) {
    cols <- rows[[i]]
    bb <- Bm[ok, cols, drop = FALSE]; vv <- Vm[ok, cols, drop = FALSE]
    ## Standardised by each estimate's OWN variance, not an average of them.
    ## The scale is arbitrary -- the simulated reference sets it -- but
    ## dividing by a common number would let the cut with the least
    ## information dominate a term whose other cuts are precisely estimated.
    cen <- sweep(bb, 2L, colMeans(bb), `-`)
    vv[!is.finite(vv) | vv <= 0] <- NA_real_
    stat[i] <- sum(cen^2 / vv, na.rm = TRUE)
    gap[i] <- max(apply(bb, 2L, function(z) diff(range(z))))
  }
  list(stat = stat, gap = gap, B = Bm)
}

#' @export
print.ilm_prop_check <- function(x, ...) {
  cat("Proportional-odds check\n")
  cat(sprintf("  reference simulated from the fitted model, B = %d\n\n",
              attr(x, "B")))
  d <- as.data.frame(x)
  d$stat <- round(d$stat, 2); d$p <- signif(d$p, 3)
  d$max_gap <- round(d$max_gap, 3)
  print(d, row.names = FALSE)
  if (any(d$status != "OK"))
    cat("\n  A term flagged here shifts the cuts by different amounts, so one\n",
        " coefficient is an average over effects that differ by category.\n",
        " family = \"multinomial\" drops the assumption; AIC() says whether\n",
        " the extra coefficients earn their place.\n", sep = "")
  invisible(x)
}
