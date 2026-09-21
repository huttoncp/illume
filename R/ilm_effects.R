## ---------------------------------------------------------------------------
## Effect sizes, on the scale the family puts them on.
##
## A coefficient from a log link is a log rate ratio and from a logit link a
## log odds ratio, so the number people want is the exponential of it. That is
## arithmetic. Three things around it are not.
##
## THE INTERVAL IS BUILT ON THE LINK SCALE AND THEN TRANSFORMED. exp() of the
## endpoints is exact and asymmetric, which is the shape a ratio's uncertainty
## actually has. A delta-method standard error on the ratio scale, plus or
## minus two of them, is symmetric, can reach below zero, and is what a
## surprising amount of software reports.
##
## A RATIO FROM A MIXED MODEL IS CONDITIONAL. exp(beta) there is the odds ratio
## for a GIVEN cluster -- two subjects in the same hospital, one exposed -- and
## not for the population. The two differ by more the larger the random
## effects, and the conditional one is routinely reported as if it were the
## marginal one. ilm_ame() has the marginal quantity; this prints both and says
## which is which.
##
## A STANDARDISED COEFFICIENT IS NOT AN EFFECT SIZE FOR A BINARY PREDICTOR.
## Dividing by the standard deviation of a 0/1 variable answers "what happens
## per standard deviation of treatment", which is not a question anyone has.
## Those are reported unstandardised, and the reason is given rather than the
## number quietly omitted.
##
## References:
##   Hosmer, D. W. and Lemeshow, S. (2013). Applied Logistic Regression, 3rd ed.
##   Gelman, A. (2008). Scaling regression inputs by dividing by two standard
##     deviations. Statistics in Medicine 27, 2865-2873.
## ---------------------------------------------------------------------------

#' What a coefficient means on this family's scale
#'
#' @keywords internal
#' @noRd
ilm_effect_kind <- function(object) {
  fam <- if (is.null(object$family)) "gaussian" else object$family$name
  if (isTRUE(object$ordinal))
    return(list(scale = "odds ratio", label = "proportional odds ratio",
                exp = TRUE,
                note = "the odds of being in a HIGHER category, and by the proportional-odds assumption the same ratio at every cut -- ilm_check_proportional() tests that"))
  switch(fam,
    binomial = list(scale = "odds ratio", label = "odds ratio", exp = TRUE,
                    note = "the odds of the modelled outcome"),
    poisson  = list(scale = "rate ratio", label = "incidence rate ratio",
                    exp = TRUE, note = "the expected count"),
    nbinom   = list(scale = "rate ratio", label = "incidence rate ratio",
                    exp = TRUE, note = "the expected count"),
    beta     = list(scale = "odds ratio", label = "odds ratio", exp = TRUE,
                    note = "the odds of the proportion, since the link is logit"),
    weibull = , lognormal = , loglogistic =
      list(scale = "time ratio", label = "time ratio", exp = TRUE,
           note = "survival time: above 1 means longer, and this is an acceleration factor rather than a hazard ratio"),
    rp = list(scale = "hazard ratio", label = "hazard ratio", exp = TRUE,
              note = "the hazard, held proportional over time by this baseline"),
    rp_odds = list(scale = "odds ratio", label = "odds ratio of failure",
                   exp = TRUE,
                   note = "the odds of failure by a given time, held proportional"),
    rp_normal = list(scale = "probit shift", label = "probit shift",
                     exp = FALSE,
                     note = "a shift on the probit scale, which does not exponentiate into a ratio"),
    list(scale = "difference", label = "difference", exp = FALSE,
         note = "the response itself"))
}

#' Effect sizes for the fixed effects
#'
#' Reports each coefficient on the scale its family puts it on -- an odds
#' ratio, an incidence rate ratio, a time ratio, a hazard ratio -- with an
#' interval built on the link scale and transformed, which is exact and
#' asymmetric rather than symmetric and possibly negative.
#'
#' For a gaussian response there is no ratio to take, so it reports the
#' coefficient, the standardised coefficient, and the partial variance each
#' term explains.
#'
#' @section Conditional against marginal:
#'
#' In a model with random effects, `exp(beta)` is the ratio for a **given
#' cluster**: two patients in the same hospital, one exposed. It is not the
#' ratio for the population, and the gap widens as the random effects grow.
#' Which one a reader wants is almost always the population one, and which one
#' they are usually given is this one. It is labelled, and [ilm_ame()] computes
#' the population-averaged quantity.
#'
#' @section Standardised coefficients:
#'
#' `standardise = TRUE` divides each numeric predictor's coefficient by that
#' predictor's standard deviation, so it reads per standard deviation. Binary
#' and categorical predictors are left alone: dividing by the standard
#' deviation of a 0/1 variable answers what happens per standard deviation of
#' treatment, which is not a question anyone asks. Gelman's suggestion of
#' dividing by **two** standard deviations makes a numeric predictor
#' comparable with a binary one; `standardise = "gelman"` does that.
#'
#' @param object A fitted [ilm_model()].
#' @param level Confidence level.
#' @param standardise `FALSE`, `TRUE` (one standard deviation) or `"gelman"`
#'   (two), for numeric predictors only.
#' @param vcov Optional covariance matrix to use instead of the model's, for
#'   instance from [ilm_vcov_cluster()] or [ilm_svy_vcov()].
#' @return A data frame of `term`, `estimate` (on the link scale), `se`,
#'   `effect` (on the family's scale), `lower`, `upper`, `z` and `p`.
#' @seealso [ilm_ame()] for the population-averaged effect,
#'   [ilm_emmeans()] for group means, [ilm_vcov_cluster()] for a robust
#'   covariance to pass in.
#' @examples
#' set.seed(1); n <- 400
#' d <- data.frame(x = rnorm(n), g = factor(sample(c("a", "b"), n, TRUE)))
#' d$y <- rbinom(n, 1, plogis(-0.4 + 0.8 * d$x + 0.5 * (d$g == "b")))
#' f <- ilm_model(y ~ x + g, data = d, family = "binomial", verbose = FALSE)
#' ilm_effects(f)
#' @export
ilm_effects <- function(object, level = 0.95, standardise = FALSE,
                        vcov = NULL) {
  if (!inherits(object, "ilm_model"))
    stop("`object` must be a fitted ilm_model, not ", class(object)[1],
         call. = FALSE)
  if (object$C > 1L)
    stop("a multinomial fit has one coefficient per category dimension and no ",
         "single ratio per predictor; use summary() or ilm_emmeans() with ",
         "ilm_contrast() for the comparisons you want.", call. = FALSE)
  if (!identical(standardise, FALSE) &&
      !identical(standardise, TRUE) && !identical(standardise, "gelman"))
    stop("`standardise` must be FALSE, TRUE or \"gelman\"", call. = FALSE)

  b <- stats::coef(object)
  V <- if (is.null(vcov)) as.matrix(suppressWarnings(stats::vcov(object)))
       else as.matrix(vcov)
  if (nrow(V) < length(b))
    stop("the covariance supplied is ", nrow(V), " by ", ncol(V),
         " against ", length(b), " coefficients.", call. = FALSE)
  se <- sqrt(pmax(diag(V)[seq_along(b)], 0))

  ## standardise before anything else, so the interval is built from the
  ## scaled coefficient rather than scaled afterwards
  sdv <- rep(1, length(b)); names(sdv) <- names(b)
  skipped <- character(0)
  if (!identical(standardise, FALSE)) {
    mult <- if (identical(standardise, "gelman")) 2 else 1
    mf <- object$model
    asg <- object$assign; tl <- object$term_labels
    for (j in seq_along(b)) {
      nm <- names(b)[j]
      if (nm == "(Intercept)") next
      k <- asg[j]
      v <- if (!is.na(k) && k <= length(tl)) tl[k] else NA_character_
      x <- if (!is.na(v) && v %in% names(mf)) mf[[v]] else NULL
      if (is.null(x) || !is.numeric(x) || length(unique(x)) <= 2L) {
        skipped <- c(skipped, nm); next
      }
      s <- stats::sd(x, na.rm = TRUE)
      if (is.finite(s) && s > 0) sdv[j] <- mult * s
    }
  }
  bs <- b * sdv; ses <- se * sdv

  kind <- ilm_effect_kind(object)
  crit <- if (isTRUE(object$exact_df))
    stats::qt(1 - (1 - level) / 2, object$resid_df) else
    stats::qnorm(1 - (1 - level) / 2)
  lo <- bs - crit * ses; hi <- bs + crit * ses
  tr <- if (kind$exp) exp else identity
  z <- b / se
  out <- data.frame(term = names(b), estimate = unname(bs), se = unname(ses),
                    effect = unname(tr(bs)), lower = unname(tr(lo)),
                    upper = unname(tr(hi)), z = unname(z),
                    p = 2 * stats::pnorm(-abs(unname(z))),
                    row.names = NULL, stringsAsFactors = FALSE)
  if (kind$exp) names(out)[names(out) == "effect"] <- gsub(" ", "_", kind$scale)
  ## without a transform the "effect" column is the estimate again, and two
  ## identical columns invite the reader to look for a difference between them
  else out$effect <- NULL

  ## for a gaussian fit a ratio is meaningless, so report what is
  extra <- NULL
  if (!kind$exp) extra <- ilm_effect_partial(object)

  structure(out, class = c("ilm_effects", "data.frame"), kind = kind,
            level = level, standardise = standardise, skipped = skipped,
            mixed = length(object$re) > 0L, partial = extra,
            robust = !is.null(vcov))
}

#' Partial variance explained by each term
#'
#' The share of the response's variance a term accounts for once the others are
#' in the model, from the drop in residual sum of squares when it is removed.
#' Reported for a gaussian fit, where a ratio would be meaningless.
#'
#' @keywords internal
#' @noRd
ilm_effect_partial <- function(object) {
  if (!identical(if (is.null(object$family)) "gaussian" else
                 object$family$name, "gaussian")) return(NULL)
  if (length(object$re)) return(NULL)
  y <- as.numeric(object$y); X <- object$X
  asg <- object$assign; tl <- object$term_labels
  if (is.null(asg) || !length(tl)) return(NULL)
  w <- if (is.null(object$weights)) rep(1, length(y)) else object$weights
  sw <- sqrt(w)
  rss <- function(cols) {
    if (!length(cols)) return(sum((y * sw - mean(y * sw))^2))
    sum(stats::lsfit(X[, cols, drop = FALSE] * sw, y * sw,
                     intercept = FALSE)$residuals^2)
  }
  full <- rss(seq_len(ncol(X)))
  tot <- rss(which(is.na(asg) | asg == 0L))
  do.call(rbind, lapply(seq_along(tl), function(k) {
    drop_cols <- which(!is.na(asg) & asg == k)
    if (!length(drop_cols)) return(NULL)
    red <- rss(setdiff(seq_len(ncol(X)), drop_cols))
    data.frame(term = tl[k],
               partial_r2 = (red - full) / red,
               semipartial_r2 = (red - full) / max(tot, .Machine$double.eps),
               row.names = NULL, stringsAsFactors = FALSE)
  }))
}

#' @export
print.ilm_effects <- function(x, ...) {
  k <- attr(x, "kind")
  cat(sprintf("Fixed effects as %ss (%.0f%% intervals)\n", k$label,
              100 * attr(x, "level")))
  cat("  ", if (k$exp)
        paste0("exp(coefficient): a change of one unit MULTIPLIES ", k$note)
      else paste0("a change of one unit SHIFTS ", k$note),
      "\n", sep = "")
  if (attr(x, "robust"))
    cat("  using the covariance supplied rather than the model's own\n")
  d <- as.data.frame(x); class(d) <- "data.frame"
  for (j in names(d)[-1]) d[[j]] <- signif(d[[j]], 4)
  print(d, row.names = FALSE)
  if (!identical(attr(x, "standardise"), FALSE)) {
    cat(sprintf("\n  Standardised by %s standard deviation of each NUMERIC predictor.\n",
                if (identical(attr(x, "standardise"), "gelman")) "two" else "one"))
    sk <- attr(x, "skipped")
    if (length(sk))
      cat("  Left alone: ", paste(sk, collapse = ", "),
          ". Dividing by the standard deviation of a 0/1 variable answers\n",
          "  what happens per standard deviation of treatment, which is not a\n",
          "  question anyone asks.\n", sep = "")
  }
  if (attr(x, "mixed") && k$exp)
    cat("\n  CONDITIONAL, not marginal: this is the ", k$scale,
        " for a given\n  cluster -- two units in the same group, one exposed -- and not for the\n",
        "  population. ilm_ame() gives the population-averaged effect.\n",
        sep = "")
  pp <- attr(x, "partial")
  if (!is.null(pp) && nrow(pp)) {
    cat("\n  Variance explained by each term\n")
    q <- pp; q$partial_r2 <- round(q$partial_r2, 4)
    q$semipartial_r2 <- round(q$semipartial_r2, 4)
    print(q, row.names = FALSE)
  }
  invisible(x)
}
