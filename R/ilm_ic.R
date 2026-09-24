## illume: information criteria and fit indices.
##
## SOUNDNESS, stated plainly because the numbers look ordinary and the caveats
## are not:
##
## 1. logLik() is the LAPLACE-APPROXIMATE marginal log-likelihood, not exact.
##    Differences between models are trustworthy when the random structure is
##    held fixed (the approximation error largely cancels) and progressively
##    less so when it is not.
## 2. df counts fixed effects AND covariance parameters -- the marginal AIC
##    convention.  Marginal AIC is biased toward SIMPLER RANDOM structures
##    (Vaida & Blanchard 2005); use it to compare fixed-effect structures, and
##    treat random-structure comparisons with suspicion.  Conditional AIC with
##    effective degrees of freedom would be the right tool there; not implemented.
## 3. BIC's n is genuinely ambiguous in a mixed model.  The effective sample size
##    for covariance parameters is closer to the number of GROUPS than the number
##    of observations.  We default to nobs to stay comparable with lme4/glmmTMB,
##    and offer n = "groups" for the Delattre-style alternative.  They can select
##    different models; that is a real disagreement, not a bug.
## 4. When a variance sits at zero the df count overstates -- the parameter is
##    not freely estimated -- so AIC/BIC are conservative there.

#' Warn before reporting a criterion from an unreliable fit
#' @param object A fitted `"ilm_model"` object.
#' @return `NULL`, invisibly.
#' @keywords internal
#' @noRd
ilm_ic_warn <- function(object) {
  ## at a boundary the optimum is the maximum, so the likelihood is fine
  if (!ilm_fixed_usable(object))
    warning("fit has a non-positive-definite Hessian; AIC/BIC are computed from ",
            "an unreliable optimum. See object$checks.", call. = FALSE)
}

#' Information criteria
#'
#' AIC and BIC computed from the Laplace-approximate marginal log-likelihood,
#' with degrees of freedom counting both fixed effects and covariance
#' parameters.
#'
#' @section Three caveats worth knowing:
#' **The likelihood is approximate.** Differences between models are most
#' trustworthy when the random structure is the same in both, so that the
#' approximation error largely cancels.
#'
#' **Marginal AIC favours simpler random structures.** It is well suited to
#' comparing fixed-effect specifications; comparisons of *random* structures are
#' biased toward the simpler model. Conditional AIC, which uses effective
#' degrees of freedom, is the right tool there and is not implemented.
#'
#' **BIC's sample size is ambiguous in a mixed model.** BIC penalises by
#' `log(n)`, but what counts as `n`? The number of observations, or the number of
#' independent groups? For covariance parameters the effective sample size is
#' closer to the group count. The default `n = "obs"` matches lme4 and glmmTMB so
#' values are comparable with those packages; `n = "groups"` gives the
#' alternative. They can select different models, and that is a genuine
#' disagreement rather than a bug.
#'
#' A variance sitting at zero also makes the degrees of freedom an overcount,
#' since that parameter is not freely estimated, making both criteria
#' conservative.
#'
#' @param object A fitted `"ilm_model"` object.
#' @param k Numeric. Penalty per parameter; 2 gives AIC.
#' @param n `"obs"` or `"groups"`; `BIC()` only.
#' @param ... Unused.
#' @return A numeric scalar.
#' @references
#' Vaida, F., & Blanchard, S. (2005). Conditional Akaike information for
#' mixed-effects models. *Biometrika*, 92(2), 351--370.
#'
#' Delattre, M., Lavielle, M., & Poursat, M.-A. (2014). A note on BIC in
#' mixed-effects models. *Electronic Journal of Statistics*, 8, 456--475.
#' @rdname AIC.ilm_model
#' @export
AIC.ilm_model <- function(object, ..., k = 2) {
  ilm_ic_warn(object)
  ll <- logLik(object)
  -2 * as.numeric(ll) + k * attr(ll, "df")
}

#' Information criteria
#'
#' AIC and BIC computed from the Laplace-approximate marginal log-likelihood,
#' with degrees of freedom counting both fixed effects and covariance
#' parameters.
#'
#' @section Three caveats worth knowing:
#' **The likelihood is approximate.** Differences between models are most
#' trustworthy when the random structure is the same in both, so that the
#' approximation error largely cancels.
#'
#' **Marginal AIC favours simpler random structures.** It is well suited to
#' comparing fixed-effect specifications; comparisons of *random* structures are
#' biased toward the simpler model. Conditional AIC, which uses effective
#' degrees of freedom, is the right tool there and is not implemented.
#'
#' **BIC's sample size is ambiguous in a mixed model.** BIC penalises by
#' `log(n)`, but what counts as `n`? The number of observations, or the number of
#' independent groups? For covariance parameters the effective sample size is
#' closer to the group count. The default `n = "obs"` matches lme4 and glmmTMB so
#' values are comparable with those packages; `n = "groups"` gives the
#' alternative. They can select different models, and that is a genuine
#' disagreement rather than a bug.
#'
#' A variance sitting at zero also makes the degrees of freedom an overcount,
#' since that parameter is not freely estimated, making both criteria
#' conservative.
#'
#' @param object A fitted `"ilm_model"` object.
#' @param k Numeric. Penalty per parameter; 2 gives AIC.
#' @param n `"obs"` or `"groups"`; `BIC()` only.
#' @param ... Unused.
#' @return A numeric scalar.
#' @references
#' Vaida, F., & Blanchard, S. (2005). Conditional Akaike information for
#' mixed-effects models. *Biometrika*, 92(2), 351--370.
#'
#' Delattre, M., Lavielle, M., & Poursat, M.-A. (2014). A note on BIC in
#' mixed-effects models. *Electronic Journal of Statistics*, 8, 456--475.
#' @rdname AIC.ilm_model
#' @export
BIC.ilm_model <- function(object, ..., n = c("obs", "groups")) {
  ilm_ic_warn(object); n <- match.arg(n)
  ll <- logLik(object)
  nn <- if (n == "obs") attr(ll, "nobs") else {
    gl <- object$nlk[vapply(object$re, function(e) e$kind != "basis", TRUE)]
    if (!length(gl)) attr(ll, "nobs") else min(gl)
  }
  -2 * as.numeric(ll) + log(nn) * attr(ll, "df")
}

#' Fitted category probabilities
#'
#' @param object A fitted `"ilm_model"` object.
#' @param conditional Logical. `TRUE` evaluates the random effects at their
#'   fitted values, giving in-sample fitted probabilities. `FALSE` sets grouping
#'   and AR terms to zero. Smooth terms are included either way, because they
#'   are mean structure rather than a population to average over.
#' @return A matrix of probabilities with one column per category.
#' @export
ilm_fitted <- function(object, conditional = TRUE) {
  Tc <- contr.sum(object$J)
  eta <- ilm_eta_hat(object, conditional)
  if (object$C > 1L) {
    P <- exp(eta %*% t(Tc)); P <- P / rowSums(P)
    colnames(P) <- object$ylevels
  } else {
    li <- if (!is.null(object$family)) object$family$linkinv else identity
    P <- matrix(li(eta[, 1]), ncol = 1L, dimnames = list(NULL, "response"))
  }
  P
}

#' Fitted linear predictor
#'
#' The linear predictor behind [ilm_fitted()], before the inverse link. Kept
#' separate because simulating a new response needs the predictor itself: the
#' family simulators are written in terms of `eta`, not of the fitted mean.
#'
#' @param object A fitted `"ilm_model"` object.
#' @param conditional Logical. Include the estimated random effects.
#' @return An N x C matrix.
#' Estimated AR/CAR latent values
#'
#' TMB names the entries of `par.random`, so the block is picked out by name
#' rather than by an offset that would silently shift if the parameter order
#' ever changed.
#'
#' @param object A fitted `"ilm_model"` object.
#' @return A latent-by-category matrix, or `NULL` when there is no such term.
#' @keywords internal
#' @noRd
ilm_Bar_hat <- function(object) {
  if (is.null(object$ar) || is.null(object$sdr$par.random)) return(NULL)
  v <- object$sdr$par.random
  i <- which(names(v) == "B_ar")
  if (length(i) != object$ar$n_cell * object$C) return(NULL)
  matrix(v[i], object$ar$n_cell, object$C)
}

#' @keywords internal
#' @noRd
ilm_eta_hat <- function(object, conditional = TRUE) {
  eta <- object$X %*% object$beta
  ## The AR/CAR latent is a random effect like any other and belongs in the
  ## conditional fitted value. Leaving it out put the whole autoregressive
  ## process into the residual, so every residual diagnostic on such a model
  ## measured the structure the model had already accounted for -- measured,
  ## residual correlation at short separations went UP after fitting CAR(1),
  ## from 0.23 to 0.49, when it should fall toward zero.
  if (conditional && !is.null(object$ar)) {
    Ba <- ilm_Bar_hat(object)
    if (!is.null(Ba)) eta <- eta + Ba[object$ar$idx, , drop = FALSE]
  }
  for (k in seq_along(object$re)) {
    e <- object$re[[k]]
    if (!conditional && e$kind != "basis") next
    B <- ilm_Bhat_term(object, k)
    isrr <- identical(object$re_struct[[k]]$type, "rr")
    if (e$kind == "basis") {
      ctb <- e$basis %*% B
      if (isrr) ctb <- ctb %*% t(object$Lambda[[k]])
      eta <- eta + ctb
    } else {
      nl <- object$nlk[k]
      for (i in seq_len(object$dk[k])) {
        Bi <- B[((i - 1L) * nl + 1L):(i * nl), , drop = FALSE]
        ctb <- Bi[e$group, , drop = FALSE]
        if (isrr) ctb <- ctb %*% t(object$Lambda[[k]])
        eta <- eta + e$Z[, i] * ctb
      }
    }
  }
  eta
}

#' Proper scoring rules for a categorical outcome
#'
#' R-squared has no agreed definition for a nominal outcome, but **proper scoring
#' rules** do. A scoring rule is "proper" if it is optimised by reporting your
#' true beliefs, which makes it a fair way to score probabilistic predictions.
#'
#' The log score is the average of `-log(p)` at the observed category: heavily
#' penalises confident mistakes. The Brier score is the average squared distance
#' between the predicted probability vector and the observed one-hot outcome,
#' ranging from 0 to 2. Lower is better for both.
#'
#' Accuracy is also returned, though it ignores how confident the predictions
#' were and is the least informative of the three.
#'
#' Defined for a multinomial fit, an ordinal one, and a binomial one with a 0/1
#' response, whose single fitted probability is scored as the two categories
#' it implies. Other families have no category probabilities to score; compare
#' those fits with `AIC()` or [ilm_anova()], and check their predictions with
#' [ilm_check_predictive()].
#'
#' An ordinal fit also gets the **ranked probability score**, which the other
#' two do not provide: it compares the cumulative distributions, so a
#' prediction that puts its weight one category away from the truth scores
#' better than one that puts it three away. The log and Brier scores treat
#' every miss alike, as they should for categories with no order.
#'
#' @param object A fitted `"ilm_model"` object.
#' @param conditional Logical, as in [ilm_fitted()].
#' @return A named numeric vector: `log_score`, `brier`, `accuracy`, and for an
#'   ordinal fit `rps`, on the 0 to 1 scale.
#' @references
#' Gneiting, T., & Raftery, A. E. (2007). Strictly proper scoring rules,
#' prediction, and estimation. *Journal of the American Statistical
#' Association*, 102(477), 359--378.
#'
#' Brier, G. W. (1950). Verification of forecasts expressed in terms of
#' probability. *Monthly Weather Review*, 78(1), 1--3.
#' @export
ilm_scores <- function(object, conditional = TRUE) {
  ## A score compares a predicted PROBABILITY for each category with the one
  ## observed, so it needs a categorical outcome. Anything else used to reach
  ## the indexing below with a continuous response and die on "subscript out
  ## of bounds" -- which is how model_performance() failed on every gaussian,
  ## count and survival fit.
  fam <- object$family$name
  y <- object$y
  bin <- identical(fam, "binomial")
  ord <- isTRUE(object$ordinal)
  if (!(identical(fam, "multinomial") || bin || ord) ||
      (bin && !all(y %in% c(0, 1))))
    stop("the scoring rules compare predicted category probabilities with ",
         "the category observed, so they need a categorical outcome: a ",
         "multinomial or ordinal fit, or a binomial one with a 0/1 response. ",
         "For ", ilm_article(fam), " ", fam, " fit, compare models with AIC() ",
         "or ilm_anova(), and check what the model predicts with ",
         "ilm_check_predictive().", call. = FALSE)
  ## An ordinal fit's fitted values are its latent predictor; the category
  ## probabilities are where the thresholds cut it
  P <- if (ord)
    ilm_ord_probs(ilm_eta_hat(object, conditional)[, 1L], object$zeta,
                  object$family$pfun)
  else ilm_fitted(object, conditional)
  ## a binary fit returns the probability of the second level only
  if (bin) { P <- cbind(1 - P[, 1], P[, 1]); y <- y + 1 }
  N <- nrow(P)
  py <- P[cbind(seq_len(N), y)]
  Y <- matrix(0, N, ncol(P)); Y[cbind(seq_len(N), y)] <- 1
  out <- c(log_score = mean(-log(pmax(py, .Machine$double.eps))),
           brier = mean(rowSums((P - Y)^2)),
           accuracy = mean(max.col(P, ties.method = "first") == y))
  if (ord) {
    ## squared distance between the cumulative distributions, over the J - 1
    ## cut points, so a near miss costs less than a far one
    J <- ncol(P)
    Fc <- t(apply(P, 1L, cumsum)); Oc <- t(apply(Y, 1L, cumsum))
    out["rps"] <- mean(rowSums((Fc - Oc)[, -J, drop = FALSE]^2) / (J - 1L))
  }
  out
}

## "a" or "an", for a word in a message
#' @keywords internal
#' @noRd
ilm_article <- function(w) if (grepl("^[aeiouAEIOU]", w)) "an" else "a"

#' Log-likelihood of an intercept-only comparison model
#'
#' Refits with the predictors removed but the grouping random effects kept, so
#' that McFadden's R-squared measures what the predictors add *over and above* a
#' random-intercept model.
#'
#' @param object A fitted `"ilm_model"` object.
#' @param restarts Integer. Optimiser restarts.
#' @return Numeric log-likelihood, or `NA` if the null model did not converge.
#' @keywords internal
#' @noRd
ilm_null_ll <- function(object, restarts = 1L) {
  ## A restricted likelihood belongs to contrasts orthogonal to the fixed
  ## design, so the intercept-only model's is the likelihood of different data
  ## and a ratio of the two means nothing.
  if (isTRUE(object$reml)) return(NA_real_)
  ## The intercept, and a flexible survival baseline's spline columns, which
  ## ARE the baseline rather than predictors -- without them the null model
  ## has no hazard at all.
  keep <- (!is.na(object$assign) & object$assign == 0) |
    seq_along(object$assign) %in% object$rp$cols
  if (!any(keep)) return(NA_real_)
  gi <- which(vapply(object$re, function(e) e$kind != "basis", TRUE))
  if (!length(gi)) return(NA_real_)
  ## Otherwise the same model, every part of it: a null model without the
  ## censoring, the dispersion model or the zero part is a different
  ## likelihood, and an R-squared built from the two compares unlike things.
  stub <- ilm_refit_stub(object)
  stub$re <- object$re[gi]; stub$re_struct <- object$re_struct[gi]
  f <- try(ilm_refit_like(stub, X = object$X[, keep, drop = FALSE],
                          keep = keep, restarts = restarts), silent = TRUE)
  if (inherits(f, "try-error") || f$opt$convergence != 0) return(NA_real_)
  -f$opt$objective
}

#' Fit indices for a multinomial mixed model
#'
#' A `performance::model_performance()` method. Returns AIC, AICc, BIC (both
#' sample-size conventions), the log-likelihood, degrees of freedom, McFadden's
#' R-squared, and the scoring rules from [ilm_scores()].
#'
#' @section What is deliberately absent:
#' Nakagawa's marginal and conditional R-squared, and the ICC, are **not**
#' reported. Both require a distribution-specific residual variance, and a
#' nominal multinomial outcome does not have one. Reporting them would mean
#' inventing a quantity rather than estimating it.
#'
#' McFadden's R-squared is relative to an intercept-only model with the same
#' grouping random effects, so it answers "how much do the predictors add?"
#' rather than "how much does the whole model explain?". Unlike R-squared in
#' linear regression, values around 0.2 to 0.4 already indicate a good fit.
#'
#' The scoring rules are computed **in sample** and conditional on the fitted
#' random effects, so both are optimistic. Use cross-validation or a held-out
#' set if you need an honest predictive comparison.
#'
#' @param model A fitted `"ilm_model"` object.
#' @param metrics Unused; present for compatibility with the generic.
#' @param verbose Logical. Message when the fit failed its checks.
#' @param ... Unused.
#' @return A one-row data frame of class `"performance_model"`.
#' @references
#' McFadden, D. (1974). Conditional logit analysis of qualitative choice
#' behavior. In P. Zarembka (Ed.), *Frontiers in Econometrics*. Academic Press.
#'
#' Nakagawa, S., & Schielzeth, H. (2013). A general and simple method for
#' obtaining R-squared from generalized linear mixed-effects models. *Methods in
#' Ecology and Evolution*, 4(2), 133--142. (The approach not used here, and why:
#' it needs a distribution-specific variance.)
#' @exportS3Method performance::model_performance
model_performance.ilm_model <- function(model, metrics = "all", ..., verbose = TRUE) {
  ll <- logLik(model); df <- attr(ll, "df"); n <- attr(ll, "nobs")
  ## the scoring rules exist for a categorical outcome only; elsewhere they
  ## are NA, with a note, rather than an error from deep inside
  sc <- tryCatch(ilm_scores(model, conditional = TRUE),
                 error = function(e) c(log_score = NA_real_, brier = NA_real_,
                                       accuracy = NA_real_))
  ll0 <- ilm_null_ll(model)
  out <- data.frame(
    AIC = AIC(model),
    AICc = -2 * as.numeric(ll) + 2 * df * n / max(1, n - df - 1),
    BIC = BIC(model),
    BIC_groups = BIC(model, n = "groups"),
    logLik = as.numeric(ll),
    df = df,
    R2_McFadden = if (is.na(ll0)) NA_real_ else 1 - as.numeric(ll) / ll0,
    Log_score = unname(sc["log_score"]),
    Brier = unname(sc["brier"]),
    Accuracy = unname(sc["accuracy"]))
  if (verbose && !ilm_fixed_usable(model))
    message("Fit did not pass all checks; these indices describe an unreliable optimum.")
  attr(out, "r2_note") <- paste0(
    "R2_McFadden is relative to an intercept-only model with the same grouping random effects. Nakagawa R2 and ICC are omitted: a nominal multinomial has no distribution-specific residual variance to define them.",
    if (isTRUE(model$reml)) " It is NA for this REML fit: restricted likelihoods of models with different fixed effects are not comparable. Refit with reml = FALSE for it." else "")
  attr(out, "score_note") <- if (all(is.na(sc)))
    paste0("Log_score, Brier and Accuracy score predicted category ",
           "probabilities, so they are NA for this ", model$family$name,
           " fit.")
  else
    "Log_score and Brier are IN-SAMPLE and conditional on the fitted random effects, so both are optimistic."
  class(out) <- c("performance_model", "data.frame")
  out
}
