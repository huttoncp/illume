## illume: extractor methods.
##
## The contract every downstream package depends on: names(coef()) matches
## colnames(vcov()) exactly, in a fixed documented order.  car::Anova,
## marginaleffects, emmeans and ggeffects all key off it.
##
## Ordering is obj$par order: the p*C fixed effects first (column-major over
## the sum-to-zero category contrasts), then the covariance parameters, then the
## AR block if present.

## `fixef` lives in nlme/lme4; define the generic only if nothing else has.
if (!exists("fixef")) fixef <- function(object, ...) UseMethod("fixef")

#' Number of fixed-effect coefficients
#'
#' `ncol(X)` times `C`: one coefficient per design column per category dimension.
#'
#' @param object A fitted `"lum_model"` object.
#' @return Integer.
#' @keywords internal
#' @noRd
lum_n_fixed <- function(object) ncol(object$X) * object$C

#' Warn when the covariance matrix is unusable
#'
#' If the Hessian is not positive definite, the optimiser has not found a proper
#' minimum and the standard errors derived from it mean nothing -- even though
#' the coefficients themselves may look entirely reasonable. Silently returning
#' those numbers is the single most dangerous failure mode of this model class,
#' so anything that hands out standard errors warns first.
#'
#' @param object A fitted `"lum_model"` object.
#' @return `NULL`, invisibly; called for the warning.
#' @keywords internal
#' @noRd
lum_warn_hess <- function(object) {
  if (!isTRUE(object$sdr$pdHess))
    warning("Hessian is not positive definite: this covariance matrix is unusable ",
            "for inference. See fit$checks.", call. = FALSE)
}

#' Coefficients and their covariance matrix
#'
#' `coef()` returns the fixed effects as a named vector and `vcov()` the matching
#' covariance matrix. Their names line up exactly, which is what `car::Anova()`,
#' `emmeans` and `marginaleffects` rely on.
#'
#' @details
#' Names take the form `category:predictor`, for example `low:x1`. Because
#' categories are coded sum-to-zero, `low:x1` is the effect of `x1` on the `low`
#' category *relative to the average across categories*, not relative to a
#' baseline category.
#'
#' `full = TRUE` appends the covariance parameters (entries of the Cholesky
#' factors) to the fixed effects. This matters for population-averaged
#' predictions, which depend on the random-effect covariance as well as the
#' fixed effects: a delta-method calculation given only the fixed effects would
#' silently report no uncertainty at all from the covariance parameters.
#'
#' @param object A fitted `"lum_model"` object.
#' @param full Logical. Include covariance parameters as well as fixed effects.
#' @param ... Unused.
#' @return A named numeric vector (`coef()`) or a named matrix (`vcov()`).
#' @seealso [lum_coef_table()] for a formatted table with tests.
#' @rdname coef.lum_model
#' @export
coef.lum_model <- function(object, full = FALSE, ...) {
  v <- setNames(object$opt$par, object$pnames)
  if (isTRUE(full)) v else v[seq_len(lum_n_fixed(object))]
}

#' Coefficients and their covariance matrix
#'
#' `coef()` returns the fixed effects as a named vector and `vcov()` the matching
#' covariance matrix. Their names line up exactly, which is what `car::Anova()`,
#' `emmeans` and `marginaleffects` rely on.
#'
#' @details
#' Names take the form `category:predictor`, for example `low:x1`. Because
#' categories are coded sum-to-zero, `low:x1` is the effect of `x1` on the `low`
#' category *relative to the average across categories*, not relative to a
#' baseline category.
#'
#' `full = TRUE` appends the covariance parameters (entries of the Cholesky
#' factors) to the fixed effects. This matters for population-averaged
#' predictions, which depend on the random-effect covariance as well as the
#' fixed effects: a delta-method calculation given only the fixed effects would
#' silently report no uncertainty at all from the covariance parameters.
#'
#' @param object A fitted `"lum_model"` object.
#' @param full Logical. Include covariance parameters as well as fixed effects.
#' @param ... Unused.
#' @return A named numeric vector (`coef()`) or a named matrix (`vcov()`).
#' @seealso [lum_coef_table()] for a formatted table with tests.
#' @rdname coef.lum_model
#' @export
vcov.lum_model <- function(object, full = FALSE, ...) {
  lum_warn_hess(object)
  V <- object$sdr$cov.fixed
  if (is.null(V)) return(NULL)
  ## Maximum likelihood divides the residual sum of squares by n; an ordinary
  ## linear model divides by n - p.  Rescaling by n / (n - p) makes the
  ## covariance, and therefore the standard errors, agree exactly with lm().
  if (isTRUE(object$exact_df)) {
    n <- nrow(object$X); k <- lum_n_fixed(object)
    V[seq_len(k), seq_len(k)] <- V[seq_len(k), seq_len(k)] * n / object$resid_df
  }
  dimnames(V) <- list(object$pnames, object$pnames)
  if (isTRUE(full)) V else {
    k <- seq_len(lum_n_fixed(object)); V[k, k, drop = FALSE]
  }
}

#' Fixed effects as a matrix
#'
#' The same numbers as [coef.lum_model()] but shaped as predictors by categories,
#' which is usually easier to read for a multinomial model. `lum_se_fixef()` returns
#' matching standard errors.
#'
#' @param object A fitted `"lum_model"` object.
#' @param ... Unused.
#' @return A numeric matrix with one row per design column and one column per
#'   category dimension.
#' @rdname fixef.lum_model
#' @export
fixef.lum_model <- function(object, ...) {
  B <- object$beta
  dimnames(B) <- list(colnames(object$X), object$ylevels[seq_len(object$C)])
  B
}

#' Fixed effects as a matrix
#'
#' The same numbers as [coef.lum_model()] but shaped as predictors by categories,
#' which is usually easier to read for a multinomial model. `lum_se_fixef()` returns
#' matching standard errors.
#'
#' @param object A fitted `"lum_model"` object.
#' @param ... Unused.
#' @return A numeric matrix with one row per design column and one column per
#'   category dimension.
#' @rdname fixef.lum_model
#' @export
lum_se_fixef <- function(object) {
  lum_warn_hess(object)
  V <- object$sdr$cov.fixed
  if (is.null(V)) return(NULL)
  k <- seq_len(lum_n_fixed(object)); d <- diag(V)[k]
  d[d < 0] <- NA_real_
  matrix(sqrt(d), ncol(object$X), object$C,
         dimnames = list(colnames(object$X), object$ylevels[seq_len(object$C)]))
}

#' Log-likelihood, and number of observations
#'
#' Returns the **Laplace-approximate** marginal log-likelihood -- the likelihood
#' with the random effects integrated out, using the approximation described in
#' [lum_fit()]. The `df` attribute counts fixed effects *and* covariance
#' parameters, which is the convention [AIC.lum_model()] and [BIC.lum_model()] use.
#'
#' Because it is an approximation, differences between models are most
#' trustworthy when the random structure is held fixed, so that the
#' approximation error largely cancels.
#'
#' @param object A fitted `"lum_model"` object.
#' @param ... Unused.
#' @return An object of class `"logLik"`.
#' @seealso [AIC.lum_model()], [BIC.lum_model()], [lum_anova()].
#' @rdname logLik.lum_model
#' @export
logLik.lum_model <- function(object, ...) {
  ## df counts fixed effects AND covariance parameters: both are estimated, and
  ## anova()-style comparisons need the total.
  val <- -object$opt$objective
  structure(val, df = length(object$opt$par), nobs = nrow(object$X), class = "logLik")
}

#' Log-likelihood, and number of observations
#'
#' Returns the **Laplace-approximate** marginal log-likelihood -- the likelihood
#' with the random effects integrated out, using the approximation described in
#' [lum_fit()]. The `df` attribute counts fixed effects *and* covariance
#' parameters, which is the convention [AIC.lum_model()] and [BIC.lum_model()] use.
#'
#' Because it is an approximation, differences between models are most
#' trustworthy when the random structure is held fixed, so that the
#' approximation error largely cancels.
#'
#' @param object A fitted `"lum_model"` object.
#' @param ... Unused.
#' @return An object of class `"logLik"`.
#' @seealso [AIC.lum_model()], [BIC.lum_model()], [lum_anova()].
#' @rdname logLik.lum_model
#' @export
nobs.lum_model <- function(object, ...) nrow(object$X)

#' Coefficient table with Wald tests
#'
#' Estimates, standard errors, z statistics and two-sided p-values, in the layout
#' `summary()` prints.
#'
#' @details
#' These are Wald tests, which are fast but have a known weakness: for very
#' strong effects, or when a category is nearly perfectly predicted, the Wald
#' statistic can *shrink* rather than grow. This is the Hauck-Donner effect. If a
#' term matters to your conclusions, confirm it with a likelihood-ratio test via
#' `lum_anova(test = "LRT")` or, for small samples, [lum_pb_lrt()].
#'
#' @param object A fitted `"lum_model"` object.
#' @return A data frame with columns `Estimate`, `Std. Error`, `z value` and
#'   `Pr(>|z|)`.
#' @references
#' Hauck, W. W., & Donner, A. (1977). Wald's test as applied to hypotheses in
#' logit analysis. *Journal of the American Statistical Association*, 72(360),
#' 851--853.
#' @export
lum_coef_table <- function(object) {
  b <- coef(object); s <- suppressWarnings(sqrt(diag(vcov(object))))
  t <- b / s
  ## With nothing integrated out the reference distribution is exactly t on
  ## n - p degrees of freedom; otherwise only the large-sample normal is
  ## available, because a multinomial or mixed model has no exact analogue.
  if (isTRUE(object$exact_df)) {
    data.frame(Estimate = b, `Std. Error` = s, `t value` = t,
               `Pr(>|t|)` = 2 * stats::pt(-abs(t), object$resid_df),
               check.names = FALSE)
  } else {
    data.frame(Estimate = b, `Std. Error` = s, `z value` = t,
               `Pr(>|z|)` = 2 * pnorm(-abs(t)), check.names = FALSE)
  }
}
