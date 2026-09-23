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
#' @param object A fitted `"ilm_model"` object.
#' @return Integer.
#' @keywords internal
#' @noRd
ilm_n_fixed <- function(object) ncol(object$X) * object$C

#' Warn when the covariance matrix is unusable
#'
#' If the Hessian is not positive definite, the optimiser has not found a proper
#' minimum and the standard errors derived from it mean nothing -- even though
#' the coefficients themselves may look entirely reasonable. Silently returning
#' those numbers is the single most dangerous failure mode of this model class,
#' so anything that hands out standard errors warns first.
#'
#' @param object A fitted `"ilm_model"` object.
#' @return `NULL`, invisibly; called for the warning.
#' @keywords internal
#' @noRd
ilm_warn_hess <- function(object, full = FALSE) {
  if (isTRUE(object$sdr$pdHess)) return(invisible())
  ## A term held at its boundary leaves the fixed effects usable, so asking for
  ## them is not an error worth a warning; asking for the covariance
  ## parameters too is, because the ones along a held direction carry no
  ## uncertainty.
  if (length(object$hessian_held)) {
    if (full)
      warning("the covariance of ", paste(object$hessian_held, collapse = ", "),
              " sits at a boundary, and the direction in which it cannot be ",
              "resolved is held at its estimate, so the rows of the ",
              "parameters along it carry no uncertainty here. The fixed ",
              "effects are unaffected; see fit$checks.", call. = FALSE)
    return(invisible())
  }
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
#' @param object A fitted `"ilm_model"` object.
#' @param full Logical. Include covariance parameters as well as fixed effects.
#' @param ... Unused.
#' @return A named numeric vector (`coef()`) or a named matrix (`vcov()`).
#' @seealso [ilm_coef_table()] for a formatted table with tests.
#' @rdname coef.ilm_model
#' @export
coef.ilm_model <- function(object, full = FALSE, ...) {
  v <- setNames(object$opt$par, object$pnames)
  if (isTRUE(full)) v else v[seq_len(ilm_n_fixed(object))]
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
#' @param object A fitted `"ilm_model"` object.
#' @param full Logical. Include covariance parameters as well as fixed effects.
#' @param ... Unused.
#' @return A named numeric vector (`coef()`) or a named matrix (`vcov()`).
#' @seealso [ilm_coef_table()] for a formatted table with tests.
#' @rdname coef.ilm_model
#' @export
vcov.ilm_model <- function(object, full = FALSE, ...) {
  ilm_warn_hess(object, full = isTRUE(full))
  V <- object$sdr$cov.fixed
  if (is.null(V)) return(NULL)
  ## A term held at its boundary is known, for everything computed from this
  ## matrix: zero uncertainty, which is what holding it means, rather than the
  ## NA that records it in the fit and would poison a delta method downstream.
  if (length(object$hessian_held)) V[!is.finite(V)] <- 0
  ## Maximum likelihood divides the residual sum of squares by n; an ordinary
  ## linear model divides by n - p.  Rescaling by n / (n - p) makes the
  ## covariance, and therefore the standard errors, agree exactly with lm().
  ## REML has already done precisely that -- restricting the likelihood to
  ## contrasts orthogonal to X is what produces the n - p divisor -- so
  ## applying the correction again would inflate every standard error twice.
  if (isTRUE(object$exact_df) && !isTRUE(object$reml)) {
    n <- nrow(object$X); k <- ilm_n_fixed(object)
    V[seq_len(k), seq_len(k)] <- V[seq_len(k), seq_len(k)] * n / object$resid_df
  }
  dimnames(V) <- list(object$pnames, object$pnames)
  if (isTRUE(full)) V else {
    k <- seq_len(ilm_n_fixed(object)); V[k, k, drop = FALSE]
  }
}

#' Fixed effects as a matrix
#'
#' The same numbers as [coef.ilm_model()] but shaped as predictors by categories,
#' which is usually easier to read for a multinomial model. `ilm_se_fixef()` returns
#' matching standard errors.
#'
#' @param object A fitted `"ilm_model"` object.
#' @param ... Unused.
#' @return A numeric matrix with one row per design column and one column per
#'   category dimension.
#' @rdname fixef.ilm_model
#' @export
fixef.ilm_model <- function(object, ...) {
  B <- object$beta
  dimnames(B) <- list(colnames(object$X), object$ylevels[seq_len(object$C)])
  B
}

#' Fixed effects as a matrix
#'
#' The same numbers as [coef.ilm_model()] but shaped as predictors by categories,
#' which is usually easier to read for a multinomial model. `ilm_se_fixef()` returns
#' matching standard errors.
#'
#' @param object A fitted `"ilm_model"` object.
#' @param ... Unused.
#' @return A numeric matrix with one row per design column and one column per
#'   category dimension.
#' @rdname fixef.ilm_model
#' @export
ilm_se_fixef <- function(object) {
  ilm_warn_hess(object)
  V <- object$sdr$cov.fixed
  if (is.null(V)) return(NULL)
  k <- seq_len(ilm_n_fixed(object)); d <- diag(V)[k]
  d[d < 0] <- NA_real_
  matrix(sqrt(d), ncol(object$X), object$C,
         dimnames = list(colnames(object$X), object$ylevels[seq_len(object$C)]))
}

#' Log-likelihood, and number of observations
#'
#' Returns the **Laplace-approximate** marginal log-likelihood -- the likelihood
#' with the random effects integrated out, using the approximation described in
#' [ilm_fit()]. The `df` attribute counts fixed effects *and* covariance
#' parameters, which is the convention [AIC.ilm_model()] and [BIC.ilm_model()] use.
#'
#' Because it is an approximation, differences between models are most
#' trustworthy when the random structure is held fixed, so that the
#' approximation error largely cancels.
#'
#' The value is on the same scale as `lme4::lmer()` and `nlme::lme()`, and
#' agrees with both to numerical tolerance for a gaussian model they can also
#' fit. It therefore includes every normalising constant, so `AIC()` and `BIC()`
#' may be compared across models with *different* random structures as well as
#' the same one.
#'
#' @param object A fitted `"ilm_model"` object.
#' @param ... Unused.
#' @return An object of class `"logLik"`.
#' @seealso [AIC.ilm_model()], [BIC.ilm_model()], [ilm_anova()].
#' @rdname logLik.ilm_model
#' @export
logLik.ilm_model <- function(object, ...) {
  ## df counts fixed effects AND covariance parameters: both are estimated, and
  ## anova()-style comparisons need the total.
  ##
  ## The correction is not cosmetic. The random-effect priors in the likelihood
  ## are written as 0.5 * u' Sigma^-1 u + 0.5 * log|Sigma|, which leaves out the
  ## (1/2) log(2*pi) each latent scalar contributes to the gaussian density.
  ## TMB's Laplace step then subtracts (q/2) log(2*pi) of its own, so the
  ## objective it reports is the true negative log-likelihood MINUS that
  ## constant -- and the constant grows with the number of latent values.
  ##
  ## Left uncorrected, every extra latent value bought a model about 0.92
  ## log-likelihood units for free, so AIC and BIC preferred whichever model
  ## had the larger random structure. Measured on a random-intercept fit with
  ## 50 groups, illume reported -595.083 where lme4 and nlme both reported
  ## -641.030: a gap of 45.947, against 25 * log(2*pi) = 45.947 predicted.
  ##
  ## Comparisons that hold the random structure fixed -- which is what
  ## ilm_anova(test = "LRT") does -- are unaffected, because the constant is
  ## then the same in both models and cancels. That is why this survived the
  ## coverage and power studies.
  q <- object$n_integrated
  if (is.null(q)) q <- 0L
  ## Under boundary = "avoid" the optimiser minimised the penalised objective;
  ## the likelihood of the data at that estimate is the objective with the
  ## penalty taken back out.
  pen <- if (is.null(object$re_penalty)) 0 else object$re_penalty
  val <- -object$opt$objective - pen - (q / 2) * log(2 * pi)
  structure(val, df = length(object$opt$par), nobs = nrow(object$X), class = "logLik")
}

#' Log-likelihood, and number of observations
#'
#' Returns the **Laplace-approximate** marginal log-likelihood -- the likelihood
#' with the random effects integrated out, using the approximation described in
#' [ilm_fit()]. The `df` attribute counts fixed effects *and* covariance
#' parameters, which is the convention [AIC.ilm_model()] and [BIC.ilm_model()] use.
#'
#' Because it is an approximation, differences between models are most
#' trustworthy when the random structure is held fixed, so that the
#' approximation error largely cancels.
#'
#' @param object A fitted `"ilm_model"` object.
#' @param ... Unused.
#' @return An object of class `"logLik"`.
#' @seealso [AIC.ilm_model()], [BIC.ilm_model()], [ilm_anova()].
#' @rdname logLik.ilm_model
#' @export
nobs.ilm_model <- function(object, ...) nrow(object$X)

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
#' `ilm_anova(test = "LRT")` or, for small samples, [ilm_pb_lrt()].
#'
#' @param object A fitted `"ilm_model"` object.
#' @return A data frame with columns `Estimate`, `Std. Error`, `z value` and
#'   `Pr(>|z|)`.
#' @references
#' Hauck, W. W., & Donner, A. (1977). Wald's test as applied to hypotheses in
#' logit analysis. *Journal of the American Statistical Association*, 72(360),
#' 851--853.
#' @export
ilm_coef_table <- function(object) {
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
