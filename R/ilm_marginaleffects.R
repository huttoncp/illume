## illume: marginaleffects glue.
##
## Four methods.  marginaleffects computes the delta method NUMERICALLY: it
## perturbs whatever get_coef() returns, calls set_coef(), and re-predicts.  Two
## consequences drive the design here.
##
## 1. get_coef() returns (beta, theta) TOGETHER, not just beta.  A
##    population-averaged marginal mean depends on the covariance parameters --
##    it integrates the softmax over the random-effect distribution -- so a
##    beta-only vector would silently report ZERO uncertainty from theta.
##    Printouts will therefore list variance parameters among the "coefficients",
##    which is cosmetically odd and statistically right.
##
## 2. set_coef() must rebuild every DERIVED quantity (Sigma, Lambda, Sigma_d,
##    rho), not just store the new vector.  If it only overwrote $opt$par, then
##    perturbing theta would change nothing downstream and the numerical
##    derivative with respect to the covariance parameters would come back as
##    exactly zero -- a silent, plausible-looking wrong answer.
##
## Determinism: predict() draws its random effects under a fixed seed, so the
## same standard normal draws are reused across perturbations and only their
## transformation by Sigma changes.  That keeps the numerical derivative smooth
## instead of drowning it in Monte Carlo noise.

#' Rebuild a fitted model from a modified parameter vector
#'
#' Replaces the parameters and recomputes **everything derived from them** --
#' covariance matrices, loadings and the AR correlation.
#'
#' Recomputing is the whole point. `marginaleffects` works by nudging parameters
#' and seeing how predictions change. If this function merely stored the new
#' numbers without rebuilding the covariances, nudging a covariance parameter
#' would change nothing downstream, and the reported uncertainty from those
#' parameters would come back as exactly zero -- a wrong answer that looks
#' perfectly reasonable.
#'
#' @param object A fitted `"ilm_model"` object.
#' @param par Numeric parameter vector in `coef(object, full = TRUE)` order.
#' @return The object with all derived quantities recomputed.
#' @keywords internal
#' @noRd
ilm_rebuild <- function(object, par) {
  pn <- object$pnames; C <- object$C; p <- ncol(object$X)
  object$opt$par <- setNames(as.numeric(par), pn)
  nb <- p * C
  object$beta <- matrix(object$opt$par[seq_len(nb)], p, C)
  thv <- object$opt$par[nb + seq_len(sum(object$npc + object$npd))]
  K <- length(object$re)
  Sig <- vector("list", K); Lam <- vector("list", K); Sd <- list()
  for (k in seq_len(K)) {
    blk <- thv[(object$toff[k] + 1L):object$toff[k + 1L]]
    th <- blk[seq_len(object$npc[k])]
    if (object$ty[k] == "rr") {
      L <- ilm_mkLam_num(th, C, object$rk[k]); Lam[[k]] <- L; Sig[[k]] <- L %*% t(L)
    } else {
      L <- if (object$ty[k] == "diag") ilm_mkD_num(th, C) else ilm_mkL_num(th, C)
      Sig[[k]] <- L %*% t(L)
    }
    if (object$dk[k] > 1L) {
      thd <- blk[object$npc[k] + seq_len(object$npd[k])]
      Ld <- ilm_mkLd_num(thd, object$dk[k], object$dcor[k])
      Sd[[names(object$re)[k]]] <- Ld %*% t(Ld)
    }
  }
  names(Sig) <- names(Lam) <- names(object$re)
  if (!is.null(object$ar)) {
    La <- ilm_mkL_num(object$opt$par[grepl("^ar:L\\[", pn)], C)   # pnames is in fill order
    Sig[["ar"]] <- La %*% t(La)
    object$rho <- tanh(object$opt$par[pn == "ar:rho_raw"])
  }
  object$Sigma <- Sig; object$Lambda <- Lam; object$Sigma_d <- Sd
  object
}

#' marginaleffects interface
#'
#' The four methods `marginaleffects` needs. It computes derivatives
#' **numerically**: it asks for the parameters, nudges them, sets them back, and
#' re-predicts.
#'
#' `get_coef()` returns fixed effects **and** covariance parameters together,
#' because a population-averaged prediction depends on both. Returning only the
#' fixed effects would silently report no uncertainty from the covariance
#' parameters. Variance parameters therefore appear among the "coefficients" in
#' printed output, which looks odd but is correct.
#'
#' @param model A fitted `"ilm_model"` object.
#' @param coefs Numeric replacement parameter vector.
#' @param vcov Passed through by `marginaleffects`.
#' @param newdata Optional data frame.
#' @param type Prediction type.
#' @param marginal Logical. Population-averaged predictions; defaults to the
#'   `ilm_model.marginal` option, which is `TRUE`.
#' @param ndraw Integer. Random-effect draws when `marginal = TRUE`.
#' @param ... Unused.
#' @return Parameters, a covariance matrix, a modified model, or a long-format
#'   data frame of predictions with `rowid`, `group` and `estimate`.
#' @references
#' Arel-Bundock, V., Greifer, N., & Heiss, A. (2024). How to interpret
#' statistical models using marginaleffects for R and Python. *Journal of
#' Statistical Software*, 111(9), 1--32.
#' @rdname marginaleffects-methods
#' @export
get_coef.ilm_model <- function(model, ...) coef(model, full = TRUE)

#' marginaleffects interface
#'
#' The four methods `marginaleffects` needs. It computes derivatives
#' **numerically**: it asks for the parameters, nudges them, sets them back, and
#' re-predicts.
#'
#' `get_coef()` returns fixed effects **and** covariance parameters together,
#' because a population-averaged prediction depends on both. Returning only the
#' fixed effects would silently report no uncertainty from the covariance
#' parameters. Variance parameters therefore appear among the "coefficients" in
#' printed output, which looks odd but is correct.
#'
#' @param model A fitted `"ilm_model"` object.
#' @param coefs Numeric replacement parameter vector.
#' @param vcov Passed through by `marginaleffects`.
#' @param newdata Optional data frame.
#' @param type Prediction type.
#' @param marginal Logical. Population-averaged predictions; defaults to the
#'   `ilm_model.marginal` option, which is `TRUE`.
#' @param ndraw Integer. Random-effect draws when `marginal = TRUE`.
#' @param ... Unused.
#' @return Parameters, a covariance matrix, a modified model, or a long-format
#'   data frame of predictions with `rowid`, `group` and `estimate`.
#' @references
#' Arel-Bundock, V., Greifer, N., & Heiss, A. (2024). How to interpret
#' statistical models using marginaleffects for R and Python. *Journal of
#' Statistical Software*, 111(9), 1--32.
#' @rdname marginaleffects-methods
#' @export
set_coef.ilm_model <- function(model, coefs, ...) ilm_rebuild(model, coefs)

#' marginaleffects interface
#'
#' The four methods `marginaleffects` needs. It computes derivatives
#' **numerically**: it asks for the parameters, nudges them, sets them back, and
#' re-predicts.
#'
#' `get_coef()` returns fixed effects **and** covariance parameters together,
#' because a population-averaged prediction depends on both. Returning only the
#' fixed effects would silently report no uncertainty from the covariance
#' parameters. Variance parameters therefore appear among the "coefficients" in
#' printed output, which looks odd but is correct.
#'
#' @param model A fitted `"ilm_model"` object.
#' @param coefs Numeric replacement parameter vector.
#' @param vcov Passed through by `marginaleffects`.
#' @param newdata Optional data frame.
#' @param type Prediction type.
#' @param marginal Logical. Population-averaged predictions; defaults to the
#'   `ilm_model.marginal` option, which is `TRUE`.
#' @param ndraw Integer. Random-effect draws when `marginal = TRUE`.
#' @param ... Unused.
#' @return Parameters, a covariance matrix, a modified model, or a long-format
#'   data frame of predictions with `rowid`, `group` and `estimate`.
#' @references
#' Arel-Bundock, V., Greifer, N., & Heiss, A. (2024). How to interpret
#' statistical models using marginaleffects for R and Python. *Journal of
#' Statistical Software*, 111(9), 1--32.
#' @rdname marginaleffects-methods
#' @export
get_vcov.ilm_model <- function(model, vcov = NULL, ...) {
  if (isFALSE(vcov)) return(NULL)
  suppressWarnings(vcov(model, full = TRUE))
}

#' marginaleffects interface
#'
#' The four methods `marginaleffects` needs. It computes derivatives
#' **numerically**: it asks for the parameters, nudges them, sets them back, and
#' re-predicts.
#'
#' `get_coef()` returns fixed effects **and** covariance parameters together,
#' because a population-averaged prediction depends on both. Returning only the
#' fixed effects would silently report no uncertainty from the covariance
#' parameters. Variance parameters therefore appear among the "coefficients" in
#' printed output, which looks odd but is correct.
#'
#' @param model A fitted `"ilm_model"` object.
#' @param coefs Numeric replacement parameter vector.
#' @param vcov Passed through by `marginaleffects`.
#' @param newdata Optional data frame.
#' @param type Prediction type.
#' @param marginal Logical. Population-averaged predictions; defaults to the
#'   `ilm_model.marginal` option, which is `TRUE`.
#' @param ndraw Integer. Random-effect draws when `marginal = TRUE`.
#' @param ... Unused.
#' @return Parameters, a covariance matrix, a modified model, or a long-format
#'   data frame of predictions with `rowid`, `group` and `estimate`.
#' @references
#' Arel-Bundock, V., Greifer, N., & Heiss, A. (2024). How to interpret
#' statistical models using marginaleffects for R and Python. *Journal of
#' Statistical Software*, 111(9), 1--32.
#' @rdname marginaleffects-methods
#' @export
get_predict.ilm_model <- function(model, newdata = NULL, type = "response",
                             marginal = NULL, ndraw = 100L, ...) {
  if (is.null(marginal))
    marginal <- isTRUE(getOption("ilm_model.marginal", TRUE))   # PA is the default
  if (is.null(newdata)) newdata <- model$model
  P <- predict(model, newdata = newdata, type = "response",
               marginal = marginal, ndraw = ndraw)
  n <- nrow(P)
  data.frame(rowid = rep(seq_len(n), times = ncol(P)),
             group = rep(colnames(P), each = n),
             estimate = as.vector(P),
             stringsAsFactors = FALSE)
}

#' Register the marginaleffects interface
#'
#' Call once per session before using `marginaleffects` with an `"ilm_model"` fit.
#'
#' Registering the S3 methods alone is not enough. `marginaleffects` checks model
#' classes against a fixed list and rejects anything unfamiliar *before* dispatch
#' happens, so the methods would never be reached. This also adds `"ilm_model"` to
#' that list through the option the package provides for the purpose.
#'
#' @return `TRUE` if `marginaleffects` is installed, `FALSE` otherwise,
#'   invisibly.
#' @export
ilm_register_marginaleffects <- function() {
  if (!requireNamespace("marginaleffects", quietly = TRUE)) return(invisible(FALSE))
  ns <- asNamespace("marginaleffects")
  for (g in c("get_coef", "set_coef", "get_vcov", "get_predict"))
    try(registerS3method(g, "ilm_model", get(paste0(g, ".ilm_model")), envir = ns), silent = TRUE)
  ## Registering the S3 methods is NOT sufficient: marginaleffects gates on a
  ## hard allowlist of model classes (sanity_model_supported_class), so an
  ## unregistered class is rejected before dispatch ever happens.  The supported
  ## extension hook is this option, which the allowlist prepends to its own list.
  cur <- getOption("marginaleffects_model_classes", default = character(0))
  if (!"ilm_model" %in% cur) options(marginaleffects_model_classes = c(cur, "ilm_model"))
  invisible(TRUE)
}
