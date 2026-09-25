## ---------------------------------------------------------------------------
## insight methods.
##
## `insight` is the layer the easystats packages read models through: register
## a handful of methods for it and `parameters`, `performance` and `report` all
## work on an ilm_model without knowing anything about it.
##
## Two reasons to do this rather than treat those packages as rivals. It costs
## a user nothing to have the option, and it gives illume an INDEPENDENT
## implementation to check against -- which, in this package's history, is the
## only thing that has ever found a real bug. [ilm_interpret()] writes prose the
## easystats stack cannot, because it can speak to illume's own diagnostics; on
## the parts that overlap, the two should agree, and the tests check that they
## do.
##
## Registration is deferred rather than declared, so insight stays in Suggests.
## ---------------------------------------------------------------------------

#' @keywords internal
#' @noRd
ilm_ins_formula <- function(x, ...) {
  f <- stats::formula(x)
  re <- names(x$re)
  out <- list(conditional = ilm_nobars(f))
  if (length(re)) {
    bars <- ilm_findbars(f)
    if (length(bars))
      out$random <- lapply(bars, function(b)
        ilm_one_sided(b, environment(f) %||% parent.frame()))
  }
  class(out) <- c("insight_formula", "list")
  out
}

#' @keywords internal
#' @noRd
ilm_ins_response <- function(x, ...) ilm_mf_name(stats::formula(x)[[2]])

#' @keywords internal
#' @noRd
ilm_ins_predictors <- function(x, effects = "fixed", flatten = FALSE, ...) {
  f <- stats::formula(x)
  cond <- all.vars(stats::delete.response(stats::terms(ilm_nobars(f))))
  re <- names(x$re)
  out <- list(conditional = cond)
  if (length(re) && effects %in% c("all", "random")) out$random <- re
  if (isTRUE(flatten)) unique(unlist(out, use.names = FALSE)) else out
}

#' @keywords internal
#' @noRd
ilm_ins_data <- function(x, ...) {
  if (!is.null(x$model)) return(ilm_data(x))
  stop("the fit did not keep its model frame", call. = FALSE)
}

#' @keywords internal
#' @noRd
ilm_ins_nobs <- function(x, ...) stats::nobs(x)

#' @keywords internal
#' @noRd
ilm_ins_parameters <- function(x, ...) {
  b <- stats::coef(x)
  data.frame(Parameter = names(b), Estimate = unname(b),
             stringsAsFactors = FALSE)
}

#' @keywords internal
#' @noRd
ilm_ins_find_parameters <- function(x, flatten = FALSE, ...) {
  out <- list(conditional = names(stats::coef(x)))
  if (isTRUE(flatten)) unlist(out, use.names = FALSE) else out
}

#' @keywords internal
#' @noRd
ilm_ins_varcov <- function(x, ...) suppressWarnings(stats::vcov(x))

#' @keywords internal
#' @noRd
ilm_ins_model_info <- function(x, ...) {
  fam <- if (!is.null(x$family)) x$family$name else "gaussian"
  ## the family object's own link, which a switch on the name got wrong for
  ## every family it did not list (beta and the ordinal links said "log"); a
  ## flexible parametric survival model keeps "log", as its hazard is on that
  ## scale
  lk <- if (fam %in% c("rp", "rp_odds", "rp_normal")) "log"
        else if (!is.null(x$family$link)) x$family$link else "identity"
  out <- list(
    is_binomial = fam == "binomial", is_count = fam %in% c("poisson", "nbinom"),
    is_poisson = fam == "poisson", is_negbin = fam == "nbinom",
    is_linear = fam == "gaussian", is_logit = fam %in% c("binomial", "multinomial"),
    is_probit = FALSE, is_ordinal = FALSE, is_multinomial = fam == "multinomial",
    is_categorical = fam == "multinomial", is_survival = fam %in%
      c("weibull", "lognormal", "loglogistic", "rp", "rp_odds", "rp_normal"),
    is_mixed = length(x$re) > 0L, is_bayesian = FALSE, is_gam = FALSE,
    is_zero_inflated = FALSE, is_dispersion = !is.null(x$disp_coef),
    is_censored = !is.null(x$censor), is_truncated = FALSE,
    is_survey = FALSE, is_trial = FALSE, is_exponential = FALSE,
    link_function = lk, family = fam,
    n_obs = tryCatch(stats::nobs(x), error = function(e) NA_integer_),
    model_terms = ilm_ins_predictors(x))
  class(out) <- c("insight_model_info", "list")
  out
}

#' Make illume models readable by the easystats packages
#'
#' Registers the `insight` methods an `ilm_model` needs, after which
#' `parameters::model_parameters()`, `performance::model_performance()` and
#' `report::report()` work on one. `insight` stays in Suggests, so this is
#' called rather than declared.
#'
#' This is not a substitute for [ilm_interpret()], which says things the
#' easystats stack cannot -- what illume's own diagnostics found, and what that
#' implies for the estimates. It is there because it costs a user nothing to
#' have the option, and because an independent implementation of the same
#' quantities is worth having to check against.
#'
#' @return `TRUE` invisibly if the methods were registered, `FALSE` if
#'   `insight` is not installed.
#' @seealso [ilm_register_marginaleffects()], [ilm_interpret()].
#' @examples
#' ilm_register_insight()
#' @export
ilm_register_insight <- function() {
  if (!requireNamespace("insight", quietly = TRUE)) return(invisible(FALSE))
  ns <- asNamespace("insight")
  reg <- list(find_formula = ilm_ins_formula, find_response = ilm_ins_response,
              find_predictors = ilm_ins_predictors, get_data = ilm_ins_data,
              n_obs = ilm_ins_nobs, get_parameters = ilm_ins_parameters,
              find_parameters = ilm_ins_find_parameters,
              get_varcov = ilm_ins_varcov, model_info = ilm_ins_model_info)
  for (g in names(reg))
    try(registerS3method(g, "ilm_model", reg[[g]], envir = ns), silent = TRUE)
  invisible(TRUE)
}
