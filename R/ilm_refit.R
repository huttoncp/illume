## ---------------------------------------------------------------------------
## Refit a model: to a new response for the rows it was fitted to, or to
## other rows altogether.
##
## The two are what code built on a fit keeps needing. A parametric bootstrap
## draws responses from the fit with ilm_simulate() and refits each one; a
## rolling forecast refits the same model to the data available at each
## origin. Both used to reach into illume's internals to do it.
## ---------------------------------------------------------------------------

#' Refit a model to a new response or to other data
#'
#' The same model -- formula, family, random effects, correlation over time,
#' dispersion and zero parts, fitting options -- fitted again, either to a new
#' response for the rows it was fitted to or to another data set.
#'
#' @details
#' **`y`** is a new response for the rows the model was fitted to, in the
#' layout [ilm_simulate()] returns: a vector with one value per fitted row, or
#' a single column of its matrix. The design is not rebuilt, so this is the
#' parametric bootstrap's refit. A censored response is censored as the data
#' were: a value at or beyond a row's censoring time is censored there, as
#' [ilm_simulate()] makes it. A categorical response is given as category
#' numbers, as [ilm_simulate()] gives it.
#'
#' **`data`** is another data set -- more rows, fewer, a later window of a
#' time series -- and the model is built from its formula again, so a factor's
#' levels, a smooth's basis and a correlation over time's cells come from the
#' new rows. A correlation over time must then have been given by name,
#' `ilm_ar1(~ time | group)`, so its columns can be found in the new data.
#'
#' With neither, the model is refitted to its own data; with both, it stops.
#'
#' @param fit A fitted `"ilm_model"` from the formula interface.
#' @param data A data frame to fit the same model to.
#' @param y A new response for the rows `fit` was fitted to.
#' @return A fitted `"ilm_model"`, with every method a fit from
#'   [ilm_model()] has.
#' @seealso [ilm_simulate()], which draws responses to refit;
#'   [ilm_apply_remedy()], which refits with one change.
#' @examples
#' set.seed(1)
#' d <- data.frame(id = factor(rep(1:12, each = 6)), x = rnorm(72))
#' d$y <- 1 + 0.5 * d$x + rnorm(12)[d$id] + rnorm(72)
#' fit <- ilm_model(y ~ x + (1 | id), data = d, family = "gaussian",
#'                  verbose = FALSE)
#' ## the parametric bootstrap's refit
#' ystar <- ilm_simulate(fit, nsim = 1, seed = 2)[, 1]
#' coef(ilm_refit(fit, y = ystar))
#' ## the same model on the first eight groups
#' coef(ilm_refit(fit, data = d[d$id %in% 1:8, ]))
#' @export
ilm_refit <- function(fit, data = NULL, y = NULL) {
  if (!inherits(fit, "ilm_model"))
    stop("`fit` must be a fitted ilm_model, not ", class(fit)[1],
         call. = FALSE)
  if (!is.null(data) && !is.null(y))
    stop("give `data` -- the same model on other rows -- or `y` -- a new ",
         "response for the rows it was fitted to -- but not both.",
         call. = FALSE)
  if (!is.null(data)) {
    dname <- substitute(data)
    return(ilm_refit_data(fit, data, dname))
  }
  ilm_refit_y(fit, if (is.null(y)) fit$y else y)
}

## A new response for the fitted rows: the engine refits the stored design,
## and the formula-level parts of the fit -- what predict(), emmeans and the
## rest rebuild a reference grid from -- are carried over.
#' @keywords internal
#' @noRd
ilm_refit_y <- function(fit, y) {
  if (is.data.frame(y)) y <- as.matrix(y)
  if (is.matrix(y)) {
    if (ncol(y) != 1L)
      stop("`y` must be one response: a vector, or a single column of what ",
           "ilm_simulate() returns -- this has ", ncol(y), " columns.",
           call. = FALSE)
    y <- y[, 1L]
  }
  N <- nrow(fit$X)
  if (length(y) != N)
    stop("`y` must have one value per row the model was fitted to (", N,
         "), not ", length(y), ".", call. = FALSE)
  if (anyNA(y))
    stop("`y` has missing values; a refit needs a response for every row ",
         "the model was fitted to.", call. = FALSE)
  ## a categorical response is held as category numbers, as ilm_simulate()
  ## draws it; labels are accepted too
  if (!is.null(fit$ylevels) && (is.factor(y) || is.character(y)))
    y <- match(as.character(y), fit$ylevels)
  r <- ilm_refit_like(fit, y = y)
  keep <- c("design", "call", "family_inferred", "formula", "fixed_formula",
            "terms", "xlev", "contrasts", "model", "data_extra", "smooths",
            "bars", "na.action", "n_dropped", "ylevels", "assign",
            "term_labels")
  for (k in intersect(keep, names(fit))) r[[k]] <- fit[[k]]
  ## the model frame's response is the new one, where it is a plain column
  mf <- r$model
  if (!is.null(mf)) {
    rn <- names(mf)[1L]; old <- mf[[1L]]
    if (is.factor(old) && is.numeric(y) && !is.null(fit$ylevels))
      mf[[rn]] <- factor(fit$ylevels[y], levels = levels(old),
                         ordered = is.ordered(old))
    else if (is.numeric(old) && is.null(dim(old)))
      mf[[rn]] <- as.numeric(y)
    r$model <- mf
  }
  r$edf <- ilm_smooth_edf(r)
  r
}

## Another data set: the model built from its formula again, as
## ilm_apply_remedy() does it, so levels, bases and cells come from the rows.
#' @keywords internal
#' @noRd
ilm_refit_data <- function(fit, data, dname) {
  cl <- fit$call
  if (is.null(cl) || !inherits(fit$formula, "formula"))
    stop("a refit to other data rebuilds the model from its formula, and ",
         "this one did not come from the formula interface. Fit it with ",
         "ilm_model() and a formula first.", call. = FALSE)
  if (!is.data.frame(data))
    stop("`data` must be a data frame.", call. = FALSE)
  ## a correlation over time given as vectors belongs to the old rows
  if (!is.null(fit$ar) && is.null(fit$ar$vars))
    stop("this model's correlation over time was given as vectors, which ",
         "belong to the rows it was fitted to. To refit it to other data, ",
         "give it by name, as ilm_ar1(~ time | group) (or ilm_car1(), ",
         "ilm_rw1()), so its columns are found in the new data.",
         call. = FALSE)
  env0 <- environment(fit$formula)
  if (is.null(env0)) env0 <- parent.frame(2L)
  new <- cl
  new[[1L]] <- as.name("ilm_model_formula")
  new$formula <- fit$formula
  new$data <- quote(.ilm_refit_data)
  env <- new.env(parent = env0)
  assign(".ilm_refit_data", data, envir = env)
  assign("ilm_model_formula", ilm_model_formula, envir = env)
  out <- eval(new, env)
  ## the call as the user would have written it, so update() and the
  ## remedies find these data
  out$call$data <- dname
  out
}
