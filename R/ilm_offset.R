## ---------------------------------------------------------------------------
## Offsets.
##
## A term of the linear predictor whose coefficient is fixed at one: for a
## count of events over an exposure, offset(log(exposure)), so the model is
## for the rate. It is written in the formula as in glm() and glmmTMB, and it
## enters the one linear predictor the conditional mean has -- not a zero
## part or a dispersion model.
## ---------------------------------------------------------------------------

## The offset() terms of a formula, as text for a model frame's formula.
#' @keywords internal
#' @noRd
ilm_offset_terms <- function(formula) {
  tt <- stats::terms(formula)
  io <- attr(tt, "offset")
  if (is.null(io)) return(character(0))
  vars <- as.list(attr(tt, "variables"))[-1L]
  vapply(vars[io], deparse1, "")
}

## New data for a model with an offset has to carry what the offset is made
## of: without it there is no exposure to predict at. Said here, by name,
## rather than as model.frame()'s "object 'e' not found".
#' @keywords internal
#' @noRd
ilm_offset_need <- function(object, newdata) {
  if (is.null(object$offset) || !is.data.frame(newdata)) return(invisible())
  ot <- ilm_offset_terms(object$terms)
  miss <- setdiff(all.vars(parse(text = paste(ot, collapse = " + "))), names(newdata))
  if (length(miss))
    stop("the model has an offset, ", paste(ot, collapse = " + "),
         ", and `newdata` has no column ", paste(sQuote(miss), collapse = " or "),
         " to make it from", call. = FALSE)
  invisible()
}

## ---- means, effects and scenarios: per unit of exposure --------------------
##
## A mean, a marginal effect or a scenario from a model with an offset is
## reported per unit of exposure -- the offset at zero, so for
## offset(log(exposure)) the rate per one unit of it, cases per person-year
## -- and says so. The argument `per` (Craig's item 166) rescales the offset,
## the exposure, to another value: per = 1e5 for a rate per 100,000. Its name
## is kept here, for the messages and notes. predict() keeps each row's own
## exposure unless told otherwise, as glm() does.
#' @keywords internal
#' @noRd
ilm_per_arg <- "per"

## The offset's value at an exposure: the offset's expression evaluated with
## each of its variables set to `per`. NULL (per unit) is an offset of 0;
## "unit" says the same inside predict(). NULL when the model has no offset.
#' @keywords internal
#' @noRd
ilm_offset_at <- function(object, per = NULL) {
  if (is.null(object$offset)) return(NULL)
  if (is.null(per) || identical(per, "unit")) return(0)
  if (!is.numeric(per) || length(per) != 1L || !is.finite(per) || per <= 0)
    stop("`", ilm_per_arg, "` must be one positive number, the exposure ",
         "to report at", call. = FALSE)
  ot <- ilm_offset_terms(object$terms)
  vv <- all.vars(parse(text = paste(ot, collapse = " + ")))
  env <- list2env(stats::setNames(rep(list(per), length(vv)), vv),
                  parent = if (is.null(environment(object$formula))) baseenv()
                           else environment(object$formula))
  val <- sum(vapply(ot, function(t) as.numeric(eval(str2lang(t), env)), 0))
  if (!is.finite(val))
    stop("the offset, ", paste(ot, collapse = " + "), ", is not finite at ",
         ilm_per_arg, " = ", format(per), call. = FALSE)
  val
}

## What a result from a model with an offset is per: said in its print.
#' @keywords internal
#' @noRd
ilm_per_note <- function(object, per = NULL) {
  if (is.null(object$offset)) return(NULL)
  ot <- ilm_offset_terms(object$terms)
  vv <- all.vars(parse(text = paste(ot, collapse = " + ")))
  if (is.null(per))
    paste0("per unit of exposure: the offset, ", paste(ot, collapse = " + "),
           ", at 0 (", paste(vv, collapse = ", "), " = 1 for a log offset); `",
           ilm_per_arg, " =` sets another")
  else
    paste0("at ", paste(vv, collapse = ", "), " = ", format(per, big.mark = ","),
           ", through the offset ", paste(ot, collapse = " + "))
}
