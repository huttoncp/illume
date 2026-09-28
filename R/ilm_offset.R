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
