## ---------------------------------------------------------------------------
## Censored responses: floors, ceilings, and detection limits.
##
## A censored observation is one you know the interval of but not the value. A
## reading of "0" from an assay with a detection limit of 0 does not mean zero,
## it means "somewhere at or below the limit", and fitting it as a zero pulls
## the mean down and the variance in. The likelihood contribution changes from a
## density at a point to the probability of the interval:
##
##   observed          f(y)
##   left  (floor)     F(lower)
##   right (ceiling)   1 - F(upper)
##
## That is the whole idea, and it is the same mechanism a parametric survival
## model uses for right-censored follow-up.
##
## WHY THE LIMITS ARE CARRIED, NOT JUST THE CODES.
##
## Every envelope diagnostic in this package simulates from the fit and refits.
## A simulated dataset has to be censored the way the real one was, or the
## envelope is built from a process that is not the one being checked. Holding
## the observed pattern fixed would understate the variability, because whether
## a particular draw lands beyond the limit is itself random. Carrying `lower`
## and `upper` means a simulated response can be censored exactly, and the codes
## recomputed from it, with no extra plumbing at the refit sites.
## ---------------------------------------------------------------------------

#' Mark censored observations
#'
#' Describes which observations are known only as an interval: values at or
#' below a floor, at or above a ceiling, or both.
#'
#' @details
#' Pass the result to `ilm_model(censor = )`. With `lower` or `upper` given, the
#' codes are derived from the response, and the limits are carried along so that
#' the simulation-based diagnostics can censor their replicates the same way.
#'
#' `status` is for the case where the limits vary by observation and only the
#' outcome is known -- right-censored follow-up in a survival study, most
#' commonly. Diagnostics then hold the censoring pattern fixed across
#' replicates rather than re-drawing it, which is stated where it matters.
#'
#' @param y The response.
#' @param lower Floor. Values at or below it are left-censored.
#' @param upper Ceiling. Values at or above it are right-censored.
#' @param status Alternatively, codes given directly: `-1` left-censored, `0`
#'   observed, `1` right-censored.
#' @return An integer vector of codes, classed `"ilm_censor"`, carrying the
#'   limits it was built from.
#' @seealso [ilm_model()], [ilm_describe()], whose `p_zero` column is often the
#'   first sign of a floor.
#' @examples
#' y <- c(0, 0, 1.4, 2.9, 5, 5)
#' ilm_censor(y, lower = 0, upper = 5)
#' table(ilm_censor(y, lower = 0))
#' @export
ilm_censor <- function(y, lower = NA, upper = NA, status = NULL) {
  if (!is.null(status)) {
    s <- suppressWarnings(as.integer(status))
    if (length(s) != length(y))
      stop("`status` has ", length(s), " values but `y` has ", length(y),
           call. = FALSE)
    if (anyNA(s) || !all(s %in% c(-1L, 0L, 1L)))
      stop("`status` must be -1 (left-censored), 0 (observed) or 1 ",
           "(right-censored). Values seen: ",
           paste(utils::head(sort(unique(status)), 5), collapse = ", "),
           ". For survival follow-up where 1 means the event was observed, ",
           "pass status = 1 - event.", call. = FALSE)
    return(structure(s, lower = NA_real_, upper = NA_real_,
                     class = "ilm_censor"))
  }
  yn <- suppressWarnings(as.numeric(y))
  if (any(is.na(yn) & !is.na(y)))
    stop("`y` must be numeric to compare against a limit, not ", class(y)[1],
         call. = FALSE)
  lower <- suppressWarnings(as.numeric(lower)[1])
  upper <- suppressWarnings(as.numeric(upper)[1])
  if (is.na(lower) && is.na(upper))
    stop("give `lower`, `upper` or `status`: ilm_censor() needs to know what ",
         "makes an observation censored.", call. = FALSE)
  if (!is.na(lower) && !is.na(upper) && lower >= upper)
    stop("`lower` (", lower, ") must be below `upper` (", upper, ")",
         call. = FALSE)
  s <- ilm_censor_codes(yn, lower, upper)
  if (!any(s != 0L))
    warning("no observation reaches the limit", if (!is.na(lower) && !is.na(upper))
            "s" else "", " given, so nothing is censored and the fit will be ",
            "the same as without `censor`.", call. = FALSE)
  structure(s, lower = lower, upper = upper, class = "ilm_censor")
}

## Derive the codes from a response and a pair of limits. Used on the observed
## data, and again on every simulated replicate, so that the two are censored by
## exactly the same rule.
#' @keywords internal
#' @noRd
ilm_censor_codes <- function(y, lower = NA_real_, upper = NA_real_) {
  s <- integer(length(y))
  if (!is.na(lower)) s[!is.na(y) & y <= lower] <- -1L
  if (!is.na(upper)) s[!is.na(y) & y >= upper] <- 1L
  s
}

## Re-derive a censoring vector for a response that may be simulated. When the
## limits are known the codes follow from the values; when only a status was
## given they are held as supplied.
#' @keywords internal
#' @noRd
ilm_censor_for <- function(spec, y) {
  if (is.null(spec)) return(NULL)
  lo <- attr(spec, "lower"); up <- attr(spec, "upper")
  if (is.null(lo)) lo <- NA_real_
  if (is.null(up)) up <- NA_real_
  if (is.na(lo) && is.na(up)) return(as.integer(spec))
  ilm_censor_codes(y, lo, up)
}

## Apply the limits to a simulated response, so a replicate is censored the way
## the data were. Whether a given draw lands beyond a limit is itself random,
## which is exactly the variability the envelope needs to include.
#' @keywords internal
#' @noRd
ilm_censor_apply <- function(spec, y) {
  if (is.null(spec)) return(y)
  lo <- attr(spec, "lower"); up <- attr(spec, "upper")
  if (!is.null(lo) && !is.na(lo)) y <- pmax(y, lo)
  if (!is.null(up) && !is.na(up)) y <- pmin(y, up)
  y
}

#' Print a censoring specification
#'
#' @param x An `"ilm_censor"` object.
#' @param ... Unused.
#' @return `x`, invisibly.
#' @export
print.ilm_censor <- function(x, ...) {
  n <- length(x); lo <- attr(x, "lower"); up <- attr(x, "upper")
  cat(sprintf("censoring: %d observed, %d left, %d right (of %d)\n",
              sum(x == 0L), sum(x == -1L), sum(x == 1L), n))
  if (!is.na(lo) || !is.na(up))
    cat(sprintf("  limits: lower %s, upper %s\n",
                if (is.na(lo)) "none" else format(lo),
                if (is.na(up)) "none" else format(up)))
  else
    cat("  limits not recorded; diagnostics hold the pattern fixed\n")
  p <- mean(x != 0L)
  if (p > 0.5)
    cat(sprintf("  %.0f%% censored; the fit rests mostly on interval probabilities\n",
                100 * p))
  invisible(x)
}

#' Right-censored follow-up, in the survival convention
#'
#' A convenience wrapper for the commonest censoring there is: a study where
#' some subjects had the event and the rest were still event-free when
#' follow-up ended.
#'
#' @details
#' The argument is `event`, in the convention `survival::Surv()` uses: `1` when
#' the event was observed, `0` when the subject was censored. That is the
#' opposite of the code [ilm_censor()] stores, which is why this wrapper exists
#' -- getting it backwards silently fits the model to the wrong subjects, and
#' the fit will not look wrong.
#'
#' @section What the diagnostics can and cannot do here:
#' Censoring times vary by subject and are not observed for anyone who had the
#' event, so unlike a detection limit they cannot be re-drawn. The
#' simulation-based diagnostics therefore hold the censoring pattern fixed
#' across replicates. That is the usual convention, and it makes the envelopes
#' slightly narrower than a fully random-censoring account would.
#'
#' @param time Follow-up time, strictly positive.
#' @param event `1` or `TRUE` if the event was observed, `0` or `FALSE` if the
#'   subject was censored.
#' @return An `"ilm_censor"` object, to pass as `ilm_model(censor = )`.
#' @seealso [ilm_censor()] for floors and ceilings, [ilm_model()] with
#'   `family = "weibull"`, `"lognormal"` or `"loglogistic"`.
#' @examples
#' t <- c(5, 9, 12, 3, 20)
#' e <- c(1, 0, 1, 1, 0)
#' ilm_surv(t, e)
#' @export
ilm_surv <- function(time, event) {
  if (missing(event))
    stop("`event` is required: without it there is no way to tell a subject ",
         "who had the event from one who was still event-free when follow-up ",
         "ended.", call. = FALSE)
  if (length(time) != length(event))
    stop("`time` has ", length(time), " values and `event` has ",
         length(event), "; they must be the same length.", call. = FALSE)
  e <- if (is.logical(event)) as.integer(event) else suppressWarnings(as.integer(event))
  if (anyNA(e) || !all(e %in% c(0L, 1L)))
    stop("`event` must be 1 (the event was observed) or 0 (censored). Values ",
         "seen: ", paste(utils::head(sort(unique(event)), 5), collapse = ", "),
         ".", call. = FALSE)
  tn <- suppressWarnings(as.numeric(time))
  if (any(is.na(tn) & !is.na(time)))
    stop("`time` must be numeric", call. = FALSE)
  if (any(tn <= 0, na.rm = TRUE))
    stop("follow-up time must be strictly positive; ", sum(tn <= 0, na.rm = TRUE),
         " value(s) are zero or negative. A time of exactly 0 usually means an ",
         "event before the first assessment.", call. = FALSE)
  if (!any(e == 1L))
    stop("no subject had the event, so there is nothing for the model to ",
         "estimate a time to.", call. = FALSE)
  if (mean(e == 1L) < 0.05)
    warning(sprintf("only %.1f%% of subjects had the event. Check that `event` ",
                    100 * mean(e == 1L)),
            "is 1 for the event and not the other way round.", call. = FALSE)
  structure(1L - e, lower = NA_real_, upper = NA_real_, class = "ilm_censor")
}
