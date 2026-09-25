## ---------------------------------------------------------------------------
## Which groups' random effects a prediction, a fitted value or a residual is
## for. Four words, shared by every function that makes the choice and by code
## built on a fit:
##
##   "fitted"      each group's own estimated effects, so the fit's groups only;
##   "new"         a group the fit has not seen, whose effects are unknown;
##   "typical"     every random effect at zero: a group exactly at the average;
##   "population"  averaged over the distribution of groups.
##
## "Conditional" and "marginal" used to carry this, and each meant two things.
## ilm_fitted(conditional = TRUE) was the fitted groups, where predict()'s help
## called the effects at zero "conditional"; and a marginal effect is a
## derivative whichever groups it is taken for. Each function takes the words
## that apply to it and says where the others are answered. The logical
## arguments the words replace still work, with a warning, for one release.
## ---------------------------------------------------------------------------

ilm_groups_words <- c("fitted", "new", "typical", "population")

## What each word is, and where it is answered when a function does not take it
#' @keywords internal
#' @noRd
ilm_groups_elsewhere <- function(word) {
  switch(word,
    fitted = paste0("each row's own group's effects are used by ",
                    "predict(groups = \"fitted\"), and for the rows the model ",
                    "was fitted to by ilm_fitted(), ilm_scores() and ",
                    "ilm_rqr(); ilm_ranef() gives the effects themselves"),
    new = paste0("a group the fit has not seen has unknown effects, so what ",
                 "it predicts is a spread rather than one value. ",
                 "groups = \"population\" gives the spread's mean, and ",
                 "groups = \"typical\" its value for a group at the average"),
    typical = "groups = \"typical\" sets every random effect to zero",
    population = paste0("predict(groups = \"population\") and ",
                        "ilm_ame(groups = \"population\") average over the ",
                        "groups"))
}

## Settles `groups` for function `fun`, which takes the words in `allowed`, the
## first its default. `given` is !missing(groups) in the caller. `old` is the
## value of the deprecated logical argument `old_name`, and `old_map` gives the
## word each of its values stands for, c(`TRUE` = ..., `FALSE` = ...). A
## logical `groups` is that argument given by position, ilm_fitted(fit, TRUE),
## and is read the same way.
#' @keywords internal
#' @noRd
ilm_groups_arg <- function(groups, allowed, given, fun, old = NULL,
                           old_name = NULL, old_map = NULL) {
  if (is.logical(groups) && !is.null(old_map)) {
    if (!is.null(old))
      stop("`", old_name, "` was given twice, once by position", call. = FALSE)
    old <- groups; given <- FALSE
  }
  if (!is.null(old)) {
    if (given)
      stop("give `groups` or `", old_name, "`, not both: `", old_name,
           "` is the old name for it", call. = FALSE)
    if (!is.logical(old) || length(old) != 1L || is.na(old))
      stop("`", old_name, "` must be TRUE or FALSE", call. = FALSE)
    word <- old_map[[if (old) "TRUE" else "FALSE"]]
    ilm_deprecated_arg(fun, old_name, old_map)
    return(word)
  }
  ## the whole default passed on by a wrapper is the default, as match.arg()
  ## takes it
  if (!given || identical(groups, allowed)) return(allowed[1L])
  if (!is.character(groups) || length(groups) != 1L || is.na(groups))
    stop("`groups` must be one word: ",
         paste0("\"", allowed, "\"", collapse = " or "), call. = FALSE)
  word <- ilm_groups_words[pmatch(groups, ilm_groups_words)]
  if (is.na(word))
    stop("`groups` must be \"fitted\", \"new\", \"typical\" or ",
         "\"population\", not \"", groups, "\"", call. = FALSE)
  if (!word %in% allowed)
    stop(fun, " takes groups = ",
         paste0("\"", allowed, "\"", collapse = " or "), ", not \"", word,
         "\": ", ilm_groups_elsewhere(word), ".", call. = FALSE)
  word
}

## Once per session for each function and argument, so a loop is not buried
## in copies of one message
ilm_deprecation_seen <- new.env(parent = emptyenv())

#' @keywords internal
#' @noRd
ilm_deprecated <- function(key, msg) {
  if (isTRUE(ilm_deprecation_seen[[key]])) return(invisible())
  assign(key, TRUE, envir = ilm_deprecation_seen)
  warning(msg, " It still works, but will be removed after the next release. ",
          "This warning is shown once per session.", call. = FALSE)
  invisible()
}

#' @keywords internal
#' @noRd
ilm_deprecated_arg <- function(fun, old_name, old_map) {
  ilm_deprecated(paste(fun, old_name), paste0(
    "`", old_name, "` in ", fun, " is deprecated: use groups = \"",
    old_map[["TRUE"]], "\" for ", old_name, " = TRUE, and groups = \"",
    old_map[["FALSE"]], "\" for ", old_name, " = FALSE."))
}
