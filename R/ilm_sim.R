## illume: whole simulated datasets, ready to analyse.
##
## ilm_simulate() returns responses -- a matrix, one column per draw --
## because that is what a bootstrap or an envelope needs. Rehearsing an
## analysis needs the data themselves: the rows the model was fitted to, or a
## scaffold's design, with the simulated response in the response's own
## column, so the analysis can be run on it exactly as it will be run on the
## real data. This lays the same draws out that way, and draws a study of
## another size or effect the way ilm_power() does, through the same code.

#' Simulated datasets from a model
#'
#' Whole datasets drawn from a fitted model or an [ilm_scaffold()]: each
#' row's own columns, with a simulated response in place of the observed one.
#' A shortcut for [ilm_simulate()], whose matrix of responses would otherwise
#' have to be bound back to the rows by hand.
#'
#' @section The rows:
#' The columns are the ones the model reads: the response, the variables its
#' terms are built from, its grouping factors and weights, and the variables
#' of an offset, a dispersion model and a zero part. The rows are the ones it
#' was fitted to, so a row it dropped for a missing value is not there. For a
#' scaffold they are its design grid.
#'
#' With neither `n` nor `effect`, the responses are [ilm_simulate()]'s, draw
#' for draw with the same `seed`. The rows and their group labels are the
#' fitted ones, and each dataset draws new random effects for those groups,
#' so it is a new sample of groups under the old labels.
#'
#' @section A study of another size or effect:
#' `n` draws a study of about that many rows, as [ilm_power()] draws one. A
#' fitted model's rows are resampled, whole clusters at a time when it has a
#' grouping factor, and a scaffold's design is laid out again at that size.
#' A cluster drawn twice is two clusters: the unit is renumbered by the copy
#' drawn, 1, 2, and so on, a grouping factor nested in it is prefixed with
#' its copy (`"2.b"`), and one crossed with it, such as items every
#' participant sees, keeps its labels. A row with a frequency weight is drawn
#' as that many rows, each with weight 1.
#'
#' `effect` draws the data with one coefficient, or one term, changed, as in
#' [ilm_power()]: for a single coefficient, its value on the link scale; for
#' a term of several coefficients, a multiple of them, so `effect = 0.5`
#' halves it. `term` names it, and defaults to the first coefficient that is
#' not an intercept.
#'
#' @section The response:
#' The simulated response goes in the response's own column, of the same
#' type: a factor stays a factor with the same levels, a logical stays
#' logical, and a whole-number count stays integer. The response has to be a
#' column the formula names, not a calculation on one such as `log(y)`; for
#' that, [ilm_simulate()] returns the response on the scale the model fits.
#'
#' A censored response is censored as [ilm_simulate()] censors it, and a
#' column beside it says what [ilm_model()]'s `censor` needs to fit it again:
#' * a survival time gets `.event`, `1` where the event was drawn and `0`
#'   where the time was cut at the subject's censoring time, as [ilm_surv()]
#'   takes it: `ilm_surv(d$time, d$.event)`;
#' * codes given to [ilm_censor()] as `status` are held fixed by the draw, and
#'   come back as `.censored`: `ilm_censor(d$y, status = d$.censored)`;
#' * a floor or ceiling needs no column, because the drawn values are held at
#'   it and [ilm_censor()] finds them again: `ilm_censor(d$y, lower = 0)`.
#'
#' @param object A fitted [ilm_model()], or an [ilm_scaffold()].
#' @param nsim Integer. Number of datasets.
#' @param seed Integer. Random seed, for reproducibility. With `NULL` the
#'   draws come from the current random stream, as [ilm_simulate()]'s do.
#' @param n The number of rows in each simulated study, as in [ilm_power()];
#'   `NULL` keeps the rows the model was fitted to.
#' @param effect The tested coefficient's value on the link scale, or a
#'   multiple of a term's coefficients, as in [ilm_power()]; `NULL` keeps the
#'   fitted values.
#' @param term The coefficient or term `effect` sets, named as in
#'   [ilm_power()].
#' @return A data frame when `nsim = 1`, and otherwise a list of `nsim` data
#'   frames.
#' @seealso [ilm_simulate()] for the responses alone, [ilm_power()] to count
#'   detections over many simulated studies, [ilm_refit()] to refit the same
#'   model, and [illumex::ilm_example_data()] for example data of a fixed
#'   design.
#' @examples
#' set.seed(1)
#' d <- data.frame(id = factor(rep(1:12, each = 6)), x = rnorm(72))
#' d$y <- 1 + 0.5 * d$x + rnorm(12)[d$id] + rnorm(72)
#' fit <- ilm_model(y ~ x + (1 | id), data = d, family = "gaussian",
#'                  verbose = FALSE)
#' head(ilm_sim(fit, seed = 2))
#'
#' ## rehearse the analysis on a study of 120 rows, with half the slope
#' d2 <- ilm_sim(fit, n = 120, effect = 0.25, term = "x", seed = 3)
#' coef(ilm_model(y ~ x + (1 | id), data = d2, family = "gaussian",
#'                verbose = FALSE))
#'
#' ## from a scaffold: three draws of a planned two-arm trial
#' s <- ilm_scaffold(y ~ arm, design = list(arm = c("control", "treatment")),
#'                   n_unit = 60, cells = c(control = 12, treatment = 14.5),
#'                   sd = 4, verbose = FALSE)
#' sims <- ilm_sim(s, nsim = 3, seed = 4)
#' sapply(sims, function(d) diff(tapply(d$y, d$arm, mean)))
#' @export
ilm_sim <- function(object, nsim = 1L, seed = NULL, n = NULL, effect = NULL,
                    term = NULL) {
  ilm_rng_restore(seed)                  # the user's random stream, put back on exit
  if (!inherits(object, "ilm_model"))
    stop("`object` must be a fitted ilm_model or an ilm_scaffold, not ",
         class(object)[1], call. = FALSE)
  nsim <- suppressWarnings(as.integer(nsim))
  if (length(nsim) != 1L || is.na(nsim) || nsim < 1L)
    stop("`nsim` must be a whole number of datasets, 1 or more.",
         call. = FALSE)
  for (a in c("n", "effect")) {
    v <- get(a)
    if (!is.null(v) && (!is.numeric(v) || length(v) != 1L || !is.finite(v)))
      stop("`", a, "` must be one number: ilm_sim() draws one study design ",
           "at a time, and ilm_power() is the one that takes several.",
           call. = FALSE)
  }
  if (!is.null(term) && is.null(effect))
    stop("`term` names the coefficient `effect` changes; give `effect` too, ",
         "or leave both out.", call. = FALSE)
  D <- ilm_sim_rows(object)
  B <- if (!is.null(effect))
    ilm_power_beta(object, ilm_power_target(object, term), effect)

  if (is.null(n)) {
    fit <- object
    if (!is.null(B)) fit$beta[] <- B
    Y <- ilm_simulate(fit, nsim, seed)
    out <- lapply(seq_len(nsim), function(s)
      ilm_sim_fill(object, D, D, Y[, s], object$censor))
  } else {
    if (!is.null(object$design))
      stop("this model was fitted to a complex sample, and a study of ",
           "another size is drawn as a simple random sample of rows or of ",
           "clusters, not by the survey's design. Leave `n` out to simulate ",
           "the sampled rows themselves.", call. = FALSE)
    n <- as.integer(round(n))
    if (n < 10L) stop("`n` below 10 is not a study", call. = FALSE)
    if (is.null(B)) B <- ilm_power_beta(object, list(idx = integer(0)), 1)
    study <- ilm_power_study(object)
    if (!is.null(seed)) set.seed(seed)
    out <- lapply(seq_len(nsim), function(s) {
      dr <- study$draw(n, B)
      rows <- if (is.null(dr$d)) ilm_sim_resample(object, D, dr, study) else dr$d
      ilm_sim_fill(object, D, rows, dr$y, dr$stub$censor)
    })
  }
  if (nsim == 1L) out[[1L]] else out
}

## The rows a fit was made from, as data: the columns the model reads, under
## their own names. A column the model frame computed -- "log(x)", a
## polynomial basis, offset(log(e)) -- is left out, and the variable it was
## computed from is there instead.
#' @keywords internal
#' @noRd
ilm_sim_rows <- function(object) {
  fo <- object$formula
  lhs <- if (inherits(fo, "formula") && length(fo) == 3L) fo[[2L]]
  if (is.null(object$model) || is.null(lhs))
    stop("ilm_sim() returns the rows a model was fitted to, and this one was ",
         "not fitted through a formula, so it has none to return. ",
         "ilm_simulate() draws its responses.", call. = FALSE)
  if (!is.name(lhs))
    stop("ilm_sim() puts the simulated response in the response's own ",
         "column, and `", deparse1(lhs), "` is a calculation on one. Make it ",
         "a column of the data and fit the model to that, or use ",
         "ilm_simulate(), which returns the response on the scale the model ",
         "fits.", call. = FALSE)
  d <- ilm_data(object)
  vars <- unique(c(as.character(lhs), all.vars(attr(object$model, "terms")),
                   names(object$data_extra)))
  d <- d[, names(d) %in% vars, drop = FALSE]
  attr(d, "terms") <- NULL
  attr(d, "na.action") <- NULL
  d
}

## A resampled study's rows, laid out as data. A cluster drawn twice has to
## be two clusters here too, labelled as ilm_power_stub() codes them: the
## unit by its copy, a grouping factor nested in it by copy and label, and
## one crossed with it unchanged.
#' @keywords internal
#' @noRd
ilm_sim_resample <- function(object, D, dr, study) {
  out <- D[dr$rows, , drop = FALSE]
  rownames(out) <- NULL
  u <- study$unit
  if (!is.null(u)) {
    gv <- unique(c(unlist(lapply(object$bars, function(b) all.vars(b[[3L]]))),
                   if (!is.null(object$ar$vars)) object$ar$vars[["group"]]))
    for (v in intersect(gv, names(D))) {
      g <- D[[v]]
      if (is.factor(g)) g <- droplevels(g)
      if (!ilm_nested_in(g, u)) next
      unit_itself <- ilm_nested_in(u, g)
      new <- if (unit_itself) dr$copy else paste(dr$copy, g[dr$rows], sep = ".")
      out[[v]] <- if (is.factor(g)) factor(new)
                  else if (is.numeric(g) && unit_itself) as.integer(new)
                  else as.character(new)
    }
  }
  ## each weighted row was drawn as that many rows, and each of them is one
  wexpr <- object$call$weights
  if (isTRUE(study$expanded) && is.name(wexpr) &&
      as.character(wexpr) %in% names(out))
    out[[as.character(wexpr)]] <- 1
  out
}

## One dataset: the model's columns of `rows`, in D's order, with the
## simulated response `y` in the response's column, and the censoring codes
## where the response is censored.
#' @keywords internal
#' @noRd
ilm_sim_fill <- function(object, D, rows, y, censor) {
  rn <- as.character(object$formula[[2L]])
  rows[[rn]] <- ilm_sim_as(object, D[[rn]], y)
  out <- rows[, intersect(names(D), names(rows)), drop = FALSE]
  if (!is.null(censor)) out <- ilm_sim_censor(out, censor, y)
  out
}

## A drawn response as the type of the column it replaces. Categories come
## from ilm_simulate() as codes, 1 to J, and a two-level binomial response as
## 0 and 1; each goes back to the labels the data used.
#' @keywords internal
#' @noRd
ilm_sim_as <- function(object, old, y) {
  lv <- object$ylevels
  if (!is.null(lv)) {
    catg <- is.null(object$family) || isTRUE(object$ordinal) ||
      identical(object$family$name, "multinomial")
    if (!catg && is.logical(old)) return(as.logical(y))
    lab <- if (catg) lv[y] else lv[y + 1L]
    if (is.factor(old))
      return(factor(lab, levels = levels(old), ordered = is.ordered(old)))
    if (is.character(old)) return(lab)
    v <- as.numeric(lab)
    return(if (is.integer(old)) as.integer(v) else v)
  }
  if (is.integer(old) && all(y == round(y))) return(as.integer(y))
  as.numeric(y)
}

## What a censored response needs beside it to be fitted again. A survival
## time drawn at or beyond its censoring time was cut there, so `.event` is
## the ilm_surv() code of each. Codes given as `status` are held fixed by
## ilm_simulate(), so they are the data's, as `.censored`. A floor or ceiling
## needs nothing: ilm_censor(y, lower, upper) finds the codes from the drawn
## values, by the rule the draw was held at them by.
#' @keywords internal
#' @noRd
ilm_sim_censor <- function(out, spec, y) {
  ct <- attr(spec, "ctime")
  if (!is.null(ct)) {
    out$.event <- as.integer(y < ct)
    return(out)
  }
  lo <- attr(spec, "lower"); up <- attr(spec, "upper")
  if ((is.null(lo) || is.na(lo)) && (is.null(up) || is.na(up)))
    out$.censored <- as.integer(unclass(spec))
  out
}
