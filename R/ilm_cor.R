## ---------------------------------------------------------------------------
## Correlation structures over time.
##
## All three describe a latent process shared by observations at the same
## place in time within the same group, and all three are written the same way
## in the likelihood: a first-order Markov chain, factorised conditionally.
##
##   AR(1)   phi is one number, the correlation one step apart, and the steps
##           are all the same size. Times must therefore sit on a common grid.
##
##   CAR(1)  phi varies by transition: phi_k = rho ^ d_k for a gap of d_k time
##           units. Equivalently exp(-d_k / range). This is the Ornstein-
##           Uhlenbeck process, and the conditional factorisation is exact, so
##           irregular spacing costs nothing in accuracy -- only in how many
##           distinct times there are, which is the real constraint.
##
##   RW(1)   A random walk: each step adds an independent change whose
##           variance is proportional to the time elapsed, so there is no level
##           to return to and the variance grows without bound. Each group's
##           walk is held at zero at its first time -- not estimated there --
##           so the level at that time is the fixed effects', and the intercept
##           stays identified. This is the local-level model of structural time
##           series, and like CAR(1) it takes irregular gaps as they are.
##
## THE CONSTRAINT WORTH UNDERSTANDING. These are latent-variable structures, not
## marginal ones. A latent value exists for every distinct (group, time) pair, so
## continuous time with no repeats means one latent per observation, which is the
## regime where the Laplace approximation stops being estimable (see the latent
## budget work: 33% convergence at 1.5 observations per latent). nlme's corCAR1
## does not have this problem because it puts the correlation on the residuals
## instead. The constructors therefore report the budget and say so out loud
## when it is thin -- rounding the time to a coarser grid is usually the fix.
## For a gaussian response none of this binds: the Laplace approximation is
## exact there, which the fit-time checks know and the constructors cannot.
## ---------------------------------------------------------------------------

#' @keywords internal
#' @noRd
ilm_cor_time <- function(time, group, what) {
  if (missing(time) || missing(group))
    stop("`time` and `group` are both required: a correlation over time needs ",
         "to know which observations belong to the same unit, and when each ",
         "was taken.", call. = FALSE)
  if (length(time) != length(group))
    stop("`time` has ", length(time), " values and `group` has ",
         length(group), "; they must be the same length.", call. = FALSE)
  if (is.factor(time))
    stop("`time` is a factor. Converting it here would use the level order ",
         "rather than the labels, and silently change what a gap means. Pass ",
         "as.integer(time) if the levels are equally spaced steps, or the real ",
         "times otherwise.", call. = FALSE)
  tnum <- suppressWarnings(as.numeric(time))
  if (any(is.na(tnum) & !is.na(time)))
    stop("`time` must be numeric, integer or Date, not ", class(time)[1],
         call. = FALSE)
  if (anyNA(tnum) || anyNA(group))
    stop("`time` and `group` cannot contain missing values; ", what,
         " has no place to put an observation whose time is unknown.",
         call. = FALSE)
  ## the scale the times came on, so they can be handed back on it: a gap is
  ## computed in days or seconds, but a person reads a date
  tz <- attr(time, "tzone")
  list(t = tnum, g = as.integer(factor(group)), glev = levels(factor(group)),
       tcls = if (inherits(time, "POSIXct")) "POSIXct"
              else if (inherits(time, "Date")) "Date" else "numeric",
       tz = if (is.null(tz)) "" else tz[1L])
}

## Numeric times back on the scale they were given on.
#' @keywords internal
#' @noRd
ilm_cor_as_time <- function(x, ar) {
  switch(if (is.null(ar$tcls)) "numeric" else ar$tcls,
         Date = as.Date(x, origin = "1970-01-01"),
         POSIXct = as.POSIXct(x, origin = "1970-01-01",
                              tz = if (is.null(ar$tz)) "" else ar$tz),
         x)
}

## Latent cells as the distinct (group, time) pairs, ordered by group then
## time. That ordering is what lets the preceding cell in a group always be the
## one before it, exactly as the evenly spaced case relies on. CAR(1) and the
## random walk share it; they differ only in what a step between cells means.
#' @keywords internal
#' @noRd
ilm_cor_cells <- function(z) {
  o <- order(z$g, z$t)
  gs <- z$g[o]; ts <- z$t[o]; n <- length(o)
  new <- c(TRUE, gs[-1L] != gs[-n] | ts[-1L] != ts[-n])
  idx <- integer(n); idx[o] <- cumsum(new)
  cg <- gs[new]; ct <- ts[new]; n_cell <- length(cg)
  first <- which(c(TRUE, cg[-1L] != cg[-n_cell]))
  rest <- if (n_cell > length(first)) setdiff(seq_len(n_cell), first) else integer(0)
  prev <- rest - 1L
  list(idx = idx, cg = cg, ct = ct, n_cell = n_cell, first = first,
       rest = rest, prev = prev, gap = ct[rest] - ct[prev], n = n)
}

## Shared budget reporting, so every structure says the same thing the same way.
#' @keywords internal
#' @noRd
ilm_cor_budget <- function(n_obs, n_cell, what, verbose) {
  r <- n_obs / n_cell
  if (verbose && r < 1.5)
    warning(what, " puts ", n_cell, " latent values under ", n_obs,
            " observations (", round(r, 2), " per latent). Below about 1.5 the ",
            "Laplace approximation frequently fails to converge, and a latent ",
            "value seen once carries no information the residual does not. ",
            "Round `time` to a coarser grid so observations share a latent.",
            call. = FALSE)
  r
}

## ---- a structure given by name -------------------------------------------
##
## `~ time | group` names two columns of the model's data instead of passing
## them. ilm_model() reads them from its model frame -- after na.action, so a
## dropped row cannot leave the structure out of step with the response, which
## vectors built beforehand from the whole data frame do -- and the structure
## remembers them, so anything that places new rows on it later knows where to
## look.

#' @keywords internal
#' @noRd
ilm_cor_named <- function(type, f, has_group, verbose) {
  what <- switch(type, ar1 = "ilm_ar1()", car1 = "ilm_car1()", rw1 = "ilm_rw1()")
  if (has_group)
    stop(what, " was given a formula and a `group` as well. Name the group ",
         "inside the formula, `~ time | group`, or pass both as vectors.",
         call. = FALSE)
  rhs <- if (length(f) == 2L) f[[2L]] else NULL
  if (!is.call(rhs) || !identical(rhs[[1L]], as.name("|")) ||
      !is.name(rhs[[2L]]) || !is.name(rhs[[3L]]))
    stop(what, " takes a one-sided formula `~ time | group`, naming one ",
         "column of the data for each, such as ~ day | id.", call. = FALSE)
  structure(list(type = type,
                 vars = c(time = as.character(rhs[[2L]]),
                          group = as.character(rhs[[3L]])),
                 verbose = verbose),
            class = c("ilm_cor_named", "ilm_cor"))
}

## The structure itself, from the columns the model frame holds.
#' @keywords internal
#' @noRd
ilm_cor_build <- function(spec, data) {
  v <- spec$vars
  ctor <- switch(spec$type, ar1 = ilm_ar1, car1 = ilm_car1, rw1 = ilm_rw1)
  out <- ctor(data[[v[["time"]]]], data[[v[["group"]]]], verbose = spec$verbose)
  out$vars <- v
  out
}

#' A first-order autoregressive structure over evenly spaced time
#'
#' Observations of the same unit at adjacent time steps share a correlated
#' latent value. Use this when the time index moves in equal steps: study
#' visits, months, waves. For irregular times use [ilm_car1()].
#'
#' @details
#' The latent grid spans every time step between the earliest and the latest
#' observation, so a unit that misses a wave is handled correctly -- the latent
#' for that step simply has no observation attached, and the chain continues
#' across it. What is *not* allowed is times that do not sit on a common grid at
#' all, because then "one step" has no meaning.
#'
#' **By name.** `ilm_ar1(~ time | group)` names two columns of the model's
#' data instead of passing them. [ilm_model()] reads them from its model
#' frame, after rows with missing values are dropped, so the structure cannot
#' come out of step with the response; vectors taken from the whole data frame
#' beforehand can. The fitted model remembers the names, so [ilm_cells()] and
#' predictions for new rows know where to find the time and the group.
#'
#' @param time Time index, one value per observation. Whole numbers, or values
#'   on a common step such as `c(0, 5, 10, 15)`. Or a one-sided formula
#'   `~ time | group` naming two columns of the model's data; see Details.
#' @param group Unit identifier, one value per observation. Omitted when
#'   `time` is a formula.
#' @param verbose Warn when the latent budget is thin.
#' @return A `"ilm_ar1"` specification, to pass as `ilm_model(ar = )`.
#' @seealso [ilm_car1()] for arbitrary gaps, [ilm_rw1()] for a level that
#'   drifts, [ilm_check_ar()] to test whether the structure is needed,
#'   [ilm_cells()] for the fitted cells.
#' @examples
#' d <- ilm_sim(n_id = 20, n_period = 8)
#' a <- ilm_ar1(as.integer(factor(d$date)), d$id)
#' a
#' @export
ilm_ar1 <- function(time, group, verbose = TRUE) {
  if (inherits(time, "formula"))
    return(ilm_cor_named("ar1", time, !missing(group), verbose))
  z <- ilm_cor_time(time, group, "AR(1)")
  u <- sort(unique(z$t))
  if (length(u) < 2L)
    stop("`time` has only one distinct value, so there is no autocorrelation ",
         "to model.", call. = FALSE)
  dif <- diff(u)
  step <- min(dif)
  ratio <- dif / step
  if (any(abs(ratio - round(ratio)) > 1e-6))
    stop("the times do not sit on a common grid -- the gaps between distinct ",
         "times are ", paste(signif(utils::head(sort(unique(dif)), 4), 4),
                             collapse = ", "),
         ", which are not whole multiples of one another. AR(1) needs equal ",
         "steps; use ilm_car1() instead, which takes the gaps as they are.",
         call. = FALSE)

  slot <- as.integer(round((z$t - min(u)) / step)) + 1L
  Tt <- max(slot); n_group <- max(z$g)
  idx <- (z$g - 1L) * Tt + slot
  n_cell <- n_group * Tt

  ## Times can sit on a common grid and still be badly served by one. Sampling
  ## at 0, 1, 2, 6, 15 is "regular" with a step of 1, but it spends sixteen
  ## latent values to describe five observations, and every empty slot is a
  ## parameter the data say nothing about. That is precisely what CAR(1) exists
  ## to avoid, so say so rather than only reporting the thin budget.
  if (verbose && Tt > 2L * length(u))
    warning("the times fit a grid of ", Tt, " steps but only ", length(u),
            " of those steps are ever observed, so most of the ", n_cell,
            " latent values carry no data. ilm_car1() takes the gaps as they ",
            "are and needs one latent per observed time instead.",
            call. = FALSE)
  r <- ilm_cor_budget(length(idx), n_cell, "AR(1)", verbose)
  structure(list(type = "ar1", idx = idx, n_group = n_group, Tt = Tt,
                 n_cell = n_cell, n_latent = n_cell, step = step,
                 origin = min(u), n_obs = length(idx),
                 obs_per_latent = r, glev = z$glev, tcls = z$tcls, tz = z$tz),
            class = c("ilm_ar1", "ilm_cor"))
}

#' A continuous-time autoregressive structure over irregular time
#'
#' The correlation between two observations of the same unit is `rho` raised to
#' the power of the gap between them, so observations close in time are similar
#' and the spacing need not be regular. This is the structure `nlme`'s
#' `corCAR1()` fits, and the AR(1) special case when every gap is one.
#'
#' @details
#' Written as a Markov chain the transition from one observation to the next is
#' `phi_k = rho ^ d_k` with innovation variance `1 - phi_k^2`, which is exact
#' for the Ornstein-Uhlenbeck process rather than an approximation to it.
#' Irregular spacing therefore costs nothing in accuracy.
#'
#' What it does cost is latent values. One exists per distinct (group, time)
#' pair, so genuinely continuous times with no repeats give one latent per
#' observation -- the regime where the Laplace approximation stops converging.
#' The constructor reports the budget and warns when it is thin; rounding `time`
#' to a coarser grid, so that observations share a latent, is the usual fix and
#' is what makes the difference between a model that fits and one that does not.
#'
#' `ilm_car1(~ time | group)` names two columns of the model's data instead of
#' passing them; see [ilm_ar1()] for why that is the safer form.
#'
#' @param time Time, one value per observation. Any numeric scale or `Date`;
#'   `rho` is the correlation one unit apart, so the units matter. Or a
#'   one-sided formula `~ time | group` naming two columns of the model's data.
#' @param group Unit identifier, one value per observation. Omitted when
#'   `time` is a formula.
#' @param verbose Warn when the latent budget is thin.
#' @return A `"ilm_car1"` specification, to pass as `ilm_model(ar = )`.
#' @seealso [ilm_ar1()] for evenly spaced time, [ilm_rw1()] for a level that
#'   drifts rather than reverting, [ilm_check_ar()] and [ilm_plot_acf()] to
#'   test whether the structure is needed, [ilm_cells()] for the fitted cells.
#' @examples
#' set.seed(1)
#' d <- data.frame(id = rep(1:20, each = 6),
#'                 day = as.vector(replicate(20, sort(sample(1:60, 6)))))
#' ilm_car1(d$day, d$id)
#' @export
ilm_car1 <- function(time, group, verbose = TRUE) {
  if (inherits(time, "formula"))
    return(ilm_cor_named("car1", time, !missing(group), verbose))
  z <- ilm_cor_time(time, group, "CAR(1)")
  if (length(unique(z$t)) < 2L)
    stop("`time` has only one distinct value, so there is no autocorrelation ",
         "to model.", call. = FALSE)
  cc <- ilm_cor_cells(z)
  if (!length(cc$rest))
    stop("every unit has only one distinct time, so there is nothing for a ",
         "correlation over time to describe. A plain random intercept, ",
         "(1 | group), is the model you want.", call. = FALSE)
  gap <- cc$gap

  ## gaps are positive by construction; a gap far below the typical one still
  ## drives phi to 1 and the innovation variance to 0
  med <- stats::median(gap)
  if (verbose && min(gap) < med * 1e-6)
    warning("the smallest gap between consecutive times is ", signif(min(gap), 3),
            " against a median of ", signif(med, 3),
            ". Nearly coincident times force the correlation to 1 and the ",
            "innovation variance to 0, which the optimiser handles badly. ",
            "Consider rounding `time`.", call. = FALSE)
  r <- ilm_cor_budget(cc$n, cc$n_cell, "CAR(1)", verbose)
  structure(list(type = "car1", idx = cc$idx, n_cell = cc$n_cell,
                 n_latent = cc$n_cell, first = cc$first, prev = cc$prev,
                 rest = cc$rest, gap = gap, ct = cc$ct, cg = cc$cg,
                 n_group = length(cc$first), n_obs = cc$n, obs_per_latent = r,
                 time_range = range(z$t), glev = z$glev, tcls = z$tcls,
                 tz = z$tz),
            class = c("ilm_car1", "ilm_cor"))
}

#' A random walk over time
#'
#' Observations of the same unit share a latent level that wanders: each step
#' adds an independent change whose variance is proportional to the time that
#' passed, so the level has no mean to return to. This is the local-level
#' model of structural time series. Use it when a unit's level drifts rather
#' than fluctuating about a fixed value -- where [ilm_car1()] would put its
#' correlation at the edge of its range.
#'
#' @details
#' **The walk starts at zero.** Each group's walk is held at zero at its first
#' time, so its level there is the one the fixed effects give, and every later
#' value is measured from it. Without that anchor the walk and the intercept
#' would describe the same thing and neither could be estimated. It also means
#' the groups are taken to start from a common level, give or take the fixed
#' effects; where they do not, add a random intercept, `(1 | group)`, and each
#' starts where it does.
#'
#' **The variance is per unit of time.** Between two times `d` units apart the
#' walk moves by a normal change with variance `d` times the fitted variance,
#' which the fitted model holds as `Sigma$ar`. So the units of `time` matter:
#' days and weeks give variances seven times apart for the same walk. Irregular
#' gaps enter as they are, and cost nothing in accuracy.
#'
#' **What is estimated.** One latent value per distinct (group, time) pair
#' after each group's first. With a gaussian response the Laplace approximation
#' is exact, so the fit is the exact likelihood of the local-level model, by
#' maximum likelihood or by REML (`reml = TRUE`). With other families it is
#' the Laplace approximation, and the fit-time checks judge the number of
#' observations per latent value as they do for [ilm_car1()].
#'
#' `ilm_rw1(~ time | group)` names two columns of the model's data instead of
#' passing them; see [ilm_ar1()] for why that is the safer form. It is also
#' what lets `predict(marginal = TRUE)` place new rows on the walk, whose
#' variance depends on how long after its group's first time each one falls.
#'
#' @param time Time, one value per observation: any numeric scale, a `Date` or
#'   a date-time. The variance is per unit of it. Or a one-sided formula
#'   `~ time | group` naming two columns of the model's data.
#' @param group Unit identifier, one value per observation. Omitted when
#'   `time` is a formula.
#' @param verbose Warn about nearly coincident times.
#' @return A `"ilm_rw1"` specification, to pass as `ilm_model(ar = )`.
#' @references
#' Durbin, J., & Koopman, S. J. (2012). *Time Series Analysis by State Space
#' Methods* (2nd ed.). Oxford University Press. (Chapter 2, the local level
#' model.)
#' @seealso [ilm_car1()] for a process that returns to its mean,
#'   [ilm_cells()] for the fitted cells.
#' @examples
#' set.seed(1)
#' d <- data.frame(t = 1:60, unit = "a")
#' d$y <- 2 + cumsum(c(0, rnorm(59, 0, 0.7))) + rnorm(60)
#' fit <- ilm_model(y ~ 1, data = d, ar = ilm_rw1(~ t | unit),
#'                  verbose = FALSE)
#' fit$Sigma$ar                 # variance of the walk per unit of time
#' head(ilm_cells(fit))
#' @export
ilm_rw1 <- function(time, group, verbose = TRUE) {
  if (inherits(time, "formula"))
    return(ilm_cor_named("rw1", time, !missing(group), verbose))
  z <- ilm_cor_time(time, group, "a random walk")
  cc <- ilm_cor_cells(z)
  if (!length(cc$rest))
    stop("every unit has only one distinct time, so there is no walk to ",
         "follow: each unit's walk is zero at its first time, and with only ",
         "one time the term is zero everywhere.", call. = FALSE)
  gap <- cc$gap
  med <- stats::median(gap)
  if (verbose && min(gap) < med * 1e-6)
    warning("the smallest gap between consecutive times is ", signif(min(gap), 3),
            " against a median of ", signif(med, 3),
            ". Across so short a gap the walk can barely move, which ties the ",
            "two latent values together and makes the fit harder for the ",
            "optimiser. Consider rounding `time`.", call. = FALSE)
  ## The anchors are held at zero, not integrated, so they are not latent
  ## values. The budget is not warned about here: one latent value per
  ## observation is the ordinary local-level model, and whether that is thin
  ## depends on the family, which only the fit knows.
  n_latent <- cc$n_cell - length(cc$first)
  structure(list(type = "rw1", idx = cc$idx, n_cell = cc$n_cell,
                 n_latent = n_latent, first = cc$first, prev = cc$prev,
                 rest = cc$rest, gap = gap, ct = cc$ct, cg = cc$cg,
                 n_group = length(cc$first), n_obs = cc$n,
                 obs_per_latent = cc$n / n_latent, time_range = range(z$t),
                 glev = z$glev, tcls = z$tcls, tz = z$tz),
            class = c("ilm_rw1", "ilm_cor"))
}

#' Print a correlation structure
#'
#' @param x An `"ilm_cor"` object.
#' @param ... Unused.
#' @return `x`, invisibly.
#' @export
print.ilm_cor <- function(x, ...) {
  lab <- switch(x$type, car1 = "CAR(1), continuous time",
                rw1 = "Random walk, continuous time",
                "AR(1), evenly spaced time")
  cat(lab, "\n")
  if (inherits(x, "ilm_cor_named")) {
    cat(sprintf("  over `%s` within `%s`, built from the data when the model is fitted\n",
                x$vars[["time"]], x$vars[["group"]]))
    return(invisible(x))
  }
  rw <- identical(x$type, "rw1")
  nl <- if (is.null(x$n_latent)) x$n_cell else x$n_latent
  cat(sprintf("  %d group%s, %d latent value%s%s, %d observation%s (%.2f per latent)\n",
              x$n_group, if (x$n_group == 1L) "" else "s",
              nl, if (nl == 1L) "" else "s",
              if (rw) sprintf(" after the %d held at zero", x$n_group) else "",
              x$n_obs, if (x$n_obs == 1L) "" else "s", x$obs_per_latent))
  if (x$type %in% c("car1", "rw1"))
    cat(sprintf("  gaps: min %s, median %s, max %s\n",
                signif(min(x$gap), 4), signif(stats::median(x$gap), 4),
                signif(max(x$gap), 4)))
  else
    cat(sprintf("  step: %s, grid of %d time slots\n", signif(x$step, 4), x$Tt))
  if (x$obs_per_latent < 1.5)
    cat(if (rw) "  whether that is thin depends on the family; the fit checks it\n"
        else "  the latent budget is thin; see ?ilm_car1\n")
  invisible(x)
}

## Accept the constructors above and the bare list the fitter took before
## they existed, so nothing written against the old contract breaks.
#' @keywords internal
#' @noRd
ilm_as_cor <- function(ar) {
  if (inherits(ar, "ilm_cor_named"))
    stop("a correlation over time given by name, `~ time | group`, is built ",
         "from the model's data, which ilm_model() reads and ilm_fit() does ",
         "not see. Fit it with ilm_model(), or give the time and the group ",
         "as vectors.", call. = FALSE)
  if (is.null(ar)) return(ar)
  if (inherits(ar, "ilm_cor")) {
    if (is.null(ar$n_latent)) ar$n_latent <- ar$n_cell
    return(ar)
  }
  if (!is.list(ar) || is.null(ar$idx))
    stop("`ar` must come from ilm_ar1(), ilm_car1() or ilm_rw1()", call. = FALSE)
  if (is.null(ar$type)) {
    if (is.null(ar$n_group) || is.null(ar$Tt))
      stop("`ar` must come from ilm_ar1(), ilm_car1() or ilm_rw1(); a bare ",
           "list needs `idx`, `n_group` and `Tt`.", call. = FALSE)
    ar$type <- "ar1"
    ar$n_cell <- ar$n_group * ar$Tt
    ar$n_obs <- length(ar$idx)
    ar$obs_per_latent <- ar$n_obs / ar$n_cell
  }
  if (is.null(ar$n_latent)) ar$n_latent <- ar$n_cell
  structure(ar, class = c(paste0("ilm_", ar$type), "ilm_cor"))
}

## What a correlation over time's parameters mean, in words: one place, so
## ilm_cells() and ilm_varcorr() say the same thing.
#' @keywords internal
#' @noRd
ilm_ar_words <- function(ar) {
  switch(ar$type,
    ar1 = paste0("AR(1): rho is the correlation one grid step apart (a step ",
                 "of ", signif(ar$step, 6), " time units), and the process is ",
                 "stationary with covariance Sigma$ar at every cell."),
    car1 = paste0("CAR(1): rho is the correlation one time unit apart, so ",
                  "cells d units apart correlate rho^d; the process is ",
                  "stationary with covariance Sigma$ar at every cell."),
    rw1 = paste0("Random walk: each group's walk is zero at its first cell; ",
                 "a step across d time units adds an independent change with ",
                 "covariance d * Sigma$ar, so Sigma$ar is the variance per unit ",
                 "of time and the variance grows with the time since the ",
                 "group's first cell."))
}

## The cells that are integrated out, in the order the fit stores them: all of
## them, except a random walk's anchors, which are held at zero.
#' @keywords internal
#' @noRd
ilm_ar_free <- function(ar)
  if (identical(ar$type, "rw1")) ar$rest else seq_len(ar$n_cell)

## How many latent values a structure puts under the data, per linear
## predictor. A fit saved before the count was stored has all its cells free.
#' @keywords internal
#' @noRd
ilm_ar_nlat <- function(ar) if (is.null(ar$n_latent)) ar$n_cell else ar$n_latent

#' The cells of a correlation over time
#'
#' One row per latent value of a fitted model's [ilm_ar1()], [ilm_car1()] or
#' [ilm_rw1()] term: which group and time it belongs to, how many observations
#' sit on it, and where the fit holds it. Code that works with the latent
#' values -- forecasts, draws, checks -- reads the layout here rather than
#' rebuilding the index arithmetic, which is how the two come apart.
#'
#' @param object A fitted `"ilm_model"`.
#' @return `NULL` when the model has no correlation over time. Otherwise a
#'   data frame with one row per cell, ordered by group and then time, and
#'   columns
#'   \describe{
#'     \item{`term`}{`"ar"`, the name the term has in `Sigma` and in the
#'       parameter names.}
#'     \item{`type`}{`"ar1"`, `"car1"` or `"rw1"`.}
#'     \item{`group`}{the unit, a factor with the fitted levels.}
#'     \item{`time`}{the cell's time, on the scale it was given: numeric,
#'       `Date` or date-time.}
#'     \item{`index`}{where the fit holds the cell's value: the `index`-th of
#'       the entries named `"B_ar"` in `object$sdr$par.random` and in the
#'       joint precision. `NA` at a random walk's anchors, which are held at
#'       zero rather than estimated. With a multinomial outcome a cell has
#'       one value per linear predictor, stored predictor after predictor, so
#'       predictor `c`'s value is at `index + (c - 1) * attr(, "n_latent")`.}
#'     \item{`n_obs`}{how many fitted rows sit on the cell; 0 at an AR(1) grid
#'       step no one was observed at.}
#'     \item{`anchor`}{`TRUE` at each group's first cell of a random walk.}
#'     \item{`last`}{`TRUE` at each group's final cell, where a forecast
#'       starts. An AR(1) grid runs to the latest time in the data for every
#'       group, so a group seen only early ends in cells without observations,
#'       whose values are the chain's own forecast.}
#'   }
#'   Attributes: `obs_cell`, each fitted row's cell (a row number of this
#'   table) in the order of the model frame; `vars`, the time and group
#'   columns when the term was given by name, `~ time | group`, and `NULL`
#'   when it was given vectors; `n_latent`, the number of estimated cells;
#'   `parameterisation`, what the fitted parameters mean, in words; and for
#'   AR(1), `step` and `origin`, the grid's spacing and first time.
#' @seealso [ilm_ar1()], [ilm_car1()], [ilm_rw1()].
#' @examples
#' set.seed(2)
#' d <- data.frame(id = rep(c("a", "b"), each = 12), day = rep(1:12, 2))
#' d$y <- rnorm(24) + rep(cumsum(rnorm(12, 0, 0.5)), 2)
#' fit <- ilm_model(y ~ 1, data = d, ar = ilm_rw1(~ day | id), verbose = FALSE)
#' cl <- ilm_cells(fit)
#' head(cl)
#' attr(cl, "parameterisation")
#' @export
ilm_cells <- function(object) {
  if (!inherits(object, "ilm_model"))
    stop("`object` must be a fitted ilm_model, not ", class(object)[1],
         call. = FALSE)
  ar <- object$ar
  if (is.null(ar)) return(NULL)
  type <- ar$type
  n_cell <- ar$n_cell
  if (identical(type, "ar1")) {
    Tt <- ar$Tt
    cg <- rep(seq_len(ar$n_group), each = Tt)
    slot <- rep(seq_len(Tt), ar$n_group)
    org <- if (is.null(ar$origin)) NA_real_ else ar$origin
    tt <- org + (slot - 1L) * ar$step
  } else {
    if (is.null(ar$ct))
      stop("this fit was made before the cells' times were stored; refit it ",
           "to read them.", call. = FALSE)
    cg <- ar$cg; tt <- ar$ct
  }
  glev <- if (is.null(ar$glev)) as.character(seq_len(ar$n_group)) else ar$glev
  free <- ilm_ar_free(ar)
  index <- rep(NA_integer_, n_cell); index[free] <- seq_along(free)
  anchor <- rep(FALSE, n_cell)
  if (identical(type, "rw1")) anchor[ar$first] <- TRUE
  last <- c(cg[-1L] != cg[-n_cell], TRUE)
  out <- data.frame(term = "ar", type = type,
                    group = factor(glev[cg], levels = glev),
                    time = ilm_cor_as_time(tt, ar),
                    index = index,
                    n_obs = tabulate(ar$idx, n_cell),
                    anchor = anchor, last = last,
                    stringsAsFactors = FALSE)
  attr(out, "obs_cell") <- ar$idx
  attr(out, "vars") <- ar$vars
  attr(out, "n_latent") <- length(free)
  attr(out, "parameterisation") <- ilm_ar_words(ar)
  if (identical(type, "ar1")) {
    attr(out, "step") <- ar$step
    attr(out, "origin") <- ilm_cor_as_time(org, ar)
  }
  out
}
