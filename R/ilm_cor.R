## ---------------------------------------------------------------------------
## Correlation structures over time.
##
## Both of these describe a latent process shared by observations at the same
## place in time within the same group, and both are written the same way in the
## likelihood: a first-order Markov chain, factorised conditionally.
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
## THE CONSTRAINT WORTH UNDERSTANDING. These are latent-variable structures, not
## marginal ones. A latent value exists for every distinct (group, time) pair, so
## continuous time with no repeats means one latent per observation, which is the
## regime where the Laplace approximation stops being estimable (see the latent
## budget work: 33% convergence at 1.5 observations per latent). nlme's corCAR1
## does not have this problem because it puts the correlation on the residuals
## instead. Both constructors therefore report the budget and say so out loud
## when it is thin -- rounding the time to a coarser grid is usually the fix.
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
  list(t = tnum, g = as.integer(factor(group)), glev = levels(factor(group)))
}

## Shared budget reporting, so both structures say the same thing the same way.
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
#' @param time Time index, one value per observation. Whole numbers, or values
#'   on a common step such as `c(0, 5, 10, 15)`.
#' @param group Unit identifier, one value per observation.
#' @param verbose Warn when the latent budget is thin.
#' @return A `"ilm_ar1"` specification, to pass as `ilm_model(ar = )`.
#' @seealso [ilm_car1()] for arbitrary gaps, [ilm_check_ar()] to test whether
#'   the structure is needed.
#' @examples
#' d <- ilm_sim(n_id = 20, n_period = 8)
#' a <- ilm_ar1(as.integer(factor(d$date)), d$id)
#' a
#' @export
ilm_ar1 <- function(time, group, verbose = TRUE) {
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
                 n_cell = n_cell, step = step, n_obs = length(idx),
                 obs_per_latent = r, glev = z$glev),
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
#' @param time Time, one value per observation. Any numeric scale or `Date`;
#'   `rho` is the correlation one unit apart, so the units matter.
#' @param group Unit identifier, one value per observation.
#' @param verbose Warn when the latent budget is thin.
#' @return A `"ilm_car1"` specification, to pass as `ilm_model(ar = )`.
#' @seealso [ilm_ar1()] for evenly spaced time, [ilm_check_ar()] and
#'   [ilm_plot_acf()] to test whether the structure is needed.
#' @examples
#' set.seed(1)
#' d <- data.frame(id = rep(1:20, each = 6),
#'                 day = as.vector(replicate(20, sort(sample(1:60, 6)))))
#' ilm_car1(d$day, d$id)
#' @export
ilm_car1 <- function(time, group, verbose = TRUE) {
  z <- ilm_cor_time(time, group, "CAR(1)")
  if (length(unique(z$t)) < 2L)
    stop("`time` has only one distinct value, so there is no autocorrelation ",
         "to model.", call. = FALSE)

  ## Latent cells are the distinct (group, time) pairs, ordered by group then
  ## time. That ordering is what lets the preceding cell in a group always be
  ## the one before it, exactly as the evenly spaced case relies on.
  o <- order(z$g, z$t)
  gs <- z$g[o]; ts <- z$t[o]; n <- length(o)
  new <- c(TRUE, gs[-1L] != gs[-n] | ts[-1L] != ts[-n])
  cell_sorted <- cumsum(new)
  idx <- integer(n); idx[o] <- cell_sorted

  cg <- gs[new]; ct <- ts[new]; n_cell <- length(cg)
  first <- which(c(TRUE, cg[-1L] != cg[-n_cell]))
  rest <- if (n_cell > length(first)) setdiff(seq_len(n_cell), first) else integer(0)
  if (!length(rest))
    stop("every unit has only one distinct time, so there is nothing for a ",
         "correlation over time to describe. A plain random intercept, ",
         "(1 | group), is the model you want.", call. = FALSE)
  prev <- rest - 1L
  gap <- ct[rest] - ct[prev]

  ## gaps are positive by construction; a gap far below the typical one still
  ## drives phi to 1 and the innovation variance to 0
  med <- stats::median(gap)
  if (verbose && min(gap) < med * 1e-6)
    warning("the smallest gap between consecutive times is ", signif(min(gap), 3),
            " against a median of ", signif(med, 3),
            ". Nearly coincident times force the correlation to 1 and the ",
            "innovation variance to 0, which the optimiser handles badly. ",
            "Consider rounding `time`.", call. = FALSE)
  r <- ilm_cor_budget(n, n_cell, "CAR(1)", verbose)
  structure(list(type = "car1", idx = idx, n_cell = n_cell,
                 first = first, prev = prev, rest = rest, gap = gap,
                 n_group = length(first), n_obs = n, obs_per_latent = r,
                 time_range = range(z$t), glev = z$glev),
            class = c("ilm_car1", "ilm_cor"))
}

#' Print a correlation structure
#'
#' @param x An `"ilm_cor"` object.
#' @param ... Unused.
#' @return `x`, invisibly.
#' @export
print.ilm_cor <- function(x, ...) {
  lab <- if (identical(x$type, "car1")) "CAR(1), continuous time" else
    "AR(1), evenly spaced time"
  cat(lab, "\n")
  cat(sprintf("  %d group%s, %d latent value%s, %d observation%s (%.2f per latent)\n",
              x$n_group, if (x$n_group == 1L) "" else "s",
              x$n_cell, if (x$n_cell == 1L) "" else "s",
              x$n_obs, if (x$n_obs == 1L) "" else "s", x$obs_per_latent))
  if (identical(x$type, "car1"))
    cat(sprintf("  gaps: min %s, median %s, max %s\n",
                signif(min(x$gap), 4), signif(stats::median(x$gap), 4),
                signif(max(x$gap), 4)))
  else
    cat(sprintf("  step: %s, grid of %d time slots\n", signif(x$step, 4), x$Tt))
  if (x$obs_per_latent < 1.5)
    cat("  the latent budget is thin; see ?ilm_car1\n")
  invisible(x)
}

## Accept both the constructors above and the bare list the fitter took before
## they existed, so nothing written against the old contract breaks.
#' @keywords internal
#' @noRd
ilm_as_cor <- function(ar) {
  if (is.null(ar) || inherits(ar, "ilm_cor")) return(ar)
  if (!is.list(ar) || is.null(ar$idx))
    stop("`ar` must come from ilm_ar1() or ilm_car1()", call. = FALSE)
  if (is.null(ar$type)) {
    if (is.null(ar$n_group) || is.null(ar$Tt))
      stop("`ar` must come from ilm_ar1() or ilm_car1(); a bare list needs ",
           "`idx`, `n_group` and `Tt`.", call. = FALSE)
    ar$type <- "ar1"
    ar$n_cell <- ar$n_group * ar$Tt
    ar$n_obs <- length(ar$idx)
    ar$obs_per_latent <- ar$n_obs / ar$n_cell
  }
  structure(ar, class = c(paste0("ilm_", ar$type), "ilm_cor"))
}
