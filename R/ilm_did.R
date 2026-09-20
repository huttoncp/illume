## ---------------------------------------------------------------------------
## Difference in differences.
##
## The estimator is an interaction in a mixed model, so this is a thin wrapper
## over ilm_model(). What is not thin is everything around it:
##
##  - Parallel trends is the identifying assumption and it is partly checkable.
##    The pre-treatment periods say whether the two groups were moving together
##    before anything happened, and that check belongs in the same object as the
##    estimate rather than in a footnote.
##  - Serial correlation is where difference in differences classically breaks.
##    Bertrand, Duflo and Mullainathan showed that ignoring it makes the
##    standard errors far too small. illume answers that natively with a unit
##    random effect and, where the panel is long enough, AR(1) -- machinery the
##    simulation studies already cover.
##  - Staggered adoption is a trap rather than a generalisation. Under
##    heterogeneous effects, two-way fixed effects uses already-treated units as
##    controls and the result is not an average treatment effect at all. This
##    detects it and says so instead of quietly returning a number.
## ---------------------------------------------------------------------------

#' @keywords internal
#' @noRd
ilm_did_binary <- function(v, nm) {
  if (is.logical(v)) return(as.integer(v))
  ## A FACTOR carries its own order and that order is the author's intent:
  ## factor(x, levels = c("before", "after")) means before comes first. Sorting
  ## the labels alphabetically instead puts "after" first and silently inverts
  ## the indicator, which turns a difference in differences upside down.
  if (is.factor(v)) {
    u <- levels(droplevels(v))
    if (length(u) != 2L)
      stop("`", nm, "` must have exactly 2 levels to mark a group; it has ",
           length(u), call. = FALSE)
    return(as.integer(as.character(v) == u[2]))
  }
  if (is.character(v)) {
    u <- sort(unique(stats::na.omit(v)))
    if (length(u) != 2L)
      stop("`", nm, "` must have exactly 2 values to mark a group; it has ",
           length(u), call. = FALSE)
    return(as.integer(v == u[2]))
  }
  u <- sort(unique(stats::na.omit(v)))
  if (length(u) != 2L || !all(u %in% c(0, 1)))
    stop("`", nm, "` must be 0/1, TRUE/FALSE or a two-level factor; it has ",
         length(u), " distinct value(s)", call. = FALSE)
  as.integer(v)
}

#' Difference in differences
#'
#' Estimates the effect of a treatment that switches on for some units at some
#' time, by comparing the change in the treated group with the change in the
#' control group. The estimate is the interaction in a mixed model, so the fit
#' is an ordinary [ilm_model()] and every diagnostic and method applies to it.
#'
#' @section What identifies the estimate:
#'
#' Difference in differences rests on **parallel trends**: without the
#' treatment, the two groups' outcomes would have moved together. That cannot be
#' checked where it matters, because it concerns a world in which the treatment
#' did not happen, but it can be checked before the treatment, and that is worth
#' doing. Two things here do it:
#'
#' * a formal test of whether the groups' pre-treatment slopes differ, reported
#'   as `parallel`;
#' * the event study, `$event`, one coefficient per period relative to
#'   treatment. Pre-treatment coefficients should sit near zero. This is the
#'   more informative of the two, because it shows *how* a violation looks, and
#'   [ilm_plot_did()] draws it.
#'
#' @section Serial correlation:
#'
#' Outcomes from the same unit over time are correlated, and difference in
#' differences with naive standard errors over-rejects badly as a result
#' (Bertrand, Duflo and Mullainathan 2004). A unit random intercept is fitted by
#' default. With enough periods, `ar = TRUE` adds an AR(1) process on top, which
#' is the part that matters when the within-unit correlation decays rather than
#' being constant.
#'
#' @section Staggered adoption:
#'
#' When units are treated at different times, the two-way fixed effects
#' regression is **not** an average treatment effect if the effect varies across
#' units or over time: already-treated units end up serving as controls for
#' later-treated ones, and those comparisons can enter with negative weight.
#' This function detects staggered timing and refuses to report a single
#' pooled number for it unless `allow_staggered = TRUE`, in which case the
#' result carries the caveat. The honest remedy is a cohort-based estimator --
#' the `did` and `fixest` packages implement them -- and illume does not yet
#' have one.
#'
#' @section What this is calibrated for:
#'
#' On 800 replicates of a 40-unit, 8-period panel with unit random effects and
#' a true effect of 0.800: the estimate was unbiased (+0.0029), its interval
#' covered 0.946 of the time, and the parallel-trends check raised a false alarm
#' 0.054 of the time under trends that really were parallel, with null p-values
#' uniform by Kolmogorov-Smirnov (p = 0.111) rather than merely correct at the
#' 5% threshold.
#'
#' Getting that check calibrated took choosing the right reference. The
#' comparison is between two groups of UNITS' pre-treatment slopes, so units
#' are the independent replicates and the large-sample normal is too generous
#' when there are few: it rejected 0.060, 0.068 and 0.050 of the time at 40, 20
#' and 80 units, where t on `units - 2` gave 0.048, 0.055 and 0.045. A check
#' that cries wolf is worse than no check here, because the remedy for a failed
#' parallel-trends test is to abandon the design.
#'
#' @param data A data frame, one row per unit-period.
#' @param y Outcome column.
#' @param unit Unit identifier.
#' @param time Period. Numeric or anything that orders.
#' @param treated Column marking units ever treated. Omit when `treatment` is
#'   given.
#' @param post Column marking periods after treatment starts. Omit when
#'   `treat_time` or `treatment` is given.
#' @param treat_time First treated period, when `post` is not supplied.
#' @param treatment Per-row treatment indicator, as an alternative to
#'   `treated` and `post`. Staggered timing is detected from it.
#' @param covariates Further columns for the mean structure. These must be
#'   things the treatment cannot have affected; a covariate the treatment
#'   changes is a mediator and adjusting for it removes part of the effect.
#' @param family Response distribution; inferred from `y` when `NULL`.
#' @param ar Add AR(1) within unit on top of the random intercept.
#' @param allow_staggered Report a pooled estimate even when adoption is
#'   staggered. Off by default, and see the section above.
#' @param verbose Narrate each step.
#' @param ... Passed to [ilm_model()].
#' @return An object of class `"ilm_did"`: `att`, `fit`, `parallel`, `event`
#'   and the settings used.
#' @references
#' Bertrand, M., Duflo, E. and Mullainathan, S. (2004). How much should we
#' trust differences-in-differences estimates? Quarterly Journal of Economics
#' 119(1).
#'
#' Goodman-Bacon, A. (2021). Difference-in-differences with variation in
#' treatment timing. Journal of Econometrics 225(2).
#' @seealso [ilm_plot_did()], [ilm_rdd()], [ilm_model()].
#' @examples
#' set.seed(1)
#' d <- expand.grid(unit = 1:40, time = 1:8)
#' d$treated <- as.integer(d$unit <= 20)
#' d$post <- as.integer(d$time >= 5)
#' d$y <- 1 + 0.3 * d$time + rnorm(40)[d$unit] +
#'        0.8 * d$treated * d$post + rnorm(nrow(d))
#' ilm_did(d, "y", "unit", "time", treated = "treated", post = "post",
#'         verbose = FALSE)
#' @export
ilm_did <- function(data, y, unit, time, treated = NULL, post = NULL,
                    treat_time = NULL, treatment = NULL, covariates = NULL,
                    family = NULL, ar = FALSE, allow_staggered = FALSE,
                    verbose = TRUE, ...) {
  if (!is.data.frame(data))
    stop("`data` must be a data frame; it is ", class(data)[1], call. = FALSE)
  need <- c(y, unit, time, treated, post, treatment, covariates)
  miss <- setdiff(need, names(data))
  if (length(miss))
    stop("column(s) not found in the data: ", paste(miss, collapse = ", "),
         ". Available: ", paste(utils::head(names(data), 12), collapse = ", "),
         call. = FALSE)
  say <- function(...) if (verbose) message(...)
  d <- data
  d[[unit]] <- factor(d[[unit]])
  tv <- d[[time]]
  tnum <- if (is.numeric(tv)) tv else as.integer(factor(tv))
  d$.time_num <- tnum

  say("== difference in differences ==")
  say("[1/5] design")

  ## ---- who is treated, and from when --------------------------------------
  staggered <- FALSE; first_treat <- NULL
  if (!is.null(treatment)) {
    tr <- ilm_did_binary(d[[treatment]], treatment)
    d$.treat <- tr
    ever <- tapply(tr, d[[unit]], function(v) as.integer(any(v == 1L)))
    d$.treated <- as.integer(ever[as.character(d[[unit]])])
    ft <- tapply(seq_len(nrow(d)), d[[unit]], function(i) {
      k <- i[tr[i] == 1L]
      if (!length(k)) NA_real_ else min(tnum[k])
    })
    ## tapply gives a 1-d ARRAY, and na.omit on one of those takes the matrix
    ## path and fails; as.vector is what makes it an ordinary vector
    first_treat <- stats::setNames(as.vector(ft), names(ft))
    ftt <- as.vector(stats::na.omit(first_treat))
    staggered <- length(unique(ftt)) > 1L
    tt0 <- if (length(ftt)) min(ftt) else NA_real_
    d$.post <- as.integer(d$.time_num >= tt0)
  } else {
    if (is.null(treated))
      stop("give either `treatment` (a per-row indicator) or `treated` ",
           "together with `post` or `treat_time`", call. = FALSE)
    d$.treated <- ilm_did_binary(d[[treated]], treated)
    if (!is.null(post)) d$.post <- ilm_did_binary(d[[post]], post)
    else if (!is.null(treat_time)) d$.post <- as.integer(d$.time_num >= treat_time)
    else stop("give `post` or `treat_time` alongside `treated`", call. = FALSE)
    tt0 <- suppressWarnings(min(d$.time_num[d$.post == 1L]))
    d$.treat <- d$.treated * d$.post
  }

  nper <- length(unique(d$.time_num))
  ntr <- length(unique(d[[unit]][d$.treated == 1L]))
  nco <- length(unique(d[[unit]][d$.treated == 0L]))
  say("  ", nlevels(d[[unit]]), " units over ", nper, " periods: ",
      ntr, " treated, ", nco, " control")
  say("  treatment starts at ", time, " = ", tt0)
  if (!ntr || !nco)
    stop("difference in differences needs both treated and control units; ",
         "found ", ntr, " treated and ", nco, " control", call. = FALSE)
  npre <- sum(unique(d$.time_num) < tt0)
  say("  ", npre, " pre-treatment period", if (npre == 1L) "" else "s",
      if (npre < 2L) " -- too few to check parallel trends" else "")

  if (staggered) {
    say("  STAGGERED ADOPTION: units are treated at ",
        length(unique(ftt)), " different times.")
    say("    Under heterogeneous effects a pooled two-way fixed effects ",
        "estimate is not an average treatment effect -- already-treated units ",
        "act as controls for later-treated ones. See the `did` or `fixest` ",
        "packages for a cohort-based estimator; illume does not have one yet.")
    if (!allow_staggered)
      stop("adoption is staggered, so a single pooled estimate would be ",
           "misleading. Pass allow_staggered = TRUE to get it anyway with the ",
           "caveat attached, restrict the data to one adoption cohort, or use ",
           "a cohort-based estimator.", call. = FALSE)
  }

  ## ---- family --------------------------------------------------------------
  say("[2/5] response")
  fam <- if (is.null(family)) {
    f <- ilm_dag_family_guess(d[[y]], y)
    say("  `", y, "` looks ", ilm_var_kind(d[[y]]), " -> family \"", f, "\"")
    f
  } else { say("  family \"", family, "\" as supplied"); family }

  ## ---- the difference in differences fit -----------------------------------
  say("[3/5] estimate")
  arspec <- NULL
  if (isTRUE(ar)) {
    arspec <- tryCatch(ilm_ar1(d$.time_num, d[[unit]], verbose = FALSE),
                       error = function(e) NULL)
    if (is.null(arspec)) say("  AR(1) could not be set up; using the random ",
                             "intercept alone")
    else say("  AR(1) within unit, on top of the random intercept")
  }
  ## period effects absorb whatever moved everyone at once; the unit random
  ## intercept absorbs the level differences the treated and control groups
  ## carry, which is what makes this a difference of differences
  rhs <- c(".treated", if (nper > 2L) "factor(.time_num)" else ".post",
           ".treat", covariates, paste0("(1 | ", unit, ")"))
  form <- stats::reformulate(rhs, response = y)
  environment(form) <- environment()
  fit <- ilm_model(form, data = d, family = fam, ar = arspec,
                   verbose = FALSE, ...)
  k <- ilm_term_cols(fit, ".treat")
  ct <- ilm_coef_table(fit)
  crit <- if (isTRUE(fit$exact_df)) stats::qt(0.975, fit$resid_df) else stats::qnorm(0.975)
  att <- data.frame(term = "ATT", estimate = ct[k, 1], se = ct[k, 2],
                    lower = ct[k, 1] - crit * ct[k, 2],
                    upper = ct[k, 1] + crit * ct[k, 2],
                    p_value = ct[k, 4], stringsAsFactors = FALSE)
  rownames(att) <- NULL
  say(sprintf("  ATT = %.4f  (%.4f, %.4f)  p = %s", att$estimate, att$lower,
              att$upper, format.pval(att$p_value, digits = 3, eps = 1e-4)))

  ## ---- parallel trends, before the treatment -------------------------------
  say("[4/5] parallel trends, before treatment")
  par_res <- NULL
  if (npre >= 2L) {
    pre <- d[d$.time_num < tt0, , drop = FALSE]
    pf <- stats::reformulate(c(".treated", ".time_num", ".treated:.time_num",
                               covariates, paste0("(1 | ", unit, ")")),
                             response = y)
    environment(pf) <- environment()
    pfit <- tryCatch(suppressWarnings(
              ilm_model(pf, data = pre, family = fam, verbose = FALSE, ...)),
            error = function(e) NULL)
    if (!is.null(pfit)) {
      pk <- ilm_term_cols(pfit, ".treated:.time_num")
      if (length(pk)) {
        pct <- ilm_coef_table(pfit)
        ## The comparison is between two groups of UNITS' slopes, so the units
        ## are the independent replicates and the large-sample normal is the
        ## wrong reference when there are few of them. Measured under true
        ## parallel trends at 400 replicates: the normal rejected 0.060, 0.068
        ## and 0.050 of the time at 40, 20 and 80 units against a nominal 0.05,
        ## where t on (units - 2) gave 0.048, 0.055 and 0.045. A check that
        ## cries wolf is worse than no check, since the remedy for a failed
        ## parallel-trends test is to abandon the design.
        nu <- nlevels(factor(pre[[unit]]))
        pv <- 2 * stats::pt(-abs(pct[pk, 1] / pct[pk, 2]), max(nu - 2L, 1L))
        par_res <- list(diff_slope = pct[pk, 1], se = pct[pk, 2], p_value = pv,
                        df = max(nu - 2L, 1L),
                        status = if (is.na(pv)) "UNTESTED"
                                 else if (pv < 0.05) "FAIL" else "OK",
                        n_pre = npre)
        say(sprintf("  difference in pre-treatment slope: %.4f (se %.4f), p = %s -- %s",
                    par_res$diff_slope, par_res$se,
                    format.pval(pv, digits = 3, eps = 1e-4), par_res$status))
        if (identical(par_res$status, "FAIL"))
          say("    The groups were already diverging before the treatment, so ",
              "the difference after it is not all treatment. Look at the event ",
              "study; consider a group-specific trend, or a different control ",
              "group.")
      }
    }
  } else {
    par_res <- list(diff_slope = NA_real_, se = NA_real_, p_value = NA_real_,
                    status = "UNTESTED", n_pre = npre)
    say("  not testable: ", npre, " pre-treatment period(s). Parallel trends ",
        "is assumed rather than checked.")
  }

  ## ---- event study ---------------------------------------------------------
  say("[5/5] event study")
  ev <- NULL
  if (nper > 2L) {
    dd <- d
    dd$.rel <- dd$.time_num - tt0
    rl <- sort(unique(dd$.rel))
    ## the reference must be a period that exists, or relevel() errors; -1 is
    ## the convention, the last pre-treatment period otherwise
    neg <- rl[rl < 0]
    ref <- if (-1 %in% rl) -1 else if (length(neg)) max(neg) else min(rl)
    dd$.relf <- stats::relevel(factor(dd$.rel), ref = as.character(ref))
    ef <- stats::reformulate(c("factor(.time_num)", ".treated",
                               ".treated:.relf", covariates,
                               paste0("(1 | ", unit, ")")), response = y)
    environment(ef) <- environment()
    efit <- tryCatch(suppressWarnings(
              ilm_model(ef, data = dd, family = fam, verbose = FALSE, ...)),
            error = function(e) NULL)
    if (!is.null(efit)) {
      ect <- ilm_coef_table(efit)
      kk <- grep("\\.treated:\\.relf", rownames(ect))
      if (length(kk)) {
        lv <- sub(".*\\.relf", "", rownames(ect)[kk])
        ecrit <- if (isTRUE(efit$exact_df)) stats::qt(0.975, efit$resid_df)
                 else stats::qnorm(0.975)
        ev <- data.frame(rel_time = as.numeric(lv),
                         estimate = ect[kk, 1], se = ect[kk, 2],
                         lower = ect[kk, 1] - ecrit * ect[kk, 2],
                         upper = ect[kk, 1] + ecrit * ect[kk, 2],
                         stringsAsFactors = FALSE)
        ev <- rbind(ev, data.frame(rel_time = ref, estimate = 0, se = 0,
                                   lower = 0, upper = 0))
        ev <- ev[order(ev$rel_time), ]; rownames(ev) <- NULL
        npre_bad <- sum(ev$rel_time < 0 & (ev$lower > 0 | ev$upper < 0))
        say("  ", nrow(ev), " periods relative to treatment (reference ", ref,
            "); ", npre_bad, " pre-treatment coefficient(s) exclude zero")
      }
    }
    if (is.null(ev)) say("  could not be fitted")
  } else say("  needs more than 2 periods")

  structure(list(att = att, fit = fit, parallel = par_res, event = ev,
                 y = y, unit = unit, time = time, family = fam,
                 treat_time = tt0, n_pre = npre, n_treated = ntr,
                 n_control = nco, staggered = staggered,
                 first_treat = first_treat, ar = !is.null(arspec),
                 data = d),
            class = "ilm_did")
}

#' @export
print.ilm_did <- function(x, ...) {
  cat("<ilm_did>", x$y, "~ treatment  |", x$n_treated, "treated,",
      x$n_control, "control units\n")
  cat("  family:", x$family, "  random intercept:", x$unit,
      if (x$ar) "  AR(1)" else "", "\n")
  if (isTRUE(x$staggered))
    cat("\n  CAVEAT: adoption is staggered. This pooled estimate is not a\n",
        "  clean average treatment effect when the effect varies across\n",
        "  units or over time.\n", sep = "")
  cat(sprintf("\n  ATT  %8.4f  (%7.4f, %7.4f)  p = %s\n",
              x$att$estimate, x$att$lower, x$att$upper,
              format.pval(x$att$p_value, digits = 3, eps = 1e-4)))
  if (!is.null(x$parallel)) {
    cat("\n  parallel trends before treatment:", x$parallel$status)
    if (!is.na(x$parallel$p_value))
      cat(sprintf("  (slope difference %.4f, p = %s)", x$parallel$diff_slope,
                  format.pval(x$parallel$p_value, digits = 3, eps = 1e-4)))
    cat("\n")
    if (identical(x$parallel$status, "FAIL"))
      cat("    the groups were already diverging; the estimate above absorbs that\n")
    if (identical(x$parallel$status, "UNTESTED"))
      cat("    ", x$parallel$n_pre, " pre-treatment period(s): assumed, not checked\n", sep = "")
  }
  if (!is.null(x$event)) {
    bad <- sum(x$event$rel_time < 0 & (x$event$lower > 0 | x$event$upper < 0))
    cat("  event study:", nrow(x$event), "periods,", bad,
        "pre-treatment coefficient(s) excluding zero  (ilm_plot_did)\n")
  }
  invisible(x)
}

#' Event-study plot for a difference in differences
#'
#' One coefficient per period relative to treatment, with intervals. The
#' pre-treatment coefficients are the diagnostic: under parallel trends they sit
#' near zero, and the shape of any departure says more than a single test does.
#'
#' @param x An [ilm_did()] result.
#' @param main,xlab,ylab Labels.
#' @param ... Passed to [graphics::plot()].
#' @return `x`, invisibly.
#' @seealso [ilm_did()].
#' @examples
#' set.seed(1)
#' d <- expand.grid(unit = 1:40, time = 1:8)
#' d$treated <- as.integer(d$unit <= 20)
#' d$post <- as.integer(d$time >= 5)
#' d$y <- 1 + 0.3 * d$time + rnorm(40)[d$unit] +
#'        0.8 * d$treated * d$post + rnorm(nrow(d))
#' fit <- ilm_did(d, "y", "unit", "time", treated = "treated", post = "post",
#'                verbose = FALSE)
#' ilm_plot_did(fit)
#' @export
ilm_plot_did <- function(x, main = "Event study", xlab = "periods from treatment",
                         ylab = "effect", ...) {
  if (!inherits(x, "ilm_did"))
    stop("`x` must be an ilm_did result, not ", class(x)[1], call. = FALSE)
  ev <- x$event
  if (is.null(ev))
    stop("no event study to plot: it needs more than two periods", call. = FALSE)
  yl <- range(c(ev$lower, ev$upper, 0), finite = TRUE)
  plot(ev$rel_time, ev$estimate, type = "n", ylim = yl, main = main,
       xlab = xlab, ylab = ylab, ...)
  ## the pre-treatment side is the part being judged, so it is marked out
  graphics::rect(min(ev$rel_time) - 1, yl[1] - 1, -0.5, yl[2] + 1,
                 col = grDevices::adjustcolor("grey85", 0.45), border = NA)
  graphics::abline(h = 0, lty = 2, col = "grey40")
  graphics::abline(v = -0.5, lty = 3, col = "grey40")
  ## the reference period is pinned at zero with no width, and an arrow of
  ## zero length warns rather than drawing
  wide <- ev$upper > ev$lower
  if (any(wide))
    graphics::arrows(ev$rel_time[wide], ev$lower[wide], ev$rel_time[wide],
                     ev$upper[wide], code = 3, angle = 90, length = 0.04,
                     col = "grey30")
  graphics::points(ev$rel_time, ev$estimate, pch = 19,
                   col = ifelse(ev$rel_time < 0, "grey30", "black"))
  graphics::mtext("shaded: before treatment, where the coefficients should be 0",
                  side = 3, line = 0.2, cex = 0.75, col = "grey30")
  invisible(x)
}
