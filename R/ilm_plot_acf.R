## ---------------------------------------------------------------------------
## Autocorrelation and partial-autocorrelation plots for a fitted model.
##
## What makes these different from acf() and pacf():
##
##   1. Correlations are computed WITHIN group, on model residuals, matching
##      pairs on exact time differences. A panel dataset stacked long is not one
##      series, and treating it as one manufactures correlation at the joins.
##
##   2. The band is not 2/sqrt(n). That band assumes independent observations
##      and an exactly-correct residual distribution, neither of which holds for
##      in-sample residuals from a mixed model. Here the band comes from
##      simulating from the fit and REFITTING, so it is the spread the model
##      itself produces.
##
##   3. The point colour IS the verdict. The band is drawn at the same critical
##      value the p-value is computed against, so a point outside the band is
##      exactly a point the check flags -- the picture and the table cannot
##      disagree.
##
## Why both panels. For an AR(p) process the partial autocorrelation cuts off
## after lag p while the ordinary one decays; for a moving-average process it is
## the other way round. One panel alone cannot tell those apart, and the order
## is what you need to know to decide what to add to the model.
## ---------------------------------------------------------------------------

ILM_STATUS_COL <- c(OK = "grey25", WARN = "#E08214", FAIL = "#D7191C",
                    INCONCLUSIVE = "grey70")

#' Draw one autocorrelation panel
#'
#' @param spec An `"ilm_ar_envelope"` object.
#' @param which `"acf"` or `"pacf"`.
#' @param colour Colour for lags the check does not flag.
#' @param fill Band fill.
#' @param alpha Band transparency.
#' @param size Point and spike scaling.
#' @param main Title.
#' @return Invisibly, the table behind the panel.
#' @keywords internal
#' @noRd
ilm_acf_panel <- function(spec, which = "acf", colour = "grey25",
                          fill = "grey85", alpha = NULL, size = 1,
                          main = NULL) {
  tab <- spec[[which]]
  L <- tab$lag
  est <- tab$estimate; lo <- tab$lo; hi <- tab$hi
  if (!is.null(alpha)) fill <- grDevices::adjustcolor(fill, alpha)
  cols <- ILM_STATUS_COL[tab$status]
  cols[tab$status == "OK"] <- colour
  ylim <- range(c(0, est, lo, hi), na.rm = TRUE)
  if (!all(is.finite(ylim))) ylim <- c(-1, 1)
  pad <- diff(ylim) * 0.18
  if (is.null(main))
    main <- if (identical(which, "acf")) "Residual autocorrelation"
            else "Residual partial autocorrelation"

  ## the annotation under the panel needs room, and mtext will not make it
  op <- par(mar = c(6.4, 4.2, 3.6, 1.2)); on.exit(par(op), add = TRUE)
  plot(L, est, type = "n", xlim = c(0.5, max(L) + 0.5),
       ylim = c(ylim[1] - pad, ylim[2] + pad), xaxt = "n", xlab = "",
       ylab = if (identical(which, "acf")) "correlation" else "partial correlation",
       main = main)
  title(xlab = "lag (time units apart, within group)", line = 2.3)
  axis(1, at = L)
  ## the band is per-lag, not a curve between lags, so it is drawn as one
  ## rectangle per lag rather than a continuous ribbon
  for (i in seq_along(L)) {
    if (!is.finite(lo[i]) || !is.finite(hi[i])) next
    polygon(L[i] + c(-0.38, 0.38, 0.38, -0.38), c(lo[i], lo[i], hi[i], hi[i]),
            col = fill, border = NA)
  }
  abline(h = 0, col = "grey40", lty = 2)
  ok <- is.finite(est)
  if (any(ok)) {
    graphics::segments(L[ok], 0, L[ok], est[ok], col = cols[ok],
                       lwd = 2 * size)
    points(L[ok], est[ok], pch = 16, cex = 1.1 * size, col = cols[ok])
  }
  if (any(!ok))
    points(L[!ok], rep(0, sum(!ok)), pch = 4, cex = 0.9 * size, col = "grey60")

  ad <- ilm_ar_advice(spec, which)
  verdict <- if (all(tab$status == "INCONCLUSIVE"))
    sprintf("INCONCLUSIVE: only %d of %d replicates refitted", spec$n_ok, spec$B)
  else if (is.null(ad))
    "OK: consistent with the fitted model"
  else sprintf("%s: lag %s outside the envelope",
               if (any(tab$status == "FAIL")) "FAIL" else "WARN",
               paste(ad$lags, collapse = ", "))
  mtext(verdict, side = 3, line = 0.2, cex = 0.72,
        col = if (is.null(ad)) "grey30" else ILM_STATUS_COL[[
          if (any(tab$status == "FAIL")) "FAIL" else "WARN"]])
  foot <- if (!is.null(ad$order_short)) ad$order_short
          else if (!is.null(ad)) ad$fix_short
          else sprintf("band: 95%% envelope from %d refits of simulated data",
                       spec$n_ok)
  ilm_mtext_wrap(foot, line = 3.6, cex = 0.62, col = "grey30")
  invisible(tab)
}

## mtext does not wrap, and a truncated sentence of advice is worse than no
## advice, so wrap to the panel width and draw a line at a time.
#' @keywords internal
#' @noRd
ilm_mtext_wrap <- function(x, line = 3.6, cex = 0.62, col = "grey30",
                           max_lines = 3L) {
  ## Character capacity of the plotting region, measured from this very string
  ## rather than from par("cin"): cin is the em width, which overstates the
  ## average character by about a third and wraps the caption short.
  w1 <- graphics::strwidth(x, units = "inches", cex = cex)
  cpi <- if (is.finite(w1) && w1 > 0) nchar(x) / w1 else 12
  wid <- max(20L, floor(par("pin")[1] * 0.96 * cpi))
  ln <- strwrap(x, width = wid)
  if (length(ln) > max_lines)
    ln <- c(ln[seq_len(max_lines - 1L)],
            paste0(substr(ln[max_lines], 1L,
                          max(1L, wid - 3L)), "..."))
  for (i in seq_along(ln))
    mtext(ln[i], side = 1, line = line + (i - 1L) * 0.85, cex = cex, col = col)
  invisible(ln)
}

#' Autocorrelation and partial-autocorrelation plots with verdicts
#'
#' Draws the residual autocorrelation of a fitted model, within group and by
#' lag, against an envelope built by simulating from the fit and refitting.
#' Lags the check flags are coloured, and the panel is annotated with what the
#' pattern means and what to do about it.
#'
#' @section Why not `acf()`:
#' `stats::acf()` treats its input as one series and draws a
#' `plus or minus 2/sqrt(n)` band. Neither is right here. Panel data stacked
#' long is many short series, and pairing across the joins manufactures
#' correlation that is not there; and the `2/sqrt(n)` band assumes independent
#' observations from an exactly-correct model, which in-sample residuals from a
#' mixed model are not. Both problems push in the same direction, toward
#' flagging autocorrelation that does not exist.
#'
#' @section Reading the two panels together:
#' For an autoregressive process of order p the partial autocorrelation cuts off
#' after lag p while the ordinary one decays geometrically; for a moving-average
#' process the ordinary one cuts off and the partial one decays. So the partial
#' panel is the one that tells you the *order*, which is what decides whether an
#' AR(1) term is enough.
#'
#' @section What the colours mean:
#' The band is drawn at the critical value the p-value is computed against, so a
#' point outside the band is exactly a lag the check flags. Flagged lags are
#' coloured by verdict -- amber for `WARN`, red for `FAIL` -- regardless of
#' `colour`, because there the colour *is* the verdict. Lags with too few
#' usable pairs to estimate are marked with a cross on the zero line rather than
#' silently omitted.
#'
#' @param object A fitted `"ilm_model"`, or the value returned by
#'   [ilm_check_ar()], which lets the plot reuse that call's refits instead of
#'   paying for them twice.
#' @param time Integer time index, one per observation. Ignored when `object`
#'   already carries an envelope.
#' @param group Grouping variable, one per observation.
#' @param maxlag Integer. Largest lag to examine.
#' @param B Integer. Simulated datasets behind the envelope. No p-value can
#'   fall below `1/(B+1)`, so roughly 100 is needed before `FAIL` is reachable.
#' @param which `"both"`, `"acf"` or `"pacf"`.
#' @param ncores Integer. Worker processes for the refits.
#' @param seed Integer. Random seed.
#' @param colour Colour for lags that are not flagged.
#' @param fill Band fill colour.
#' @param alpha Band transparency, in the ggplot2 sense.
#' @param size Scaling for points and spikes.
#' @param main Title. With `which = "both"` it is used for the first panel.
#' @param verbose Logical. Also print the per-lag table.
#' @param ... Unused, present so the `plot()` method can pass through.
#' @return Invisibly, the `"ilm_ar_envelope"` object behind the plot, which can
#'   be passed straight back in to redraw without refitting.
#' @seealso [ilm_check_ar()] for the table and the p-values,
#'   [ilm_check_variance()] for the other within-group assumption.
#' @examples
#' set.seed(1)
#' d <- ilm_sim(n_id = 20, n_period = 8)
#' f <- ilm_model(score ~ income + (1 | id), data = d, family = "gaussian",
#'                verbose = FALSE)
#' ilm_plot_acf(f, time = as.integer(factor(d$date)), group = d$id,
#'              maxlag = 3, B = 12)
#' @export
ilm_plot_acf <- function(object, time, group, maxlag = 8L, B = 100L,
                         which = "both", ncores = 1L, seed = 1L,
                         colour = "grey25", fill = "grey85", alpha = NULL,
                         size = 1, main = NULL, verbose = FALSE, ...) {
  if (length(which) != 1L || !which %in% c("both", "acf", "pacf"))
    stop("unknown `which`: ", paste(sQuote(which), collapse = ", "),
         ". Options are 'both', 'acf' and 'pacf'.", call. = FALSE)

  spec <- if (inherits(object, "ilm_ar_envelope")) object
          else if (is.list(object) && inherits(object$envelope, "ilm_ar_envelope"))
            object$envelope
          else if (inherits(object, "ilm_model")) {
            if (1 / (B + 1) > 0.01)
              warning("B = ", B, " puts the smallest achievable p-value at ",
                      signif(1 / (B + 1), 2),
                      ", so a FAIL verdict is unreachable. Use B >= 100.",
                      call. = FALSE)
            ilm_ar_envelope(object, time, group, maxlag, B, ncores, seed)
          }
          else stop("`object` must be a fitted ilm_model object or the result ",
                    "of ilm_check_ar(), not ", class(object)[1], call. = FALSE)

  panels <- if (identical(which, "both")) c("acf", "pacf") else which
  if (length(panels) > 1L) {
    op <- par(mfrow = c(1, 2)); on.exit(par(op), add = TRUE)
  }
  for (i in seq_along(panels))
    ilm_acf_panel(spec, panels[i], colour, fill, alpha, size,
                  if (i == 1L) main else NULL)
  if (verbose) for (w in panels) ilm_ar_report(spec, w)
  invisible(spec)
}
