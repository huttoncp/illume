## ---------------------------------------------------------------------------
## Seeing an anomaly scan, rather than only reading its table.
##
## ilm_anomaly() returns a ranked table, and a table cannot say the one thing
## that decides what happens next: whether the flagged rows STAND APART from
## the rest, or are only the top few percent of a smooth continuum that the
## threshold happened to cut. Five genuine outliers and the upper tail of a
## heavy-tailed column can produce the same list. So four views, one function:
##
##   "scores"   sorted score against rank, against the reference the scan
##              simulated. The view that answers the question above.
##   "drivers"  which column drives the flags. One column driving nearly all
##              of them is a problem IN that column, not a multivariate one.
##   "map"      the rows on the first two dimensions of the data, so it is
##              visible whether the flagged ones are alike.
##   "row"      one row, column by column: how unusual each value is on its
##              own, beside how far it sits off the shared structure.
##
## THE REFERENCE BAND RESTARTS AT THE LINE. The band is the per-rank null the
## scan kept (`null_curve`): where the k-th highest score falls in data with
## no anomalies at all. Drawn as it stands, that is the wrong reference for
## the rows after the flagged ones. Setting m rows aside moves every remaining
## row up m ranks, so a perfectly clean remainder sits ABOVE the band for a
## long stretch, and five genuine outliers look like a continuum -- the one
## picture this plot exists to rule out. Over 60 scans of 400 rows each, the
## share of the ranks after the line above the band as simulated averaged
## 0.69 with 2% planted anomalies and 0.65 with 5%: not far short of the 0.82
## to 1.00 of noise that genuinely has heavy tails.
##
## Right of the line the band is therefore the null for the unflagged rows ON
## THEIR OWN: their j-th highest score against the null at the rank that
## matches it in quantile, j * n / (n - m). Left of the line it is the null's
## own top ranks, which is what the flagged rows are being set against.
##
## THE VERDICT is a reading of that band, not a test. The first thing tried --
## "a continuum when half the ranks after the line are above it" -- called
## data with NO anomalies a continuum in 21 scans of 100, because the band
## inherits the reference's slightly light upper tail (?ilm_anomaly reports
## the same thing as the raw p running hot). So nothing is read off the band
## unless something was flagged, and the reading has three outcomes rather
## than two. The measured rates are in the help page.
##
## References:
##   Atkinson, A. C. (1981). Two graphical displays for outlying and
##     influential observations in regression. Biometrika 68, 13-20.
## ---------------------------------------------------------------------------

#' See an anomaly scan
#'
#' Four views of an [ilm_anomaly()] result. The default, `"scores"`, is the one
#' that decides what to do next: whether the flagged rows stand apart from the
#' rest, or are only the top few percent of a smooth continuum that the
#' threshold happened to cut. A ranked table cannot tell those apart; the
#' shape of the sorted scores can.
#'
#' @section The four views:
#'
#' * `"scores"` -- every row's score against its rank, highest first, with a
#'   line at the number flagged. For the default reconstruction method a band
#'   shows where 95% of the datasets the scan simulated, which have the same
#'   structure and no anomalies, put the score at each rank. Genuine outliers
#'   are a few points well above the band, after which the curve drops back
#'   into it. A continuum has no break: the rows after the line keep scoring
#'   above what clean data produce. When that happens the plot says so, and
#'   the choice is to say plainly that you are reporting a fixed share of the
#'   rows, or to treat the tail as structure rather than as outliers --
#'   [ilm_profile()] looks for the groups in it.
#' * `"drivers"` -- how many flagged rows each column drives, with a tick at
#'   the count the drivers of the unflagged rows would predict. When one
#'   column drives nearly every flag the problem is usually in that column --
#'   a unit, a sentinel code, a misplaced decimal -- rather than in how the
#'   columns combine. Look at it on its own with [ilm_outliers()] or
#'   [ilm_describe()], and recode sentinels with [ilm_recode_errors()].
#' * `"map"` -- every scored row on the first two dimensions of [ilm_reduce()],
#'   in grey, with the flagged rows in red and sized by score. Flagged rows
#'   that sit together may be one problem repeated or one subpopulation, and
#'   `ilm_profile(x)` takes the flagged rows and describes them. A
#'   reconstruction anomaly is odd *off* the main dimensions, so it need not
#'   sit at the edge of this map: the view shows whether the flagged rows are
#'   alike, not whether they are extreme. Needs the PCAmixdata package, as
#'   [ilm_reduce()] does.
#' * `"row"` -- one row, column by column: its z-score, which is how unusual
#'   each value is in its own column, beside its residual, which is how far the
#'   value sits from where the structure shared by the other rows puts it. The
#'   residual is divided by that column's typical residual (its median absolute
#'   deviation across rows), so that 3 means unusual on both bars; unscaled, a
#'   residual of 0.8 is large in a column the structure explains well and
#'   ordinary in one it barely explains. A small |z| beside a large residual is
#'   a row that is odd only as a combination -- the case the reconstruction
#'   exists for, and the one a column-at-a-time scan cannot see. The columns
#'   are ordered by that scaled residual, so the top bar need not be the
#'   `driver` in the table, which is the column contributing most to the score.
#'
#' @section Why the band restarts at the line:
#'
#' The band the scan simulates is for data with no anomalies at all. Set the
#' flagged rows aside and every remaining row moves up by that many ranks, so
#' against the band as it stands a perfectly clean remainder scores too high
#' for a long stretch and genuine outliers look like a continuum. Over 60
#' scans of 400 rows each, the share of the ranks after the line sitting above
#' the band as simulated averaged 0.69 with 2% planted anomalies and 0.65 with
#' 5% -- not far short of the 0.82 to 1.00 of noise that genuinely has heavy
#' tails.
#'
#' So right of the line the band is redrawn for the unflagged rows on their
#' own: their j-th highest score against the reference at the rank that
#' matches it in quantile, `j * n / (n - m)` for `m` rows flagged out of `n`.
#' The band steps up at the line because of this, and the step is the point:
#' the remainder should begin where a clean dataset's highest scores do.
#'
#' @section How far to trust the verdict:
#'
#' It is a reading of the band, not a test. Over the first `max(20, 2m)`
#' ranks after the line it takes the share still scoring above the band:
#' under a quarter, the flagged rows stand clear; three quarters or more,
#' there is no break, and the message names the remedy; in between, the break
#' is not clean and the plot says only that.
#'
#' Nothing is read off the band when nothing is flagged. The band inherits
#' the reference's slightly light upper tail -- [ilm_anomaly()] reports the
#' same thing as the raw p-value running hot -- and in data with no anomalies
#' at all the top twenty scores sat at least half above it in 21 scans of 100.
#'
#' Measured on fresh scans of 400 rows by 8 columns with a rank-2 structure,
#' 100 per condition, counting the scans that flagged anything:
#'
#' ```
#'                              flagged   stand   no clean     no
#'                              anything  clear    break      break
#'   2% planted anomalies          100      89       11          0
#'   5% planted anomalies          100      97        3          0
#'   a single planted anomaly       90      65       17          8
#'   noise with t tails, 10 df      65      14       19         32
#'   noise with t tails, 5 df       98       7       26         65
#'   noise with t tails, 3 df      100      11       31         58
#' ```
#'
#' Planted anomalies are almost never called a continuum; heavy tails are
#' called one in half to two thirds of the scans that flag anything, and said
#' to stand clear in 7 to 22% of them. Size matters: at 200 rows by 6 columns the t tails were called
#' a continuum in 21 of 53 and 16 of 60 scans, and at 1000 rows by 12 in 58
#' and 52 of 60, with planted anomalies called one in at most 4 of 60 at
#' either size. When the verdict and the picture seem to disagree, trust the
#' picture and look at how far above the band the points after the line sit.
#'
#' @section Isolation forests:
#'
#' `method = "iforest"` has no null behind its score, so `"scores"` draws no
#' band and labels the line as what it is: the share of rows you chose to call
#' anomalous, not a test. `"row"` needs the residuals from a fitted structure
#' and stops. `"drivers"` and `"map"` work for both methods.
#'
#' @param x An [ilm_anomaly()] result.
#' @param type Which view: `"scores"` (the default), `"drivers"`, `"map"` or
#'   `"row"`.
#' @param row For `type = "row"`, the row to show, numbered as in the scanned
#'   data -- the `row` column of the result. Defaults to the highest-scoring
#'   row.
#' @param top_n How much to show: ranks for `"scores"` (by default the larger
#'   of 40 and three times the number flagged), columns for `"drivers"` (15)
#'   and `"row"` (20). Ignored by `"map"`.
#' @param data The data frame that was scanned, for `"map"` and `"row"`. Only
#'   needed if the result did not keep it -- see `keep_data` in
#'   [ilm_anomaly()].
#' @param main Plot title.
#' @param ... Passed to [tinyplot::tinyplot()].
#' @param pch Plotting character for the points in `"scores"` and `"map"`.
#'   Takes a name as well as a number: `"filled circle"` is 16, and every code
#'   from 0 to 25 has one.
#' @return `NULL`, invisibly. The verdict is written on the plot, and where it
#'   calls for a remedy it is also given as a message.
#' @seealso [ilm_anomaly()], [ilm_anomalous()] for the flagged rows as a data
#'   frame, [ilm_profile()], [ilm_outliers()].
#' @references Atkinson, A. C. (1981). Two graphical displays for outlying and
#'   influential observations in regression. *Biometrika* 68, 13-20.
#' @examples
#' set.seed(1)
#' n <- 300
#' f <- rnorm(n)
#' d <- data.frame(north = f + rnorm(n, 0, .3), central = 2 * f + rnorm(n, 0, .3),
#'                 south = -f + rnorm(n, 0, .3))
#' ## three rows that are ordinary in every column but break the pattern
#' d[1:3, ] <- rbind(c(1.2, -2.4, 1.2), c(-1, 2, 1), c(0.8, 1.6, 0.8))
#' a <- ilm_anomaly(d, progress = FALSE)
#' ilm_plot_anomaly(a)                  # do the flagged rows stand apart?
#' ilm_plot_anomaly(a, "row")           # the top row, column by column
#' ilm_plot_anomaly(a, "drivers")
#' \donttest{
#' if (requireNamespace("PCAmixdata", quietly = TRUE))
#'   ilm_plot_anomaly(a, "map")
#' }
#' @export
ilm_plot_anomaly <- function(x, type = c("scores", "drivers", "map", "row"),
                             row = NULL, top_n = NULL, data = NULL,
                             main = NULL, ..., pch = NULL) {
  if (!inherits(x, "ilm_anomaly"))
    stop("`x` must be an ilm_anomaly() result; it is ", class(x)[1],
         call. = FALSE)
  type <- match.arg(type)
  if (!is.null(top_n) && (!is.numeric(top_n) || length(top_n) != 1L ||
                          !is.finite(top_n) || top_n < 1))
    stop("`top_n` must be a single number of at least 1", call. = FALSE)
  d <- as.data.frame(x)
  class(d) <- "data.frame"
  if (!nrow(d)) stop("the result has no rows to plot", call. = FALSE)
  d <- d[order(-d$score), , drop = FALSE]
  row.names(d) <- NULL
  pch <- ilm_pch(if (is.null(pch)) 16L else pch)
  switch(type,
         scores  = ilm_plot_anom_scores(x, d, top_n, main, pch, ...),
         drivers = ilm_plot_anom_drivers(x, d, top_n, main, ...),
         map     = ilm_plot_anom_map(x, d, data, main, pch, ...),
         row     = ilm_plot_anom_row(x, d, row, top_n, data, main, ...))
  invisible(NULL)
}

## ---- "scores" ---------------------------------------------------------------

## The reference at fractional ranks. Left of the line (and everywhere when
## nothing is flagged) it is the null as simulated; right of it, the null for
## the n - nf unflagged rows on their own, matched in quantile. rule = 2 holds
## the ends, so a rank above the first reads as the first.
#' @keywords internal
#' @noRd
ilm_anom_band <- function(nq, nf, x, shifted) {
  n <- nrow(nq)
  r <- if (shifted) (x - nf) * n / (n - nf) else x
  out <- vapply(c("lo", "mid", "hi"), function(k)
    stats::approx(seq_len(n), nq[, k], xout = r, rule = 2)$y,
    numeric(length(x)))
  if (length(x) == 1L)
    out <- matrix(out, 1L, dimnames = list(NULL, c("lo", "mid", "hi")))
  out
}

## Stand apart, or a continuum? Read over the first max(20, 2 * nf) ranks
## after the line, against the restarted band, as the share of them that
## still score above it: under a quarter, the flagged rows stand clear; three
## quarters or more, there is no break; between, the break is not clean.
##
## Only read when something was flagged. The band inherits the reference's
## slightly light upper tail -- the same thing ?ilm_anomaly reports as the
## raw p running hot -- and in data with no anomalies at all the top twenty
## scores sat half above it in 21 scans of 100. When nothing is flagged
## there is no line to ask about, so nothing is read off the band.
#' @keywords internal
#' @noRd
ilm_anom_shape <- function(score, nq, nf) {
  s <- sort(score, decreasing = TRUE)
  n <- length(s); m <- n - nf
  none <- list(shape = NA_character_, above = NA_integer_, J = 0L,
               share = NA_real_)
  if (is.null(nq) || nrow(nq) != n || nf < 1L || m < 1L) return(none)
  J <- as.integer(min(max(20L, 2L * nf), m))
  i <- nf + seq_len(J)
  hi <- ilm_anom_band(nq, nf, i, shifted = TRUE)[, "hi"]
  above <- sum(s[i] > hi)
  share <- above / J
  list(shape = if (share >= 0.75) "continuum" else if (share < 0.25)
         "separate" else "unclear",
       above = as.integer(above), J = J, share = share)
}

## What the "scores" view says, on the plot and -- where there is a remedy to
## name -- as a message. The subtitles are kept to about 65 characters, which
## is what fits a panel of a 2 x 2 layout; the message carries the rest.
#' @keywords internal
#' @noRd
ilm_anom_scores_text <- function(iforest, nf, shape, alpha) {
  if (iforest)
    return(list(sub = "No null for a forest, so no band: the line is a share, not a test",
                message = NULL))
  if (!nf)
    return(list(sub = sprintf("Nothing flagged at a false-discovery rate of %s",
                              format(alpha)), message = NULL))
  if (is.null(shape) || is.na(shape$shape))
    return(list(sub = sprintf("%d flagged at a false-discovery rate of %s",
                              nf, format(alpha)), message = NULL))
  rows <- if (nf == 1L) "row" else "rows"
  switch(shape$shape,
    separate = list(sub = sprintf(
      "%d flagged %s stand%s clear; the rest fall back into the band",
      nf, rows, if (nf == 1L) "s" else ""), message = NULL),
    unclear = list(sub = sprintf(
      "No clean break: %d of the next %d rows still score above the band",
      shape$above, shape$J), message = NULL),
    continuum = list(
      sub = sprintf("No break: %d of the next %d rows still score above the band",
                    shape$above, shape$J),
      message = paste0(
        "ilm_plot_anomaly(): the scores do not stop at the line. With the ",
        nf, " flagged ", rows, " set aside, ", shape$above, " of the next ",
        shape$J, " still score above what data with no anomalies produce at ",
        "those ranks, so this looks like the upper tail of a continuum rather ",
        "than a few rows that stand apart. Either say plainly that you are ",
        "reporting a fixed share of the rows, or treat the tail as structure ",
        "rather than as outliers: ilm_profile() looks for the groups in it.")))
}

## The band as drawing code, with every number inlined. tinyplot evaluates
## `draw` itself, after the axes and before the points, which is what puts the
## band BEHIND them; inlining means it does not matter which environment it
## evaluates in.
#' @keywords internal
#' @noRd
ilm_anom_band_draw <- function(pieces) {
  fill <- grDevices::adjustcolor("grey60", 0.35)
  calls <- lapply(pieces, function(p) bquote({
    graphics::polygon(c(.(p$x), rev(.(p$x))), c(.(p$lo), rev(.(p$hi))),
                      col = .(fill), border = NA)
    graphics::lines(.(p$x), .(p$mid), lty = 2, col = "grey40")
  }))
  as.call(c(as.name("{"), calls))
}

#' @keywords internal
#' @noRd
ilm_plot_anom_scores <- function(x, d, top_n, main, pch, ...) {
  n <- nrow(d); nf <- sum(d$flag)
  iforest <- identical(attr(x, "method"), "iforest")
  alpha <- attr(x, "alpha") %||% 0.05
  nq <- if (iforest) NULL else attr(x, "null_curve", exact = TRUE)
  if (!is.null(nq) && nrow(nq) != n)
    stop("this result has ", n, " rows but its reference was simulated for ",
         nrow(nq), ", so it looks like part of a scan. A rank only means ",
         "something against the whole of one: plot the full ilm_anomaly() ",
         "result, and take the flagged rows afterwards with ilm_anomalous().",
         call. = FALSE)
  if (!iforest && is.null(nq))
    message("ilm_plot_anomaly(): this result carries no reference curve, so ",
            "no band is drawn. Re-running ilm_anomaly() adds one.")
  k <- as.integer(if (is.null(top_n)) min(n, max(40L, 3L * nf))
                  else min(n, top_n))
  rk <- seq_len(k)
  s <- d$score
  shape <- if (is.null(nq)) NULL else ilm_anom_shape(s, nq, nf)
  txt <- ilm_anom_scores_text(iforest, nf, shape, alpha)

  ylim <- range(s[rk][is.finite(s[rk])])
  draw <- NULL
  if (!is.null(nq)) {
    pieces <- list()
    if (nf > 0L) {
      xl <- c(seq_len(min(nf, k)), if (nf < k) nf + 0.5)
      if (length(xl) > 1L) {
        b <- ilm_anom_band(nq, nf, xl, shifted = FALSE)
        pieces$left <- list(x = xl, lo = b[, "lo"], mid = b[, "mid"],
                            hi = b[, "hi"])
      }
    }
    if (nf < k) {
      xr <- if (nf > 0L) c(nf + 0.5, (nf + 1L):k) else rk
      if (length(xr) > 1L) {
        b <- ilm_anom_band(nq, nf, xr, shifted = nf > 0L)
        pieces$right <- list(x = xr, lo = b[, "lo"], mid = b[, "mid"],
                             hi = b[, "hi"])
      }
    }
    if (length(pieces)) {
      ylim <- range(c(ylim, unlist(lapply(pieces, function(p) c(p$lo, p$hi)))))
      draw <- ilm_anom_band_draw(pieces)
    }
  }
  args <- list(x = rk, y = s[rk], type = "p", pch = pch, col = "grey40",
               xlab = "rank, highest score first",
               ylab = if (iforest) "isolation score" else "score (log scale)",
               main = if (is.null(main)) "Anomaly scores against rank" else main,
               sub = txt$sub, ylim = ylim)
  if (!iforest) args$log <- "y"
  if (!is.null(draw)) args$draw <- draw
  do.call(tinyplot::tinyplot, utils::modifyList(args, list(...)))
  fl <- d$flag[rk]
  if (any(fl))
    graphics::points(rk[fl], s[rk][fl], pch = pch, col = "firebrick")
  if (nf > 0L && nf < k) {
    graphics::abline(v = nf + 0.5, lty = 2, col = "grey20")
    ilm_anom_line_label(nf + 0.5, if (iforest)
      paste0("top ", format(100 * alpha), "%, by choice") else
      sprintf("%d flagged", nf))
  }
  if (!is.null(txt$message)) message(txt$message)
  invisible(NULL)
}

## A label beside the line, just inside the top of the plotting region, in
## data units whichever way the axis is scaled.
#' @keywords internal
#' @noRd
ilm_anom_line_label <- function(x, label) {
  u <- graphics::par("usr")
  y <- u[4] - 0.04 * (u[4] - u[3])
  if (graphics::par("ylog")) y <- 10^y
  graphics::text(x, y, label, pos = 4, cex = 0.8, col = "grey20")
}

## ---- "drivers" ---------------------------------------------------------------

## Counts of each column as the driver among the flagged rows, against the
## count the unflagged rows' drivers predict. The prediction is smoothed by
## half a row per column, so a column that never drives an unflagged row does
## not make one flag look infinitely surprising.
##
## "One column carries the flags" needs three things at once: at least five
## flags, at least 80% of them on one column, and more than that column's
## share among the other rows explains (binomial tail below 0.01). The last
## matters because columns are not equally likely to drive even in clean
## data -- one the shared structure explains poorly keeps more of its variance
## off it, and drives more rows for that reason alone.
#' @keywords internal
#' @noRd
ilm_anom_drivers_table <- function(d, cols) {
  lev <- unique(c(cols, d$driver[!is.na(d$driver)]))
  nf <- sum(d$flag)
  cnt <- as.numeric(table(factor(d$driver[d$flag], levels = lev)))
  oth <- as.numeric(table(factor(d$driver[!d$flag], levels = lev)))
  q <- (oth + 0.5) / (sum(oth) + 0.5 * length(lev))
  tab <- data.frame(column = lev, flagged = cnt, expected = nf * q,
                    stringsAsFactors = FALSE)
  tab <- tab[order(-tab$flagged, -tab$expected), , drop = FALSE]
  row.names(tab) <- NULL
  top <- tab[1L, ]
  share <- if (nf) top$flagged / nf else 0
  p <- if (nf) stats::pbinom(top$flagged - 1, nf, top$expected / nf,
                             lower.tail = FALSE) else 1
  one <- nf >= 5L && share >= 0.8 && p < 0.01
  sub <- if (one)
    sprintf("`%s` drives %d of %d flags; the other rows predict %s",
            top$column, top$flagged, nf, format(round(top$expected, 1)))
  else if (nf < 5L)
    sprintf("%d flag%s: too few to read a pattern into", nf,
            if (nf == 1L) "" else "s")
  else sprintf("No single column carries the flags (at most %d of %d)",
               top$flagged, nf)
  msg <- if (one) paste0(
    "ilm_plot_anomaly(): one column, `", top$column, "`, drives ", top$flagged,
    " of the ", nf, " flags, where the drivers of the other rows predict ",
    "about ", format(round(top$expected, 1)), ". A multivariate anomaly ",
    "rarely looks like that; a problem IN that column does -- a unit, a ",
    "sentinel code, a misplaced decimal. Look at the column on its own with ",
    "ilm_outliers() or ilm_describe(), and recode sentinel values with ",
    "ilm_recode_errors().")
  list(table = tab, one_column = one, sub = sub, message = msg)
}

#' @keywords internal
#' @noRd
ilm_plot_anom_drivers <- function(x, d, top_n, main, ...) {
  if (!any(d$flag))
    stop("no rows were flagged, so there are no drivers to count. ",
         "type = \"scores\" shows the whole scan against its reference.",
         call. = FALSE)
  v <- ilm_anom_drivers_table(d, attr(x, "columns"))
  k <- min(nrow(v$table), if (is.null(top_n)) 15L else as.integer(top_n))
  tab <- v$table[rev(seq_len(k)), , drop = FALSE]  # flipped: the first level is drawn at the bottom
  f <- factor(tab$column, levels = tab$column)
  op <- ilm_anom_label_margin(tab$column)
  on.exit(graphics::par(op), add = TRUE)
  args <- list(x = f, y = tab$flagged, type = "barplot", flip = TRUE,
               xlab = "", xaxt = "n",
               ylab = "flagged rows it drives (red tick: what the other rows predict)",
               main = if (is.null(main)) "Which column drives each flag" else main,
               sub = v$sub)
  do.call(tinyplot::tinyplot, utils::modifyList(args, list(...)))
  graphics::axis(2, at = seq_len(k), labels = tab$column, las = 1, tick = FALSE)
  graphics::segments(tab$expected, seq_len(k) - 0.35, tab$expected,
                     seq_len(k) + 0.35, col = "firebrick", lwd = 2)
  if (!is.null(v$message)) message(v$message)
  invisible(NULL)
}

## Room for the column names, written horizontally. Left to tinyplot, the
## names on a flipped bar chart run along the axis and whichever would overlap
## are dropped without a word -- usually including the one on the top bar,
## which is the one that matters. So the left margin is sized to the longest
## name here, the caller draws the axis itself, and the caller restores `mar`
## on exit.
#' @keywords internal
#' @noRd
ilm_anom_label_margin <- function(labels) {
  w <- max(graphics::strwidth(labels, units = "inches",
                              cex = graphics::par("cex.axis")))
  mar <- graphics::par("mar")
  mar[2L] <- max(mar[2L], w / graphics::par("csi") + 1.5)
  graphics::par(mar = mar)
}

## ---- shared: the data that was scanned ----------------------------------------

#' @keywords internal
#' @noRd
ilm_anom_data <- function(x, d, data, type) {
  dat <- if (!is.null(data)) data else attr(x, "data", exact = TRUE)
  if (is.null(dat))
    stop("type = \"", type, "\" needs the data that was scanned, and this ",
         "result did not keep it. Pass `data =` -- the frame ilm_anomaly() ",
         "was given -- or re-run it with keep_data = TRUE.", call. = FALSE)
  if (!is.data.frame(dat)) dat <- as.data.frame(dat)
  cols <- attr(x, "columns")
  miss <- setdiff(cols, names(dat))
  if (length(miss))
    stop("`data` has no column ", paste(miss, collapse = ", "), ", which the ",
         "scan used. Pass the frame that was scanned.", call. = FALSE)
  if (max(d$row) > nrow(dat))
    stop("the result refers to row ", max(d$row), " but `data` has ",
         nrow(dat), " rows. Pass the frame that was scanned.", call. = FALSE)
  dat
}

## ---- "map" -------------------------------------------------------------------

#' @keywords internal
#' @noRd
ilm_plot_anom_map <- function(x, d, data, main, pch, ...) {
  if (!requireNamespace("PCAmixdata", quietly = TRUE))
    stop("type = \"map\" places the rows with ilm_reduce(), which needs the ",
         "PCAmixdata package. Install it with install.packages(\"PCAmixdata\"), ",
         "or use type = \"scores\" or \"drivers\", which need nothing further.",
         call. = FALSE)
  dat <- ilm_anom_data(x, d, data, "map")
  cols <- attr(x, "columns")
  rows <- sort(d$row)
  ## Only the rows that were scored. For the reconstruction that leaves out
  ## rows with gaps, which PCAmix would otherwise fill with column means and
  ## draw at the centre, where nothing about them was measured.
  r <- ilm_reduce(dat[rows, cols, drop = FALSE], ndim = 2L)
  co <- r$ind_coord
  i <- match(rows, d$row)
  fl <- d$flag[i]
  sc <- d$score[i]
  pv <- r$eig$pct_var
  iforest <- identical(attr(x, "method"), "iforest")
  sub <- if (!any(fl)) "Nothing flagged" else if (iforest)
    "Red: flagged, by score. Together? ilm_profile(x) describes them"
  else "Red: flagged, by score. Odd OFF these axes, so may sit mid-cloud"
  args <- list(x = co$dim1, y = co$dim2, type = "p", pch = pch, col = "grey70",
               xlab = sprintf("Dim 1 (%.1f%%)", pv[1]),
               ylab = sprintf("Dim 2 (%.1f%%)", pv[2]),
               main = if (is.null(main)) "Where the flagged rows sit" else main,
               sub = sub)
  do.call(tinyplot::tinyplot, utils::modifyList(args, list(...)))
  graphics::abline(h = 0, v = 0, lty = 3, col = "grey70")
  if (any(fl)) {
    s <- sc[fl]
    cx <- if (length(s) > 1L && diff(range(s)) > 0)
      1 + 2 * (s - min(s)) / diff(range(s)) else rep(2, length(s))
    o <- order(s)                   # the highest score is drawn last, on top
    graphics::points(co$dim1[fl][o], co$dim2[fl][o], pch = pch, cex = cx[o],
                     col = grDevices::adjustcolor("firebrick", 0.8))
    ## labelled beside each point and clear of it: at a fixed offset the
    ## label of a large point sits underneath the point
    if (sum(fl) <= 12L)
      for (j in seq_along(s))
        graphics::text(co$dim1[fl][j], co$dim2[fl][j], labels = rows[fl][j],
                       pos = 4, offset = 0.3 + 0.35 * cx[j], cex = 0.7,
                       col = "firebrick4")
  }
  invisible(NULL)
}

## ---- "row" -------------------------------------------------------------------

## One row's values against its own columns (z) and against the structure
## (the residual, in units of that column's typical residual). Both are in
## "typical deviations", which is what lets them share an axis.
#' @keywords internal
#' @noRd
ilm_anom_row_values <- function(x, d, row, dat) {
  R <- attr(x, "residual", exact = TRUE)
  if (is.null(R))
    stop("this result carries no residuals; re-running ilm_anomaly() adds them.",
         call. = FALSE)
  cols <- attr(x, "columns")
  X <- as.matrix(as.data.frame(dat)[cols])
  scored <- which(stats::complete.cases(X))
  if (length(scored) != nrow(R))
    stop("`data` does not line up with the scan: it has ", length(scored),
         " complete rows across the scanned columns, and the scan scored ",
         nrow(R), ". Pass the frame that was scanned.", call. = FALSE)
  if (is.null(row)) row <- d$row[1]
  if (!is.numeric(row) || length(row) != 1L || !is.finite(row) ||
      row != round(row))
    stop("`row` must be a single row number, as in the `row` column of the ",
         "result", call. = FALSE)
  i <- match(row, scored)
  if (is.na(i))
    stop(if (row < 1 || row > nrow(X))
      paste0("there is no row ", row, " in the data, which has ", nrow(X))
      else paste0("row ", row, " was not scored: it has a missing value in at ",
                  "least one of the scanned columns"), call. = FALSE)
  Xc <- X[scored, , drop = FALSE]
  scl <- apply(Xc, 2L, stats::sd)
  scl[!is.finite(scl) | scl <= 0] <- 1
  z <- (Xc[i, ] - colMeans(Xc)) / scl
  ## the typical residual: robust, because the rows being looked for are in
  ## the same columns and would inflate a standard deviation
  rs <- apply(R, 2L, stats::mad)
  bad <- !is.finite(rs) | rs <= 0
  if (any(bad)) rs[bad] <- apply(R[, bad, drop = FALSE], 2L, stats::sd)
  rs[!is.finite(rs) | rs <= 0] <- 1
  e <- R[i, ] / rs
  vals <- data.frame(column = cols, z = as.numeric(z), residual = as.numeric(e),
                     stringsAsFactors = FALSE)
  vals <- vals[order(-abs(vals$residual)), , drop = FALSE]
  row.names(vals) <- NULL
  j <- match(row, d$row)
  jz <- which.max(abs(vals$z)); je <- which.max(abs(vals$residual))
  sub <- if (abs(vals$z[jz]) >= 3)
    sprintf("`%s` is extreme on its own (z = %.1f): ilm_outliers() sees it too",
            vals$column[jz], vals$z[jz])
  else if (abs(vals$residual[je]) >= 3)
    sprintf("Odd only as a combination: largest |z| %.1f, largest residual %.1f",
            abs(vals$z[jz]), abs(vals$residual[je]))
  else "Nothing unusual here, on its own or as a combination"
  list(values = vals, row = as.integer(row), score = d$score[j],
       flag = d$flag[j], sub = sub)
}

#' @keywords internal
#' @noRd
ilm_plot_anom_row <- function(x, d, row, top_n, data, main, ...) {
  if (identical(attr(x, "method"), "iforest"))
    stop("type = \"row\" sets each value against where a fitted structure puts ",
         "it, and an isolation forest fits none. For this view re-run ",
         "ilm_anomaly() with method = \"reconstruction\", which uses the ",
         "numeric columns; type = \"drivers\" works for a forest.", call. = FALSE)
  dat <- ilm_anom_data(x, d, data, "row")
  v <- ilm_anom_row_values(x, d, row, dat)
  k <- min(nrow(v$values), if (is.null(top_n)) 20L else as.integer(top_n))
  vals <- v$values[rev(seq_len(k)), , drop = FALSE]  # flipped: largest on top
  lv <- vals$column
  meas <- c("in its own column (z)", "off the structure (residual)")
  op <- ilm_anom_label_margin(lv)
  on.exit(graphics::par(op), add = TRUE)
  args <- list(x = factor(rep(lv, 2L), levels = lv),
               y = c(vals$z, vals$residual),
               by = factor(rep(meas, each = k), levels = meas),
               type = tinyplot::type_barplot(beside = TRUE), flip = TRUE,
               xlab = "", xaxt = "n", ylab = "typical deviations for that column",
               main = if (is.null(main))
                 sprintf("Row %d: score %s, %s", v$row, format(signif(v$score, 3)),
                         if (isTRUE(v$flag)) "flagged" else "not flagged")
                 else main,
               sub = v$sub, palette = c("grey55", "firebrick"),
               ## Given explicitly: tinyplot titles a legend by deparsing
               ## `by`, and through do.call() that is the evaluated factor --
               ## "invalid graphics state" as soon as the labels are long.
               legend = list(title = "how far out"))
  dots <- list(...)
  ## tinyplot reserves legend room by changing the device layout, which
  ## blanks a figure that already has one (par(mfrow), layout()).
  if (!identical(graphics::par("mfrow"), c(1L, 1L)) && is.null(dots$legend)) {
    args$legend <- FALSE
    message("ilm_plot_anomaly(): legend dropped inside a multi-panel layout; ",
            "grey is the z-score, red the residual.")
  }
  do.call(tinyplot::tinyplot, utils::modifyList(args, dots))
  graphics::axis(2, at = seq_len(k), labels = lv, las = 1, tick = FALSE)
  graphics::abline(v = 0, col = "grey30")
  graphics::abline(v = c(-3, 3), lty = 3, col = "grey50")
  invisible(NULL)
}
