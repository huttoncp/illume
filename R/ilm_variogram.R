## ---------------------------------------------------------------------------
## The residual variogram: correlation as a function of separation.
##
## WHY THIS HAD TO EXIST BEFORE ilm_car1() WAS USEFUL.
##
## ilm_check_ar() matches pairs on EXACT lags, which is right for evenly spaced
## data and nearly useless without it. Measured on a 60 x 6 panel with times
## drawn from 1..30, exact matching found 53 pairs at lag 1, 5 at lag 2 and none
## at all beyond: five of six lags came back INCONCLUSIVE. Fitting a CAR(1) term
## and then having no way to check it is exactly the gap this package is meant
## not to have.
##
## Binning by DISTANCE instead uses every within-group pair. The same 60 x 6
## panel yields 900 pairs spread across the whole range of separations.
##
## The reference is the same as everywhere else: simulate from the fit, refit,
## and compare. That matters more here than usual, because residual correlation
## within a group is biased downward at every separation -- the residuals sum to
## roughly zero within a group -- so a band centred on zero is the wrong
## reference and would flag correct models.
## ---------------------------------------------------------------------------

#' Within-group pairs, and how far apart they are
#'
#' Every pair of observations in the same group, with the separation between
#' them. Computed once and reused across every refit, since it depends only on
#' the design.
#'
#' @param grp Integer group codes.
#' @param tim Numeric times.
#' @param max_pairs Cap on the number of pairs kept. A group of `m`
#'   observations contributes `m(m-1)/2` pairs, so a few long series can
#'   produce millions; beyond the cap a random subset is taken, which changes
#'   the precision but not what is being estimated.
#' @param seed Random seed for that subset.
#' @return A list with `i`, `j`, `d` and `n_total`.
#' @keywords internal
#' @noRd
ilm_pair_index <- function(grp, tim, max_pairs = 2e5, seed = 1L) {
  o <- order(grp, tim)
  g <- grp[o]
  runs <- rle(g)
  starts <- cumsum(c(0L, utils::head(runs$lengths, -1L)))
  ii <- jj <- vector("list", length(runs$lengths))
  for (k in seq_along(runs$lengths)) {
    m <- runs$lengths[k]
    if (m < 2L) next
    cb <- utils::combn(m, 2L)
    ii[[k]] <- o[starts[k] + cb[1, ]]
    jj[[k]] <- o[starts[k] + cb[2, ]]
  }
  i <- unlist(ii); j <- unlist(jj)
  n_total <- length(i)
  if (!n_total)
    stop("no group has two or more observations, so there are no pairs to ",
         "measure a correlation over", call. = FALSE)
  if (n_total > max_pairs) {
    set.seed(seed)
    keep <- sample.int(n_total, max_pairs)
    i <- i[keep]; j <- j[keep]
  }
  list(i = i, j = j, d = abs(tim[j] - tim[i]), n_total = n_total)
}

#' Correlation of residual pairs, by separation bin
#'
#' @param z Numeric residual vector.
#' @param pr A pair index from `ilm_pair_index()`.
#' @param bin Integer bin membership, one per pair.
#' @param nb Number of bins.
#' @param min_pairs Minimum pairs needed to report a bin.
#' @return A numeric vector of correlations, one per bin.
#' @keywords internal
#' @noRd
ilm_bin_cor <- function(z, pr, bin, nb, min_pairs = 30L) {
  a <- z[pr$i]; b <- z[pr$j]
  out <- rep(NA_real_, nb)
  for (k in seq_len(nb)) {
    s <- bin == k
    if (sum(s) < min_pairs) next
    x <- a[s]; y <- b[s]
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) >= min_pairs && stats::sd(x[ok]) > 0 && stats::sd(y[ok]) > 0)
      out[k] <- stats::cor(x[ok], y[ok])
  }
  out
}

#' Residual correlation as a function of separation in time
#'
#' Bins every pair of observations within a group by how far apart they are, and
#' reports the residual correlation in each bin against an envelope built by
#' simulating from the fitted model and refitting. This is the check for a
#' [ilm_car1()] term, and the one to use whenever observation times are
#' irregular.
#'
#' @section Why not [ilm_check_ar()]:
#' That function matches pairs at exact lags, which is right for evenly spaced
#' data and close to useless without it. On a 60-unit panel with six
#' observations each, drawn from thirty possible times, exact matching found 53
#' pairs at lag 1, five at lag 2 and none beyond -- five of six lags returned no
#' verdict. Binning by distance uses all 900 pairs.
#'
#' @section Reading it:
#' Correlation that starts high and decays toward zero as separation grows is
#' what an autoregressive process looks like, and is what [ilm_car1()] fits.
#' Correlation that is flat and positive at every separation is a group effect
#' the random intercept has not absorbed. Correlation that rises again at some
#' separation is a cycle; see [ilm_fourier()].
#'
#' The envelope is not centred on zero, and should not be. Residuals within a
#' group sum to roughly zero, so pairs of them correlate negatively even when
#' the model is exactly right, and the further apart they are the more strongly
#' this bites.
#'
#' @section Relation to the semivariogram:
#' For standardised residuals the semivariance at separation `d` is `1 - r(d)`,
#' so `type = "semivariance"` plots the same information the other way up, in
#' the form `nlme::Variogram()` uses.
#'
#' @param object A fitted `"ilm_model"` object.
#' @param time Time, one value per observation.
#' @param group Grouping variable, one value per observation.
#' @param breaks Number of separation bins, or explicit break points. Bins are
#'   chosen to hold roughly equal numbers of pairs, so a bin far out is not
#'   estimated from a handful of them.
#' @param B Simulated datasets behind the envelope. A p-value cannot fall below
#'   `1/(B+1)`, so about 100 is needed before `FAIL` is reachable.
#' @param type `"correlation"` or `"semivariance"`.
#' @param ncores Worker processes for the refits.
#' @param seed Random seed.
#' @param max_pairs Cap on the number of within-group pairs used. A group of
#'   `m` observations contributes `m(m-1)/2` of them, so a few long series can
#'   run to millions; beyond the cap a random subset is taken, which changes the
#'   precision but not what is being estimated.
#' @param plot Draw the variogram.
#' @param verbose Print the table.
#' @return Invisibly, a list with the per-bin `table`, the simulated null, and
#'   the number of replicates that refitted.
#' @references
#' Pinheiro, J. C., & Bates, D. M. (2000). *Mixed-Effects Models in S and
#' S-PLUS*. Springer. (Chapter 5 covers the residual variogram.)
#' @seealso [ilm_car1()] for the remedy, [ilm_check_ar()] for evenly spaced
#'   data, [ilm_plot_acf()].
#' @examples
#' set.seed(1)
#' d <- data.frame(id = factor(rep(1:20, each = 5)),
#'                 t = as.vector(replicate(20, sort(sample(1:20, 5)))),
#'                 x = rnorm(100))
#' d$y <- d$x + rnorm(100)
#' f <- ilm_model(y ~ x + (1 | id), data = d, family = "gaussian",
#'                verbose = FALSE)
#' ilm_variogram(f, d$t, d$id, breaks = 4, B = 15, plot = FALSE)
#' @export
ilm_variogram <- function(object, time, group, breaks = 8L, B = 100L,
                          type = "correlation", ncores = 1L, seed = 1L,
                          max_pairs = 2e5, plot = TRUE, verbose = TRUE) {
  if (!inherits(object, "ilm_model"))
    stop("`object` must be a fitted ilm_model object, not ", class(object)[1],
         call. = FALSE)
  if (length(type) != 1L || !type %in% c("correlation", "semivariance"))
    stop("unknown `type`: ", paste(sQuote(type), collapse = ", "),
         ". Options are 'correlation' and 'semivariance'.", call. = FALSE)
  N <- nrow(object$X)
  if (missing(time) || missing(group))
    stop("`time` and `group` are both required: a variogram measures how ",
         "correlation falls away with separation within a unit. Supply one ",
         "value of each per observation (the model has ", N, " rows).",
         call. = FALSE)
  if (length(time) != N || length(group) != N)
    stop("`time` has ", length(time), " values and `group` has ", length(group),
         ", but the model has ", N, " rows.", call. = FALSE)
  if (is.factor(time))
    stop("`time` is a factor. Converting it here would use the level order ",
         "rather than the labels, and silently change what a separation means. ",
         "Pass as.integer(time) if the levels are equally spaced steps.",
         call. = FALSE)
  tim <- suppressWarnings(as.numeric(time))
  if (anyNA(tim))
    stop("`time` must be numeric, integer or Date, not ", class(time)[1],
         call. = FALSE)
  B <- suppressWarnings(as.integer(B)[1])
  if (is.na(B) || B < 2L)
    stop("`B` must be at least 2; the envelope is the spread across ",
         "simulated replicates", call. = FALSE)
  if (1 / (B + 1) > 0.01)
    warning("B = ", B, " puts the smallest achievable p-value at ",
            signif(1 / (B + 1), 2),
            ", so a FAIL verdict is unreachable. Use B >= 100.", call. = FALSE)

  grp <- as.integer(factor(group))
  pr <- ilm_pair_index(grp, tim, max_pairs, seed)

  ## Equal-count bins, so a bin out at the far separations is not estimated
  ## from six pairs while the near ones rest on thousands.
  if (length(breaks) == 1L) {
    nb <- as.integer(breaks)
    if (is.na(nb) || nb < 2L)
      stop("`breaks` must be at least 2 bins, or a vector of break points",
           call. = FALSE)
    br <- unique(stats::quantile(pr$d, seq(0, 1, length.out = nb + 1L),
                                 na.rm = TRUE))
    if (length(br) < 3L)
      stop("the separations take only ", length(br),
           " distinct values, which is too few to bin. Use ilm_check_ar() ",
           "instead, which is built for a small set of exact lags.",
           call. = FALSE)
  } else {
    br <- sort(unique(as.numeric(breaks)))
  }
  bin <- as.integer(cut(pr$d, br, include.lowest = TRUE))
  nb <- length(br) - 1L
  npair <- tabulate(bin, nb)
  dmid <- vapply(seq_len(nb), function(k) {
    s <- bin == k
    if (any(s)) mean(pr$d[s]) else NA_real_
  }, 0)

  stat <- function(fit) {
    R <- ilm_pearson_ovr(fit)
    m <- matrix(NA_real_, nb, ncol(R))
    for (j in seq_len(ncol(R))) m[, j] <- ilm_bin_cor(R[, j], pr, bin, nb)
    m
  }
  obs <- stat(object)
  C <- ncol(obs)
  labels <- if (C == 1L) "residual" else object$ylevels

  env <- ilm_refit_stat(object, stat, c(nb, C), B, ncores, seed,
                        exports = c("pr", "bin", "nb"), where = environment())
  v <- ilm_lag_verdict(obs, env$null, labels, npair)
  tab <- v$table
  names(tab)[names(tab) == "lag"] <- "bin"
  names(tab)[names(tab) == "n_pairs"] <- "n_pairs"
  tab <- cbind(tab[1L], separation = round(dmid, 3), tab[-1L])
  if (env$n_ok < 10L) tab$status <- "INCONCLUSIVE"

  ## The semivariance of standardised residuals is 1 - correlation, so this is
  ## the same information the other way up. Subtracting from 1 reverses the
  ## order of the interval, so lo and hi swap; p and status are unchanged,
  ## because the test is the same test.
  if (identical(type, "semivariance")) {
    lo <- round(1 - tab$hi, 4); hi <- round(1 - tab$lo, 4)
    tab$estimate <- round(1 - tab$estimate, 4)
    tab$null_mean <- round(1 - tab$null_mean, 4)
    tab$lo <- lo; tab$hi <- hi
  }

  res <- list(table = tab, null = env$null, n_ok = env$n_ok, B = B,
              type = type, n_pairs_total = pr$n_total, breaks = br,
              rho = object$rho)
  if (verbose) ilm_variogram_report(res)
  if (plot) ilm_plot_variogram(res)
  invisible(res)
}

#' @keywords internal
#' @noRd
ilm_variogram_report <- function(res) {
  cat(sprintf("\nresidual %s by separation (%d refits, %d pairs)\n",
              res$type, res$n_ok, res$n_pairs_total))
  print(res$table, row.names = FALSE)
  bad <- which(res$table$status %in% c("WARN", "FAIL"))
  if (!length(bad)) {
    if (res$n_ok >= 10L)
      cat("\n>> correlation falls away with separation as the fitted model says it should.\n")
    else
      cat(sprintf("\n>> only %d of %d replicates refitted, so no verdict.\n",
                  res$n_ok, res$B))
    return(invisible(NULL))
  }
  cat(sprintf("\n>> %d of %d separation bins outside the envelope (bin %s).\n",
              length(bad), nrow(res$table), paste(bad, collapse = ", ")))
  cat("   why: residual correlation the model does not account for\n")
  cat("   try: ", ilm_variogram_advice(res), "\n", sep = "")
  cat("   note: bins are not independent, so treat the count as indicative\n")
  invisible(NULL)
}

## The shape of what is left over says which remedy to reach for, and the three
## shapes need three different ones.
##
## Judged on the DEPARTURE from the simulated null, not on the raw correlation.
## The null is not zero here -- residuals within a group sum to roughly zero, so
## every bin is biased negative -- and a rule written against zero reads that
## bias as signal.
##
##   alternating sign   a cycle. An autoregressive decay crosses the null at
##                      most once; an oscillation crosses it repeatedly.
##   flat and positive  a group effect the random structure has not absorbed.
##   decaying           autoregression, which is what ilm_car1() fits.
#' @keywords internal
#' @noRd
ilm_variogram_advice <- function(res) {
  tb <- res$table
  e <- tb$estimate; m <- tb$null_mean
  if (identical(res$type, "semivariance")) { e <- 1 - e; m <- 1 - m }
  dev <- e - m
  ok <- is.finite(dev)
  if (sum(ok) < 3L)
    return("add a correlation over time, e.g. ilm_car1(time, group)")
  dv <- dev[ok]
  sgn <- sign(dv); sgn <- sgn[sgn != 0]
  flips <- if (length(sgn) > 1L) sum(diff(sgn) != 0) else 0L
  if (flips >= 2L)
    return(paste0("the correlation alternates in sign across separations ",
                  "rather than decaying away, which is a cycle and not ",
                  "autoregression: see ilm_fourier() and ilm_cyclic()"))
  if (all(dv > 0) && diff(range(dv)) < 0.4 * mean(dv))
    return(paste0("the correlation sits at about the same level at every ",
                  "separation rather than decaying, which is a group effect: ",
                  "add or widen a random effect for the unit before reaching ",
                  "for ilm_car1()"))
  if (dv[1] > 0)
    return(paste0("the correlation is strongest at the shortest separations ",
                  "and decays from there, which is what ilm_car1(time, group) ",
                  "fits; ilm_ar1() if the times are evenly spaced"))
  paste0("residual correlation departs from what the model implies, but not ",
         "in a shape any one correlation structure fixes; check the mean ",
         "structure first with ilm_check_omitted()")
}

#' Plot a residual variogram
#'
#' @param res The value returned by [ilm_variogram()].
#' @param colour Colour for bins that are not flagged.
#' @param fill Envelope fill.
#' @param alpha Envelope transparency.
#' @param size Point scaling.
#' @param main Title.
#' @return Invisibly, the table behind the plot.
#' @seealso [ilm_variogram()].
#' @examples
#' set.seed(1)
#' d <- data.frame(id = factor(rep(1:20, each = 5)),
#'                 t = as.vector(replicate(20, sort(sample(1:20, 5)))),
#'                 x = rnorm(100))
#' d$y <- d$x + rnorm(100)
#' f <- ilm_model(y ~ x + (1 | id), data = d, family = "gaussian",
#'                verbose = FALSE)
#' v <- ilm_variogram(f, d$t, d$id, breaks = 4, B = 15, plot = FALSE,
#'                    verbose = FALSE)
#' ilm_plot_variogram(v)
#' @export
ilm_plot_variogram <- function(res, colour = "grey25", fill = "grey85",
                               alpha = NULL, size = 1, main = NULL) {
  if (!is.list(res) || is.null(res$table) || is.null(res$type))
    stop("`res` must be the value returned by ilm_variogram(), not ",
         class(res)[1], call. = FALSE)
  tb <- res$table
  x <- tb$separation; est <- tb$estimate; lo <- tb$lo; hi <- tb$hi
  if (!is.null(alpha)) fill <- grDevices::adjustcolor(fill, alpha)
  cols <- ILM_STATUS_COL[tb$status]; cols[tb$status == "OK"] <- colour
  ylim <- range(c(0, est, lo, hi), na.rm = TRUE)
  if (!all(is.finite(ylim))) ylim <- c(-1, 1)
  pad <- diff(ylim) * 0.18
  if (is.null(main))
    main <- if (identical(res$type, "semivariance")) "Residual semivariogram"
            else "Residual correlation by separation"

  op <- par(mar = c(6.4, 4.2, 3.6, 1.2)); on.exit(par(op), add = TRUE)
  plot(x, est, type = "n", ylim = c(ylim[1] - pad, ylim[2] + pad), xlab = "",
       ylab = if (identical(res$type, "semivariance")) "semivariance"
              else "correlation", main = main)
  title(xlab = "separation in time (within group)", line = 2.3)
  w <- diff(range(x, na.rm = TRUE)) / (2.6 * max(nrow(tb), 1L))
  for (k in seq_len(nrow(tb))) {
    if (!is.finite(lo[k]) || !is.finite(hi[k]) || !is.finite(x[k])) next
    polygon(x[k] + c(-w, w, w, -w), c(lo[k], lo[k], hi[k], hi[k]),
            col = fill, border = NA)
  }
  abline(h = if (identical(res$type, "semivariance")) 1 else 0,
         col = "grey40", lty = 2)
  ok <- is.finite(est) & is.finite(x)
  if (any(ok)) {
    lines(x[ok], est[ok], col = "grey55")
    points(x[ok], est[ok], pch = 16, cex = 1.1 * size, col = cols[ok])
  }
  ## An earlier version drew the fitted rho^d curve here for comparison. It is
  ## the MARGINAL correlation of the latent process, while the points are
  ## CONDITIONAL residual correlations -- shrunk by having conditioned on that
  ## process, and biased negative by the within-group centring. The two are not
  ## on the same scale and the curve sat well above the points on a model that
  ## fitted perfectly well. The grey envelope is what the model implies for
  ## this statistic, on this statistic's own scale, and is the honest
  ## comparison.

  bad <- which(tb$status %in% c("WARN", "FAIL"))
  verdict <- if (all(tb$status == "INCONCLUSIVE"))
    sprintf("INCONCLUSIVE: only %d of %d replicates refitted", res$n_ok, res$B)
  else if (!length(bad)) "OK: consistent with the fitted model"
  else sprintf("%s: bin %s outside the envelope",
               if (any(tb$status == "FAIL")) "FAIL" else "WARN",
               paste(bad, collapse = ", "))
  mtext(verdict, side = 3, line = 0.2, cex = 0.72,
        col = if (!length(bad)) "grey30" else
          ILM_STATUS_COL[[if (any(tb$status == "FAIL")) "FAIL" else "WARN"]])
  foot <- if (length(bad)) ilm_variogram_advice(res)
          else sprintf("band: 95%% envelope from %d refits of simulated data",
                       res$n_ok)
  ilm_mtext_wrap(foot, line = 3.6, cex = 0.62, col = "grey30")
  invisible(tb)
}
