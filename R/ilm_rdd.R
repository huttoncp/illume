## ---------------------------------------------------------------------------
## Regression discontinuity.
##
## Treatment is assigned by whether a running variable crosses a cutoff, so
## units just either side are comparable and the jump at the cutoff is the
## effect. Estimation is a local linear fit on each side; the interesting part
## is whether the design holds up:
##
##  - did anyone move themselves across the cutoff? (a density jump)
##  - do covariates jump too? (they should not -- only treatment changes there)
##  - is there a jump at a cutoff where nothing happened? (a placebo)
##  - how much does the answer depend on the bandwidth?
##
## On inference: `rdrobust` implements the bias-corrected intervals of Calonico,
## Cattaneo and Titiunik, which are the state of the art, and this does not
## reimplement them. What it reports is an honest local linear fit with the
## bandwidth sensitivity laid out, plus the design checks above, inside the same
## object and workflow as everything else in the package.
## ---------------------------------------------------------------------------

## Triangular kernel weights: linear down from 1 at the cutoff to 0 at the
## bandwidth. Standard for local linear regression discontinuity, and it is what
## makes the fit LOCAL rather than a global polynomial in disguise.
#' @keywords internal
#' @noRd
ilm_rdd_weights <- function(r, h, kernel = "triangular") {
  u <- abs(r) / h
  switch(kernel,
    triangular = pmax(0, 1 - u),
    uniform = as.numeric(u <= 1),
    epanechnikov = pmax(0, 0.75 * (1 - u^2)),
    stop("unknown `kernel`: ", kernel, call. = FALSE))
}

## A rule-of-thumb bandwidth. Deliberately simple and deliberately labelled as
## such: it is a starting point for the sensitivity sweep rather than a claim to
## be the optimal one, which is what rdrobust computes properly.
#' @keywords internal
#' @noRd
ilm_rdd_bw <- function(r, y, cutoff) {
  s <- stats::sd(r, na.rm = TRUE); n <- sum(!is.na(r))
  h <- 1.84 * s * n^(-1 / 5)
  ## must leave enough points on both sides to fit a line
  rep_ok <- function(hh) min(sum(r >= cutoff - hh & r < cutoff),
                             sum(r >= cutoff & r <= cutoff + hh))
  while (rep_ok(h) < 10L && h < 4 * s) h <- h * 1.25
  h
}

#' @keywords internal
#' @noRd
ilm_rdd_fit1 <- function(d, y, run, cutoff, h, kernel, covariates, poly,
                         family, ...) {
  dd <- d[abs(d[[run]] - cutoff) <= h, , drop = FALSE]
  dd$.r <- dd[[run]] - cutoff
  dd$.above <- as.integer(dd$.r >= 0)
  dd$.w <- ilm_rdd_weights(dd$.r, h, kernel)
  dd <- dd[dd$.w > 0 & !is.na(dd[[y]]), , drop = FALSE]
  nb <- sum(dd$.above == 0L); na <- sum(dd$.above == 1L)
  if (nb < 5L || na < 5L)
    return(list(fit = NULL, n_below = nb, n_above = na,
                msg = sprintf("only %d below and %d above the cutoff within h", nb, na)))
  ## separate slopes either side: a common slope would force the two sides to
  ## share a trend and push any curvature into the jump
  trend <- if (poly == 1L) c(".r", ".above:.r")
           else c(sprintf("I(.r^%d)", seq_len(poly)),
                  sprintf(".above:I(.r^%d)", seq_len(poly)))
  form <- stats::reformulate(c(".above", trend, covariates), response = y)
  environment(form) <- environment()
  fit <- tryCatch(suppressWarnings(
           ilm_model(form, data = dd, family = family, weights = dd$.w,
                     verbose = FALSE, ...)),
         error = function(e) NULL)
  list(fit = fit, n_below = nb, n_above = na, data = dd,
       msg = if (is.null(fit)) "the model would not fit" else "")
}

## Does the density of the running variable jump at the cutoff? If people can
## move themselves across it, the ones just above are not comparable to the ones
## just below, and no bandwidth repairs that.
##
## Comparing the raw counts either side against a 50/50 split is NOT this test,
## though it looks like it: that asks whether the density is SYMMETRIC about the
## cutoff, where the question is whether it is CONTINUOUS there. A density with
## an ordinary slope fails the first while passing the second, and slopes are
## the common case. Measured over 500 unmanipulated data sets:
##
##                                 count split   local linear
##   uniform                          0.044         0.046
##   normal centred at the cutoff     0.038         0.050
##   normal centred away from it      0.760         0.056
##   exponential                      1.000         0.060
##
## so the count split flags three quarters of clean sloped data and all of a
## clean exponential one. It is not less powerful either way round: with 30% of
## the units just below the cutoff moved above it, the count split caught 0.808
## and the local linear fit 0.984.
##
## So: a local linear density on each side, compared at the cutoff, which is
## McCrary's test in its binned form.
#' @keywords internal
#' @noRd
ilm_rdd_density <- function(r, cutoff, h, nbins = 20L) {
  win <- r[abs(r - cutoff) <= h & !is.na(r)]
  nb <- sum(win < cutoff); na <- sum(win >= cutoff)
  out <- function(p, note = "") list(n_below = nb, n_above = na, p_value = p,
    status = if (is.na(p)) "UNTESTED" else if (p < 0.05) "FAIL" else "OK",
    note = note)
  if (nb < 20L || na < 20L)
    return(out(NA_real_, "too few observations within the bandwidth"))

  side_fit <- function(v, lo, hi) {
    br <- seq(lo, hi, length.out = nbins + 1L)
    mid <- (br[-1] + br[-length(br)]) / 2
    cnt <- as.integer(table(cut(v, br, include.lowest = TRUE)))
    ## a Poisson log-density, which is what a count per bin is, extrapolated to
    ## the cutoff; the width of each bin is common so it drops out
    g <- tryCatch(suppressWarnings(
           stats::glm(cnt ~ mid, family = stats::poisson())),
         error = function(e) NULL)
    if (is.null(g)) return(NULL)
    pr <- tryCatch(stats::predict(g, newdata = data.frame(mid = cutoff),
                                  se.fit = TRUE), error = function(e) NULL)
    if (is.null(pr)) NULL else c(fit = unname(pr$fit), se = unname(pr$se.fit))
  }
  a <- side_fit(win[win < cutoff], cutoff - h, cutoff)
  b <- side_fit(win[win >= cutoff], cutoff, cutoff + h)
  if (is.null(a) || is.null(b) || !all(is.finite(c(a, b))))
    return(out(NA_real_, "the density model would not fit"))
  se <- sqrt(a[["se"]]^2 + b[["se"]]^2)
  if (!is.finite(se) || se <= 0) return(out(NA_real_, "degenerate density fit"))
  z <- (b[["fit"]] - a[["fit"]]) / se
  out(2 * stats::pnorm(-abs(z)),
      sprintf("log-density jump %+.3f (se %.3f)", b[["fit"]] - a[["fit"]], se))
}

#' @keywords internal
#' @noRd
ilm_rdd_jump <- function(fit) {
  if (is.null(fit)) return(NULL)
  k <- ilm_term_cols(fit, ".above")
  if (!length(k)) return(NULL)
  ct <- ilm_coef_table(fit)
  crit <- if (isTRUE(fit$exact_df)) stats::qt(0.975, fit$resid_df) else stats::qnorm(0.975)
  data.frame(estimate = ct[k, 1], se = ct[k, 2],
             lower = ct[k, 1] - crit * ct[k, 2],
             upper = ct[k, 1] + crit * ct[k, 2],
             p_value = ct[k, 4], stringsAsFactors = FALSE)
}

#' Regression discontinuity
#'
#' Estimates the jump in an outcome at the point where treatment switches on,
#' by fitting a straight line on each side of the cutoff within a bandwidth and
#' taking the difference at the cutoff.
#'
#' @section What this reports, and what it does not:
#'
#' A local linear fit with a triangular kernel, the bandwidth sensitivity laid
#' out, and the design checks below. It does **not** implement the
#' bias-corrected robust intervals of Calonico, Cattaneo and Titiunik; the
#' `rdrobust` package does that and does it well. The interval here is the
#' ordinary one for a weighted local fit, which understates uncertainty when the
#' bandwidth is large enough for curvature to matter -- which is exactly what
#' `$bandwidth` is there to show you.
#'
#' The default is local **linear**. A high-order global polynomial produces
#' estimates driven by points far from the cutoff and is a well-documented way
#' to manufacture a discontinuity (Gelman and Imbens 2019). `poly` is available
#' and going above 2 warns.
#'
#' @section The design checks:
#'
#' A regression discontinuity is only as good as the claim that crossing the
#' cutoff is the *only* thing that changes there.
#'
#' * **Density** -- if people can move themselves across the cutoff, the ones
#'   just above are not comparable to the ones just below. The test is whether
#'   the density is *discontinuous* at the cutoff, not whether it is symmetric
#'   about it: a local linear density is fitted on each side and compared there.
#'   Over 500 unmanipulated data sets it flagged 0.046 to 0.060 of them across
#'   uniform, normal and exponential running variables, and caught 0.984 of the
#'   cases where 30% of the units just below the cutoff had been moved above
#'   it.
#' * **Covariate balance** -- anything measured before treatment should not jump
#'   at the cutoff. If it does, something other than treatment changes there.
#' * **Placebo cutoffs** -- the same estimate at points where nothing happened.
#'   These should be null; a comparable jump elsewhere means the method is
#'   finding jumps in noise.
#' * **Bandwidth sensitivity** -- the estimate across a range of bandwidths. An
#'   effect that appears only in a narrow window is not an effect.
#'
#' @section Fuzzy assignment:
#'
#' This estimates a **sharp** design, where crossing the cutoff determines
#' treatment. If a treatment column is supplied and assignment turns out fuzzy
#' -- the probability of treatment jumps but not from 0 to 1 -- this says so and
#' declines to report the sharp estimate as if it were the effect of treatment.
#' A fuzzy design needs an instrumental-variables estimator, which illume does
#' not have.
#'
#' @param data A data frame.
#' @param y Outcome column.
#' @param running Running (forcing) variable.
#' @param cutoff Value of `running` at which treatment switches on.
#' @param treatment Optional column of actual treatment received, used to check
#'   whether assignment is sharp.
#' @param h Bandwidth. `NULL` uses a rule-of-thumb starting point.
#' @param kernel `"triangular"` (default), `"uniform"` or `"epanechnikov"`.
#' @param poly Polynomial order on each side; 1 is local linear.
#' @param covariates Further columns for the mean structure.
#' @param family Response distribution; inferred from `y` when `NULL`.
#' @param bw_range Multipliers of `h` for the sensitivity sweep.
#' @param placebo Placebo cutoffs, as quantiles of the running variable either
#'   side of the real cutoff. `NULL` picks a few.
#' @param verbose Narrate each step.
#' @param ... Passed to [ilm_model()].
#' @return An object of class `"ilm_rdd"`: `jump`, `fit`, `bandwidth`,
#'   `density`, `balance`, `placebo` and the settings used.
#' @references
#' Gelman, A. and Imbens, G. (2019). Why high-order polynomials should not be
#' used in regression discontinuity designs. Journal of Business and Economic
#' Statistics 37(3).
#'
#' Calonico, S., Cattaneo, M. D. and Titiunik, R. (2014). Robust nonparametric
#' confidence intervals for regression-discontinuity designs. Econometrica
#' 82(6).
#' @seealso [ilm_plot_rdd()], [ilm_did()], [ilm_model()].
#' @examples
#' set.seed(1); n <- 2000
#' r <- runif(n, -1, 1)
#' y <- 0.5 * r + 0.8 * (r >= 0) + rnorm(n, 0, 0.5)
#' ilm_rdd(data.frame(r = r, y = y), "y", "r", cutoff = 0, verbose = FALSE)
#' @export
ilm_rdd <- function(data, y, running, cutoff = 0, treatment = NULL, h = NULL,
                    kernel = "triangular", poly = 1L, covariates = NULL,
                    family = NULL, bw_range = c(0.5, 0.75, 1, 1.5, 2),
                    placebo = NULL, verbose = TRUE, ...) {
  if (!is.data.frame(data))
    stop("`data` must be a data frame; it is ", class(data)[1], call. = FALSE)
  miss <- setdiff(c(y, running, treatment, covariates), names(data))
  if (length(miss))
    stop("column(s) not found in the data: ", paste(miss, collapse = ", "),
         ". Available: ", paste(utils::head(names(data), 12), collapse = ", "),
         call. = FALSE)
  if (!is.numeric(data[[running]]))
    stop("`running` (", running, ") must be numeric; it is ",
         class(data[[running]])[1], call. = FALSE)
  poly <- as.integer(poly)
  if (poly < 1L) stop("`poly` must be at least 1", call. = FALSE)
  if (poly > 2L)
    warning("a polynomial of order ", poly, " lets points far from the cutoff ",
            "drive the estimate, which is a known way to manufacture a ",
            "discontinuity. Prefer poly = 1 with a smaller bandwidth.",
            call. = FALSE)
  say <- function(...) if (verbose) message(...)
  d <- data[!is.na(data[[running]]), , drop = FALSE]
  r <- d[[running]]

  say("== regression discontinuity ==")
  say("[1/6] design")
  nb <- sum(r < cutoff); na <- sum(r >= cutoff)
  say("  ", nrow(d), " rows: ", nb, " below the cutoff, ", na, " at or above")
  if (nb < 10L || na < 10L)
    stop("too few observations on one side of the cutoff (", nb, " below, ",
         na, " above)", call. = FALSE)

  ## ---- sharp or fuzzy ------------------------------------------------------
  sharp <- TRUE; comply <- NULL
  if (!is.null(treatment)) {
    tr <- ilm_did_binary(d[[treatment]], treatment)
    pb <- mean(tr[r < cutoff] == 1L); pa <- mean(tr[r >= cutoff] == 1L)
    comply <- c(below = pb, above = pa)
    sharp <- pb < 0.02 && pa > 0.98
    say(sprintf("  treatment received: %.1f%% below, %.1f%% above",
                100 * pb, 100 * pa))
    if (!sharp) {
      say("  FUZZY assignment: crossing the cutoff changes the probability of ",
          "treatment but does not determine it. The jump below is then the ",
          "effect of BEING ELIGIBLE, not of being treated; scaling it to the ",
          "latter needs an instrumental-variables estimator, which illume ",
          "does not have.")
    } else say("  sharp assignment")
  }

  ## ---- family and bandwidth ------------------------------------------------
  say("[2/6] response and bandwidth")
  fam <- if (is.null(family)) {
    f <- ilm_dag_family_guess(d[[y]], y)
    say("  `", y, "` looks ", ilm_var_kind(d[[y]]), " -> family \"", f, "\"")
    f
  } else { say("  family \"", family, "\" as supplied"); family }
  hh <- if (is.null(h)) {
    hb <- ilm_rdd_bw(r, d[[y]], cutoff)
    say(sprintf("  bandwidth %.4g (rule of thumb -- a starting point, not an optimum;",
                hb)); say("    see $bandwidth, and rdrobust for a chosen one)")
    hb
  } else { say(sprintf("  bandwidth %.4g as supplied", h)); h }

  ## ---- the estimate --------------------------------------------------------
  say("[3/6] estimate")
  m <- ilm_rdd_fit1(d, y, running, cutoff, hh, kernel, covariates, poly, fam, ...)
  if (is.null(m$fit))
    stop("could not fit at the cutoff: ", m$msg, call. = FALSE)
  jump <- ilm_rdd_jump(m$fit)
  say("  ", m$n_below, " below and ", m$n_above, " above within the bandwidth")
  say(sprintf("  jump = %.4f  (%.4f, %.4f)  p = %s", jump$estimate,
              jump$lower, jump$upper,
              format.pval(jump$p_value, digits = 3, eps = 1e-4)))

  ## ---- bandwidth sensitivity ----------------------------------------------
  say("[4/6] bandwidth sensitivity")
  bw <- do.call(rbind, lapply(sort(unique(c(bw_range, 1))), function(k) {
    mk <- ilm_rdd_fit1(d, y, running, cutoff, hh * k, kernel, covariates,
                       poly, fam, ...)
    jk <- ilm_rdd_jump(mk$fit)
    if (is.null(jk)) return(NULL)
    data.frame(multiplier = k, h = hh * k, n = mk$n_below + mk$n_above,
               estimate = jk$estimate, se = jk$se, lower = jk$lower,
               upper = jk$upper, p_value = jk$p_value)
  }))
  if (!is.null(bw)) {
    rownames(bw) <- NULL
    say(sprintf("  estimate ranges %.4f to %.4f across %.2gx to %.2gx the bandwidth",
                min(bw$estimate), max(bw$estimate), min(bw$multiplier),
                max(bw$multiplier)))
    sgn <- all(bw$lower > 0) || all(bw$upper < 0)
    if (!sgn)
      say("    the interval covers zero at some bandwidths: the result depends ",
          "on how wide a window is taken, which is worth reporting")
  }

  ## ---- density and balance -------------------------------------------------
  say("[5/6] design checks")
  dens <- ilm_rdd_density(r, cutoff, hh)
  k_below <- dens$n_below; k_above <- dens$n_above
  dens_p <- dens$p_value
  say(sprintf("  density at the cutoff: %d below, %d above, p = %s -- %s",
              k_below, k_above, format.pval(dens_p, digits = 3, eps = 1e-4),
              dens$status))
  if (identical(dens$status, "FAIL"))
    say("    observations pile up on one side, which is what manipulating the ",
        "running variable looks like. The units either side may not be ",
        "comparable, and no bandwidth fixes that.")

  bal <- NULL
  if (length(covariates)) {
    bal <- do.call(rbind, lapply(covariates, function(cv) {
      if (!is.numeric(d[[cv]])) return(NULL)
      mk <- ilm_rdd_fit1(d, cv, running, cutoff, hh, kernel, NULL, poly,
                         "gaussian", ...)
      jk <- ilm_rdd_jump(mk$fit)
      if (is.null(jk)) return(NULL)
      data.frame(covariate = cv, estimate = jk$estimate, se = jk$se,
                 p_value = jk$p_value,
                 status = if (jk$p_value < 0.05) "FAIL" else "OK",
                 stringsAsFactors = FALSE)
    }))
    if (!is.null(bal)) {
      rownames(bal) <- NULL
      nbad <- sum(bal$status == "FAIL")
      say("  covariate balance: ", nbad, " of ", nrow(bal),
          " jump at the cutoff")
      if (nbad)
        say("    a covariate that jumps means something besides treatment ",
            "changes at the cutoff, and the estimate absorbs it")
    }
  } else say("  covariate balance: no covariates given")

  ## ---- placebo cutoffs -----------------------------------------------------
  say("[6/6] placebo cutoffs")
  pc <- if (is.null(placebo)) {
    lo <- stats::quantile(r[r < cutoff], c(0.33, 0.66), na.rm = TRUE)
    hi <- stats::quantile(r[r >= cutoff], c(0.33, 0.66), na.rm = TRUE)
    unname(c(lo, hi))
  } else placebo
  plc <- do.call(rbind, lapply(pc, function(cc) {
    ## only the side the placebo sits on, so the real cutoff is not inside
    sub <- if (cc < cutoff) d[r < cutoff, , drop = FALSE]
           else d[r >= cutoff, , drop = FALSE]
    mk <- ilm_rdd_fit1(sub, y, running, cc, hh, kernel, covariates, poly, fam, ...)
    jk <- ilm_rdd_jump(mk$fit)
    if (is.null(jk)) return(NULL)
    data.frame(cutoff = cc, estimate = jk$estimate, se = jk$se,
               p_value = jk$p_value,
               status = if (jk$p_value < 0.05) "FAIL" else "OK",
               stringsAsFactors = FALSE)
  }))
  if (!is.null(plc)) {
    rownames(plc) <- NULL
    nbad <- sum(plc$status == "FAIL")
    say("  ", nbad, " of ", nrow(plc), " placebo cutoffs show a jump")
    if (nbad)
      say("    a jump where nothing happens means the method is finding jumps ",
          "in noise; treat the estimate at the real cutoff with caution")
  } else say("  none could be fitted")

  structure(list(jump = jump, fit = m$fit, bandwidth = bw, density = dens,
                 balance = bal, placebo = plc, h = hh, cutoff = cutoff,
                 kernel = kernel, poly = poly, family = fam, y = y,
                 running = running, sharp = sharp, compliance = comply,
                 n_below = m$n_below, n_above = m$n_above, data = d),
            class = "ilm_rdd")
}

#' @export
print.ilm_rdd <- function(x, ...) {
  cat("<ilm_rdd>", x$y, "at", x$running, "=", x$cutoff, "\n")
  cat("  local", if (x$poly == 1L) "linear" else paste0("poly-", x$poly),
      "fit,", x$kernel, "kernel, h =", signif(x$h, 4),
      sprintf("(%d below, %d above)\n", x$n_below, x$n_above))
  cat("  family:", x$family, "\n")
  if (!isTRUE(x$sharp))
    cat("\n  FUZZY assignment: the jump below is the effect of being ELIGIBLE,\n",
        "  not of being treated.\n", sep = "")
  cat(sprintf("\n  jump %8.4f  (%7.4f, %7.4f)  p = %s\n", x$jump$estimate,
              x$jump$lower, x$jump$upper,
              format.pval(x$jump$p_value, digits = 3, eps = 1e-4)))
  if (!is.null(x$bandwidth))
    cat(sprintf("  across bandwidths %.2gx-%.2gx: %.4f to %.4f\n",
                min(x$bandwidth$multiplier), max(x$bandwidth$multiplier),
                min(x$bandwidth$estimate), max(x$bandwidth$estimate)))
  cat("\n  design checks\n")
  cat("    density at the cutoff:", x$density$status,
      sprintf("(%d below, %d above)\n", x$density$n_below, x$density$n_above))
  if (!is.null(x$balance))
    cat("    covariate balance:   ", sum(x$balance$status == "FAIL"), "of",
        nrow(x$balance), "jump\n")
  if (!is.null(x$placebo))
    cat("    placebo cutoffs:     ", sum(x$placebo$status == "FAIL"), "of",
        nrow(x$placebo), "show a jump\n")
  cat("\n  Interval is the ordinary one for a weighted local fit. For\n",
      "  bias-corrected robust intervals see the rdrobust package.\n", sep = "")
  invisible(x)
}

#' Regression discontinuity plot
#'
#' Binned means of the outcome against the running variable, with the fitted
#' line either side of the cutoff. The canonical figure: it shows the jump, and
#' it shows whether the jump is the only thing going on.
#'
#' @param x An [ilm_rdd()] result.
#' @param nbins Bins on each side.
#' @param main,xlab,ylab Labels.
#' @param ... Passed to [graphics::plot()].
#' @return `x`, invisibly.
#' @seealso [ilm_rdd()].
#' @examples
#' set.seed(1); n <- 2000
#' r <- runif(n, -1, 1)
#' y <- 0.5 * r + 0.8 * (r >= 0) + rnorm(n, 0, 0.5)
#' fit <- ilm_rdd(data.frame(r = r, y = y), "y", "r", cutoff = 0,
#'                verbose = FALSE)
#' ilm_plot_rdd(fit)
#' @export
ilm_plot_rdd <- function(x, nbins = 20L, main = "Regression discontinuity",
                         xlab = NULL, ylab = NULL, ...) {
  if (!inherits(x, "ilm_rdd"))
    stop("`x` must be an ilm_rdd result, not ", class(x)[1], call. = FALSE)
  d <- x$data; r <- d[[x$running]]; yv <- d[[x$y]]
  keep <- abs(r - x$cutoff) <= x$h & !is.na(yv)
  r <- r[keep]; yv <- yv[keep]
  if (is.null(xlab)) xlab <- x$running
  if (is.null(ylab)) ylab <- x$y

  side <- r >= x$cutoff
  bin_stats <- function(rr, yy) {
    if (length(rr) < 2L) return(NULL)
    br <- seq(min(rr), max(rr), length.out = nbins + 1L)
    k <- cut(rr, br, include.lowest = TRUE)
    data.frame(x = tapply(rr, k, mean), y = tapply(yy, k, mean),
               n = tapply(yy, k, length))
  }
  b0 <- bin_stats(r[!side], yv[!side]); b1 <- bin_stats(r[side], yv[side])
  bs <- rbind(b0, b1); bs <- bs[!is.na(bs$x), , drop = FALSE]

  plot(bs$x, bs$y, type = "n", main = main, xlab = xlab, ylab = ylab,
       ylim = range(c(bs$y, yv), finite = TRUE), ...)
  graphics::points(r, yv, pch = 16, cex = 0.25,
                   col = grDevices::adjustcolor("grey60", 0.35))
  graphics::points(bs$x, bs$y, pch = 21, bg = "white", col = "black", cex = 0.9)
  graphics::abline(v = x$cutoff, lty = 2, col = "grey30")

  ## the fitted line either side, from the fit that produced the estimate
  for (s in c(FALSE, TRUE)) {
    rr <- r[side == s]
    if (length(rr) < 3L) next
    g <- seq(min(rr), max(rr), length.out = 60L)
    nd <- data.frame(.r = g - x$cutoff, .above = as.integer(s))
    nd$.w <- ilm_rdd_weights(nd$.r, x$h, x$kernel)
    pr <- tryCatch(as.numeric(stats::predict(x$fit, newdata = nd,
                                             type = "response")),
                   error = function(e) NULL)
    if (!is.null(pr) && length(pr) == length(g))
      graphics::lines(g, pr, lwd = 2, col = if (s) "black" else "grey25")
  }
  graphics::mtext(sprintf("jump %.3f (%.3f, %.3f); h = %.3g",
                          x$jump$estimate, x$jump$lower, x$jump$upper, x$h),
                  side = 3, line = 0.2, cex = 0.75, col = "grey30")
  invisible(x)
}
