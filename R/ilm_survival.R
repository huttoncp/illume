## ---------------------------------------------------------------------------
## Survival curves, and the check that the family was the right one.
##
## A parametric survival model buys you a smooth curve, extrapolation past the
## last event, and coefficients that mean something. What it costs is an
## assumption about the SHAPE of the baseline, and that assumption is the one
## worth checking. The Kaplan-Meier estimator makes no such assumption, so
## laying the fitted curve over it is the check: where the two part company, the
## family is wrong.
##
## The envelope around the comparison is built the way every other one in this
## package is -- simulate from the fit, refit, recompute -- because the
## Kaplan-Meier of a finite sample wobbles, and a fitted curve that tracks it to
## within that wobble is not evidence of anything wrong.
## ---------------------------------------------------------------------------

#' Kaplan-Meier estimate
#'
#' The product-limit estimator, with Greenwood standard errors. Written out
#' rather than taken from `survival`, which is a suggested package and so may
#' not be present when a diagnostic wants to run.
#'
#' @param time Follow-up time.
#' @param event `1` if the event was observed, `0` if censored.
#' @return A list with `time`, `surv`, `se` and `n_risk`.
#' @keywords internal
#' @noRd
ilm_km <- function(time, event) {
  o <- order(time)
  tt <- time[o]; ev <- event[o]
  ut <- sort(unique(tt[ev == 1L]))
  if (!length(ut))
    return(list(time = numeric(0), surv = numeric(0), se = numeric(0),
                n_risk = numeric(0)))
  n_risk <- vapply(ut, function(z) sum(tt >= z), 0)
  n_ev <- vapply(ut, function(z) sum(tt == z & ev == 1L), 0)
  surv <- cumprod(1 - n_ev / n_risk)
  ## Greenwood, guarding the last event time where everyone still at risk has
  ## the event and the term is undefined
  den <- n_risk * (n_risk - n_ev)
  v <- cumsum(ifelse(den > 0, n_ev / den, 0))
  list(time = ut, surv = surv, se = surv * sqrt(v), n_risk = n_risk)
}

## Step the estimate onto a fixed grid, so a statistic computed on it has the
## same shape for every replicate.
#' @keywords internal
#' @noRd
ilm_km_at <- function(km, grid) {
  if (!length(km$time)) return(rep(NA_real_, length(grid)))
  vapply(grid, function(g) {
    i <- sum(km$time <= g)
    if (i == 0L) 1 else km$surv[i]
  }, 0)
}

#' Predicted survival curve from a fitted model
#'
#' The probability of surviving past each time, for one or more covariate
#' patterns, with confidence intervals.
#'
#' @details
#' All three accelerated failure time families have a closed-form survivor
#' function in `z = (log t - eta) / scale`, so the curve and its interval come
#' straight from the fitted linear predictor. Intervals are computed on the
#' complementary log-log scale, `log(-log S)`, which is linear in `eta` and so
#' keeps the interval inside `(0, 1)` without truncation.
#'
#' Uncertainty in the scale parameter is not included: the interval reflects
#' uncertainty in the coefficients only. It is therefore slightly narrow, most
#' noticeably in the tail beyond the last observed event.
#'
#' @param object A fitted `"ilm_model"` with an accelerated failure time family.
#' @param newdata Covariate patterns, one row each. With none, the curve is
#'   drawn at the median of each numeric predictor and the commonest level of
#'   each factor.
#' @param times Times at which to evaluate. With none, a grid spanning the
#'   observed follow-up.
#' @param conf Confidence level.
#' @return A data frame with `row`, `time`, `surv`, `lower` and `upper`.
#' @seealso [ilm_plot_survival()], which draws it against the Kaplan-Meier
#'   estimate, [ilm_surv()].
#' @examples
#' set.seed(1)
#' d <- data.frame(x = rnorm(200))
#' tt <- exp(1.5 + 0.8 * d$x + 0.7 * log(rexp(200)))
#' ct <- rexp(200, rate = 1 / (2 * median(tt)))
#' d$time <- pmin(tt, ct); d$event <- as.integer(tt <= ct)
#' f <- ilm_model(time ~ x, data = d, family = "weibull",
#'                censor = ilm_surv(d$time, d$event), verbose = FALSE)
#' head(ilm_survival(f, newdata = data.frame(x = c(-1, 1))))
#' @export
ilm_survival <- function(object, newdata = NULL, times = NULL, conf = 0.95) {
  if (!inherits(object, "ilm_model"))
    stop("`object` must be a fitted ilm_model object, not ", class(object)[1],
         call. = FALSE)
  if (!isTRUE(object$family$aft))
    stop("a survival curve needs an accelerated failure time family. This ",
         "model is ", object$family$name, "; use family = \"weibull\", ",
         "\"lognormal\" or \"loglogistic\".", call. = FALSE)

  if (is.null(newdata)) newdata <- ilm_typical_row(object)
  pr <- stats::predict(object, newdata = newdata, se.fit = TRUE)
  eta <- log(as.matrix(if (is.list(pr)) pr$fit else pr)[, 1])
  se_eta <- if (is.list(pr) && !is.null(pr$se.fit)) {
    ## predict() returns the standard error on the response scale, and the
    ## curve needs it on the linear-predictor scale it was computed on
    as.matrix(pr$se.fit)[, 1] / as.matrix(if (is.list(pr)) pr$fit else pr)[, 1]
  } else rep(0, length(eta))
  sc <- unname(object$dispersion[[1]])

  if (is.null(times)) {
    y <- as.numeric(object$y)
    times <- seq(min(y[y > 0]), max(y), length.out = 100L)
  }
  times <- sort(unique(times[times > 0]))
  crit <- stats::qnorm(1 - (1 - conf) / 2)
  nm <- object$family$name

  out <- lapply(seq_along(eta), function(i) {
    z <- (log(times) - eta[i]) / sc
    S <- switch(nm,
      lognormal   = stats::pnorm(-z),
      loglogistic = 1 / (1 + exp(z)),
      weibull     = exp(-exp(z)))
    ## on the complementary log-log scale the curve is linear in eta, with
    ## slope -1/scale, so the delta method is exact there
    cll <- log(-log(pmin(pmax(S, 1e-12), 1 - 1e-12)))
    half <- crit * se_eta[i] / sc
    data.frame(row = i, time = times, surv = S,
               lower = exp(-exp(cll + half)),
               upper = exp(-exp(cll - half)))
  })
  do.call(rbind, out)
}

## The covariate pattern an effect plot would use: median for a number, the
## commonest level for a factor.
#' @keywords internal
#' @noRd
ilm_typical_row <- function(object) {
  mf <- object$model
  if (is.null(mf))
    stop("this model does not carry its data, so there is no typical ",
         "covariate pattern to draw at. Pass `newdata`.", call. = FALSE)
  mf <- as.data.frame(mf)
  resp <- if (!is.null(object$formula)) all.vars(object$formula[[2]]) else character(0)
  nd <- mf[1L, , drop = FALSE]
  for (cn in setdiff(names(mf), resp)) {
    v <- mf[[cn]]
    nd[[cn]] <- if (is.numeric(v) && !is.logical(v)) stats::median(v, na.rm = TRUE)
                else {
                  tb <- sort(table(v), decreasing = TRUE)
                  if (is.factor(v)) factor(names(tb)[1], levels = levels(v))
                  else names(tb)[1]
                }
  }
  nd
}

#' Survival curve against the Kaplan-Meier estimate
#'
#' Draws the fitted survival curve over the non-parametric estimate, which is
#' the check that the family was the right choice: the Kaplan-Meier makes no
#' assumption about the shape of the baseline, so where the two part company,
#' the parametric assumption is what is wrong.
#'
#' @section What the envelope is:
#' The Kaplan-Meier of a finite sample wobbles, so a fitted curve that tracks it
#' to within that wobble is not evidence of anything. The grey band is the
#' spread of the Kaplan-Meier across datasets simulated from the fit and
#' refitted, which is the same reference every other check in this package uses.
#' With `B = 0` the band is omitted and only the two curves are drawn.
#'
#' @param object A fitted `"ilm_model"` with an accelerated failure time family.
#' @param time,event The follow-up used to fit it. With `event` missing it is
#'   taken from the model's censoring specification.
#' @param by Optional grouping variable, one value per observation: one fitted
#'   curve and one Kaplan-Meier per level.
#' @param B Simulated datasets behind the envelope; `0` to skip it.
#' @param conf Confidence level for the Kaplan-Meier band when `B = 0`.
#' @param ncores Worker processes for the refits.
#' @param seed Random seed.
#' @param colour,fill,alpha,size Appearance.
#' @param main Title.
#' @param verbose Print the verdict.
#' @return Invisibly, a list with the fitted curve, the Kaplan-Meier, the
#'   largest gap between them and a `status`.
#' @seealso [ilm_survival()], [ilm_surv()].
#' @examples
#' set.seed(1)
#' d <- data.frame(x = rnorm(300))
#' tt <- exp(1.5 + 0.8 * d$x + 0.7 * log(rexp(300)))
#' ct <- rexp(300, rate = 1 / (2 * median(tt)))
#' d$time <- pmin(tt, ct); d$event <- as.integer(tt <= ct)
#' f <- ilm_model(time ~ x, data = d, family = "weibull",
#'                censor = ilm_surv(d$time, d$event), verbose = FALSE)
#' ilm_plot_survival(f, d$time, d$event, B = 0)
#' @export
ilm_plot_survival <- function(object, time, event, by = NULL, B = 60L,
                              conf = 0.95, ncores = 1L, seed = 1L,
                              colour = "#2C7FB8", fill = "grey85",
                              alpha = NULL, size = 1, main = NULL,
                              verbose = TRUE) {
  if (!inherits(object, "ilm_model"))
    stop("`object` must be a fitted ilm_model object, not ", class(object)[1],
         call. = FALSE)
  if (!isTRUE(object$family$aft))
    stop("a survival curve needs an accelerated failure time family. This ",
         "model is ", object$family$name, "; use family = \"weibull\", ",
         "\"lognormal\" or \"loglogistic\".", call. = FALSE)
  N <- nrow(object$X)
  if (missing(time)) time <- as.numeric(object$y)
  if (missing(event) || is.null(event)) {
    cs <- object$censor
    if (is.null(cs))
      stop("`event` is required when the model carries no censoring ",
           "specification. Pass the event indicator used to fit it.",
           call. = FALSE)
    event <- as.integer(ilm_censor_for(cs, as.numeric(object$y)) == 0L)
  }
  if (length(time) != N || length(event) != N)
    stop("`time` has ", length(time), " values and `event` has ", length(event),
         ", but the model has ", N, " rows.", call. = FALSE)
  if (!is.null(by) && length(by) != N)
    stop("`by` has ", length(by), " values but the model has ", N, " rows.",
         call. = FALSE)

  time <- as.numeric(time); event <- as.integer(event)
  grid <- seq(min(time[time > 0]), max(time), length.out = 100L)
  gvec <- if (is.null(by)) rep("all", N) else as.character(by)
  lev <- unique(gvec)

  ## observed Kaplan-Meier per level, on the shared grid
  obs <- vapply(lev, function(L) {
    s <- gvec == L
    ilm_km_at(ilm_km(time[s], event[s]), grid)
  }, numeric(length(grid)))
  obs <- matrix(obs, nrow = length(grid))

  ## the fitted curve, at the typical covariate pattern within each level
  mf <- object$model
  fitc <- vapply(lev, function(L) {
    nd <- if (is.null(by) || is.null(mf)) ilm_typical_row(object) else {
      s <- gvec == L
      r <- ilm_typical_row(object)
      ## honour the level being drawn when it is itself a predictor
      for (cn in names(r)) if (cn %in% names(mf) && all(mf[[cn]][s][1] == mf[[cn]][s]))
        r[[cn]] <- mf[[cn]][s][1]
      r
    }
    ilm_survival(object, newdata = nd, times = grid)$surv
  }, numeric(length(grid)))
  fitc <- matrix(fitc, nrow = length(grid))

  ## the envelope: how much the Kaplan-Meier moves across refits of data
  ## simulated from this model
  nullarr <- NULL; n_ok <- 0L
  if (B > 0L) {
    stat <- function(fit) {
      yv <- as.numeric(fit$y)
      ev <- as.integer(ilm_censor_for(object$censor, yv) == 0L)
      if (is.null(object$censor)) ev <- rep(1L, length(yv))
      m <- vapply(lev, function(L) {
        s <- gvec == L
        ilm_km_at(ilm_km(yv[s], ev[s]), grid)
      }, numeric(length(grid)))
      matrix(m, nrow = length(grid))
    }
    env <- ilm_refit_stat(object, stat, c(length(grid), length(lev)), B,
                          ncores, seed,
                          exports = c("grid", "gvec", "lev"),
                          where = environment())
    nullarr <- env$null; n_ok <- env$n_ok
  }

  lo <- hi <- NULL
  if (!is.null(nullarr) && n_ok >= 10L) {
    lo <- apply(nullarr, 1:2, stats::quantile, (1 - conf) / 2, na.rm = TRUE)
    hi <- apply(nullarr, 1:2, stats::quantile, 1 - (1 - conf) / 2, na.rm = TRUE)
  }

  gap <- max(abs(obs - fitc), na.rm = TRUE)
  outside <- if (!is.null(lo)) mean(obs < lo | obs > hi, na.rm = TRUE) else NA_real_
  status <- if (is.na(outside)) "INCONCLUSIVE" else
            if (outside > 0.20) "FAIL" else
            if (outside > 0.05) "WARN" else "OK"

  ## ---- draw ----------------------------------------------------------------
  if (!is.null(alpha)) fill <- grDevices::adjustcolor(fill, alpha)
  cols <- if (length(lev) == 1L) colour else ilm_curve_colours(length(lev))
  op <- par(mar = c(6.4, 4.2, 3.6, 1.2)); on.exit(par(op), add = TRUE)
  plot(range(grid), c(0, 1), type = "n", xlab = "", ylab = "surviving",
       main = main %||% "Survival: fitted against Kaplan-Meier")
  title(xlab = "time", line = 2.3)
  ## With one group the empirical curve is grey and the fitted one coloured.
  ## With several, both are coloured by group -- drawing every empirical curve
  ## in the same grey makes them impossible to pair with their fit, which is
  ## the whole point of splitting by group.
  one <- length(lev) == 1L
  kmcol <- if (one) "grey30" else grDevices::adjustcolor(cols, 0.55)
  for (k in seq_along(lev)) {
    if (!is.null(lo))
      polygon(c(grid, rev(grid)), c(lo[, k], rev(hi[, k])), col = fill,
              border = NA)
    graphics::lines(grid, obs[, k], type = "s",
                    col = if (one) kmcol else kmcol[k], lwd = 1.5)
    graphics::lines(grid, fitc[, k], col = cols[k], lwd = 2.2 * size)
  }
  if (one)
    legend("topright", legend = c("Kaplan-Meier", "fitted"),
           col = c("grey30", cols[1]), lwd = c(1.5, 2.2), bty = "n", cex = 0.75)
  else {
    legend("topright", legend = lev, col = cols, lwd = 2.2, bty = "n",
           cex = 0.7, title = "fitted")
    legend("right", legend = lev, col = kmcol, lwd = 1.5, lty = 1, bty = "n",
           cex = 0.7, title = "Kaplan-Meier")
  }
  mtext(sprintf("%s: largest gap %.3f%s", status, gap,
                if (is.na(outside)) "" else
                  sprintf(", %.0f%% of the curve outside the envelope",
                          100 * outside)),
        side = 3, line = 0.2, cex = 0.72,
        col = if (status == "OK") "grey30" else
          ILM_STATUS_COL[[if (status == "FAIL") "FAIL" else "WARN"]])
  foot <- if (status %in% c("WARN", "FAIL"))
    paste0("the fitted shape does not track the data: try another family ",
           "(weibull, lognormal, loglogistic) and compare by AIC")
  else if (is.na(outside))
    "no envelope; raise B for a verdict"
  else sprintf("band: %.0f%% envelope from %d refits of simulated data",
               100 * conf, n_ok)
  ilm_mtext_wrap(foot, line = 3.6, cex = 0.62, col = "grey30")

  if (verbose) {
    cat(sprintf("\nsurvival curve against Kaplan-Meier (%d refits)\n", n_ok))
    cat(sprintf("largest vertical gap: %.4f\n", gap))
    if (!is.na(outside))
      cat(sprintf("share of the curve outside the envelope: %.1f%%\n",
                  100 * outside))
    cat(status, if (status %in% c("WARN", "FAIL"))
      " -- try another family and compare by AIC" else "", "\n", sep = "")
  }
  invisible(list(grid = grid, km = obs, fitted = fitc, lo = lo, hi = hi,
                 gap = gap, outside = outside, status = status, n_ok = n_ok,
                 levels = lev))
}
