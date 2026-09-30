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
#' function in `z = (log t - eta) / scale`, and a flexible (Royston-Parmar)
#' model one in its baseline spline plus `eta`, so the curve comes straight
#' from the linear predictor.
#'
#' @section Which groups:
#' In a model with random effects a curve is for some group, and `groups` says
#' which, as it does for [predict.ilm_model()]:
#' * `"fitted"` (the default) draws each row at **its own group's** estimated
#'   effects. A row whose group the fit has not seen, or `newdata` without the
#'   grouping columns, is an error naming the two choices below.
#' * `"typical"` draws each row with every random effect at zero: the curve of
#'   a **typical group**. With no `newdata` the curve is at a typical covariate
#'   pattern, which belongs to no group, and this is what it gives.
#' * `"population"` averages the **survival curve itself** over the
#'   distribution of random effects, by quadrature at each time: the curve for
#'   the population of groups as a whole. It is flatter than the typical
#'   group's curve, not the curve at zero, since survival is not linear in the
#'   linear predictor.
#'
#' Intervals for a typical group are computed on the complementary log-log
#' scale for the accelerated failure time families, which is linear in `eta`,
#' and on the linear predictor for a flexible model. For `"fitted"` and
#' `"population"` they are percentile intervals from `nsim` joint draws of the
#' coefficients -- and, for `"fitted"`, of the groups' own effects, as
#' [predict.ilm_model()] draws them. Uncertainty in an accelerated failure time
#' family's scale parameter is not included, so its intervals are slightly
#' narrow, most in the tail beyond the last observed event.
#'
#' @param object A fitted `"ilm_model"` with a time-to-event family.
#' @param newdata Covariate patterns, one row each. With none, the curve is
#'   drawn at the median of each numeric predictor and the commonest level of
#'   each factor, for a typical group.
#' @param times Times at which to evaluate. With none, a grid spanning the
#'   observed follow-up.
#' @param conf Confidence level.
#' @param groups `"fitted"` (the default), `"typical"` or `"population"`; see
#'   "Which groups". Only a model with random effects or a correlation over
#'   time has a choice to make.
#' @param nsim Draws behind the intervals of `"fitted"` and `"population"`.
#' @param seed Random seed for those draws.
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
ilm_survival <- function(object, newdata = NULL, times = NULL, conf = 0.95,
                         groups = c("fitted", "typical", "population"),
                         nsim = 200L, seed = 1L) {
  ilm_rng_restore(seed)                  # the user's random stream, put back on exit
  if (!inherits(object, "ilm_model"))
    stop("`object` must be a fitted ilm_model object, not ", class(object)[1],
         call. = FALSE)
  if (!isTRUE(object$family$aft) && is.null(object$rp))
    stop("a survival curve needs a time-to-event family. This model is ",
         object$family$name, "; use family = \"weibull\", \"lognormal\", ",
         "\"loglogistic\", or \"rp\" for a flexible baseline.",
         call. = FALSE)
  given <- !missing(groups)
  groups <- ilm_groups_arg(groups, c("fitted", "typical", "population"), given,
                           "ilm_survival()")
  if (is.null(newdata)) {
    ## a typical covariate pattern belongs to no group
    if (given && identical(groups, "fitted"))
      stop("with no `newdata` the curve is at a typical covariate pattern, which ",
           "belongs to no group. Give `newdata` with each row's group, or use ",
           "groups = \"typical\" or groups = \"population\".", call. = FALSE)
    if (!given) groups <- "typical"
    newdata <- ilm_typical_row(object)
  }
  if (is.null(times)) {
    y <- as.numeric(object$y)
    times <- seq(min(y[y > 0]), max(y), length.out = 100L)
  }
  times <- sort(unique(times[times > 0]))
  gk <- which(vapply(object$re, function(e) !identical(e$kind, "basis"), TRUE))
  ## with nothing of a group's own, every choice is the typical group's curve
  if (!length(gk) && is.null(object$ar)) groups <- "typical"
  if (identical(groups, "typical")) {
    if (!is.null(object$rp)) return(ilm_rp_survival(object, newdata, times, conf))
    return(ilm_aft_survival(object, newdata, times, conf))
  }
  ilm_survival_groups(object, newdata, times, conf, groups, gk, nsim, seed)
}

## A typical group's curve for an accelerated failure time family, with its
## interval on the complementary log-log scale.
#' @keywords internal
#' @noRd
ilm_aft_survival <- function(object, newdata, times, conf) {
  pr <- stats::predict(object, newdata = newdata, se.fit = TRUE, groups = "typical")
  eta <- log(as.matrix(if (is.list(pr)) pr$fit else pr)[, 1])
  se_eta <- if (is.list(pr) && !is.null(pr$se.fit)) {
    ## predict() returns the standard error on the response scale, and the
    ## curve needs it on the linear-predictor scale it was computed on
    as.matrix(pr$se.fit)[, 1] / as.matrix(if (is.list(pr)) pr$fit else pr)[, 1]
  } else rep(0, length(eta))
  sc <- unname(object$dispersion[[1]])
  crit <- stats::qnorm(1 - (1 - conf) / 2)
  nm <- object$family$name
  out <- lapply(seq_along(eta), function(i) {
    S <- ilm_aft_surv(nm, (log(times) - eta[i]) / sc)
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

## an accelerated failure time family's survivor function in z
#' @keywords internal
#' @noRd
ilm_aft_surv <- function(nm, z)
  switch(nm, lognormal = stats::pnorm(-z), loglogistic = 1 / (1 + exp(z)),
         weibull = exp(-exp(z)))

## A curve at each row's own group, or averaged over the groups (item 213).
## The survival at time t for a linear predictor e is S(t | e): for an
## accelerated failure time family f((log t - e) / scale), for a flexible model
## Sfun(B(t) gamma + e). "fitted" adds each row's own group's effects to e;
## "population" averages S(t | e + b) over b ~ N(0, the row's latent
## variance) by quadrature -- the average of the curve, never the curve of an
## averaged mean. Intervals: percentiles over joint draws.
#' @keywords internal
#' @noRd
ilm_survival_groups <- function(object, newdata, times, conf, groups, gk, nsim, seed) {
  rp <- object$rp
  nd <- ilm_newX(object, newdata)
  Xn <- nd$X; n <- nrow(Xn)
  off <- if (is.null(nd$offset)) 0 else as.numeric(nd$offset)
  keep <- if (is.null(rp)) seq_len(ncol(object$X))
          else setdiff(seq_len(ncol(object$X)), rp$cols)
  if (ncol(Xn) != length(keep))
    stop("the new data give ", ncol(Xn), " covariate columns where the fit ",
         "has ", length(keep), call. = FALSE)
  Bs <- if (!is.null(rp)) ilm_rcs(log(times), rp$knots)
  sc <- if (is.null(rp)) unname(object$dispersion[[1]])
  nm <- object$family$name
  ## S(t | e) at every time, one column per linear predictor
  surv_at <- function(beta, e) {
    if (is.null(rp))
      vapply(e, function(ei) ilm_aft_surv(nm, (log(times) - ei) / sc),
             numeric(length(times)))
    else {
      base <- as.vector(Bs %*% beta[rp$cols])
      vapply(e, function(ei) object$family$surv(base + ei), numeric(length(times)))
    }
  }
  pl <- if (identical(groups, "fitted")) ilm_fitted_place(object, newdata)
  sdrow <- NULL
  if (identical(groups, "population")) {
    has_ar <- !is.null(object$ar) && !is.null(object$Sigma[["ar"]])
    ar_el <- if (has_ar && identical(object$ar$type, "rw1")) ilm_rw_elapsed(object, newdata)
    tms <- ilm_latent_terms(object, gk, newdata)
    sdrow <- sqrt(ilm_latent_var(object, tms, has_ar, ar_el, n))
  }
  ## the curves, one column per row, for one set of parameters
  curves <- function(beta, bvec = NULL, bar = NULL) {
    e <- as.vector(Xn %*% beta[keep]) + off
    if (!is.null(pl)) e <- e + ilm_fitted_shift(object, pl, n, bvec, bar)[, 1L]
    if (is.null(sdrow)) return(surv_at(beta, e))
    ## the average of the curve over each row's latent normal, time by time
    vapply(seq_len(n), function(i)
      vapply(seq_along(times), function(k)
        ilm_normal_expect(function(x) surv_at(beta, x)[k, ], e[i], sdrow[i]), 0),
      numeric(length(times)))
  }
  est <- curves(as.vector(object$beta[, 1L]))
  ## joint draws: the coefficients, and for "fitted" the groups' own effects
  blk <- intersect(c("beta", if (!is.null(pl)) c("bvec", "B_ar")),
                   names(object$obj$env$par))
  dr <- tryCatch(ilm_draws(object, nsim = nsim, seed = seed, blocks = blk, natural = FALSE),
                 error = function(e) e)
  lo <- hi <- matrix(NA_real_, nrow(est), ncol(est))
  if (inherits(dr, "error")) {
    warning("no intervals: ", conditionMessage(dr), call. = FALSE)
  } else {
    w <- rownames(dr$draws)
    acc <- array(NA_real_, c(nrow(est), ncol(est), nsim))
    for (s in seq_len(nsim))
      acc[, , s] <- curves(dr$draws[w == "beta", s],
                           if (any(w == "bvec")) dr$draws[w == "bvec", s],
                           if (any(w == "B_ar")) dr$draws[w == "B_ar", s])
    a <- (1 - conf) / 2
    lo <- apply(acc, 1:2, stats::quantile, probs = a, names = FALSE, na.rm = TRUE)
    hi <- apply(acc, 1:2, stats::quantile, probs = 1 - a, names = FALSE, na.rm = TRUE)
  }
  do.call(rbind, lapply(seq_len(n), function(i)
    data.frame(row = i, time = times, surv = est[, i], lower = lo[, i], upper = hi[, i])))
}

## The covariate pattern an effect plot would use: median for a number, the
## commonest level for a factor.
#' @keywords internal
#' @noRd
ilm_typical_row <- function(object) {
  mf <- ilm_data(object)
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
#' @param progress Show a progress bar. Defaults to [interactive()], so a
#'   bar appears when someone is watching and nothing is written in a
#'   script or a knitted document. See [illumex::ilm_progress_arg].
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
                              progress = NULL,
                              colour = "#2C7FB8", fill = "grey85",
                              alpha = NULL, size = 1, main = NULL,
                              verbose = TRUE) {
  ilm_rng_restore(seed)                  # the user's random stream, put back on exit
  if (!inherits(object, "ilm_model"))
    stop("`object` must be a fitted ilm_model object, not ", class(object)[1],
         call. = FALSE)
  if (!isTRUE(object$family$aft) && is.null(object$rp))
    stop("a survival curve needs a time-to-event family. This model is ",
         object$family$name, "; use family = \"weibull\", \"lognormal\", ",
         "\"loglogistic\", or \"rp\" for a flexible baseline.",
         call. = FALSE)
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
  mf <- ilm_data(object)
  fitc <- vapply(lev, function(L) {
    nd <- if (is.null(by) || is.null(mf)) ilm_typical_row(object) else {
      s <- gvec == L
      r <- ilm_typical_row(object)
      ## honour the level being drawn when it is itself a predictor
      for (cn in names(r)) if (cn %in% names(mf) && all(mf[[cn]][s][1] == mf[[cn]][s]))
        r[[cn]] <- mf[[cn]][s][1]
      r
    }
    ilm_survival(object, newdata = nd, times = grid, groups = "typical")$surv
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
                          where = environment(), progress = progress)
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
