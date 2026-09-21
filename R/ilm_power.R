## ---------------------------------------------------------------------------
## A-priori power, by simulation.
##
## Take a model as the account of how the data arise, generate datasets from
## it at a given size and effect, refit each one, and count how often the term
## of interest is detected. That works for any family this package fits --
## zero-inflated counts, ordinal outcomes, mixed models -- because nothing in
## it depends on a closed-form variance existing.
##
## Two things about a simulated power estimate are routinely dropped, and both
## change what the number means.
##
## IT IS ITSELF AN ESTIMATE. A power of 0.80 from 200 replicates has a standard
## error of 2.8%, so it is 0.74 to 0.86. Reading a required sample size off it
## to three figures is false precision, and reporting the interval is the
## difference between "you need 340 people" and "somewhere between 300 and
## 400, and here is how to narrow that".
##
## SOME FITS DO NOT CONVERGE. Dropping those and dividing by the rest gives
## power CONDITIONAL ON CONVERGENCE, which is not power: a study that fails to
## fit has not detected anything. Both numbers are reported, along with the
## failure rate, and the default counts a failure as a non-detection because
## that is what it is.
##
## References:
##   Arnold, B. F., Hogan, D. R., Colford, J. M. and Hubbard, A. E. (2011).
##     Simulation methods to estimate design power. BMC Medical Research
##     Methodology 11, 94.
##   Green, P. and MacLeod, C. J. (2016). SIMR: an R package for power analysis
##     of generalised linear mixed models by simulation. Methods in Ecology and
##     Evolution 7, 493-498.
## ---------------------------------------------------------------------------

#' Resample a design to a given size, keeping clusters whole
#'
#' Rows are drawn with replacement. When the model has a grouping factor the
#' CLUSTERS are drawn instead, because a mixed design's power depends far more
#' on how many groups there are than on how many rows sit inside each, and
#' resampling rows would quietly hold the number of groups fixed.
#'
#' @keywords internal
#' @noRd
ilm_power_resample <- function(mf, n, group = NULL) {
  if (is.null(group) || !group %in% names(mf)) {
    i <- sample.int(nrow(mf), n, replace = TRUE)
    out <- mf[i, , drop = FALSE]
    rownames(out) <- NULL
    return(out)
  }
  g <- factor(mf[[group]])
  lv <- levels(g)
  idx <- split(seq_len(nrow(mf)), g)
  per <- nrow(mf) / length(lv)
  ng <- max(2L, round(n / per))
  pick <- sample(lv, ng, replace = TRUE)
  rows <- unlist(lapply(pick, function(l) idx[[l]]), use.names = FALSE)
  out <- mf[rows, , drop = FALSE]
  ## a cluster drawn twice has to become two clusters, or the design has
  ## fewer independent groups than it appears to
  out[[group]] <- factor(rep(seq_along(pick),
                             times = vapply(pick, function(l)
                               length(idx[[l]]), 0L)))
  rownames(out) <- NULL
  out
}

#' Power for a term, by simulating from a fitted model
#'
#' Treats a fitted model as the truth, generates datasets of the requested
#' sizes, refits each, and counts how often the named term is detected. It
#' works for any family [ilm_model()] fits, because it never needs a
#' closed-form variance.
#'
#' @section The estimate has a standard error:
#'
#' Power from `sims` replicates is a proportion, so it carries a standard error
#' of `sqrt(p (1 - p) / sims)` -- 2.8% at 0.80 from 200 replicates. The
#' returned interval is a Wilson interval on that, and the sample size implied
#' by a target power is a range rather than a number. Raising `sims` narrows
#' it, and nothing else does.
#'
#' @section Fits that do not converge:
#'
#' A replicate that fails to fit has not detected anything, so `power` counts
#' it as a non-detection. `power_converged` divides by the replicates that
#' worked, which is the number most software reports and is **power
#' conditional on convergence** -- a different and more flattering quantity.
#' When `converged` is below one the two differ, and the gap is the size of the
#' problem rather than something to smooth over.
#'
#' @param object A fitted [ilm_model()] to treat as the truth.
#' @param n Sample sizes to try. Defaults to a spread around the fitted size.
#'   For a model with a grouping factor this is the number of ROWS, and the
#'   number of clusters moves with it.
#' @param term The fixed-effect coefficient to test. Defaults to the first
#'   that is not an intercept.
#' @param effect Values for that coefficient on the LINK scale. Defaults to
#'   the fitted value. Supplying several traces power across effect sizes.
#' @param sims Replicates per cell.
#' @param alpha Two-sided level.
#' @param seed Random seed.
#' @param progress Show a progress bar; see [ilm_progress_arg].
#' @return An object of class `"ilm_power"`: one row per `n` by `effect` cell
#'   with `power`, its Monte Carlo interval, `power_converged` and `converged`.
#' @seealso [ilm_power_n()] to read off the size for a target power,
#'   [plot.ilm_power()] for the curve, [ilm_simulate()] for the generator.
#' @references Arnold, B. F., Hogan, D. R., Colford, J. M. and Hubbard, A. E.
#'   (2011). Simulation methods to estimate design power. *BMC Medical Research
#'   Methodology* 11, 94.
#' @examples
#' set.seed(1); n <- 200
#' d <- data.frame(x = rnorm(n))
#' d$y <- rbinom(n, 1, plogis(-0.5 + 0.5 * d$x))
#' f <- ilm_model(y ~ x, data = d, family = "binomial", verbose = FALSE)
#' ilm_power(f, n = c(200, 400), sims = 50)
#' @export
ilm_power <- function(object, n = NULL, term = NULL, effect = NULL,
                      sims = 200L, alpha = 0.05, seed = 1L, progress = NULL) {
  if (!inherits(object, "ilm_model"))
    stop("`object` must be a fitted ilm_model, not ", class(object)[1],
         call. = FALSE)
  if (is.null(object$formula) || is.null(object$model))
    stop("this fit did not keep its formula and data, so nothing can be ",
         "simulated from it; fit through the formula interface.",
         call. = FALSE)
  b <- stats::coef(object)
  cand <- names(b)[names(b) != "(Intercept)"]
  if (!length(cand))
    stop("the model has no term to find: it is an intercept only.",
         call. = FALSE)
  if (is.null(term)) term <- cand[1L]
  if (!term %in% names(b))
    stop("`term` (", term, ") is not a coefficient. Available: ",
         paste(cand, collapse = ", "), call. = FALSE)
  n0 <- nrow(object$model)
  if (is.null(n)) n <- unique(round(n0 * c(0.5, 1, 2)))
  n <- sort(unique(as.integer(n)))
  if (any(n < 10L)) stop("`n` below 10 is not a study", call. = FALSE)
  if (is.null(effect)) effect <- unname(b[term])
  grid <- expand.grid(n = n, effect = effect, KEEP.OUT.ATTRS = FALSE)

  ## the grouping factor, if any, so clusters are resampled whole
  group <- if (length(object$re)) names(object$re)[1L] else NULL
  fam <- if (is.null(object$family)) "gaussian" else object$family$name
  fo <- object$formula
  resp <- all.vars(fo)[1L]
  wnm <- NULL

  set.seed(seed)
  pb <- ilm_progress(nrow(grid) * sims, progress, "simulating studies")
  tick <- 0L
  res <- vector("list", nrow(grid))
  for (g in seq_len(nrow(grid))) {
    rej <- ok <- logical(sims)
    for (s in seq_len(sims)) {
      d <- ilm_power_resample(object$model, grid$n[g], group)
      bb <- b; bb[term] <- grid$effect[g]
      d[[resp]] <- ilm_power_draw(object, d, bb)
      fit <- try(suppressWarnings(suppressMessages(
        ilm_model(fo, data = d, family = fam, verbose = FALSE,
                  restarts = 1L))), silent = TRUE)
      tick <- tick + 1L; pb$tick(tick)
      if (inherits(fit, "try-error") || !isTRUE(fit$ok)) next
      se <- try(sqrt(diag(as.matrix(suppressWarnings(stats::vcov(fit))))),
                silent = TRUE)
      cf <- stats::coef(fit)
      if (inherits(se, "try-error") || !term %in% names(cf) ||
          !is.finite(se[term]) || se[term] <= 0) next
      ok[s] <- TRUE
      rej[s] <- 2 * stats::pnorm(-abs(unname(cf[term] / se[term]))) < alpha
    }
    nok <- sum(ok)
    p <- mean(rej)                       # failures count as non-detections
    ci <- ilm_wilson(sum(rej), sims, alpha = 0.05)
    res[[g]] <- data.frame(
      n = grid$n[g], effect = grid$effect[g], power = p,
      mc_lower = ci[1L], mc_upper = ci[2L],
      power_converged = if (nok) sum(rej) / nok else NA_real_,
      converged = nok / sims, sims = sims, row.names = NULL)
  }
  pb$done()
  out <- do.call(rbind, res)
  structure(out, class = c("ilm_power", "data.frame"), term = term,
            alpha = alpha, family = fam, formula = fo, grouped = !is.null(group),
            fitted_n = n0)
}

#' Wilson interval for a proportion
#'
#' Rather than the Wald one, which at a power of 0.98 from 200 replicates puts
#' the upper end above 1 and is worst exactly where power estimates live.
#'
#' @keywords internal
#' @noRd
ilm_wilson <- function(x, n, alpha = 0.05) {
  z <- stats::qnorm(1 - alpha / 2)
  p <- x / n
  d <- 1 + z^2 / n
  ctr <- (p + z^2 / (2 * n)) / d
  hw <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / d
  c(max(ctr - hw, 0), min(ctr + hw, 1))
}

#' Draw a response for a simulated design
#'
#' @keywords internal
#' @noRd
ilm_power_draw <- function(object, d, beta) {
  X <- ilm_newX(object, d)$X
  eta <- as.numeric(X %*% beta)
  ## the random effects are redrawn for each simulated study, since a new
  ## study means new clusters rather than the ones that were observed
  if (length(object$re)) {
    for (k in seq_along(object$re)) {
      nm <- names(object$re)[k]
      if (!nm %in% names(d)) next
      gg <- factor(d[[nm]])
      S <- object$Sigma[[k]]
      u <- stats::rnorm(nlevels(gg), 0, sqrt(max(S[1L, 1L], 0)))
      eta <- eta + u[as.integer(gg)]
    }
  }
  fam <- object$family
  if (is.null(fam)) return(eta + stats::rnorm(length(eta)))
  dv <- ilm_disp_vec(object)
  dsp <- if (is.null(dv)) numeric(0) else log(unname(dv[1L]))
  w <- rep(1, length(eta))
  if (isTRUE(object$ordinal)) {
    z <- eta + fam$qfun(stats::runif(length(eta)))
    yi <- as.integer(rowSums(outer(z, as.numeric(object$zeta), `>`))) + 1L
    return(factor(object$ylevels[yi], levels = object$ylevels, ordered = TRUE))
  }
  y <- as.numeric(fam$sim(matrix(eta, ncol = 1L), w, dsp))
  if (!is.null(object$Zzi)) {
    Z <- ilm_zi_design(object$zi_formula, d, colnames(object$Zzi))
    pz <- ilm_zi_p(object, Z)
    y <- ilm_zi_rng(object, exp(eta), if (is.null(dv)) 1 else unname(dv[1L]),
                    pz, y)
  }
  y
}

#' The sample size a target power implies
#'
#' Interpolates the simulated curve, and reports the range the Monte Carlo
#' error allows rather than a single number the simulation cannot support.
#'
#' @param object An [ilm_power()] result.
#' @param target Power to reach.
#' @return A data frame with one row per effect size: the interpolated `n`, and
#'   the `n_lower`/`n_upper` implied by the Monte Carlo interval.
#' @seealso [ilm_power()].
#' @export
ilm_power_n <- function(object, target = 0.8) {
  if (!inherits(object, "ilm_power"))
    stop("`object` must be an ilm_power() result, not ", class(object)[1],
         call. = FALSE)
  d <- as.data.frame(object); class(d) <- "data.frame"
  if (length(unique(d$n)) < 2L)
    stop("at least two sample sizes are needed to interpolate one; rerun ",
         "ilm_power() with several values of `n`.", call. = FALSE)
  cross <- function(x, y) {
    o <- order(x); x <- x[o]; y <- y[o]
    if (all(y < target) || all(y > target)) return(NA_real_)
    stats::approx(y, x, xout = target, ties = "ordered")$y
  }
  do.call(rbind, lapply(split(d, d$effect), function(z)
    data.frame(effect = z$effect[1L], n = cross(z$n, z$power),
               n_lower = cross(z$n, z$mc_upper),
               n_upper = cross(z$n, z$mc_lower),
               row.names = NULL)))
}

#' @export
print.ilm_power <- function(x, ...) {
  cat(sprintf("Simulated power for %s (%s family, alpha = %.3g)\n",
              attr(x, "term"), attr(x, "family"), attr(x, "alpha")))
  cat(sprintf("  %d replicates per cell; the fitted study had %d rows\n",
              x$sims[1L], attr(x, "fitted_n")))
  if (attr(x, "grouped"))
    cat("  clusters are resampled whole, so the number of GROUPS moves with n\n")
  d <- as.data.frame(x); class(d) <- "data.frame"
  d$sims <- NULL
  for (j in c("power", "mc_lower", "mc_upper", "power_converged", "converged"))
    d[[j]] <- round(d[[j]], 3)
  d$effect <- signif(d$effect, 4)
  print(d, row.names = FALSE)
  cat("\n  mc_lower/mc_upper is a Wilson interval on the power ESTIMATE: at\n",
      "  0.80 from ", x$sims[1L], " replicates the standard error is ",
      sprintf("%.3f", sqrt(0.8 * 0.2 / x$sims[1L])),
      ", so a\n  sample size read off this curve is a range. More replicates ",
      "narrow it;\n  nothing else does.\n", sep = "")
  if (any(x$converged < 1, na.rm = TRUE))
    cat("\n  Some replicates did not converge. `power` counts those as ",
        "non-detections,\n  because a study that will not fit has not detected ",
        "anything.\n  `power_converged` divides by the ones that worked, which ",
        "is power\n  CONDITIONAL ON CONVERGENCE and is the more flattering ",
        "number.\n", sep = "")
  invisible(x)
}

#' Power curve
#'
#' @param x An [ilm_power()] result.
#' @param target Draw a line at this power.
#' @param ... Passed to the plotting function.
#' @return `x`, invisibly.
#' @seealso [ilm_power()].
#' @export
plot.ilm_power <- function(x, target = 0.8, ...) {
  d <- as.data.frame(x); class(d) <- "data.frame"
  eff <- unique(d$effect)
  op <- graphics::par(no.readonly = TRUE); on.exit(graphics::par(op))
  graphics::plot(range(d$n), c(0, 1), type = "n", xlab = "sample size",
                 ylab = "power",
                 main = sprintf("power for %s", attr(x, "term")), ...)
  graphics::abline(h = target, lty = 2, col = "grey50")
  for (i in seq_along(eff)) {
    z <- d[d$effect == eff[i], ]
    z <- z[order(z$n), ]
    ## the Monte Carlo band, because the curve is itself estimated
    graphics::polygon(c(z$n, rev(z$n)), c(z$mc_lower, rev(z$mc_upper)),
                      col = grDevices::adjustcolor(i + 1, alpha.f = 0.18),
                      border = NA)
    graphics::lines(z$n, z$power, col = i + 1, lwd = 2)
    graphics::points(z$n, z$power, col = i + 1, pch = 19)
  }
  if (length(eff) > 1L)
    graphics::legend("bottomright", legend = sprintf("effect %.3g", eff),
                     col = seq_along(eff) + 1, lwd = 2, bty = "n")
  invisible(x)
}
