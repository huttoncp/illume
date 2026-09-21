## ---------------------------------------------------------------------------
## Causal mediation.
##
## How much of a treatment's effect on an outcome runs THROUGH a mediator, and
## how much goes around it. The quantities are counterfactual:
##
##   ACME(t)  = E[ Y(t, M(1)) - Y(t, M(0)) ]   the indirect effect: what
##              changes when the mediator moves as it would under treatment,
##              while the treatment itself is held at t
##   ADE(t)   = E[ Y(1, M(t)) - Y(0, M(t)) ]   the direct effect: what changes
##              when the treatment moves and the mediator is held where it
##              would have been under t
##   total    = ACME(1) + ADE(0) = ACME(0) + ADE(1)
##
## The Baron-Kenny product of coefficients is the special case of ACME when the
## outcome model is linear and has no treatment-by-mediator interaction. It is
## the case people fit and it is often not the case they have: with an
## interaction the indirect effect differs by treatment arm, and ACME(0) and
## ACME(1) are two different numbers that a single product cannot be. With a
## binary outcome the product is not the effect at all.
##
## This estimates the counterfactual definitions by simulation, so it works for
## any family ilm_model() fits, and reproduces the product exactly where the
## product is right -- which is the test.
##
## THE ASSUMPTION. Decomposing an effect needs no unmeasured confounding of
## treatment-outcome, of treatment-mediator, AND of MEDIATOR-OUTCOME. The last
## is the one that is rarely plausible and never testable: the mediator was not
## randomised. Nothing in the data can check it, so ilm_mediate_sens() asks the
## only useful question instead -- how strong would such confounding have to be
## before the conclusion changed.
##
## References:
##   Imai, K., Keele, L. and Tingley, D. (2010). A general approach to causal
##     mediation analysis. Psychological Methods 15, 309-334.
##   VanderWeele, T. J. (2015). Explanation in Causal Inference. OUP.
## ---------------------------------------------------------------------------

#' Draw a parameter vector from a fit's sampling distribution
#'
#' @keywords internal
#' @noRd
ilm_med_draw <- function(fit) {
  b <- stats::coef(fit)
  V <- as.matrix(suppressWarnings(stats::vcov(fit)))
  V <- V[seq_along(b), seq_along(b), drop = FALSE]
  as.numeric(b + ilm_msqrt(V) %*% stats::rnorm(length(b)))
}

#' Predicted mean from a design and a coefficient vector
#'
#' @keywords internal
#' @noRd
ilm_med_mu <- function(fit, nd, beta) {
  X <- ilm_newX(fit, nd)$X
  eta <- as.numeric(X %*% beta)
  li <- if (!is.null(fit$family)) fit$family$linkinv else identity
  li(eta)
}

#' How much of an effect runs through a mediator
#'
#' Splits the effect of a treatment on an outcome into the part that operates
#' through a mediator and the part that does not, using the counterfactual
#' definitions rather than a product of coefficients.
#'
#' @section What it estimates:
#'
#' For each simulation draw, parameters are drawn from both models' sampling
#' distributions, the mediator is predicted under each treatment value, and the
#' outcome is predicted under each of the four combinations of treatment and
#' mediator. Averaging over the observed units gives
#'
#' * `ACME` -- the average causal mediation (indirect) effect, reported at each
#'   treatment value because they differ whenever the outcome model has a
#'   treatment-by-mediator interaction;
#' * `ADE` -- the average direct effect, likewise;
#' * `total` -- their sum, which is the same whichever pairing is used;
#' * `prop_mediated` -- ACME over total, which is a ratio of estimates and
#'   behaves badly when the total is near zero. It is reported and should be
#'   read with that in mind.
#'
#' With a linear outcome model and no interaction this reduces exactly to the
#' product of the treatment-to-mediator and mediator-to-outcome coefficients.
#' The tests check that it does, to 1e-10.
#'
#' @section The assumption that cannot be checked:
#'
#' A decomposition needs no unmeasured confounding of the treatment-outcome
#' relationship, of the treatment-mediator relationship, and of the
#' **mediator-outcome** relationship. The first two can be addressed by design;
#' the third rarely can, because the mediator was not randomised -- whatever
#' makes someone's mediator high may also make their outcome high for reasons
#' that have nothing to do with the treatment.
#'
#' Nothing in the data tests this. [ilm_mediate_sens()] asks the answerable
#' question instead: how strong would such confounding have to be before the
#' indirect effect went away.
#'
#' @param model_m A fitted [ilm_model()] for the mediator.
#' @param model_y A fitted [ilm_model()] for the outcome, including both the
#'   treatment and the mediator among its predictors.
#' @param treat Name of the treatment variable, in both models.
#' @param mediator Name of the mediator: the response of `model_m` and a
#'   predictor in `model_y`.
#' @param control_value,treat_value The two treatment values to contrast.
#'   Defaults are 0 and 1, or the first two levels of a factor.
#' @param sims Simulation draws.
#' @param level Confidence level for the percentile intervals.
#' @param seed Random seed.
#' @param progress Show a progress bar; see [ilm_progress_arg].
#' @return An object of class `"ilm_mediate"`: a data frame of estimates and
#'   intervals, with the draws attached.
#' @seealso [ilm_mediate_sens()] for the untestable assumption,
#'   [ilm_dag_model()] for whether the adjustment sets are available at all.
#' @references Imai, K., Keele, L. and Tingley, D. (2010). A general approach
#'   to causal mediation analysis. *Psychological Methods* 15, 309-334.
#' @examples
#' set.seed(1); n <- 400
#' d <- data.frame(x = rbinom(n, 1, 0.5), c = rnorm(n))
#' d$m <- 0.3 + 0.7 * d$x + 0.2 * d$c + rnorm(n)
#' d$y <- 1 + 0.4 * d$x + 0.6 * d$m + 0.1 * d$c + rnorm(n)
#' fm <- ilm_model(m ~ x + c, data = d, family = "gaussian", verbose = FALSE)
#' fy <- ilm_model(y ~ x + m + c, data = d, family = "gaussian",
#'                 verbose = FALSE)
#' ilm_mediate(fm, fy, treat = "x", mediator = "m", sims = 200)
#' @export
ilm_mediate <- function(model_m, model_y, treat, mediator,
                        control_value = NULL, treat_value = NULL,
                        sims = 1000L, level = 0.95, seed = 1L,
                        progress = NULL) {
  for (nm in c("model_m", "model_y")) {
    o <- get(nm)
    if (!inherits(o, "ilm_model"))
      stop("`", nm, "` must be a fitted ilm_model, not ", class(o)[1],
           call. = FALSE)
    if (length(o$re))
      stop("`", nm, "` has random effects. A counterfactual mediator has to ",
           "be predicted for each unit, and with a random effect that is a ",
           "conditional quantity while the effects here are population ",
           "ones; fit both models without the random effect and cluster the ",
           "standard errors instead, or treat the groups as fixed.",
           call. = FALSE)
  }
  dm <- model_m$model; dy <- model_y$model
  if (nrow(dm) != nrow(dy))
    stop("the two models were fitted to ", nrow(dm), " and ", nrow(dy),
         " rows. A mediation decomposition compares counterfactuals for the ",
         "same units, so both have to see the same ones -- drop the ",
         "incomplete rows before fitting, or give both models the same ",
         "na.action and variables.", call. = FALSE)
  if (!treat %in% names(dm) || !treat %in% names(dy))
    stop("`treat` (", treat, ") must be a predictor in BOTH models; it is ",
         "missing from ", if (!treat %in% names(dm)) "the mediator model"
         else "the outcome model", ".", call. = FALSE)
  ry <- names(dy)[1L]
  if (!mediator %in% names(dy))
    stop("`mediator` (", mediator, ") is not a predictor in the outcome ",
         "model, so nothing can run through it.", call. = FALSE)
  rm_ <- names(dm)[1L]
  if (!identical(rm_, mediator))
    stop("`mediator` (", mediator, ") is not the response of the mediator ",
         "model, which is ", rm_, ".", call. = FALSE)
  if (mediator %in% names(dm)[-1L])
    stop("the mediator appears among its own predictors.", call. = FALSE)

  tv <- dm[[treat]]
  if (is.null(control_value) || is.null(treat_value)) {
    if (is.factor(tv)) {
      if (nlevels(tv) < 2L) stop("`treat` has one level", call. = FALSE)
      control_value <- levels(tv)[1L]; treat_value <- levels(tv)[2L]
    } else { control_value <- 0; treat_value <- 1 }
  }
  ## An interaction is what makes ACME differ by arm. Its absence is an
  ## assumption, not a simplification, and it should be visible.
  has_int <- any(grepl(paste0("(^|:)", treat, "(:|$)"),
                       model_y$term_labels) &
                 grepl(paste0("(^|:)", mediator, "(:|$)"),
                       model_y$term_labels))

  set.seed(seed)
  n <- nrow(dm)
  mk <- function(d, var, val) { d[[var]] <- ilm_med_set(d[[var]], val); d }
  d0 <- mk(dm, treat, control_value); d1 <- mk(dm, treat, treat_value)
  y0 <- mk(dy, treat, control_value); y1 <- mk(dy, treat, treat_value)

  out <- matrix(NA_real_, sims, 5L,
                dimnames = list(NULL, c("acme0", "acme1", "ade0", "ade1",
                                        "total")))
  pb <- ilm_progress(sims, progress, "simulating counterfactuals")
  for (s in seq_len(sims)) {
    bm <- ilm_med_draw(model_m); by <- ilm_med_draw(model_y)
    ## the counterfactual mediator is a DRAW, not a fitted mean: Y(t, M(t'))
    ## averages over the mediator's own variation, and using its mean instead
    ## would understate the indirect effect wherever the outcome model is not
    ## linear in it
    m0 <- ilm_med_rng(model_m, d0, bm)
    m1 <- ilm_med_rng(model_m, d1, bm)
    pred <- function(base, mv) {
      nd <- base; nd[[mediator]] <- mv
      mean(ilm_med_mu(model_y, nd, by))
    }
    y00 <- pred(y0, m0); y01 <- pred(y0, m1)
    y10 <- pred(y1, m0); y11 <- pred(y1, m1)
    out[s, ] <- c(y01 - y00, y11 - y10, y10 - y00, y11 - y01,
                  (y01 - y00) + (y11 - y01))
    pb$tick(s)
  }
  pb$done()

  a <- (1 - level) / 2
  est <- colMeans(out)
  ci <- apply(out, 2L, stats::quantile, c(a, 1 - a), na.rm = TRUE)
  pm <- out[, "acme0"] / out[, "total"]
  res <- data.frame(
    effect = c("ACME (control)", "ACME (treated)", "ADE (control)",
               "ADE (treated)", "Total effect", "Proportion mediated"),
    estimate = c(est, mean(pm, na.rm = TRUE)),
    lower = c(ci[1L, ], stats::quantile(pm, a, na.rm = TRUE)),
    upper = c(ci[2L, ], stats::quantile(pm, 1 - a, na.rm = TRUE)),
    p = c(vapply(seq_len(5L), function(j)
      2 * min(mean(out[, j] <= 0), mean(out[, j] >= 0)), 0), NA_real_),
    row.names = NULL, stringsAsFactors = FALSE)
  structure(res, class = c("ilm_mediate", "data.frame"), draws = out,
            sims = sims, level = level, interaction = has_int,
            treat = treat, mediator = mediator,
            values = c(control_value, treat_value), n = n,
            model_m = model_m, model_y = model_y)
}

#' Set a variable to a counterfactual value, keeping its type
#'
#' @keywords internal
#' @noRd
ilm_med_set <- function(x, value) {
  if (is.factor(x)) factor(rep(as.character(value), length(x)),
                           levels = levels(x))
  else rep(as.numeric(value), length(x))
}

#' Draw a counterfactual mediator for every unit
#'
#' @keywords internal
#' @noRd
ilm_med_rng <- function(fit, nd, beta) {
  X <- ilm_newX(fit, nd)$X
  eta <- matrix(as.numeric(X %*% beta), ncol = 1L)
  fam <- fit$family
  if (is.null(fam)) return(as.numeric(eta))
  d <- ilm_disp_vec(fit)
  dsp <- if (is.null(d)) numeric(0) else log(unname(d[1L]))
  w <- if (is.null(fit$weights)) rep(1, nrow(eta)) else fit$weights
  as.numeric(fam$sim(eta, w, dsp))
}

#' @export
print.ilm_mediate <- function(x, ...) {
  cat(sprintf("Causal mediation: %s -> %s -> outcome\n",
              attr(x, "treat"), attr(x, "mediator")))
  cat(sprintf("  %s = %s against %s, %d units, %d simulation draws\n",
              attr(x, "treat"), attr(x, "values")[2L],
              attr(x, "values")[1L], attr(x, "n"), attr(x, "sims")))
  d <- as.data.frame(x); class(d) <- "data.frame"
  d$estimate <- signif(d$estimate, 4); d$lower <- signif(d$lower, 4)
  d$upper <- signif(d$upper, 4); d$p <- signif(d$p, 3)
  print(d, row.names = FALSE)
  if (!attr(x, "interaction"))
    cat("\n  The outcome model has no ", attr(x, "treat"), ":",
        attr(x, "mediator"), " interaction, so the two ACMEs are equal by\n",
        "  construction rather than by evidence. Adding one lets the indirect\n",
        "  effect differ by arm, and is worth doing before concluding it does not.\n",
        sep = "")
  cat("\n  All of this rests on there being no unmeasured confounding of the\n",
      "  MEDIATOR and the outcome, which the data cannot check because the\n",
      "  mediator was not randomised. ilm_mediate_sens() asks how strong such\n",
      "  confounding would have to be to remove the indirect effect.\n", sep = "")
  invisible(x)
}

## ---------------------------------------------------------------------------
## How strong would the confounding have to be?
##
## Posit an unmeasured U, standardised, entering both models:
##
##     M = a X + lambda_m U + e_M
##     Y = c X + b M + lambda_y U + e_Y
##
## Omitting U biases the mediator-to-outcome coefficient by the usual
## omitted-variable amount. The regression of Y on M already controls for X, so
## the relevant variation in M is its residual variation, and
##
##     b_observed = b_true + lambda_y * lambda_m / var(M | X, U)
##
## where var(M | X) = lambda_m^2 + var(e_M) is what the fitted mediator model
## reports as its residual variance. Since ACME = a * b in the linear case and
## `a` is unaffected -- the treatment IS randomised or adjusted, which is the
## assumption not in question -- the corrected indirect effect is
##
##     ACME(lambda_m, lambda_y) = ACME_observed
##                                - a * lambda_y * lambda_m / s2_m
##
## with s2_m the mediator model's residual variance. That is an exact
## consequence of the posited model rather than an approximation, and the test
## file checks it by generating data WITH a U, hiding it, and confirming that
## the correction at the true (lambda_m, lambda_y) recovers the true ACME.
## ---------------------------------------------------------------------------

#' How strong would unmeasured confounding have to be?
#'
#' The one assumption a mediation analysis cannot check is that nothing
#' unmeasured causes both the mediator and the outcome. This does not test it,
#' because nothing can. It asks the answerable question instead: given an
#' unmeasured confounder of stated strength, what would the indirect effect
#' have been, and how strong would it have to be for that effect to vanish.
#'
#' @section Reading the result:
#'
#' `lambda_m` and `lambda_y` are the confounder's coefficients in the mediator
#' and outcome models, with the confounder standardised, so they are on the
#' scale of those models' own coefficients: a `lambda_y` of 0.5 means an
#' unmeasured variable with the same pull on the outcome as a predictor whose
#' coefficient is 0.5.
#'
#' The `threshold` is the product `lambda_m * lambda_y` at which the corrected
#' indirect effect reaches zero. Compare it against the products of the
#' coefficients you DID measure: if a confounder half as strong as your
#' strongest covariate would overturn the finding, the finding is fragile, and
#' if it would take one twice as strong as anything observed, it is not.
#'
#' @param object An [ilm_mediate()] result on linear models.
#' @param lambda_m,lambda_y Grids of confounder strengths.
#' @return An object of class `"ilm_mediate_sens"`: the grid of corrected
#'   indirect effects, with the zero-crossing threshold attached.
#' @seealso [ilm_mediate()].
#' @references VanderWeele, T. J. (2015). *Explanation in Causal Inference*.
#'   Oxford University Press, chapter 3.
#' @examples
#' set.seed(1); n <- 400
#' d <- data.frame(x = rbinom(n, 1, 0.5))
#' d$m <- 0.3 + 0.7 * d$x + rnorm(n)
#' d$y <- 1 + 0.4 * d$x + 0.6 * d$m + rnorm(n)
#' fm <- ilm_model(m ~ x, data = d, family = "gaussian", verbose = FALSE)
#' fy <- ilm_model(y ~ x + m, data = d, family = "gaussian", verbose = FALSE)
#' md <- ilm_mediate(fm, fy, treat = "x", mediator = "m", sims = 200)
#' ilm_mediate_sens(md)
#' @export
ilm_mediate_sens <- function(object,
                             lambda_m = seq(0, 1.5, length.out = 16L),
                             lambda_y = seq(0, 1.5, length.out = 16L)) {
  if (!inherits(object, "ilm_mediate"))
    stop("`object` must be an ilm_mediate() result, not ", class(object)[1],
         call. = FALSE)
  mm <- attr(object, "model_m"); my <- attr(object, "model_y")
  fam <- function(f) if (is.null(f$family)) "gaussian" else f$family$name
  if (fam(mm) != "gaussian" || fam(my) != "gaussian")
    stop("this correction is exact for linear models and is not derived for ",
         "others; the mediator model is ", fam(mm), " and the outcome model ",
         fam(my), ". Report ilm_mediate() with the assumption stated rather ",
         "than a sensitivity analysis that does not apply to it.",
         call. = FALSE)

  tr <- attr(object, "treat")
  bm <- stats::coef(mm)
  a <- bm[grep(paste0("^", tr), names(bm))[1L]]
  if (is.na(a))
    stop("could not find the treatment's coefficient in the mediator model.",
         call. = FALSE)
  s2m <- unname(mm$dispersion[[1L]])^2
  acme <- object$estimate[object$effect == "ACME (control)"]

  g <- expand.grid(lambda_m = lambda_m, lambda_y = lambda_y)
  g$acme <- acme - unname(a) * g$lambda_m * g$lambda_y / s2m
  ## the product at which it reaches zero
  thr <- if (abs(unname(a)) > 0) acme * s2m / unname(a) else NA_real_
  structure(g, class = c("ilm_mediate_sens", "data.frame"),
            observed = acme, threshold = thr, a = unname(a), s2m = s2m,
            model_y = my)
}

#' @export
print.ilm_mediate_sens <- function(x, ...) {
  thr <- attr(x, "threshold")
  cat("Sensitivity of the indirect effect to unmeasured mediator-outcome confounding\n")
  cat(sprintf("  observed ACME: %.4f\n", attr(x, "observed")))
  cat(sprintf("  it reaches zero when lambda_m * lambda_y = %.4f\n", thr))
  cat(sprintf("  so a confounder equally strong in both models would need\n  lambda = %.3f\n",
              sqrt(abs(thr))))
  ## the honest comparison: against the coefficients that WERE measured
  b <- stats::coef(attr(x, "model_y"))
  b <- b[names(b) != "(Intercept)"]
  if (length(b)) {
    cat(sprintf("\n  For scale, the outcome model's own coefficients run from %.3f to %.3f.\n",
                min(abs(b)), max(abs(b))))
    cat("  If a confounder that strong is plausible, the indirect effect is not\n",
        "  established; if it is larger than anything you measured, it is.\n",
        sep = "")
  }
  invisible(x)
}
