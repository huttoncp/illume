## ---------------------------------------------------------------------------
## Scenario projection.
##
## "If we rolled this out at the higher dose, what would happen?" is a
## prediction question, and the obvious way to answer it -- build a row with
## the predictors set to the scenario's values and predict -- answers a
## DIFFERENT question from the one being asked.
##
## That row is a representative unit: what happens to someone who looks like
## THIS. The question behind a decision to scale something is what happens to
## OUR POPULATION: set the intervention for every observed unit, leave each
## one's own covariates alone, predict, average. The two coincide under an
## identity link and part company under every other one, because the mean of a
## non-linear function is not that function of the mean, and the gap grows with
## the spread of the covariates.
##
## The second is standardisation, or g-computation. When the adjustment set is
## valid it is also the causal estimand, which is why ilm_dag_model() and this
## belong to the same workflow.
##
## EXTRAPOLATION is the other thing a scenario quietly does. A model will
## return a number for a dose nobody received, and the interval around it will
## look no different from one inside the data, because the interval reflects
## uncertainty about the COEFFICIENTS and not about whether the functional form
## survives out there. So each scenario is checked against the observed support
## -- each value against its own range, and the COMBINATION against how far it
## sits from the nearest real unit.
##
## References:
##   Robins, J. (1986). A new approach to causal inference in mortality studies.
##     Mathematical Modelling 7, 1393-1512.
##   Hernan, M. A. and Robins, J. M. (2020). Causal Inference: What If. CRC.
## ---------------------------------------------------------------------------

#' A representative row: numeric predictors at their mean, factors at the mode
#'
#' @keywords internal
#' @noRd
ilm_scen_reference <- function(mf) {
  as.data.frame(lapply(mf, function(x) {
    if (is.factor(x)) factor(names(sort(table(x), decreasing = TRUE))[1L],
                             levels = levels(x))
    else if (is.logical(x)) mean(x, na.rm = TRUE) > 0.5
    else if (is.numeric(x)) mean(x, na.rm = TRUE)
    else x[1L]
  }), stringsAsFactors = FALSE)
}

#' Is a scenario inside the data, as values and as a combination?
#'
#' @keywords internal
#' @noRd
ilm_scen_support <- function(mf, scen, vars) {
  out <- list(outside = character(0), distance = NA_real_,
              typical = NA_real_, far = FALSE)
  for (v in vars) {
    x <- mf[[v]]; s <- scen[[v]]
    if (is.factor(x)) {
      if (!as.character(s) %in% levels(x))
        out$outside <- c(out$outside, sprintf("%s = %s is not a level seen",
                                              v, s))
    } else if (is.numeric(x)) {
      r <- range(x, na.rm = TRUE)
      ## A round number just past the edge -- dose = 0 where the smallest
      ## observed is 0.013 -- is not extrapolation in any sense a reader
      ## cares about, and flagging it every time is how a warning gets
      ## ignored. The tolerance is 2% of the observed spread.
      tol <- 0.02 * diff(r)
      if (s < r[1L] - tol || s > r[2L] + tol)
        out$outside <- c(out$outside,
                         sprintf("%s = %.4g is outside the observed %.4g to %.4g",
                                 v, s, r[1L], r[2L]))
    }
  }
  ## The combination. Every value can sit inside its own range while the point
  ## sits in a corner no unit occupies -- a 25-year-old with 40 years of
  ## service. Standardised nearest-neighbour distance answers that, compared
  ## against how far apart real units are from each other.
  num <- vars[vapply(vars, function(v) is.numeric(mf[[v]]) &&
                       !is.factor(mf[[v]]), TRUE)]
  if (length(num) >= 2L) {
    M <- as.matrix(mf[num])
    ctr <- colMeans(M, na.rm = TRUE)
    sdv <- apply(M, 2L, stats::sd, na.rm = TRUE)
    sdv[!is.finite(sdv) | sdv <= 0] <- 1
    Z <- sweep(sweep(M, 2L, ctr, "-"), 2L, sdv, "/")
    p <- (unlist(scen[num]) - ctr) / sdv
    d <- sqrt(rowSums(sweep(Z, 2L, p, "-")^2))
    out$distance <- min(d, na.rm = TRUE)
    ## how far a real unit typically sits from its own nearest neighbour,
    ## on a subsample so this stays cheap
    idx <- if (nrow(Z) > 400L) sample.int(nrow(Z), 400L) else seq_len(nrow(Z))
    nn <- vapply(idx, function(i) {
      dd <- sqrt(rowSums(sweep(Z, 2L, Z[i, ], "-")^2))
      min(dd[-i], na.rm = TRUE)
    }, 0)
    out$typical <- stats::quantile(nn, 0.99, na.rm = TRUE)
    out$far <- is.finite(out$distance) && is.finite(out$typical) &&
      out$distance > out$typical
  }
  out
}

#' What the model says would happen under specified scenarios
#'
#' Sets named predictors to chosen values and reports the outcome the model
#' implies, with an interval, for each combination. Built for the question a
#' client asks before scaling something: what would this do to our population.
#'
#' @section Which population:
#'
#' `over = "sample"` (the default) sets the named predictors for **every
#' observed unit**, leaves each unit's other covariates as they are, predicts,
#' and averages. That is standardisation, and it answers "what would happen to
#' this population". When the adjustment set is valid it is also the causal
#' estimand.
#'
#' `over = "reference"` builds a single row with the named predictors at their
#' scenario values and everything else at its mean or modal value, and predicts
#' once. That answers "what would happen to a unit that looks like this".
#'
#' They are different questions and they give different answers under any
#' non-linear link, because the average of a prediction is not the prediction
#' at the average. The printed output says which was used.
#'
#' @section Extrapolation:
#'
#' A model returns a number for a dose nobody received, and the interval around
#' it looks like any other, because the interval carries uncertainty about the
#' coefficients rather than about whether the functional form holds out there.
#' Each scenario is therefore checked twice: each value against that
#' predictor's observed range, and the whole combination against how far it
#' sits from the nearest real unit, in standardised units, compared with how
#' far real units sit from each other. A flagged scenario is not refused -- it
#' is reported, because sometimes extrapolating is the point and the reader
#' should know they are.
#'
#' @param object A fitted [ilm_model()].
#' @param ... Named predictors and the values to set them to, for instance
#'   `dose = c(0, 10, 20), age = 40`. Every combination is used.
#' @param over `"sample"` or `"reference"`; see above.
#' @param sims Parameter draws for the interval.
#' @param level Confidence level.
#' @param contrast Compare scenarios against one another: `FALSE`, `"first"`
#'   (each against the first) or `"pairwise"`.
#' @param seed Random seed.
#' @param progress Show a progress bar; see [illumex::ilm_progress_arg].
#' @return An object of class `"ilm_scenario"`: a data frame of the grid with
#'   `estimate`, `lower` and `upper`, plus any contrasts.
#' @seealso [ilm_ame()] for the effect of a one-unit change rather than a named
#'   scenario, [ilm_emmeans()] for group means, [ilm_dag_model()] for whether
#'   the adjustment set licenses a causal reading.
#' @references Hernán, M. A. and Robins, J. M. (2020). *Causal Inference: What
#'   If*. Chapman & Hall/CRC.
#' @examples
#' set.seed(1); n <- 400
#' d <- data.frame(dose = runif(n, 0, 20), age = rnorm(n, 50, 10))
#' d$y <- rbinom(n, 1, plogis(-3 + 0.12 * d$dose + 0.03 * d$age))
#' f <- ilm_model(y ~ dose + age, data = d, family = "binomial",
#'                verbose = FALSE)
#' ilm_scenario(f, dose = c(0, 10, 20), sims = 200, contrast = "first")
#' @export
ilm_scenario <- function(object, ..., over = c("sample", "reference"),
                         sims = 1000L, level = 0.95, contrast = FALSE,
                         seed = 1L, progress = NULL) {
  over <- match.arg(over)
  if (!inherits(object, "ilm_model"))
    stop("`object` must be a fitted ilm_model, not ", class(object)[1],
         call. = FALSE)
  if (object$C > 1L)
    stop("a multinomial fit predicts a vector of probabilities per row, not ",
         "one number, so a scenario has no single value to report; use ",
         "predict() with newdata.", call. = FALSE)
  spec <- list(...)
  if (!length(spec) || is.null(names(spec)) || any(names(spec) == ""))
    stop("name the predictors to set, for instance dose = c(0, 10, 20).",
         call. = FALSE)
  mf <- object$model
  miss <- setdiff(names(spec), names(mf))
  if (length(miss))
    stop("not in the model: ", paste(miss, collapse = ", "),
         ". The model has ", paste(setdiff(names(mf), names(mf)[1L]),
                                   collapse = ", "), ".", call. = FALSE)
  grid <- expand.grid(spec, stringsAsFactors = FALSE,
                      KEEP.OUT.ATTRS = FALSE)

  ## what is being held fixed, and at what
  others <- setdiff(names(mf)[-1L], names(spec))
  ref <- ilm_scen_reference(mf)

  set.seed(seed)
  mixed <- length(object$re) > 0L
  nsc <- nrow(grid)
  draws <- matrix(NA_real_, sims, nsc)
  ## the base data each scenario is applied to
  base <- if (over == "sample") mf else ref
  pb <- ilm_progress(sims, progress, "projecting scenarios")
  ## random-effect draws, so a mixed model's answer is population-averaged
  ## rather than the value for a cluster whose effect happens to be zero
  reb <- if (mixed) ilm_scen_re_draws(object, 200L) else NULL
  for (s in seq_len(sims)) {
    b <- ilm_med_draw(object)
    for (j in seq_len(nsc)) {
      nd <- base
      for (v in names(spec)) nd[[v]] <- ilm_med_set(nd[[v]], grid[j, v])
      draws[s, j] <- ilm_scen_mean(object, nd, b, reb)
    }
    pb$tick(s)
  }
  pb$done()

  a <- (1 - level) / 2
  res <- grid
  res$estimate <- colMeans(draws)
  res$lower <- apply(draws, 2L, stats::quantile, a, na.rm = TRUE)
  res$upper <- apply(draws, 2L, stats::quantile, 1 - a, na.rm = TRUE)

  ## support, one scenario at a time
  sup <- lapply(seq_len(nsc), function(j)
    ilm_scen_support(mf, grid[j, , drop = FALSE], names(spec)))
  res$outside <- vapply(sup, function(z)
    length(z$outside) > 0L || isTRUE(z$far), TRUE)

  ct <- NULL
  if (!identical(contrast, FALSE)) {
    contrast <- match.arg(as.character(contrast), c("first", "pairwise"))
    prs <- if (contrast == "first")
      cbind(1L, seq_len(nsc)[-1L]) else t(utils::combn(nsc, 2L))
    if (nrow(prs)) {
      lab <- function(j) paste(sprintf("%s=%s", names(spec),
                                       unlist(grid[j, names(spec)])),
                               collapse = ", ")
      dd <- draws[, prs[, 2L], drop = FALSE] - draws[, prs[, 1L], drop = FALSE]
      ct <- data.frame(
        contrast = vapply(seq_len(nrow(prs)), function(i)
          paste(lab(prs[i, 2L]), "-", lab(prs[i, 1L])), ""),
        estimate = colMeans(dd),
        lower = apply(dd, 2L, stats::quantile, a, na.rm = TRUE),
        upper = apply(dd, 2L, stats::quantile, 1 - a, na.rm = TRUE),
        p = apply(dd, 2L, function(z)
          2 * min(mean(z <= 0), mean(z >= 0))),
        row.names = NULL, stringsAsFactors = FALSE)
    }
  }

  structure(res, class = c("ilm_scenario", "data.frame"), draws = draws,
            contrasts = ct, over = over, held = others, reference = ref,
            support = sup, level = level, mixed = mixed,
            vars = names(spec), sims = sims)
}

#' Draws of the total random-effect shift, for a population-averaged mean
#'
#' @keywords internal
#' @noRd
ilm_scen_re_draws <- function(object, ndraw) {
  gk <- which(vapply(object$re, function(e) e$kind != "basis", TRUE))
  if (!length(gk)) return(NULL)
  rowSums(vapply(gk, function(k) {
    S <- object$Sigma[[k]]
    as.numeric(matrix(stats::rnorm(ndraw), ndraw, 1L) *
                 sqrt(max(S[1L, 1L], 0)))
  }, numeric(ndraw)))
}

#' The average predicted outcome under one scenario
#'
#' @keywords internal
#' @noRd
ilm_scen_mean <- function(object, nd, beta, reb) {
  X <- ilm_newX(object, nd)$X
  eta <- as.numeric(X %*% beta)
  li <- if (is.null(object$family)) identity else object$family$linkinv
  mu <- if (is.null(reb)) li(eta) else
    ## average the inverse link over the random-effect distribution, which is
    ## what makes this a population mean rather than the value for a cluster
    ## whose effect is exactly zero
    rowMeans(vapply(reb, function(u) li(eta + u), numeric(length(eta))))
  if (!is.null(object$Zzi)) {
    Z <- ilm_zi_design(object$zi_formula, nd, colnames(object$Zzi))
    mu <- ilm_zi_mean(object, mu, Z)
  }
  mean(mu)
}

#' @export
print.ilm_scenario <- function(x, ...) {
  cat(sprintf("Scenario projection (%s)\n",
              if (attr(x, "over") == "sample")
                "averaged over the observed units" else
                "one representative unit"))
  if (attr(x, "over") == "sample")
    cat("  every unit set to the scenario's values, its own covariates kept\n")
  else
    cat("  a single row: the scenario's values, everything else at its mean\n",
        "  or modal value. This is a unit that looks like this, NOT the\n",
        "  population -- use over = \"sample\" for that.\n", sep = "")
  held <- attr(x, "held")
  if (length(held))
    cat("  held as observed: ", paste(utils::head(held, 8L), collapse = ", "),
        if (length(held) > 8L) sprintf(" and %d more", length(held) - 8L) else "",
        "\n", sep = "")
  if (attr(x, "mixed"))
    cat("  averaged over the random effects, so this is a population mean\n")
  d <- as.data.frame(x); class(d) <- "data.frame"
  for (j in c("estimate", "lower", "upper")) d[[j]] <- signif(d[[j]], 4)
  names(d)[names(d) == "outside"] <- "extrapolating"
  print(d, row.names = FALSE)
  ct <- attr(x, "contrasts")
  if (!is.null(ct)) {
    cat("\n  Differences between scenarios\n")
    q <- ct
    for (j in c("estimate", "lower", "upper")) q[[j]] <- signif(q[[j]], 4)
    q$p <- signif(q$p, 3)
    print(q, row.names = FALSE)
  }
  sup <- attr(x, "support")
  bad <- which(vapply(sup, function(z)
    length(z$outside) > 0L || isTRUE(z$far), TRUE))
  if (length(bad)) {
    cat("\n  Outside the data\n")
    for (i in bad) {
      for (m in sup[[i]]$outside) cat("    row ", i, ": ", m, "\n", sep = "")
      if (isTRUE(sup[[i]]$far))
        cat(sprintf("    row %d: the COMBINATION is %.2f standardised units from the\n               nearest observed unit, against %.2f for a typical one\n",
                    i, sup[[i]]$distance, sup[[i]]$typical))
    }
    cat("  The interval there carries uncertainty about the coefficients, not\n",
        "  about whether the model's shape survives outside the data. Nothing\n",
        "  in the number tells the reader that; this does.\n", sep = "")
  }
  invisible(x)
}
