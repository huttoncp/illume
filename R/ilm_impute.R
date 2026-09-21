## ---------------------------------------------------------------------------
## Multiple imputation by chained equations, using ilm_model() as the engine.
##
## One variable at a time: regress it on the others among the rows where it was
## observed, then DRAW replacement values for the rows where it was not. Cycle
## until the fills stop moving, and do the whole thing m times.
##
## The word that matters is DRAW. Filling in the fitted mean is not imputation,
## it is smoothing: it puts every imputed row exactly on the regression surface,
## which shrinks the variance of the completed data and narrows every interval
## computed from it afterwards. A proper draw carries two sources of noise --
## the parameters are not known (so beta is drawn from its sampling
## distribution) and the value is not determined by them (so residual noise is
## added on top). Leaving either out is how imputation quietly breaks coverage.
##
## The m completed data sets are then analysed separately and combined by
## Rubin's rules, which is what turns "we guessed" into an honest interval: the
## spread of the estimates ACROSS imputations is the cost of not having observed
## the values, and it is added to the usual within-imputation variance.
## ---------------------------------------------------------------------------

#' @keywords internal
#' @noRd
ilm_imp_family <- function(v) {
  k <- ilm_var_kind(v)
  switch(k, continuous = "gaussian", binary = "binomial", count = "poisson",
         nominal = "multinomial", NA_character_)
}

## Draw from the posterior predictive of one variable, for the rows where it is
## missing. Parameter uncertainty first, then the response's own randomness.
#' @keywords internal
#' @noRd
ilm_imp_draw <- function(fit, newdata, fam, yobs) {
  b <- suppressWarnings(stats::coef(fit, full = TRUE))
  V <- tryCatch(suppressWarnings(stats::vcov(fit, full = TRUE)),
                error = function(e) NULL)
  ## beta* from its sampling distribution: without this every imputation would
  ## use the same fitted surface and the between-imputation variance would be
  ## only residual noise, understating it
  f2 <- fit
  if (!is.null(V) && length(b) == nrow(V)) {
    ev <- eigen((V + t(V)) / 2, symmetric = TRUE)
    lam <- pmax(ev$values, 0)
    z <- stats::rnorm(length(b))
    bs <- as.numeric(b + ev$vectors %*% (sqrt(lam) * z))
    f2 <- tryCatch(ilm_rebuild(fit, bs), error = function(e) fit)
  }
  mu <- tryCatch(suppressWarnings(
          stats::predict(f2, newdata = newdata, type = "response")),
        error = function(e) NULL)
  if (is.null(mu)) return(NULL)
  nmiss <- nrow(newdata)

  if (fam == "gaussian") {
    mu <- as.numeric(mu)
    ## The RESIDUAL scale, not the marginal one. `sigma()` and `residuals()`
    ## are not available for an ilm_model, and sd(y) is the wrong number --
    ## it is the spread before the predictors explain any of it, so using it
    ## would inject far more noise into each imputation than the model says is
    ## there. $dispersion is the residual standard deviation and agrees with
    ## lm()'s sigma to the printed digits.
    s <- fit$dispersion
    if (is.null(s) || !is.finite(s) || s <= 0)
      s <- stats::sd(yobs, na.rm = TRUE)
    ## and draw it too, rather than fixing it at its estimate: sigma is no more
    ## known than beta is, and a proper imputation carries both
    nu <- max(length(yobs) - length(b), 1L)
    s <- s * sqrt(nu / stats::rchisq(1L, nu))
    return(stats::rnorm(nmiss, mu, s))
  }
  if (fam == "binomial") {
    p <- if (is.matrix(mu)) mu[, ncol(mu)] else as.numeric(mu)
    dr <- stats::rbinom(nmiss, 1L, pmin(pmax(p, 0), 1))
    lv <- levels(factor(yobs))
    return(if (length(lv) == 2L) factor(lv[dr + 1L], levels = lv) else dr)
  }
  if (fam == "poisson")
    return(stats::rpois(nmiss, pmax(as.numeric(mu), 0)))
  if (fam == "multinomial") {
    P <- if (is.matrix(mu)) mu else matrix(mu, nrow = nmiss)
    lv <- colnames(P)
    if (is.null(lv)) lv <- levels(factor(yobs))
    pick <- apply(P, 1L, function(pr) {
      pr[!is.finite(pr) | pr < 0] <- 0
      if (sum(pr) <= 0) return(sample.int(length(pr), 1L))
      sample.int(length(pr), 1L, prob = pr)
    })
    return(factor(lv[pick], levels = levels(factor(yobs))))
  }
  NULL
}

#' Multiple imputation by chained equations
#'
#' Fills in missing values `m` times, each fill drawn rather than fitted, so
#' that the analysis afterwards can price in not having observed them. Each
#' incomplete variable is regressed on the others with [ilm_model()], using the
#' family its own type calls for, and the missing entries are drawn from the
#' resulting predictive distribution.
#'
#' @section Why several, and why drawn:
#'
#' Filling in one set of values and analysing it as though those values had been
#' observed treats a guess as data. The estimates may be fine; the standard
#' errors are not, because nothing in them knows that part of the data set was
#' invented. Intervals come out too narrow by an amount that grows with the
#' proportion missing.
#'
#' Two things fix that. Each imputed value is **drawn** from the predictive
#' distribution -- the coefficients are drawn from their sampling distribution
#' and the response's own randomness is added on top -- so the fills differ from
#' one imputation to the next in the way the unknown values really could. And
#' the `m` completed data sets are analysed separately and pooled by Rubin's
#' rules in [ilm_mi_pool()], where the spread of the estimates across imputations
#' becomes part of the reported uncertainty.
#'
#' `single = TRUE` returns one completed data set. It exists because it is
#' occasionally what is wanted -- a plot, a rough look -- and it warns, because
#' inference computed from it will be overconfident.
#'
#' @section What it assumes:
#'
#' That the data are missing at random given the variables supplied: the chance
#' a value is missing may depend on what is observed, but not on the missing
#' value itself once the observed variables are accounted for. That assumption
#' is **not testable** ([ilm_check_missing()] explains why), and imputation does
#' not make it true. Including variables that predict both the missingness and
#' the missing value makes it more plausible.
#'
#' Imputing an outcome from nothing but the predictors of the model you intend
#' to fit adds no information at all -- the imputations lie on the very surface
#' you are about to estimate. Auxiliary variables are what make imputing an
#' outcome worthwhile.
#'
#' @section Which method, and what it costs:
#'
#' `method = "auto"` uses chained equations and falls back to the low-rank
#' route only when the per-variable regressions cannot be fitted -- when there
#' are more predictors than rows on which a variable was observed. That
#' ordering is not a preference, it is what the ablation says. Hiding known
#' cells and scoring the imputations against them, with a noise floor of 0.500:
#'
#' ```
#'   design                          mean-fill    fcs   lowrank
#'   n=200 p=8  rank 3                   2.289  1.177     1.963
#'   n=200 p=8  full rank                2.739  2.622     3.581
#'   n=400 p=12 rank 4, 30% missing      2.010  1.054     1.683
#'   n=60  p=80 rank 3  (p > n)          1.920      -     0.887
#' ```
#'
#' Where chained equations can be fitted it is clearly better: the per-variable
#' regressions use everything, where a rank-k fit discards whatever falls
#' outside those k components. Where `p > n` it cannot be fitted at all, and
#' the low-rank route gets within twice the noise floor.
#'
#' The middle row is the warning. With no low-rank structure to find, the
#' reconstruction imposes one, and does **worse than filling in column means**.
#' `ilm_impute()` warns when cross-validation picks the largest rank it was
#' offered, which is the signature of that case.
#'
#' Reconstruction accuracy is not the whole story either. It says how close the
#' filled values are, not whether inference computed afterwards is calibrated --
#' single imputation scores well on the first and fails the second. The
#' coverage table below is the one that matters for inference.
#'
#' @section What this is calibrated for:
#'
#' 200 replicates, 400 rows, `y = 0.5x + 0.3z + e`, with `x` made missing three
#' ways. Coverage of a nominal 95% interval for the coefficient on `x`, Monte
#' Carlo error 0.015:
#'
#' ```
#'                        full   complete   multiple    single
#'                        data      cases  imputation  imputation
#'   MCAR, 30% missing
#'     bias            -0.0028    -0.0012     -0.0035    -0.0069
#'     coverage          0.955      0.960       0.970      0.890
#'   MAR on a covariate, 40% missing
#'     bias            -0.0028     0.0001     -0.0017     0.0002
#'     coverage          0.955      0.980       0.970      0.835
#'   MAR on the OUTCOME, 41% missing
#'     bias            -0.0028    -0.0989     -0.0090    -0.0058
#'     coverage          0.955      0.615       0.955      0.790
#' ```
#'
#' Three things to read off it. Complete cases are **unbiased** when missingness
#' depends on a covariate, even at 40% missing -- which is why
#' [ilm_check_missing()] distinguishes that case and tells you not to bother
#' imputing. Complete cases **fail badly** when missingness depends on the
#' outcome: a bias of -0.099 is a fifth of the effect, and coverage collapses to
#' 0.615. Multiple imputation repairs exactly that case, 0.955.
#'
#' And single imputation is the cautionary column. Its point estimates are
#' respectable throughout -- the bias is no worse than multiple imputation's --
#' but its coverage runs 0.790 to 0.890, because nothing in its standard errors
#' knows that part of the data was invented. That is the whole reason `m`
#' defaults to more than one.
#'
#' Pooled coverage sits slightly high, 0.970 where 0.950 is nominal, in the two
#' cells where complete cases were already valid. Rubin's rules are known to be
#' mildly conservative, and erring wide is the right direction for a method
#' whose job is to stop intervals being too narrow.
#'
#' @param data A data frame.
#' @param m Number of imputations. 20 or more is cheap insurance; the classic
#'   advice of 5 dates from when it was not.
#' @param predictors Columns to use as predictors in the imputation models.
#'   Default is every column except those being imputed in that step.
#' @param exclude Columns never to impute and never to use, such as an
#'   identifier.
#' @param maxit Cycles through the variables per imputation.
#' @param method `"auto"` uses chained equations and falls back to the
#'   `"glrm"` fits a generalized low rank model instead, which uses a loss
#'   suited to each column's type and so can impute CATEGORICAL columns,
#'   which `"lowrank"` leaves alone -- see [ilm_glrm()].
#'   low-rank route when they cannot be fitted; `"fcs"` and `"lowrank"` force
#'   one. See the section below for what each costs.
#' @param ncp Rank for the low-rank route. `NULL` chooses it by
#'   cross-validation over held-out observed cells.
#' @param single Return a single completed data frame instead. Warns.
#' @param seed Random seed.
#' @param verbose Narrate progress.
#' @param progress Show a progress bar. Defaults to [interactive()], so a
#'   bar appears when someone is watching and nothing is written in a
#'   script or a knitted document. See [ilm_progress_arg].
#' @return An object of class `"ilm_mids"` holding the `m` completed data sets,
#'   or a data frame when `single = TRUE`.
#' @seealso [ilm_mi_pool()] to analyse them, [ilm_check_missing()] to decide
#'   whether you need to.
#' @references
#' Rubin, D. B. (1987). Multiple Imputation for Nonresponse in Surveys. Wiley.
#'
#' van Buuren, S. and Groothuis-Oudshoorn, K. (2011). mice: Multivariate
#' Imputation by Chained Equations in R. Journal of Statistical Software 45(3).
#'
#' Josse, J. and Husson, F. (2016). missMDA: A Package for Handling Missing
#' Values in Multivariate Data Analysis. Journal of Statistical Software 70(1).
#' The low-rank route follows their multiple imputation PCA: a regularised
#' iterative fit, bootstrapped so the imputations differ, then pooled.
#' @examples
#' set.seed(1); n <- 200
#' d <- data.frame(x = rnorm(n), z = rnorm(n))
#' d$y <- 0.5 * d$x + 0.3 * d$z + rnorm(n)
#' d$x[sample(n, 40)] <- NA
#' imp <- ilm_impute(d, m = 5, seed = 1, verbose = FALSE)
#' imp
#' @export
ilm_impute <- function(data, m = 20L, predictors = NULL, exclude = NULL,
                       maxit = 5L,
                       method = c("auto", "fcs", "lowrank", "glrm"),
                       ncp = NULL, single = FALSE, seed = NULL,
                       verbose = TRUE, progress = NULL) {
  method <- match.arg(method)
  if (!is.data.frame(data))
    stop("`data` must be a data frame; it is ", class(data)[1], call. = FALSE)
  m <- as.integer(m); maxit <- as.integer(maxit)
  if (m < 1L) stop("`m` must be at least 1", call. = FALSE)
  if (!is.null(seed)) set.seed(seed)
  say <- function(...) if (verbose) message(...)

  usable <- setdiff(names(data), exclude)
  inc <- usable[vapply(data[usable], function(v) anyNA(v), TRUE)]
  if (!length(inc)) {
    say("ilm_impute: no missing values to impute")
    return(if (single) data else structure(
      list(imputations = list(data), m = 1L, incomplete = character(),
           data = data, families = character()), class = "ilm_mids"))
  }
  ## a variable illume has no family for cannot be drawn from, so it is left
  ## alone and said so rather than filled with something invented
  fams <- vapply(inc, function(v) ilm_imp_family(data[[v]]), "")
  bad <- inc[is.na(fams)]
  if (length(bad)) {
    say("  cannot impute (no family for their type, left as NA): ",
        paste(bad, collapse = ", "))
    inc <- setdiff(inc, bad); fams <- fams[inc]
  }
  if (!length(inc)) {
    warning("nothing could be imputed", call. = FALSE)
    return(if (single) data else structure(
      list(imputations = list(data), m = 1L, incomplete = character(),
           data = data, families = character()), class = "ilm_mids"))
  }
  ## Chained equations regresses each incomplete variable on all the others,
  ## which needs more complete rows than predictors. Past that point there is
  ## no regression to fit and the old behaviour was to skip the variable
  ## entirely, leaving it NA. A low-rank fit has no such limit: it never fits a
  ## p-predictor regression, only a rank-ncp approximation.
  pred_pool <- if (is.null(predictors)) usable else intersect(predictors, usable)
  num_inc <- inc[vapply(inc, function(v) is.numeric(data[[v]]), TRUE)]
  ## What limits chained equations is the rows on which each variable was
  ## OBSERVED, since the others are filled in before it is regressed on them --
  ## not the rows complete across every column. Counting complete cases instead
  ## makes the constraint look far tighter than it is: 8 columns at 20% missing
  ## leaves only 17% of rows complete while every one of them is usable here.
  n_fit <- if (length(inc))
    min(vapply(inc, function(v) sum(!is.na(data[[v]])), 1L)) else nrow(data)
  n_pred <- length(pred_pool) - 1L
  too_wide <- n_pred >= max(n_fit - 2L, 1L)
  use_lowrank <- method == "lowrank" ||
    (method == "auto" && too_wide && length(num_inc))
  if (method == "fcs" && too_wide)
    warning("there are ", n_pred, " predictors and about ", n_fit,
            " usable rows, so the per-variable regressions cannot be fitted ",
            "and those variables will be left as they are. method = ",
            "\"lowrank\" imputes them from a low-rank approximation instead.",
            call. = FALSE)
  if (isTRUE(single) && m > 1L) m <- 1L
  if (isTRUE(single))
    warning("a single imputation treats the filled-in values as if they had ",
            "been observed, so anything computed from it will have standard ",
            "errors that are too small. Use the default and ilm_mi_pool() for ",
            "inference.", call. = FALSE)

  say("== multiple imputation ==")
  say("  ", length(inc), " variable(s) to impute: ",
      paste(sprintf("%s (%s)", inc, fams[inc]), collapse = ", "))
  say("  ", m, " imputation(s), ", maxit, " cycle(s) each")

  miss_idx <- lapply(inc, function(v) which(is.na(data[[v]])))
  names(miss_idx) <- inc

  one_imputation <- function(i) {
    cur <- data
    ## start from a random observed value, so the first cycle has something to
    ## condition on without that something being a constant
    for (v in inc) {
      obs <- data[[v]][!is.na(data[[v]])]
      cur[[v]][miss_idx[[v]]] <- sample(obs, length(miss_idx[[v]]), replace = TRUE)
    }
    for (it in seq_len(maxit)) {
      for (v in inc) {
        mi <- miss_idx[[v]]
        if (!length(mi)) next
        rhs <- setdiff(pred_pool, v)
        ## drop predictors that are constant or would blow the design up
        rhs <- rhs[vapply(rhs, function(w) {
          z <- cur[[w]]
          if (is.numeric(z)) stats::sd(z, na.rm = TRUE) > 0
          else nlevels(droplevels(factor(z))) > 1L &&
               nlevels(droplevels(factor(z))) <= 20L
        }, TRUE)]
        if (!length(rhs)) next
        tr <- cur[-mi, c(v, rhs), drop = FALSE]
        tr <- tr[stats::complete.cases(tr), , drop = FALSE]
        if (nrow(tr) < 10L) next
        form <- stats::reformulate(rhs, response = v)
        environment(form) <- environment()
        fit <- tryCatch(suppressWarnings(
                 ilm_model(form, data = tr, family = fams[[v]],
                           verbose = FALSE)), error = function(e) NULL)
        if (is.null(fit)) next
        nd <- cur[mi, rhs, drop = FALSE]
        dr <- tryCatch(ilm_imp_draw(fit, nd, fams[[v]], tr[[v]]),
                       error = function(e) NULL)
        if (is.null(dr) || length(dr) != length(mi) || anyNA(dr)) next
        if (is.factor(cur[[v]]) && !is.factor(dr))
          dr <- factor(levels(cur[[v]])[dr + 1L], levels = levels(cur[[v]]))
        cur[[v]][mi] <- dr
      }
    }
    cur
  }

  if (method == "glrm") {
    if (verbose)
      say("  using a generalized low-rank route: a loss per column type, so ",
          "categorical columns are imputed too")
    k <- if (is.null(ncp)) max(1L, min(2L, length(inc))) else as.integer(ncp)
    gcols <- intersect(pred_pool, names(data))
    ## the rank is fitted once on the full data so every replicate answers to
    ## the same structure; only the row weights differ between them
    base <- ilm_glrm(data[gcols], rank = k, progress = FALSE)
    if (verbose)
      say("  rank ", base$rank, ", ridge ", signif(base$lambda, 3),
          if (is.null(ncp)) " (chosen by cross-validation)" else "")
    pb <- ilm_progress(m, progress)
    imps <- vector("list", m)
    for (i in seq_len(m)) {
      w <- tabulate(sample.int(nrow(data), nrow(data), replace = TRUE),
                    nbins = nrow(data))
      fit <- tryCatch(ilm_glrm(data[gcols], rank = k, lambda = base$lambda,
                               weights = w, progress = FALSE),
                      error = function(e) base)
      imps[[i]] <- ilm_glrm_draw(fit, data, inc)
      pb$tick(i)
    }
    pb$done()
    if (isTRUE(single)) return(imps[[1]])
    return(structure(list(imputations = imps, m = m, incomplete = inc,
                          data = data, families = stats::setNames(
                            base$loss[inc], inc),
                          maxit = maxit, method = "glrm", ncp = base$rank),
                     class = "ilm_mids"))
  }

  if (use_lowrank) {
    if (verbose)
      say("  using a low-rank (MIPCA) route: ", n_pred, " predictors against ",
          "about ", n_fit, " usable rows")
    cats <- setdiff(inc, num_inc)
    if (length(cats) && verbose)
      say("  NOT imputed, since a low-rank reconstruction is defined for ",
          "numeric columns only: ", paste(cats, collapse = ", "))
    X <- as.matrix(data[intersect(pred_pool, names(data))[
      vapply(intersect(pred_pool, names(data)), function(v)
        is.numeric(data[[v]]), TRUE)]])
    miss <- is.na(X)
    k <- if (is.null(ncp)) ilm_lowrank_ncp(X, miss, seed = seed) else
      as.integer(ncp)
    if (verbose) say("  rank ", k,
                     if (is.null(ncp)) " (chosen by cross-validation)" else "")
    ## A rank at the ceiling of the search means cross-validation never found a
    ## point where more components stopped helping -- which is what happens when
    ## there is no low-rank structure to find. That is the regime where this
    ## method does real harm: on full-rank data it scored 3.581 against 2.739
    ## for simply filling in column means.
    if (is.null(ncp) && isFALSE(attr(k, "beats_mean")))
      warning("a rank-", as.integer(k), " reconstruction predicts held-out ",
              "cells no better than the column means do (",
              signif(attr(k, "cv_error"), 3), " against ",
              signif(attr(k, "mean_error"), 3),
              "), so these columns have no low-rank structure to exploit. ",
              "Imposing one is worse than not imputing at all; prefer ",
              "method = \"fcs\" wherever it can be fitted.", call. = FALSE)
    k <- as.integer(k)
    pb <- ilm_progress(m, progress)
    imps <- vector("list", m)
    for (i in seq_len(m)) {
      ## the bootstrap, as row weights: a row drawn twice counts twice
      w <- tabulate(sample.int(nrow(X), nrow(X), replace = TRUE),
                    nbins = nrow(X))
      fit <- ilm_lowrank_fit(X, miss, k, w = w)
      cur <- data
      for (v in intersect(num_inc, colnames(X))) {
        mi <- miss[, v]
        ## the reconstruction is the fitted value; the imputation is a DRAW
        ## from around it, which is what keeps the m data sets apart
        cur[[v]][mi] <- fit$filled[mi, v] +
          stats::rnorm(sum(mi), 0, fit$sigma[[v]])
      }
      imps[[i]] <- cur; pb$tick(i)
    }
    pb$done()
    if (isTRUE(single)) return(imps[[1]])
    return(structure(list(imputations = imps, m = m, incomplete = num_inc,
                          data = data, families = stats::setNames(
                            rep("lowrank", length(num_inc)), num_inc),
                          maxit = maxit, method = "lowrank", ncp = k),
                     class = "ilm_mids"))
  }

  pb <- ilm_progress(m, progress)
  imps <- vector("list", m)
  for (i in seq_len(m)) { imps[[i]] <- one_imputation(i); pb$tick(i) }
  pb$done()
  if (isTRUE(single)) return(imps[[1]])
  structure(list(imputations = imps, m = m, incomplete = inc, data = data,
                 families = fams[inc], maxit = maxit, method = "fcs"),
            class = "ilm_mids")
}

#' @export
print.ilm_mids <- function(x, ...) {
  cat("<ilm_mids>", x$m, "imputation(s) of", nrow(x$data), "rows\n")
  if (!length(x$incomplete)) { cat("  nothing was missing\n"); return(invisible(x)) }
  cat("  imputed:\n")
  for (v in x$incomplete)
    cat(sprintf("    %-18s %-12s %d value(s)\n", v, x$families[[v]],
                sum(is.na(x$data[[v]]))))
  cat("\n  Analyse with ilm_mi_pool(), which combines the", x$m,
      "fits by Rubin's rules.\n")
  cat("  Using one of these on its own would understate the uncertainty.\n")
  invisible(x)
}

## ---- Rubin's rules ---------------------------------------------------------

## Named ilm_mi_pool rather than ilm_pool: the latter already exists in
## ilm_parallel.R, where it builds a cluster of worker processes, and six
## diagnostics call it.
#' Fit a model across imputations and pool the results
#'
#' Fits the model to each completed data set and combines them by Rubin's rules.
#' The point estimate is the average; the variance is the average
#' within-imputation variance plus the variance of the estimates *between*
#' imputations, inflated by `1 + 1/m`. That second term is the price of not
#' having observed the values, and leaving it out is exactly what makes single
#' imputation overconfident.
#'
#' Degrees of freedom follow Barnard and Rubin (1999), which corrects the
#' original formula when the complete-data degrees of freedom are small -- the
#' uncorrected version can return more degrees of freedom than the data could
#' possibly supply.
#'
#' @param object An [ilm_impute()] result, or a list of fitted models.
#' @param formula Model formula, when `object` holds data sets.
#' @param ... Passed to [ilm_model()].
#' @return An object of class `"ilm_pooled"`: a coefficient table with
#'   `estimate`, `se`, `df`, `lower`, `upper`, `p_value`, plus `fmi` (the
#'   fraction of information lost to missingness) per coefficient.
#' @references
#' Barnard, J. and Rubin, D. B. (1999). Small-sample degrees of freedom with
#' multiple imputation. Biometrika 86(4).
#' @seealso [ilm_impute()].
#' @examples
#' set.seed(1); n <- 200
#' d <- data.frame(x = rnorm(n), z = rnorm(n))
#' d$y <- 0.5 * d$x + 0.3 * d$z + rnorm(n)
#' d$x[sample(n, 40)] <- NA
#' imp <- ilm_impute(d, m = 5, seed = 1, verbose = FALSE)
#' ilm_mi_pool(imp, y ~ x + z, family = "gaussian")
#' @export
ilm_mi_pool <- function(object, formula = NULL, ...) {
  fits <- if (inherits(object, "ilm_mids")) {
    if (is.null(formula))
      stop("`formula` is needed to fit a model to the imputations",
           call. = FALSE)
    lapply(object$imputations, function(d)
      tryCatch(suppressWarnings(ilm_model(formula, data = d, verbose = FALSE, ...)),
               error = function(e) NULL))
  } else if (is.list(object)) object else
    stop("`object` must be an ilm_impute() result or a list of fits",
         call. = FALSE)
  fits <- Filter(Negate(is.null), fits)
  if (!length(fits)) stop("no imputation could be fitted", call. = FALSE)
  m <- length(fits)
  if (m < 2L) {
    warning("pooling ", m, " fit: the between-imputation variance cannot be ",
            "estimated from one imputation, so this interval is the ",
            "complete-data one and is too narrow.", call. = FALSE)
  }

  B <- lapply(fits, function(f) stats::coef(f))
  nm <- names(B[[1]])
  if (!all(vapply(B, function(b) identical(names(b), nm), TRUE)))
    stop("the imputations did not produce the same coefficients; a factor ",
         "level is probably absent from some of them", call. = FALSE)
  Bm <- do.call(rbind, B)
  Um <- do.call(rbind, lapply(fits, function(f)
    suppressWarnings(diag(as.matrix(stats::vcov(f))))))

  qbar <- colMeans(Bm)                       # the pooled estimate
  ubar <- colMeans(Um)                       # average within-imputation variance
  bvar <- if (m > 1L) apply(Bm, 2L, stats::var) else rep(0, length(qbar))
  tvar <- ubar + (1 + 1 / m) * bvar          # total variance
  se <- sqrt(tvar)

  ## fraction of missing information, and the relative increase in variance
  riv <- ifelse(ubar > 0, (1 + 1 / m) * bvar / ubar, 0)
  lam <- ifelse(tvar > 0, (1 + 1 / m) * bvar / tvar, 0)
  df_old <- ifelse(lam > 0, (m - 1) / lam^2, Inf)
  ## Barnard-Rubin: without this the degrees of freedom can exceed what the
  ## complete data could ever supply
  f1 <- fits[[1]]
  dcom <- if (isTRUE(f1$exact_df) && !is.null(f1$resid_df)) f1$resid_df else
          max(stats::nobs(f1) - length(qbar), 1)
  df_obs <- ((dcom + 1) / (dcom + 3)) * dcom * (1 - lam)
  df <- ifelse(is.finite(df_old), 1 / (1 / df_old + 1 / df_obs), df_obs)
  df <- pmax(df, 1)
  fmi <- (riv + 2 / (df + 3)) / (riv + 1)

  crit <- stats::qt(0.975, df)
  out <- data.frame(term = nm, estimate = unname(qbar), se = unname(se),
                    df = unname(df),
                    lower = unname(qbar - crit * se),
                    upper = unname(qbar + crit * se),
                    p_value = unname(2 * stats::pt(-abs(qbar / se), df)),
                    fmi = unname(fmi), stringsAsFactors = FALSE)
  rownames(out) <- NULL
  structure(out, class = c("ilm_pooled", "data.frame"), m = m, fits = fits)
}

#' @export
print.ilm_pooled <- function(x, ...) {
  m <- attr(x, "m")
  cat("<ilm_pooled>", m, "imputations, combined by Rubin's rules\n\n")
  cat(sprintf("  %-18s %9s %9s %8s %9s %9s %7s\n", "term", "estimate", "se",
              "df", "lower", "upper", "fmi"))
  for (i in seq_len(nrow(x)))
    cat(sprintf("  %-18s %9.4f %9.4f %8.1f %9.4f %9.4f %7.3f\n",
                x$term[i], x$estimate[i], x$se[i], x$df[i], x$lower[i],
                x$upper[i], x$fmi[i]))
  cat("\n  fmi is the fraction of information lost to missingness.\n")
  cat("  Above about 0.5, more imputations are worth having.\n")
  invisible(x)
}

## ---- low-rank imputation, for when a regression per variable cannot be fitted

## Regularised iterative PCA on a matrix with holes (Josse and Husson 2016).
##
## Fill the gaps with column means, take the rank-`ncp` reconstruction, put the
## reconstruction back into the gaps only, and repeat until the fill stops
## moving. The singular values are SHRUNK by the noise variance before
## reconstructing, which is what the "regularised" means and what stops the
## fit chasing the very values it just invented: without it the procedure
## overfits the observed cells and the imputations are too confident.
##
## Row weights carry the bootstrap. A row drawn twice counts twice, which
## perturbs the estimated subspace exactly as refitting on the resample would,
## without having to re-index the missing cells.
#' @keywords internal
#' @noRd
ilm_lowrank_fit <- function(X, miss, ncp, w = NULL, maxit = 200L, tol = 1e-6) {
  n <- nrow(X); p <- ncol(X)
  if (is.null(w)) w <- rep(1, n)
  sw <- sqrt(w / mean(w))
  ctr <- vapply(seq_len(p), function(j) {
    v <- X[!miss[, j], j]; if (length(v)) mean(v) else 0 }, 1)
  sc <- vapply(seq_len(p), function(j) {
    v <- X[!miss[, j], j]; s <- if (length(v) > 1L) stats::sd(v) else 1
    if (!is.finite(s) || s <= 0) 1 else s }, 1)
  Z <- sweep(sweep(X, 2L, ctr, "-"), 2L, sc, "/")
  Z[miss] <- 0                                   # the mean, once centred
  ncp <- max(1L, min(ncp, min(n, p) - 1L))

  old <- Z[miss]; sig2 <- 0
  for (it in seq_len(maxit)) {
    sv <- svd(Z * sw, nu = ncp, nv = ncp)
    d <- sv$d[seq_len(ncp)]
    ## the noise variance, from the discarded directions
    tot <- sum(sv$d^2); kept <- sum(d^2)
    dfres <- max(n * p - length(miss[miss]) - (n + p) * ncp, 1)
    sig2 <- max((tot - kept) / dfres, 0)
    ## Shrink each component by how much of it is signal rather than noise.
    ## As a RATIO in [0, 1], applied to an unweighted projection onto the
    ## weighted subspace -- not by rescaling a weighted reconstruction back.
    ## A bootstrap draw leaves some rows with weight zero, so dividing by the
    ## weights produces Inf for exactly those rows and the next svd() fails
    ## with "infinite or missing values in 'x'".
    f <- pmax(d^2 - sig2, 0) / pmax(d^2, .Machine$double.eps)
    V <- sv$v[, seq_len(ncp), drop = FALSE]
    rec <- Z %*% V %*% diag(f, ncp, ncp) %*% t(V)
    Z[miss] <- rec[miss]
    if (it > 1L && mean(abs(Z[miss] - old)) < tol) break
    old <- Z[miss]
  }
  ## sigma is indexed by column name downstream, so it carries them
  sg <- sqrt(sig2) * sc
  names(sg) <- colnames(X)
  list(filled = sweep(sweep(Z, 2L, sc, "*"), 2L, ctr, "+"),
       sigma = sg, ncp = ncp)
}

## How many components? Too few over-smooths, too many fits the noise, and the
## number is not free -- so it is chosen by holding out observed cells and
## seeing which rank predicts them best, rather than fixed at a guess.
#' @keywords internal
#' @noRd
ilm_lowrank_ncp <- function(X, miss, ncp_max = NULL, folds = 3L, seed = NULL) {
  n <- nrow(X); p <- ncol(X)
  hi <- min(ncp_max %||% (min(n, p) - 1L), min(n, p) - 1L, 8L)
  if (hi < 2L) return(1L)
  obs <- which(!miss)
  if (length(obs) < 40L) return(min(2L, hi))
  if (!is.null(seed)) set.seed(seed)
  fold <- sample(rep_len(seq_len(folds), length(obs)))
  err <- numeric(hi); base <- 0; nb <- 0L
  for (k in seq_len(hi)) {
    e <- 0; m <- 0L
    for (f in seq_len(folds)) {
      hold <- obs[fold == f]
      Xh <- X; mh <- miss; Xh[hold] <- NA; mh[hold] <- TRUE
      ## the column means WITHOUT the held-out cells: what a rank-0 fit, i.e.
      ## no structure at all, would predict. Scored once, on the first pass.
      if (k == 1L) {
        cm <- colMeans(Xh, na.rm = TRUE)
        pred <- cm[((hold - 1L) %/% nrow(X)) + 1L]
        base <- base + sum((pred - X[hold])^2, na.rm = TRUE)
        nb <- nb + length(hold)
      }
      fit <- tryCatch(ilm_lowrank_fit(Xh, mh, k, maxit = 30L),
                      error = function(z) NULL)
      if (is.null(fit)) next
      e <- e + sum((fit$filled[hold] - X[hold])^2); m <- m + length(hold)
    }
    err[k] <- if (m) e / m else Inf
  }
  if (all(!is.finite(err))) return(min(2L, hi))
  k <- which.min(err)
  ## Does a low-rank fit beat no structure at all? If it does not, these
  ## columns have none to find, and imposing one is worse than the column
  ## means -- measured at 3.581 against 2.739 on full-rank data.
  structure(k, beats_mean = nb > 0 && err[k] < base / nb,
            cv_error = err[k], mean_error = if (nb) base / nb else NA_real_)
}
