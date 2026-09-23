## ---------------------------------------------------------------------------
## The DAG-guided workflow: graph in, fitted models and an effect table out,
## narrating each step.
##
## One rule governs the whole thing and is worth stating before any code. The
## GRAPH fixes the mean structure -- the exposure and its adjustment set -- and
## nothing in here searches over it. Only the ERROR structure is adjusted when a
## diagnostic fails. That split is what keeps the reported intervals meaning
## what they say: searching over the mean structure and then reporting the
## winner's p-values is post-selection inference, and it would quietly undo the
## nominal coverage the simulation studies establish. A DAG is a
## pre-registration device, and this function treats it as one.
##
## Where several adjustment sets are admissible, all of them are fitted and
## reported side by side. If the graph is right they should agree, so the spread
## across them is a sensitivity analysis that costs nothing but compute.
## ---------------------------------------------------------------------------

## Columns that look like a grouping factor: repeated values, enough distinct
## levels to estimate a variance, and not so many that every row is its own
## group.
##
## Only variables the GRAPH DOES NOT MENTION are eligible. A variable in the
## graph has a causal role, and a confounder belongs in the mean structure as a
## fixed effect rather than in the error structure; anything outside the graph
## can only be structure in how the data were collected, which is exactly what
## a random effect is for.
#' @keywords internal
#' @noRd
ilm_dag_clusters <- function(data, exclude, min_levels = 5L, min_per = 2) {
  cand <- setdiff(names(data), exclude)
  out <- list()
  for (nm in cand) {
    v <- data[[nm]]
    if (is.numeric(v) && !is.integer(v) && length(unique(v)) > 50L) next
    if (inherits(v, "Date") || inherits(v, "POSIXt")) next
    u <- unique(stats::na.omit(v))
    nl <- length(u)
    if (nl < min_levels) next
    per <- sum(!is.na(v)) / nl
    if (per < min_per) next
    out[[nm]] <- list(levels = nl, per = per)
  }
  ## the most repeated structure first: that is the one most likely to matter
  out[order(-vapply(out, `[[`, 1, "per"))]
}

## Fit, diagnose, and if the spread is wrong, refit with a model for it.
## Deliberately narrow: the only remediation is to the error structure, and each
## step taken is recorded so the path is auditable.
#' @keywords internal
#' @noRd
ilm_dag_remediate <- function(fit, form, data, family, verbose, ...) {
  steps <- character()
  cv <- tryCatch(suppressWarnings(
          ilm_check_variance(fit, plot = FALSE, verbose = FALSE)),
        error = function(e) NULL)
  if (!is.null(cv) && cv$status %in% c("WARN", "FAIL") &&
      family %in% c("gaussian", "nbinom", "weibull", "lognormal", "loglogistic")) {
    if (verbose)
      message("      spread check: ", cv$status,
              sprintf(" (rho = %.3f, p = %s)", cv$trend, format(cv$p_trend)),
              " -- refitting with dispformula = ~ mu")
    ## Match the estimator of the fit being remediated. The two differ only in
    ## the variance structure, not in X, so the AIC comparison below stays
    ## valid under REML -- which it would not be if the fixed effects moved.
    alt <- tryCatch(suppressWarnings(
             ilm_model(form, data = data, family = family,
                       dispformula = ~ mu, verbose = FALSE,
                       reml = isTRUE(fit$reml), ...)),
           error = function(e) NULL)
    if (!is.null(alt)) {
      better <- tryCatch(stats::AIC(alt) < stats::AIC(fit), error = function(e) FALSE)
      if (isTRUE(better)) {
        steps <- c(steps, sprintf(
          "modelled the dispersion (~ mu) after a %s spread check; AIC %.1f -> %.1f",
          cv$status, stats::AIC(fit), stats::AIC(alt)))
        if (verbose) message("      kept it: AIC ", round(stats::AIC(fit), 1),
                             " -> ", round(stats::AIC(alt), 1))
        fit <- alt
      } else {
        steps <- c(steps, sprintf(
          "a dispersion model did not improve on the constant-variance fit (AIC %.1f vs %.1f), so it was not kept",
          stats::AIC(alt), stats::AIC(fit)))
        if (verbose) message("      did not help; kept the simpler fit")
      }
    }
  } else if (!is.null(cv)) {
    steps <- c(steps, paste0("spread check: ", cv$status))
    if (verbose) message("      spread check: ", cv$status)
  }
  list(fit = fit, steps = steps)
}

#' Fit the model a causal graph implies
#'
#' Takes a DAG and a data frame and does what an analyst would: checks the graph
#' against the data, finds what has to be adjusted for, picks a response
#' distribution, looks for grouping structure, fits, runs the diagnostics, and
#' fixes the error structure if an assumption fails. It narrates each step.
#'
#' @section What is fixed and what is searched:
#'
#' **The graph fixes the mean structure.** The exposure and its adjustment set
#' are determined by the causal question, and nothing here searches over them.
#' That is the point of supplying a DAG: it commits to the specification before
#' the data are looked at, so the intervals reported at the end mean what they
#' say. Searching over which covariates to include and then reporting the
#' winner's p-values is post-selection inference, and no amount of care
#' elsewhere repairs it.
#'
#' **Only the error structure is adjusted**, and only when a diagnostic asks
#' for it -- a model for the dispersion when the spread is not constant, a
#' random effect for grouping the graph does not speak to. These change what the
#' standard errors are, not which effect is being estimated. Every step taken is
#' recorded in `$steps` and printed as it happens.
#'
#' @section Several adjustment sets:
#'
#' If more than one minimal set is admissible, every one is fitted and reported.
#' They estimate the same quantity, so if the graph is right they should agree;
#' the spread across them is a sensitivity analysis that comes free. A set that
#' disagrees markedly is worth more attention than any single point estimate.
#'
#' @section When the effect is not identified:
#'
#' If no set of measured variables suffices, that is the result, and this
#' function returns it rather than fitting something that cannot answer the
#' question. It is the most useful thing a DAG can tell you, and it can only be
#' said before the modelling, not after.
#'
#' @param dag An [ilm_dag()], or anything [ilm_dag()] accepts.
#' @param data A data frame.
#' @param exposure,outcome Variable names; taken from the graph if declared.
#' @param family Response distribution. Inferred from the outcome when `NULL`
#'   (or `"auto"`), by the same rules as `ilm_model(family = "auto")`, and the
#'   inference is stated rather than assumed.
#' @param cluster Grouping variable(s) for random intercepts. `NULL` looks for
#'   them among variables the graph does not mention; `character(0)` fits none.
#' @param auto_error Adjust the error structure when a diagnostic asks.
#' @param test_dag Test the graph's implied conditional independencies first.
#' @param max_sets Stop after this many adjustment sets.
#' @param verbose Narrate each step.
#' @param reml Estimate the variance components by restricted maximum
#'   likelihood. Defaults to `TRUE` here, unlike [ilm_model()], because the
#'   graph fixed the adjustment set before any data were looked at: the fixed
#'   effects are not being selected, so the one thing REML forbids -- comparing
#'   likelihoods across different fixed structures -- never arises, and its
#'   unbiased variance components are simply better. Applied to gaussian
#'   responses only: for any other family the restricted likelihood is only
#'   approximate (see [ilm_model()]), so this route stays with maximum
#'   likelihood there.
#' @param ... Passed to [ilm_model()].
#' @return An object of class `"ilm_dag_model"`: the graph, the sets, the fits,
#'   an `effects` table with one row per set, the DAG test, and `steps`.
#' @seealso [ilm_adjust_sets()], [ilm_dag_test()], [ilm_model()].
#' @examples
#' set.seed(1); n <- 300
#' z <- rnorm(n); x <- 0.5 * z + rnorm(n); y <- 0.4 * x + 0.6 * z + rnorm(n)
#' d <- data.frame(x = x, y = y, z = z)
#' g <- ilm_dag("dag { x [exposure] ; y [outcome] ; z -> x -> y ; z -> y }")
#' ilm_dag_model(g, d, verbose = FALSE)
#' @export
ilm_dag_model <- function(dag, data, exposure = NULL, outcome = NULL,
                          family = NULL, cluster = NULL, auto_error = TRUE,
                          test_dag = TRUE, max_sets = 8L, verbose = TRUE,
                          reml = TRUE, ...) {
  g <- ilm_dag(dag)
  if (!is.data.frame(data))
    stop("`data` must be a data frame; it is ", class(data)[1], call. = FALSE)
  x <- if (is.null(exposure)) g$exposure else as.character(exposure)
  y <- if (is.null(outcome))  g$outcome  else as.character(outcome)
  if (length(x) != 1L || length(y) != 1L)
    stop("one `exposure` and one `outcome` are needed; the graph declares ",
         length(g$exposure), " and ", length(g$outcome),
         ". Pass them to ilm_dag() or to this function.", call. = FALSE)
  steps <- character()
  say <- function(...) if (verbose) message(...)

  ## ---- 1. the graph against the data ---------------------------------------
  say("== DAG-guided analysis ==")
  say("[1/6] graph and data")
  declared <- setdiff(g$nodes, g$latent)
  present <- intersect(declared, names(data))
  absent <- setdiff(declared, names(data))
  say("  ", length(g$nodes), " variables in the graph; ", length(present),
      " of the ", length(declared), " observed ones found in the data")
  if (length(g$latent))
    say("  unobserved by declaration: ", paste(g$latent, collapse = ", "))
  if (length(absent)) {
    say("  NOT FOUND in the data, so treated as unobserved: ",
        paste(absent, collapse = ", "))
    steps <- c(steps, paste0("graph variables absent from the data: ",
                             paste(absent, collapse = ", ")))
  }
  for (v in c(x, y))
    if (!v %in% names(data))
      stop("the ", if (v == x) "exposure" else "outcome", " `", v,
           "` is not a column of `data`. Columns are: ",
           paste(utils::head(names(data), 12), collapse = ", "), call. = FALSE)

  ## ---- 2. does the data agree with the graph? ------------------------------
  dt <- NULL
  say("[2/6] does the data agree with the graph?")
  if (test_dag) {
    dt <- ilm_dag_test(g, data, verbose = verbose)
    steps <- c(steps, paste0("graph vs data: ", attr(dt, "verdict")))
  } else say("  skipped (test_dag = FALSE)")

  ## ---- 3. identification ---------------------------------------------------
  say("[3/6] identification")
  say("  effect of `", x, "` on `", y, "`")
  sets <- ilm_adjust_sets(g, x, y, observed = intersect(names(data), g$nodes))
  if (!length(sets)) {
    say("  NO adjustment set: with the variables measured here, no regression ",
        "identifies this effect.")
    say("  This is a result, not a failure -- the graph says the data cannot ",
        "answer the question. Measure a confounder, or change the graph.")
    steps <- c(steps, "not identifiable from the measured variables")
    return(structure(list(dag = g, exposure = x, outcome = y, family = NULL,
                          sets = list(), fits = list(), effects = NULL,
                          dag_test = dt, steps = steps, identified = FALSE),
                     class = "ilm_dag_model"))
  }
  if (length(sets) > max_sets) {
    say("  ", length(sets), " minimal adjustment sets; keeping the ", max_sets,
        " smallest (max_sets)")
    sets <- sets[seq_len(max_sets)]
  }
  say("  ", length(sets), " minimal adjustment set",
      if (length(sets) == 1L) "" else "s", ":")
  for (s in sets)
    say("    {", if (length(s)) paste(s, collapse = ", ") else "empty", "}")
  steps <- c(steps, paste0(length(sets), " minimal adjustment set(s)"))

  ## ---- 4. response and error structure -------------------------------------
  say("[4/6] response and error structure")
  wf <- ilm_workflow_family(family, data[[y]], y, say)
  fam <- wf$family
  if (wf$inferred) steps <- c(steps, paste0("family inferred as ", fam))

  re <- character()
  if (is.null(cluster)) {
    cl <- ilm_dag_clusters(data, exclude = union(g$nodes, c(x, y)))
    if (length(cl)) {
      re <- names(cl)[1]
      say("  grouping structure the graph does not mention: `", re, "` (",
          cl[[1]]$levels, " levels, ", round(cl[[1]]$per, 1),
          " observations each) -> (1 | ", re, ")")
      if (!fam %in% c("gaussian"))
        say("    NOTE: with a non-identity link a random effect changes what ",
            "the coefficients mean -- they become conditional on the group ",
            "rather than population-averaged. Pass cluster = character(0) if ",
            "you want the marginal quantity.")
      steps <- c(steps, paste0("added (1 | ", re, ") for grouping outside the graph"))
      if (length(cl) > 1L)
        say("    (also candidates: ", paste(utils::head(names(cl)[-1], 4),
                                            collapse = ", "), ")")
    } else say("  no grouping structure found outside the graph")
  } else {
    re <- as.character(cluster)
    if (length(re)) say("  random intercepts as supplied: ",
                        paste(re, collapse = ", "))
  }

  ## ---- 5. fitting ----------------------------------------------------------
  say("[5/6] fitting ", length(sets), " model",
      if (length(sets) == 1L) "" else "s",
      " (the mean structure is fixed by the graph)")
  fits <- vector("list", length(sets))
  rows <- vector("list", length(sets))
  for (i in seq_along(sets)) {
    z <- sets[[i]]
    ## column names become formula text, so each is quoted (see ilm_names.R)
    rhs <- c(ilm_bq(c(x, z)), if (length(re)) sprintf("(1 | %s)", ilm_bq(re)))
    form <- stats::reformulate(rhs, response = as.name(y))
    environment(form) <- environment()
    say("  set ", i, " of ", length(sets), ": ",
        paste(deparse(form), collapse = " "))
    ## The DAG fixed the adjustment set before any data were looked at, so the
    ## fixed effects are not up for selection here and REML is the right
    ## estimator: unbiased variance components, and no likelihood comparison
    ## across fixed structures that it would invalidate. Outside a linear model
    ## the restricted likelihood is only approximate, so it applies only where
    ## it is exact.
    use_reml <- isTRUE(reml) && identical(fam, "gaussian")
    fit <- tryCatch(suppressWarnings(
             ilm_model(form, data = data, family = fam, verbose = FALSE,
                       reml = use_reml, ...)),
           error = function(e) structure(list(msg = conditionMessage(e)),
                                         class = "ilm_failed"))
    if (inherits(fit, "ilm_failed")) {
      say("    did not fit: ", fit$msg)
      steps <- c(steps, paste0("set ", i, " did not fit: ", fit$msg))
      fits[[i]] <- NULL; next
    }
    if (auto_error) {
      rr <- ilm_dag_remediate(fit, form, data, fam, verbose, ...)
      fit <- rr$fit
      steps <- c(steps, paste0("set ", i, ": ", rr$steps))
    }
    fits[[i]] <- fit

    ct <- ilm_coef_table(fit)
    k <- ilm_term_cols(fit, x)
    crit <- if (isTRUE(fit$exact_df))
              stats::qt(0.975, fit$resid_df) else stats::qnorm(0.975)
    rows[[i]] <- if (!length(k)) NULL else data.frame(
      set = i,
      adjusted_for = if (length(z)) paste(z, collapse = " + ") else "(nothing)",
      term = rownames(ct)[k],
      estimate = ct[k, 1], se = ct[k, 2],
      lower = ct[k, 1] - crit * ct[k, 2],
      upper = ct[k, 1] + crit * ct[k, 2],
      p_value = ct[k, 4], stringsAsFactors = FALSE)
  }

  ## ---- 6. the answer -------------------------------------------------------
  say("[6/6] results")
  eff <- do.call(rbind, rows)
  if (!is.null(eff)) rownames(eff) <- NULL
  if (!is.null(eff) && length(unique(eff$term)) >= 1L && length(sets) > 1L) {
    sp <- tapply(eff$estimate, eff$term, function(v) diff(range(v)))
    say("  spread of the exposure estimate across adjustment sets: ",
        paste(sprintf("%s %.4f", names(sp), sp), collapse = ", "))
    steps <- c(steps, paste0("spread across sets: ",
                             paste(sprintf("%s %.4f", names(sp), sp),
                                   collapse = ", ")))
  }
  structure(list(dag = g, exposure = x, outcome = y, family = fam,
                 sets = sets, fits = fits, effects = eff, dag_test = dt,
                 random = re, steps = steps, identified = TRUE),
            class = "ilm_dag_model")
}

## Which coefficients belong to one term. Matching on the name prefix is wrong
## -- an exposure `x` would also collect a covariate called `xray` -- so this
## uses the model matrix's `assign` attribute, which maps each column to the
## term that produced it and so handles a factor's levels exactly.
#' @keywords internal
#' @noRd
ilm_term_cols <- function(fit, term) {
  tryCatch({
    mm <- stats::model.matrix(fit)
    tl <- attr(stats::terms(fit), "term.labels")
    j <- match(ilm_as_label(term, tl), tl)
    if (is.na(j)) integer() else which(attr(mm, "assign") == j)
  }, error = function(e) which(names(stats::coef(fit)) == term))
}

#' @export
print.ilm_dag_model <- function(x, ...) {
  cat("<ilm_dag_model>", x$exposure, "->", x$outcome, "\n")
  if (!isTRUE(x$identified)) {
    cat("\n  NOT IDENTIFIABLE from the measured variables.\n")
    cat("  No set of them closes the back-door paths, so no regression on\n")
    cat("  these data estimates this effect. Measure a confounder, or revise\n")
    cat("  the graph.\n")
    return(invisible(x))
  }
  cat("  family:", x$family)
  if (length(x$random)) cat("   random intercepts:", paste(x$random, collapse = ", "))
  cat("\n")
  if (!is.null(x$dag_test))
    cat("  graph vs data:", attr(x$dag_test, "verdict"),
        sprintf("(%d claim%s tested)", nrow(x$dag_test),
                if (nrow(x$dag_test) == 1L) "" else "s"), "\n")
  cat("\n  effect of", x$exposure, "by adjustment set:\n")
  e <- x$effects
  if (is.null(e)) { cat("    (nothing fitted)\n"); return(invisible(x)) }
  for (i in seq_len(nrow(e)))
    cat(sprintf("    [%d] %-28s %-14s %8.4f  (%7.4f, %7.4f)  p = %s\n",
                e$set[i], substr(e$adjusted_for[i], 1, 28), e$term[i],
                e$estimate[i], e$lower[i], e$upper[i],
                format.pval(e$p_value[i], digits = 3, eps = 1e-4)))
  if (length(unique(e$set)) > 1L)
    cat("\n  These estimate the same quantity. Agreement is evidence the graph\n",
        "  holds up; a set that disagrees is worth more attention than any\n",
        "  single estimate.\n", sep = "")
  invisible(x)
}
