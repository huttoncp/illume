## ---------------------------------------------------------------------------
## Missing data: what the pattern says, and whether it matters.
##
## Three facts govern everything here, and two of them are routinely got wrong.
##
## 1. MAR and MNAR CANNOT be told apart from the observed data. That is a
##    theorem, not a limitation of effort: the data that would distinguish them
##    are the ones that are missing. Any function claiming to test it is lying.
##    What IS testable is MCAR -- whether missingness is unrelated to anything
##    observed -- and that is what is tested here.
##
## 2. Complete-case analysis is unbiased under conditions much weaker than
##    MCAR. For a regression of y on X, dropping incomplete rows is unbiased
##    whenever missingness is independent of Y GIVEN X -- even when it depends
##    on X strongly. So "not MCAR" does not imply "must impute", and the
##    distinction this makes is between missingness tied to the covariates
##    (complete cases are fine, imputation buys precision at best) and
##    missingness tied to the outcome (complete cases are biased).
##
## 3. Single imputation destroys the standard errors. Filling a value in and
##    then analysing as though it had been observed treats a guess as data:
##    intervals come out too narrow, by an amount that grows with how much is
##    missing. In a package whose whole claim is nominal coverage that would be
##    self-defeating, so multiple imputation is the default and single
##    imputation has to be asked for by name.
## ---------------------------------------------------------------------------

## Association between a missingness indicator and an observed variable, on one
## scale whatever the variable's type: a correlation for numeric, Cramer's V for
## a factor. Both sit in [0, 1] and are comparable enough to sort by.
#' @keywords internal
#' @noRd
ilm_miss_assoc <- function(r, w) {
  keep <- !is.na(w)
  if (sum(keep) < 10L || length(unique(r[keep])) < 2L)
    return(c(effect = NA_real_, p = NA_real_))
  r <- r[keep]; w <- w[keep]
  if (is.numeric(w)) {
    if (stats::sd(w) == 0) return(c(effect = NA_real_, p = NA_real_))
    ct <- suppressWarnings(stats::cor.test(as.numeric(r), w))
    return(c(effect = unname(abs(ct$estimate)), p = unname(ct$p.value)))
  }
  w <- factor(w)
  if (nlevels(w) < 2L) return(c(effect = NA_real_, p = NA_real_))
  tb <- table(r, w)
  if (any(dim(tb) < 2L)) return(c(effect = NA_real_, p = NA_real_))
  cs <- suppressWarnings(stats::chisq.test(tb))
  v <- sqrt(unname(cs$statistic) / (sum(tb) * (min(dim(tb)) - 1L)))
  c(effect = unname(v), p = unname(cs$p.value))
}

## Does missingness in `v` depend on the outcome BEYOND what the covariates
## already explain? That conditional question -- not the marginal one -- is what
## decides whether complete cases are biased.
##
## The distinction is not academic and is easy to get wrong. If missingness
## depends on a covariate z, and the outcome also depends on z, then missingness
## and the outcome are MARGINALLY associated while carrying no information about
## each other once z is held fixed. Complete cases are unbiased there. A
## marginal test says otherwise and sends the user off to impute for nothing.
#' @keywords internal
#' @noRd
ilm_miss_outcome_test <- function(data, v, y, covars) {
  r <- as.integer(is.na(data[[v]]))
  covars <- setdiff(covars, c(v, y))
  d <- data[, unique(c(y, covars)), drop = FALSE]
  d$.r <- r
  d <- d[stats::complete.cases(d), , drop = FALSE]
  if (nrow(d) < 20L || length(unique(d$.r)) < 2L)
    return(list(effect = NA_real_, p_value = NA_real_))
  ## keep only covariates that can actually enter a model
  covars <- covars[vapply(covars, function(w) {
    z <- d[[w]]
    if (is.numeric(z)) stats::sd(z, na.rm = TRUE) > 0
    else nlevels(droplevels(factor(z))) > 1L &&
         nlevels(droplevels(factor(z))) <= 20L
  }, TRUE)]
  f0 <- stats::reformulate(if (length(covars)) covars else "1", response = ".r")
  f1 <- stats::reformulate(c(covars, y), response = ".r")
  g0 <- tryCatch(suppressWarnings(stats::glm(f0, data = d,
                   family = stats::binomial())), error = function(e) NULL)
  g1 <- tryCatch(suppressWarnings(stats::glm(f1, data = d,
                   family = stats::binomial())), error = function(e) NULL)
  if (is.null(g0) || is.null(g1)) return(list(effect = NA_real_, p_value = NA_real_))
  an <- tryCatch(suppressWarnings(stats::anova(g0, g1, test = "Chisq")),
                 error = function(e) NULL)
  pc <- grep("^Pr", names(an), value = TRUE)
  pv <- if (is.null(an) || !length(pc)) NA_real_ else an[[pc[1]]][2]
  ## magnitude on one scale: the correlation left between missingness and the
  ## outcome once the covariates are taken out of both
  eff <- tryCatch({
    ra <- stats::residuals(g0, type = "pearson")
    rb <- if (!length(covars)) d[[y]] else
      stats::residuals(stats::lm(stats::reformulate(covars, response = y), data = d))
    abs(unname(stats::cor(ra, rb, use = "complete.obs")))
  }, error = function(e) NA_real_)
  list(effect = eff, p_value = pv)
}

#' What the missing values look like, and whether they matter
#'
#' Reports how much is missing, which variables go missing together, and what
#' missingness is related to -- then says what that implies for the analysis.
#'
#' @section What can and cannot be established:
#'
#' **MCAR** -- missingness unrelated to anything -- is testable, and is tested
#' here by asking whether each variable's missingness is associated with the
#' observed values of the others.
#'
#' **MAR versus MNAR is not testable.** The data that would separate them are
#' precisely the data that are missing. No amount of pattern analysis settles
#' it, and this function does not pretend otherwise: it reports what
#' missingness is associated with among the things you *can* see, and says
#' plainly that a dependence on the unseen values themselves cannot be ruled
#' out.
#'
#' @section Why "not MCAR" does not mean "impute":
#'
#' For a regression of `y` on covariates, dropping incomplete rows is unbiased
#' whenever missingness is independent of `y` **given** the covariates -- even
#' when it depends on the covariates strongly. That is a far weaker condition
#' than MCAR, and it is why naming `y` changes the advice:
#'
#' * missingness related to the covariates but not to `y`: complete cases are
#'   unbiased. Imputation can recover precision but is not needed for
#'   correctness, and it adds assumptions.
#' * missingness related to `y` itself: complete cases are biased, and
#'   [ilm_impute()] is the remedy.
#'
#' @param data A data frame.
#' @param covariates Variables to hold fixed when asking whether missingness
#'   depends on the outcome. Default is every other column; pass the model's
#'   predictors when they are a subset.
#' @param y Optional outcome column. Naming it is what separates "complete
#'   cases are fine" from "complete cases are biased", so it is worth naming.
#' @param min_effect Smallest association worth reporting, as a correlation or
#'   Cramer's V. With several thousand rows an association of 0.02 is
#'   significant and means nothing.
#' @param alpha Level for the adjusted p-values.
#' @param adjust Multiplicity adjustment across pairs, passed to
#'   [stats::p.adjust()].
#' @param verbose Narrate the findings.
#' @return An object of class `"ilm_missing"`: `variables`, `patterns`,
#'   `associations`, `monotone` and a `verdict`.
#' @seealso [ilm_impute()], [ilm_describe_na_all()], [ilm_plot_missing()].
#' @examples
#' d <- ilm_sim()
#' ilm_check_missing(d, y = "score", verbose = FALSE)
#' @export
ilm_check_missing <- function(data, y = NULL, covariates = NULL,
                              min_effect = 0.1, alpha = 0.05,
                              adjust = "holm", verbose = TRUE) {
  if (!is.data.frame(data))
    stop("`data` must be a data frame; it is ", class(data)[1], call. = FALSE)
  if (!is.null(y) && !y %in% names(data))
    stop("`y` (", y, ") is not a column of `data`", call. = FALSE)
  say <- function(...) if (verbose) message(...)
  n <- nrow(data)

  ## ---- how much, and where -------------------------------------------------
  pna <- vapply(data, function(v) mean(is.na(v)), 1)
  vars <- data.frame(variable = names(pna), n_missing = vapply(data,
                       function(v) sum(is.na(v)), 1L),
                     p_missing = unname(pna), stringsAsFactors = FALSE)
  vars <- vars[order(-vars$p_missing), , drop = FALSE]
  rownames(vars) <- NULL
  inc <- names(pna)[pna > 0]
  cc <- sum(stats::complete.cases(data))

  say("== missing data ==")
  if (!length(inc)) {
    say("  no missing values in any of the ", ncol(data), " columns")
    return(structure(list(variables = vars, patterns = NULL,
                          associations = NULL, monotone = NA,
                          n = n, n_complete = n, y = y, verdict = "NONE"),
                     class = "ilm_missing"))
  }
  say("  ", length(inc), " of ", ncol(data), " columns have missing values")
  say("  ", cc, " of ", n, " rows are complete (",
      sprintf("%.1f%%", 100 * cc / n), ")")

  ## ---- which variables go missing together ---------------------------------
  M <- is.na(data[inc])
  key <- apply(M, 1L, function(z) paste0(as.integer(z), collapse = ""))
  tb <- sort(table(key), decreasing = TRUE)
  patterns <- data.frame(pattern = names(tb), n = as.integer(tb),
                         n_missing = vapply(names(tb), function(k)
                           sum(strsplit(k, "")[[1]] == "1"), 1L),
                         stringsAsFactors = FALSE)
  rownames(patterns) <- NULL
  attr(patterns, "variables") <- inc
  say("  ", nrow(patterns), " distinct missingness pattern",
      if (nrow(patterns) == 1L) "" else "s",
      " across those ", length(inc), " column",
      if (length(inc) == 1L) "" else "s")

  ## A monotone pattern means the variables can be ordered so that once one is
  ## missing every later one is too -- which is what a drop-out study produces,
  ## and it makes imputation a chain of ordinary regressions rather than a loop.
  ord <- inc[order(colSums(M))]
  Mo <- M[, ord, drop = FALSE]
  monotone <- all(vapply(seq_len(ncol(Mo) - 1L), function(j)
    all(Mo[, j] <= Mo[, j + 1L]), TRUE))
  if (ncol(Mo) < 2L) monotone <- TRUE
  say("  pattern is ", if (monotone) "MONOTONE" else "non-monotone",
      if (monotone) " (missingness accumulates in one order, as drop-out does)"
      else " (variables go missing independently of each other)")

  ## ---- what is missingness related to? -------------------------------------
  rows <- list()
  for (v in inc) {
    r <- is.na(data[[v]])
    for (w in setdiff(names(data), v)) {
      a <- ilm_miss_assoc(r, data[[w]])
      if (is.na(a[["effect"]])) next
      rows[[length(rows) + 1L]] <- data.frame(
        missing_in = v, related_to = w, effect = a[["effect"]],
        p_value = a[["p"]], is_outcome = !is.null(y) && identical(w, y),
        stringsAsFactors = FALSE)
    }
  }
  assoc <- if (length(rows)) do.call(rbind, rows) else NULL
  if (!is.null(assoc)) {
    assoc$p_adj <- stats::p.adjust(assoc$p_value, method = adjust)
    assoc$flag <- assoc$p_adj < alpha & assoc$effect >= min_effect
    assoc <- assoc[order(-assoc$effect), , drop = FALSE]
    rownames(assoc) <- NULL
  }

  ## ---- does missingness depend on the OUTCOME, given the covariates? -------
  ## The marginal table above is description. This is the question that decides
  ## whether complete cases are biased, and it has to be asked conditionally:
  ## missingness driven by a covariate the outcome also depends on is
  ## marginally associated with the outcome while telling you nothing about it
  ## once that covariate is held fixed.
  cvs <- if (is.null(covariates)) setdiff(names(data), y)
         else intersect(as.character(covariates), names(data))
  outcome_test <- NULL
  if (!is.null(y)) {
    ot <- lapply(inc, function(v) {
      r <- ilm_miss_outcome_test(data, v, y, cvs)
      data.frame(missing_in = v, effect = r$effect, p_value = r$p_value,
                 stringsAsFactors = FALSE)
    })
    outcome_test <- do.call(rbind, ot)
    outcome_test$p_adj <- stats::p.adjust(outcome_test$p_value, method = adjust)
    outcome_test$flag <- !is.na(outcome_test$p_adj) &
      outcome_test$p_adj < alpha & !is.na(outcome_test$effect) &
      outcome_test$effect >= min_effect
    rownames(outcome_test) <- NULL
  }

  ## ---- the verdict ---------------------------------------------------------
  flagged <- if (is.null(assoc)) assoc else assoc[assoc$flag, , drop = FALSE]
  n_flag <- if (is.null(flagged)) 0L else nrow(flagged)
  y_flag <- if (is.null(outcome_test)) 0L else sum(outcome_test$flag)

  verdict <- if (y_flag) "RELATED_TO_OUTCOME"
             else if (n_flag) "RELATED_TO_COVARIATES" else "MCAR_NOT_REJECTED"

  if (verbose) {
    if (verdict == "MCAR_NOT_REJECTED") {
      say("  missingness is not associated with any observed variable above ",
          "min_effect = ", min_effect, ": consistent with MCAR.")
      say("    Complete-case analysis is unbiased under MCAR. Imputation would ",
          "recover the ", n - cc, " dropped row(s) of precision and nothing else.")
    } else {
      say("  missingness IS associated with observed variables:")
      for (i in seq_len(min(n_flag, 6L)))
        say(sprintf("    %s missing <- %s (effect %.3f, p = %s)%s",
                    flagged$missing_in[i], flagged$related_to[i],
                    flagged$effect[i],
                    format.pval(flagged$p_adj[i], digits = 2, eps = 1e-4),
                    if (flagged$is_outcome[i]) "   <- THE OUTCOME" else ""))
      if (n_flag > 6L) say("    ... and ", n_flag - 6L, " more")
      if (verdict == "RELATED_TO_OUTCOME") {
        say("  Missingness depends on the OUTCOME even after the other ",
            "variables are held fixed. Complete-case analysis of a model for `",
            y, "` is biased here, because the rows that survive are not a fair ",
            "sample of the rows that matter.")
        say("    Remedy: ilm_impute(), then ilm_mi_pool().")
      } else {
        say("  Missingness is related to covariates but NOT to ",
            if (is.null(y)) "any named outcome"
            else paste0("the outcome `", y, "`"), ".")
        say("    Held against the other variables, missingness carries no ",
            "further information about the outcome, so complete-case analysis ",
            "stays UNBIASED: it needs missingness independent of the outcome ",
            "GIVEN the covariates, not MCAR. Imputation would buy back the ",
            n - cc, " dropped row(s) of precision and nothing else.")
        if (is.null(y))
          say("    Name `y` to have this checked against your outcome rather ",
              "than assumed.")
      }
    }
    say("  MNAR cannot be ruled out by any of this. If the chance of a value ",
        "going missing depends on the value itself, nothing observed reveals ",
        "it; the remedy is a sensitivity analysis, not a test.")
  }

  structure(list(variables = vars, patterns = patterns, associations = assoc,
                 outcome_test = outcome_test, monotone = monotone, n = n,
                 n_complete = cc, y = y, min_effect = min_effect,
                 verdict = verdict),
            class = "ilm_missing")
}

#' @export
print.ilm_missing <- function(x, ...) {
  cat("<ilm_missing>", x$n_complete, "of", x$n, "rows complete",
      sprintf("(%.1f%%)\n", 100 * x$n_complete / x$n))
  if (identical(x$verdict, "NONE")) {
    cat("  no missing values\n"); return(invisible(x))
  }
  v <- x$variables[x$variables$n_missing > 0, , drop = FALSE]
  cat("\n  missing by variable\n")
  for (i in seq_len(nrow(v)))
    cat(sprintf("    %-20s %6d  %5.1f%%\n", v$variable[i], v$n_missing[i],
                100 * v$p_missing[i]))
  cat("\n  pattern:", if (isTRUE(x$monotone)) "monotone" else "non-monotone",
      sprintf("(%d distinct)\n", nrow(x$patterns)))
  cat("  verdict:", x$verdict, "\n")
  msg <- switch(x$verdict,
    MCAR_NOT_REJECTED = "    consistent with MCAR; complete cases are unbiased",
    RELATED_TO_COVARIATES = "    related to covariates, not to the outcome: complete cases stay unbiased",
    RELATED_TO_OUTCOME = "    related to the OUTCOME: complete cases are biased -- see ilm_impute()",
    "")
  if (nzchar(msg)) cat(msg, "\n")
  if (!is.null(x$outcome_test)) {
    cat(sprintf("
  does missingness depend on `%s` once the others are held fixed?
",
                x$y))
    ot <- x$outcome_test
    for (i in seq_len(nrow(ot)))
      cat(sprintf("    %-16s partial r = %s, p = %s%s
", ot$missing_in[i],
                  ifelse(is.na(ot$effect[i]), "-", sprintf("%.3f", ot$effect[i])),
                  format.pval(ot$p_adj[i], digits = 2, eps = 1e-4),
                  if (isTRUE(ot$flag[i])) "   <- yes" else ""))
  }
  if (!is.null(x$associations)) {
    f <- x$associations[x$associations$flag, , drop = FALSE]
    if (nrow(f)) {
      cat("\n  strongest associations with missingness\n")
      for (i in seq_len(min(nrow(f), 8L)))
        cat(sprintf("    %-16s <- %-16s %.3f%s\n", f$missing_in[i],
                    f$related_to[i], f$effect[i],
                    if (f$is_outcome[i]) "  (outcome)" else ""))
    }
  }
  cat("\n  MNAR cannot be ruled out by any test; that needs a sensitivity analysis.\n")
  invisible(x)
}
