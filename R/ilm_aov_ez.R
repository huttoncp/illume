## ---------------------------------------------------------------------------
## The user-facing entry point, the table assembly, and the follow-ups.
## See ilm_aov_internals.R for the marshalling and the sums of squares.
## ---------------------------------------------------------------------------

#' Factorial and repeated-measures ANOVA
#'
#' Specify the design by naming columns -- the participant, the outcome, which
#' factors vary between participants and which vary within them -- rather than
#' by writing a formula with an error term in it. Reports F on a denominator
#' degrees of freedom, mean squared error, generalized eta squared, and a
#' sphericity correction where one is needed.
#'
#' @section What it reports, and why not just a mixed model:
#' On a balanced complete design a mixed model and a classical
#' repeated-measures ANOVA answer the same question. They report it
#' differently, and one of them reports things the other cannot:
#' [ilm_model()] gives a Wald chi-square, because once random effects are
#' integrated out there is no exact residual degrees of freedom to divide by.
#' An F, a mean squared error, a generalized eta squared and a
#' Greenhouse-Geisser correction do not fall out of that fit. So the omnibus
#' table here is computed classically, and a mixed model is fitted alongside
#' for the things the classical route cannot do -- marginal means, contrasts,
#' simple slopes, diagnostics, and an analysis that keeps participants with
#' incomplete data rather than dropping them.
#'
#' @section When the classical answer is the wrong one:
#' Three situations make it the wrong tool, and each is checked and reported
#' with the alternative named rather than left to be noticed:
#'
#' * **Sphericity fails.** The within-participant variances and covariances
#'   are not what the F test assumes. Greenhouse-Geisser and Huynh-Feldt
#'   corrections are reported, and so is the other option: an unstructured
#'   mixed model does not assume sphericity, so it needs no correction at all.
#' * **Participants have missing cells.** Classical repeated measures has to
#'   drop them completely, because a participant missing one occasion
#'   contributes to no within-participant contrast. How many were dropped is
#'   reported. A mixed model uses them.
#' * **A covariate varies within participant.** Refused, because a single
#'   coefficient for a time-varying covariate blends two effects that can
#'   point in opposite directions -- how participants who score higher on
#'   average differ, and what happens when a participant scores higher than
#'   usual.
#'
#' @section Generalized eta squared:
#' `ges` is comparable across designs in a way that partial eta squared is
#' not, which is why it is the default here and in `afex`. Its denominator
#' includes the sums of squares of any **measured** rather than manipulated
#' variable. Covariates are measured by definition, so they are added to
#' `observed` automatically; name any measured factors there too.
#'
#' @param id Column identifying the participant.
#' @param dv Column holding the numeric outcome.
#' @param data A long data frame: one row per participant per occasion.
#' @param between Factors that vary BETWEEN participants, as a character
#'   vector.
#' @param within Factors that vary WITHIN participants.
#' @param covariate Numeric covariates. Must be constant within participant.
#' @param observed Variables measured rather than manipulated, affecting the
#'   `ges` denominator. Covariates are included automatically.
#' @param engine `"aov"` (the default) computes the omnibus table classically
#'   and fits a mixed model alongside. `"mixed"` makes the mixed model the
#'   primary analysis, translating this specification into the corresponding
#'   [ilm_model()] formula.
#' @param reml Passed to [ilm_model()] for the companion fit. `TRUE` here,
#'   because the design fixed the fixed effects before any data were seen.
#' @param fun_aggregate Applied when a participant has several rows per cell.
#'   Defaults to the mean, with a message saying it happened.
#' @param posthoc Follow a significant effect with pairwise comparisons.
#'   `"auto"` does it for significant effects involving a factor with three or
#'   more levels; `FALSE` never does; a character vector names effects.
#' @param alpha Threshold for deciding what `posthoc` follows up, and for
#'   flagging sphericity.
#' @param correction Sphericity correction applied to the reported p-value:
#'   `"auto"` (Greenhouse-Geisser when Mauchly's test rejects), `"GG"`,
#'   `"HF"`, or `"none"`.
#' @param verbose Narrate the checks.
#' @return An object of class `"ilm_aov_ez"`: `anova` (the table),
#'   `sphericity`, `fit` (the companion [ilm_model()]), `posthoc`, `notes`,
#'   and the design.
#' @seealso [ilm_model()] for the general engine, [ilm_emmeans()] and
#'   [ilm_contrast()] for follow-ups, [ilm_trends()] when the interaction is
#'   with a continuous predictor.
#' @examples
#' set.seed(1)
#' n <- 24
#' d <- expand.grid(id = factor(seq_len(n)), time = factor(1:3))
#' d$grp <- factor(rep(c("ctl", "trt"), each = n / 2))[as.integer(d$id)]
#' d$score <- 10 + 2 * (d$grp == "trt") * as.integer(d$time) +
#'   rnorm(nrow(d)) + rep(rnorm(n), 3)
#' a <- ilm_aov_ez("id", "score", d, between = "grp", within = "time",
#'                 verbose = FALSE)
#' a
#' @export
ilm_aov_ez <- function(id, dv, data, between = NULL, within = NULL,
                       covariate = NULL, observed = NULL,
                       engine = c("aov", "mixed"), reml = TRUE,
                       fun_aggregate = NULL, posthoc = "auto", alpha = 0.05,
                       correction = c("auto", "GG", "HF", "none"),
                       verbose = TRUE) {
  engine <- match.arg(engine); correction <- match.arg(correction)
  if (!is.data.frame(data))
    stop("`data` must be a data frame; it is ", class(data)[1], call. = FALSE)
  id <- as.character(id)[1]; dv <- as.character(dv)[1]
  between <- as.character(between); within <- as.character(within)
  covariate <- as.character(covariate); observed <- as.character(observed)
  if (!length(between) && !length(within))
    stop("name at least one `between` or `within` factor; with neither there ",
         "is no effect to test", call. = FALSE)
  dup <- intersect(between, within)
  if (length(dup))
    stop("`", paste(dup, collapse = "`, `"), "` is named as both between- and ",
         "within-participant. A factor varies one way or the other: it is ",
         "`within` if each participant is observed at more than one of its ",
         "levels.", call. = FALSE)
  ## a covariate is measured, not manipulated, so it belongs in the ges
  ## denominator whether or not anyone remembered to say so
  observed <- unique(c(observed, covariate))

  m <- ilm_aov_marshal(data, id, dv, between, within, covariate, fun_aggregate)
  notes <- character(0)
  say <- function(...) {
    txt <- paste0(...)
    notes <<- c(notes, txt)
    if (verbose) message("  ", txt)
    invisible(NULL)
  }
  if (verbose)
    message("ilm_aov_ez: ", m$nsub, " participants, ", length(m$cells),
            " within-participant cell(s)")
  if (m$n_na) say(m$n_na, " row(s) dropped for missing values.")
  if (m$aggregated)
    say("More than one row per participant per cell: averaged within cells. ",
        "If those repeats are meaningful they are a `within` factor, or use ",
        "engine = \"mixed\", which can model them.")
  if (m$dropped)
    say(m$dropped, " participant(s) dropped for incomplete cells -- classical ",
        "repeated measures cannot use a participant missing an occasion. ",
        "engine = \"mixed\" keeps them.")
  if (m$nsub < 3L)
    stop("only ", m$nsub, " participant(s) have complete data; there is ",
         "nothing to test.", call. = FALSE)

  ## ---- the between-participants design -----------------------------------
  ## contr.sum throughout, which is what makes these Type III sums of squares
  ## the average over the other factor rather than its effect at whichever
  ## level happens to be first. afex does the same and says so.
  ## Names are quoted on their way into the formula, so a column called
  ## `dose group` is one variable and not two symbols; the effect labels are
  ## turned back into plain names below, where they are matched and printed.
  bterms <- c(between, covariate)
  bform <- if (length(bterms)) {
    parts <- c(if (length(between)) paste(ilm_bq(between), collapse = " * "),
               ilm_bq(covariate))
    stats::as.formula(paste("~", paste(parts, collapse = " + ")))
  } else ~ 1
  ctr <- if (length(between))
    stats::setNames(rep(list("contr.sum"), length(between)), between) else NULL
  X <- stats::model.matrix(bform, m$bdat, contrasts.arg = ctr)
  basg <- attr(X, "assign")
  blabs <- c("(Intercept)", attr(stats::terms(bform), "term.labels"))

  Ps <- ilm_aov_within_P(m$idata, within)
  wlabs <- names(Ps)

  ## ---- every effect: each between term crossed with each within term ------
  rows <- list(); sph <- list(); ss_intercept <- NULL
  for (wi in seq_along(wlabs)) for (bi in seq_along(blabs)) {
    wl <- wlabs[wi]; bl <- blabs[bi]
    if (wl == "(Intercept)" && bl == "(Intercept)") {
      ## The grand mean is never reported -- nobody wants a test that the
      ## outcome differs from zero -- but its ERROR is the between-participant
      ## stratum, and generalized eta squared needs every stratum's error in
      ## its denominator. In a design with a between-participants factor that
      ## stratum arrives anyway through that factor's row; in a purely
      ## within-participants design it would otherwise be missing, and every
      ## ges would come out more than twice too large.
      si <- ilm_aov_ss(m$Y, X, which(basg == 0L), Ps[[wi]])
      if (!is.null(si)) ss_intercept <- si$SSE
      next
    }
    lab <- if (wl == "(Intercept)") bl else
           if (bl == "(Intercept)") wl else paste0(bl, ":", wl)
    lab <- paste(ilm_unbq(strsplit(lab, ":", fixed = TRUE)[[1L]]),
                 collapse = ":")
    idx <- which(basg == bi - 1L)
    if (!length(idx)) next
    s <- ilm_aov_ss(m$Y, X, idx, Ps[[wi]])
    if (is.null(s)) next
    sp <- ilm_aov_sphericity(s$E, s$P, s$df_e_b)
    sph[[lab]] <- sp
    rows[[lab]] <- data.frame(
      effect = lab, num_df = s$df_h, den_df = s$df_e,
      MSE = s$SSE / s$df_e, SS = s$SS, SSE = s$SSE,
      F_value = (s$SS / s$df_h) / (s$SSE / s$df_e),
      gg = sp$gg, hf = sp$hf, stringsAsFactors = FALSE)
  }
  if (!length(rows))
    stop("no testable effect was produced from this specification",
         call. = FALSE)
  tab <- do.call(rbind, rows); rownames(tab) <- NULL

  ## ---- generalized eta squared -------------------------------------------
  ## SS / (SS + sum of the DISTINCT error sums of squares + the SS of any
  ## measured variable, less its own when it is the effect being sized).
  ## Olejnik and Algina (2003); the form is afex's.
  ## compared as names, not as a regular expression: a factor called
  ## `time (h)` is not a pattern, and matching it as one finds nothing
  obs <- vapply(tab$effect, function(e)
    any(observed %in% strsplit(e, ":", fixed = TRUE)[[1L]]), TRUE)
  obs_all <- sum(tab$SS * obs)
  err_pool <- sum(unique(c(tab$SSE, ss_intercept)))
  tab$ges <- tab$SS / (tab$SS + err_pool + obs_all - tab$SS * obs)
  tab$pes <- tab$SS / (tab$SS + tab$SSE)

  ## ---- sphericity correction ---------------------------------------------
  tab$correction <- "none"
  tab$p_value <- stats::pf(tab$F_value, tab$num_df, tab$den_df,
                           lower.tail = FALSE)
  for (i in seq_len(nrow(tab))) {
    sp <- sph[[tab$effect[i]]]
    if (is.na(sp$gg)) next
    use <- switch(correction,
      none = "none",
      GG = "GG", HF = "HF",
      auto = if (!is.na(sp$p_mauchly) && sp$p_mauchly < alpha) "GG" else "none")
    if (use == "none") next
    eps <- if (use == "GG") sp$gg else sp$hf
    eps <- max(min(eps, 1), 1e-6)
    tab$correction[i] <- use
    tab$num_df[i] <- tab$num_df[i] * eps
    tab$den_df[i] <- tab$den_df[i] * eps
    tab$p_value[i] <- stats::pf(tab$F_value[i], tab$num_df[i], tab$den_df[i],
                                lower.tail = FALSE)
  }

  ## ---- the companion mixed model ------------------------------------------
  fit <- NULL; fit_err <- NULL
  fo <- ilm_aov_formula(dv, between, within, covariate, id)
  fit <- tryCatch(suppressWarnings(
           ilm_model(fo, data = data, family = "gaussian", reml = reml,
                     verbose = FALSE)),
         error = function(e) { fit_err <<- conditionMessage(e); NULL })
  if (is.null(fit))
    say("the companion mixed model did not fit (", fit_err, "), so marginal ",
        "means and contrasts are unavailable.")

  ## ---- sphericity narration ------------------------------------------------
  bad <- names(sph)[vapply(sph, function(z)
    !is.na(z$p_mauchly) && z$p_mauchly < alpha, TRUE)]
  if (length(bad))
    say("sphericity rejected for ", paste(bad, collapse = ", "),
        ". The reported p-value uses ",
        if (correction == "none") "NO correction, which is anti-conservative here"
        else paste0("a ", if (correction == "auto") "Greenhouse-Geisser" else
                    correction, " correction"),
        ". A correction adjusts a test whose assumption has failed; the ",
        "alternative is a model that never made it -- ilm_model(..., ",
        "re_struct = \"us\") leaves the within-participant covariance free.")

  ph <- ilm_aov_posthoc(fit, tab, m, between, within, covariate, posthoc,
                        alpha, say)

  structure(list(anova = tab, sphericity = sph, fit = fit, posthoc = ph,
                 notes = notes, formula = fo, engine = engine,
                 design = list(id = id, dv = dv, between = between,
                               within = within, covariate = covariate,
                               observed = observed),
                 nsub = m$nsub, dropped = m$dropped, alpha = alpha,
                 correction = correction),
            class = "ilm_aov_ez")
}

#' The ilm_model() formula this specification describes
#'
#' The translation people otherwise do by hand and get wrong: between and
#' within factors are crossed, covariates are added, and the participant is a
#' random intercept.
#'
#' @param dv,between,within,covariate,id Column names.
#' @return A formula.
#' @keywords internal
#' @noRd
ilm_aov_formula <- function(dv, between, within, covariate, id) {
  fx <- ilm_bq(c(between, within))
  rhs <- if (length(fx)) paste(fx, collapse = " * ") else "1"
  if (length(covariate))
    rhs <- paste(rhs, "+", paste(ilm_bq(covariate), collapse = " + "))
  rhs <- paste0(rhs, " + (1 | ", ilm_bq(id), ")")
  stats::as.formula(paste(ilm_bq(dv), "~", rhs))
}

#' Pairwise follow-ups for the effects that earned them
#'
#' An omnibus F with three or more levels says the means are not all equal and
#' nothing about which differ, so the comparison is where the finding actually
#' is. Only significant effects are followed, and only those with something to
#' compare.
#'
#' @param fit The companion mixed model, or `NULL`.
#' @param tab The ANOVA table.
#' @param m The marshalled data.
#' @param between,within Factor names.
#' @param posthoc `"auto"`, `FALSE`, or effect names.
#' @param alpha Significance threshold.
#' @param say Note-recording function.
#' @return A named list of [ilm_contrast()] results.
#' @keywords internal
#' @noRd
ilm_aov_posthoc <- function(fit, tab, m, between, within, covariate, posthoc,
                            alpha, say) {
  if (isFALSE(posthoc) || is.null(fit)) return(list())
  facs <- c(between, within)
  nlev <- vapply(facs, function(v) {
    z <- if (v %in% names(m$bdat)) m$bdat[[v]] else m$idata[[v]]
    length(levels(droplevels(factor(z))))
  }, integer(1))
  want <- if (identical(posthoc, "auto")) {
    sig <- tab$effect[tab$p_value < alpha]
    Filter(function(e) {
      parts <- strsplit(e, ":", fixed = TRUE)[[1]]
      ## a covariate interaction is always worth following: the question is
      ## which slopes differ, and there is no "number of levels" to qualify on
      if (any(parts %in% covariate)) return(TRUE)
      ## otherwise only where there is something to compare -- a two-level
      ## factor's omnibus test IS the comparison, and repeating it adds nothing
      any(nlev[parts] >= 3, na.rm = TRUE)
    }, sig)
  } else as.character(posthoc)

  out <- list()
  for (e in want) {
    parts <- strsplit(e, ":", fixed = TRUE)[[1]]
    cv <- parts[parts %in% covariate]
    fp <- parts[parts %in% facs]

    ## ---- a factor by covariate interaction: whose slope, and is it zero? ---
    if (length(cv) == 1L && length(fp) >= 1L) {
      r <- tryCatch({
        tr <- ilm_trends(fit, fp, var = cv)
        list(trends = tr, contrasts = ilm_contrast(tr))
      }, error = function(err) NULL)
      if (!is.null(r)) { out[[e]] <- r; next }
    }
    if (!length(fp)) next

    ## ---- a two-way interaction: simple effects, not every cell pair -------
    ## Comparing all nine cells of a 3x3 gives thirty-six differences, most of
    ## which cross both factors at once and answer nothing anybody asked. What
    ## takes an interaction apart is the simple effect: hold one factor, compare
    ## the other, in both directions.
    if (length(fp) == 2L) {
      r <- tryCatch(ilm_aov_simple(fit, fp, m), error = function(err) NULL)
      if (!is.null(r)) { out[[e]] <- r; next }
    }
    r <- tryCatch(ilm_contrast(ilm_emmeans(fit, fp)),
                  error = function(err) NULL)
    if (!is.null(r)) out[[e]] <- r
  }
  if (length(out))
    say("pairwise comparisons computed for: ", paste(names(out), collapse = ", "),
        ". These come from the mixed model, so they use all the data and a ",
        "single pooled error; they are not the classical stratum-specific ",
        "error terms.")
  out
}

#' @export
print.ilm_aov_ez <- function(x, digits = 3, ...) {
  d <- x$design
  cat("Analysis of Variance\n")
  cat("  ", d$dv, " ~ ", sep = "")
  bits <- c(if (length(d$between)) paste0(paste(d$between, collapse = " * "),
                                          " (between)"),
            if (length(d$within)) paste0(paste(d$within, collapse = " * "),
                                         " (within)"),
            if (length(d$covariate)) paste0(paste(d$covariate, collapse = " + "),
                                            " (covariate)"))
  cat(paste(bits, collapse = " + "), "\n", sep = "")
  cat("  ", x$nsub, " participants",
      if (x$dropped) paste0(" (", x$dropped, " dropped for incomplete cells)")
      else "", "\n\n", sep = "")

  t2 <- x$anova
  sig <- stats::symnum(t2$p_value, corr = FALSE, na = FALSE,
                       cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                       symbols = c("***", "**", "*", ".", " "))
  out <- data.frame(
    Effect = t2$effect,
    df = paste0(format(round(t2$num_df, 2), trim = TRUE), ", ",
                format(round(t2$den_df, 2), trim = TRUE)),
    MSE = round(t2$MSE, digits),
    F = round(t2$F_value, digits),
    ges = round(t2$ges, 3),
    p = format.pval(t2$p_value, digits = 3, eps = 1e-4),
    ` ` = as.character(sig), check.names = FALSE, stringsAsFactors = FALSE)
  if (any(t2$correction != "none"))
    out$corr <- t2$correction
  print(out, row.names = FALSE)
  cat("---\nSignif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")

  sp <- x$sphericity
  keep <- names(sp)[vapply(sp, function(z) !is.na(z$W), TRUE)]
  if (length(keep)) {
    cat("\nSphericity (Mauchly):\n")
    ms <- data.frame(Effect = keep,
                     W = round(vapply(sp[keep], `[[`, 1, "W"), 3),
                     p = format.pval(vapply(sp[keep], `[[`, 1, "p_mauchly"),
                                     digits = 3, eps = 1e-4),
                     GG = round(vapply(sp[keep], `[[`, 1, "gg"), 3),
                     HF = round(vapply(sp[keep], `[[`, 1, "hf"), 3),
                     stringsAsFactors = FALSE)
    print(ms, row.names = FALSE)
  }
  if (length(x$posthoc)) {
    cat("\nPairwise comparisons:\n")
    for (nm in names(x$posthoc)) {
      cat("\n-- ", nm, " --\n", sep = "")
      print(x$posthoc[[nm]])
    }
  }
  if (length(x$notes)) {
    cat("\nNotes:\n")
    for (n in x$notes) cat("  * ", strwrap(n, width = 76, exdent = 4), "\n",
                           sep = "")
  }
  invisible(x)
}
