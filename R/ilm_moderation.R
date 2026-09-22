## ---------------------------------------------------------------------------
## Does the effect of X on Y depend on something else?
##
## WHY THIS EXISTS, and why it is not just "fit the interaction": an analyst
## who tries several moderators and reports the strongest is running a search,
## and the p-value from that search is not the p-value from a test. Measured
## with six candidates in a model with no moderation at all, the naive
## pick-the-best-and-report-it procedure rejected 29.3% of the time at a
## nominal 5%. That is the number this function exists to prevent.
##
## WHAT FIXES IT, measured rather than assumed. Three options were compared on
## data with known moderation, each candidate tested by the JOINT test of its
## interaction block:
##
##                              naive   holm   BH    bonf   honest split
##   null (want 0.05)           0.293   0.027  0.027 0.027  0.053
##   moderation by numeric      0.807   0.547  0.553 0.547  0.320
##   moderation by 3-level fac  0.707   0.380  0.380 0.380  0.200
##
## A multiplicity adjustment wins, and it is not close. Sample splitting -- fit
## on half, test the winner on the other half -- is better CALIBRATED (0.053
## against Holm's conservative 0.027) and costs roughly half the power. Worse,
## the cost compounds with degrees of freedom: going from a 1 df numeric
## moderator to a 2 df three-level factor takes Holm down 31% and splitting
## down 38%, so splitting is at its worst exactly in the three-level-factor
## case that experimental work is full of.
##
## So Holm is the default and `split = TRUE` is offered for anyone who wants
## the calibration and can pay for it. Splitting is only strictly NECESSARY
## when the hypotheses cannot be enumerated -- an open-ended tree search over
## covariates and split points, where no adjustment is available to apply.
## That is not implemented here; this function searches a named, countable set.
##
## THE TEST IS THE JOINT ONE. A single interaction coefficient is a test of
## moderation only when the moderator has one degree of freedom. On a 4-level
## factor the three coefficients gave p = 0.62, 0.0018 and 0.00002 for the
## same moderator; ilm_anova() gives the one answer, F on 3 df.
## ---------------------------------------------------------------------------

#' Does the effect of X depend on something else?
#'
#' Searches a set of candidate moderators for evidence that the effect of a
#' treatment or exposure varies, testing each by the joint test of its
#' interaction block and adjusting for having looked at several.
#'
#' @section What this is for, and what it is not:
#'
#' This is the **exploratory** tool. A moderator you hypothesised a priori
#' belongs in the model, where [ilm_anova()] tests it and [ilm_emmeans()] and
#' [ilm_contrast()] describe it -- and it owes no multiplicity penalty, because
#' you did not go looking. Any such interaction already present in `object` is
#' therefore excluded from the search and said so, rather than being charged
#' for a search it was not part of.
#'
#' @section Why the adjustment is not optional:
#'
#' Trying several moderators and reporting the strongest is a search. With six
#' candidates and no moderation present at all, that procedure rejected 29.3%
#' of the time at a nominal 5%. Every adjustment method brings it back: Holm,
#' BH and Bonferroni all landed at 0.027 in the same simulation.
#'
#' @section Holm, or an honest split:
#'
#' `split = TRUE` picks the strongest moderator on half the data and tests it
#' on the other half, which needs no adjustment because the test half was never
#' used to choose. It is better calibrated -- 0.053 against Holm's conservative
#' 0.027 -- and costs about half the power (0.320 against 0.547 for a numeric
#' moderator, 0.200 against 0.380 for a three-level factor). The gap widens
#' with degrees of freedom, so splitting is weakest exactly where experimental
#' designs live. Holm is the default for that reason.
#'
#' Splitting is genuinely *required* only when the hypotheses cannot be
#' counted, as in an open-ended tree search over covariates and split points.
#' This function searches a named set, so they can be.
#'
#' @param object A fitted [ilm_model()], or an [ilm_dag_model()].
#' @param x The treatment, exposure or intervention whose effect might be
#'   moderated. Required, **except** for an [ilm_dag_model()], where the graph
#'   has already named the exposure. There is no sensible default otherwise:
#'   `x:m` is the same term whichever of the two you call the moderator, and
#'   only the interpretation distinguishes them.
#' @param moderators Candidate moderators. Defaults to the model's other
#'   fixed-effect terms -- variables you already judged worth adjusting for --
#'   which keeps the multiplicity burden honest. Naming them explicitly can
#'   reach anything in the model frame, or in `data`.
#' @param adjust Passed to [stats::p.adjust()]; any method it supports,
#'   including `"holm"` (the default), `"BH"`, `"bonferroni"`, `"BY"`,
#'   `"hochberg"`, `"hommel"` and `"none"`.
#' @param split Use an honest sample split instead of an adjustment. See above:
#'   better calibrated, roughly half the power.
#' @param data Needed only for a moderator that is not in the model frame.
#' @param n_keep Refits to retain on the result, strongest first, so
#'   [ilm_plot_moderation()] can draw them without refitting. The rest are
#'   refitted on demand.
#' @param seed Random seed, used only when `split = TRUE`.
#' @param progress Show a progress bar. Defaults to [interactive()].
#' @return An object of class `"ilm_moderation"`: a data frame with one row per
#'   candidate (`moderator`, `df`, `statistic`, `p`, `p_adj`), carrying the
#'   retained refits and the exposure.
#' @seealso [ilm_plot_moderation()] to see one, [ilm_anova()] for a moderator
#'   you specified in advance, [ilm_trends()] for slopes within levels.
#' @examples
#' \donttest{
#' set.seed(1)
#' d <- data.frame(tx = rbinom(300, 1, 0.5), age = rnorm(300),
#'                 site = factor(sample(c("a", "b", "c"), 300, TRUE)))
#' d$y <- 0.4 * d$tx + 0.3 * d$age + 0.8 * d$tx * (d$site == "c") + rnorm(300)
#' fit <- ilm_model(y ~ tx + age + site, data = d, verbose = FALSE)
#' ilm_moderation(fit, x = "tx")
#' }
#' @export
ilm_moderation <- function(object, x = NULL, moderators = NULL,
                           adjust = "holm", split = FALSE, data = NULL,
                           n_keep = 3L, seed = 1L, progress = NULL) {
  dag <- inherits(object, "ilm_dag_model")
  fit <- if (dag) object$fit else object
  if (!inherits(fit, "ilm_model"))
    stop("`object` must be a fitted ilm_model() or an ilm_dag_model(); it is ",
         class(object)[1], call. = FALSE)
  if (is.null(fit$call) || is.null(fit$model))
    stop("this fit did not keep its call and model frame, so the moderation ",
         "models cannot be built from it; fit through the formula interface.",
         call. = FALSE)
  if (!adjust %in% stats::p.adjust.methods)
    stop("`adjust` must be one of stats::p.adjust.methods: ",
         paste(stats::p.adjust.methods, collapse = ", "), call. = FALSE)

  ## ---- which variable's effect is being moderated -------------------------
  if (is.null(x) && dag) {
    x <- object$exposure
    message("ilm_moderation(): taking the exposure from the DAG: '", x, "'.")
  }
  ## Two spellings of every name from here on. Internally a term is its LABEL,
  ## as R writes it -- `x 1` with its backticks -- because that is what goes
  ## into a formula and what ilm_anova() names its rows by. What comes back to
  ## the user, and what ilm_plot_moderation() looks up in the model frame, is
  ## the plain name, x 1.
  tl <- fit$term_labels
  main <- tl[!grepl(":", tl)]
  if (is.null(x))
    stop("`x` must name the treatment, exposure or intervention whose effect ",
         "might be moderated. There is no default: `x:m` is the same term ",
         "whichever of the two you call the moderator, so only you know which ",
         "is which. This model's fixed effects are: ",
         paste(ilm_unbq(main), collapse = ", "), ".", call. = FALSE)
  x <- ilm_as_label(x, main)
  if (!x %in% main)
    stop("`x` (", x, ") is not a fixed effect of this model. It has: ",
         paste(ilm_unbq(main), collapse = ", "), ".", call. = FALSE)

  ## ---- which candidates ----------------------------------------------------
  dat <- if (is.null(data)) fit$model else data
  apriori <- ilm_mod_apriori(tl, x)
  if (is.null(moderators)) {
    moderators <- setdiff(main, c(x, apriori))
  } else {
    ## named as the data knows them, or as R writes them; checked in the
    ## first form and carried on in the second
    plain <- ilm_unbq(as.character(moderators))
    bad <- setdiff(plain, names(dat))
    if (length(bad))
      stop("moderator(s) not found: ", paste(bad, collapse = ", "),
           ". The model frame has: ", paste(names(dat), collapse = ", "),
           ". Pass `data =` for a moderator that is not in the model.",
           call. = FALSE)
    moderators <- ilm_bq(plain)
    if (x %in% moderators) {
      moderators <- setdiff(moderators, x)
      message("ilm_moderation(): `x` cannot moderate itself; '", ilm_unbq(x),
              "' dropped from the candidates.")
    }
    hit <- intersect(moderators, apriori)
    if (length(hit)) moderators <- setdiff(moderators, hit)
  }
  ## An interaction already in the model was specified A PRIORI, so it owes no
  ## multiplicity penalty -- charging it one would make a pre-registered
  ## moderator weaker for having been tested beside exploratory ones.
  if (length(apriori))
    message("ilm_moderation(): ",
            paste(sprintf("`%s:%s`", ilm_unbq(x), ilm_unbq(apriori)),
                  collapse = ", "),
            if (length(apriori) == 1L) " is" else " are",
            " already in the model, so ",
            if (length(apriori) == 1L) "it was" else "they were",
            " specified a priori and ",
            if (length(apriori) == 1L) "is" else "are",
            " not charged a multiplicity penalty here. Test ",
            if (length(apriori) == 1L) "it" else "them",
            " with ilm_anova().")
  if (!length(moderators))
    stop("no candidate moderators left to test. The model's other fixed ",
         "effects are ", paste(ilm_unbq(setdiff(main, x)), collapse = ", "),
         if (length(apriori))
           paste0(", and ", paste(ilm_unbq(apriori), collapse = ", "),
                  " already interact with `", ilm_unbq(x), "`") else "",
         ".", call. = FALSE)

  ## ---- test each candidate -------------------------------------------------
  pb <- ilm_progress(length(moderators) * (if (split) 2L else 1L), progress,
                     "testing moderators")
  tick <- 0L
  if (split) set.seed(seed)
  idx <- if (split) sample(nrow(dat), floor(nrow(dat) / 2)) else NULL
  fit_dat <- if (split) dat[idx, , drop = FALSE] else dat
  res <- lapply(moderators, function(m) {
    tick <<- tick + 1L; pb$tick(tick)
    ilm_mod_one(fit, x, m, fit_dat, apriori, keep = TRUE)
  })
  names(res) <- moderators
  pv <- vapply(res, function(z) z$p, 0)
  dfv <- vapply(res, function(z) z$df, 0)
  stat <- vapply(res, function(z) z$statistic, 0)

  if (split) {
    ## the winner is chosen on the first half and tested on the second, which
    ## has not been used for anything yet
    pick <- names(which.min(pv))
    tick <- tick + 1L; pb$tick(tick)
    ho <- ilm_mod_one(fit, x, pick, dat[-idx, , drop = FALSE], apriori,
                      keep = TRUE)
    pb$done()
    out <- data.frame(moderator = ilm_unbq(pick), df = ho$df,
                      statistic = ho$statistic,
                      p = ho$p, p_adj = ho$p, row.names = NULL,
                      stringsAsFactors = FALSE)
    return(structure(out, class = c("ilm_moderation", "data.frame"),
                     x = ilm_unbq(x), adjust = "none (honest split)", split = TRUE,
                     n_split = c(length(idx), nrow(dat) - length(idx)),
                     apriori = ilm_unbq(apriori), fits = list(ho$fit),
                     formulas = ho$formula,
                     n_candidates = length(moderators)))
  }
  pb$done()
  out <- data.frame(moderator = ilm_unbq(moderators), df = dfv, statistic = stat,
                    p = pv, p_adj = stats::p.adjust(pv, adjust),
                    row.names = NULL, stringsAsFactors = FALSE)
  o <- order(out$p_adj, out$p)
  out <- out[o, , drop = FALSE]; row.names(out) <- NULL
  res <- res[o]
  keep <- seq_len(min(as.integer(n_keep), length(res)))
  structure(out, class = c("ilm_moderation", "data.frame"),
            x = ilm_unbq(x), adjust = adjust, split = FALSE,
            apriori = ilm_unbq(apriori),
            fits = lapply(res[keep], `[[`, "fit"),
            formulas = vapply(res, function(z) z$formula, ""),
            n_candidates = length(moderators))
}

## Interactions with `x` that are ALREADY in the model: specified a priori.
#' @keywords internal
#' @noRd
ilm_mod_apriori <- function(tl, x) {
  it <- grep(":", tl, value = TRUE)
  if (!length(it)) return(character(0))
  parts <- strsplit(it, ":", fixed = TRUE)
  keep <- vapply(parts, function(p) x %in% p, TRUE)
  unique(unlist(lapply(parts[keep], function(p) setdiff(p, x))))
}

## One candidate: refit with x * m and take the JOINT test of the block.
##
## The refit goes through the fit's own CALL with only the formula changed, so
## family, ziformula, dispformula, ar, weights, design, contrasts and reml all
## come along. Rebuilding ilm_model(f, data = ...) by hand drops every one of
## them silently -- a zero-inflated model would be tested without its zero
## part and nothing would say so.
#' @keywords internal
#' @noRd
ilm_mod_one <- function(fit, x, m, dat, apriori, keep = FALSE) {
  ## Each candidate is tested on equal footing: the base model is the user's,
  ## minus any OTHER candidate's interaction, plus this one's.
  tl <- fit$term_labels
  drop_it <- grep(":", tl, value = TRUE)
  drop_it <- drop_it[vapply(strsplit(drop_it, ":", fixed = TRUE),
                            function(p) x %in% p, TRUE)]
  base <- setdiff(tl, c(drop_it, m, x))
  rhs <- c(paste0(x, " * ", m), base)
  ## The response as it was WRITTEN, not its first variable: all.vars() would
  ## turn log(y) into y and refit the moderation models on a different
  ## outcome, and would drop the backticks from `my y`.
  resp <- ilm_term_text(fit$formula[[2L]])
  ## random-effect bars and smooths are carried across verbatim, each on one
  ## line with its backticks -- a bare deparse() splits a long bar in two
  bars <- if (length(fit$bars))
    vapply(fit$bars, function(b) paste0("(", ilm_term_text(b), ")"), "") else
      character(0)
  f <- stats::as.formula(paste(resp, "~", paste(c(rhs, bars), collapse = " + ")),
                         env = environment(fit$formula))
  cl <- fit$call
  cl$formula <- f
  cl$data <- quote(.ilm_mod_data)
  cl$verbose <- FALSE
  env <- new.env(parent = parent.frame())
  assign(".ilm_mod_data", dat, envir = env)
  nf <- try(suppressMessages(suppressWarnings(eval(cl, env))), silent = TRUE)
  fail <- list(p = NA_real_, df = NA_real_, statistic = NA_real_,
               formula = paste(deparse(f), collapse = " "), fit = NULL)
  if (inherits(nf, "try-error")) return(fail)
  a <- try(suppressWarnings(ilm_anova(nf)), silent = TRUE)
  if (inherits(a, "try-error")) return(fail)
  rn <- paste0(x, ":", m)
  if (!rn %in% rownames(a)) rn <- paste0(m, ":", x)
  if (!rn %in% rownames(a)) return(fail)
  list(p = a[rn, ncol(a)], df = a[rn, "Df"],
       statistic = a[rn, 2L],
       formula = paste(deparse(f), collapse = " "),
       fit = if (keep) nf else NULL)
}

#' @export
print.ilm_moderation <- function(x, ...) {
  d <- as.data.frame(x); class(d) <- "data.frame"
  ex <- attr(x, "x")
  if (isTRUE(attr(x, "split"))) {
    ns <- attr(x, "n_split")
    cat(sprintf("Moderation of `%s`, honest split (%d chose, %d tested)\n",
                ex, ns[1], ns[2]))
    cat(sprintf("  strongest of %d candidates on the first half, tested on the second\n",
                attr(x, "n_candidates")))
  } else {
    cat(sprintf("Moderation of `%s`: %d candidate%s, %s-adjusted\n", ex,
                attr(x, "n_candidates"),
                if (attr(x, "n_candidates") == 1L) "" else "s",
                attr(x, "adjust")))
  }
  d$statistic <- round(d$statistic, 3)
  d$p <- signif(d$p, 3); d$p_adj <- signif(d$p_adj, 3)
  print(d, row.names = FALSE)
  ap <- attr(x, "apriori")
  if (length(ap))
    cat("\n  ", paste(sprintf("%s:%s", ex, ap), collapse = ", "),
        if (length(ap) == 1L) " is" else " are",
        " in the model already, so ", if (length(ap) == 1L) "it was" else
          "they were", " specified a priori\n",
        "  and ", if (length(ap) == 1L) "is" else "are",
        " not tested here. ilm_anova() is where ",
        if (length(ap) == 1L) "it belongs" else "they belong", ".\n", sep = "")
  hit <- d$moderator[d$p_adj < 0.05]
  cat("\n  ", if (length(hit))
    paste0("`", paste(hit, collapse = "`, `"), "`: the effect of `", ex,
           "` is not the same for everyone.")
    else paste0("No candidate survives adjustment. That is not evidence the ",
                "effect is constant --"),
    "\n", sep = "")
  if (!length(hit))
    cat("  a moderation search is underpowered by construction, and a ",
        ifelse(max(d$df, na.rm = TRUE) > 1, "multi-df candidate costs more still.\n",
               "null here is a weak null.\n"), sep = "")
  cat("  ilm_plot_moderation() to see one. Each test is the joint test of the\n",
      "  interaction block, so `df` is what it cost.\n", sep = "")
  invisible(x)
}
