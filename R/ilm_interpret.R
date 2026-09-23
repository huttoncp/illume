## ---------------------------------------------------------------------------
## Turning a fit into sentences.
##
## Three rules hold this together, and they are the difference between a useful
## interpreter and a confident liar:
##
##  1. The prose is TEMPLATED, never generated. The same fit gives the same
##     words, and each sentence is testable.
##  2. The effects are reported on the RESPONSE scale. A log-odds is not
##     something to report to anybody, and an odds ratio is misread as a risk
##     ratio by almost everyone, so the change in probability is computed and
##     given as well.
##  3. Nothing is said about bias that a diagnostic did not measure. Every such
##     sentence traces to a check illume actually ran, with its verdict and the
##     remedy that check names. Where no diagnostic speaks, the interpreter is
##     silent rather than reassuring.
##
## And causal language is LICENSED, not assumed. "increases" requires that the
## model came from a design that identifies an effect -- a DAG with a valid
## adjustment set, a difference in differences, a regression discontinuity.
## Otherwise it is "associated with", every time.
## ---------------------------------------------------------------------------

#' @keywords internal
#' @noRd
ilm_wrap <- function(x, width = 76L, indent = "") {
  paste0(indent, strwrap(x, width = width - nchar(indent)), collapse = "\n")
}

#' @keywords internal
#' @noRd
ilm_fmt <- function(v, digits = 3) formatC(v, format = "f", digits = digits)

## A p-value is reported as a number and as a strength of evidence. The
## dichotomy on its own is what makes people write "no effect" for a wide
## interval, so the wording is graded and the interval always travels with it.
#' @keywords internal
#' @noRd
ilm_evidence <- function(p) {
  if (is.na(p)) return("of unknown strength")
  if (p < 0.001) "very strong evidence"
  else if (p < 0.01) "strong evidence"
  else if (p < 0.05) "moderate evidence"
  else if (p < 0.1) "weak evidence, not conventionally significant,"
  else "little evidence"
}

## ---- average marginal effects ----------------------------------------------

#' Average marginal effect on the response scale
#'
#' What a coefficient means where the response lives. For a gaussian model that
#' is the coefficient itself; for anything with a link it is not, and the
#' difference matters: an odds ratio of 2 can be a 3-point change in probability
#' or a 20-point one depending on where the data sit.
#'
#' For a numeric predictor this is the average of the derivative of the fitted
#' mean with respect to it, taken over the rows actually observed. For a factor
#' it is the average change in fitted mean from moving every row to that level
#' from the reference, which is a contrast rather than a derivative.
#'
#' Standard errors come from the delta method: the effect is differentiated
#' numerically with respect to the parameters and combined with their covariance.
#' The parameter vector includes the covariance parameters, because a
#' population-averaged prediction depends on them and pretending otherwise would
#' understate the uncertainty.
#'
#' @param object An [ilm_model()].
#' @param terms Which predictors. Default is every fixed-effect term.
#' @param eps Relative step for the numerical derivatives.
#' @return A data frame with `term`, `level`, `estimate`, `se`, `lower`,
#'   `upper`, and `kind` (`"slope"` or `"contrast"`). For an outcome with
#'   categories -- multinomial or ordinal -- there is also a `category` column
#'   and one row per category: the effect on the probability of each. Those
#'   rows sum to zero for every term, because the probabilities sum to one.
#' @seealso [ilm_interpret()], and the `marginaleffects` package, which does
#'   this and a great deal more once [ilm_register_marginaleffects()] is called.
#' @examples
#' set.seed(1); n <- 400
#' d <- data.frame(x = rnorm(n), g = factor(sample(c("a", "b"), n, TRUE)))
#' d$y <- rbinom(n, 1, plogis(0.4 * d$x + 0.6 * (d$g == "b")))
#' fit <- ilm_model(y ~ x + g, data = d, family = "binomial", verbose = FALSE)
#' ilm_ame(fit)
#' @export
ilm_ame <- function(object, terms = NULL, eps = 1e-4) {
  if (!inherits(object, "ilm_model"))
    stop("`object` must be a fitted ilm_model, not ", class(object)[1],
         call. = FALSE)
  mf <- object$model
  if (is.null(mf))
    stop("the fit did not keep its model frame, so marginal effects cannot ",
         "be computed", call. = FALSE)
  tl <- attr(stats::terms(object), "term.labels")
  ## only terms that are a bare variable: a marginal effect for `poly(x, 3)` or
  ## an interaction is a different question and is not answered by pretending.
  ## Compared as PLAIN names: a label keeps the backticks a column name does
  ## not have, so `x 1` used to match nothing and quietly got no effect.
  vars <- ilm_unbq(tl)[ilm_unbq(tl) %in% names(mf)]
  if (!is.null(terms)) vars <- intersect(vars, ilm_unbq(terms))
  if (!length(vars)) return(NULL)

  ## An outcome with categories -- multinomial or ordinal -- has a probability
  ## for each one, and so an effect on each one; they sum to zero, because the
  ## probabilities sum to one. This used to take the LAST column of those
  ## probabilities and report it as THE effect, with nothing to say which
  ## category it was, and ilm_interpret() then quoted that one number beside
  ## every category's coefficient.
  p1 <- suppressWarnings(stats::predict(object, newdata = mf[1L, , drop = FALSE],
                                        type = "response"))
  cats <- if (is.matrix(p1) && ncol(p1) > 1L) colnames(p1) else NULL
  mu <- function(dd) {
    p <- suppressWarnings(stats::predict(object, newdata = dd,
                                         type = "response"))
    if (!is.null(cats)) as.matrix(p)
    else if (is.matrix(p)) p[, ncol(p)] else as.numeric(p)
  }
  avg <- function(z) if (is.null(cats)) mean(z) else unname(colMeans(z))
  mk <- function(v, lev, kind, val) {
    if (is.null(cats))
      data.frame(term = v, level = lev, kind = kind, value = val,
                 stringsAsFactors = FALSE)
    else
      data.frame(term = v, level = lev, category = cats, kind = kind,
                 value = val, stringsAsFactors = FALSE)
  }
  ## one number per (term, level), as a function of the parameter vector
  eff_of <- function(obj) {
    old <- object
    object <<- obj
    on.exit(object <<- old)
    out <- list()
    for (v in vars) {
      x <- mf[[v]]
      if (is.numeric(x)) {
        h <- eps * max(stats::sd(x, na.rm = TRUE), 1e-8)
        d1 <- mf; d1[[v]] <- x + h
        d0 <- mf; d0[[v]] <- x - h
        out[[length(out) + 1L]] <- mk(v, NA_character_, "slope",
                                      avg((mu(d1) - mu(d0)) / (2 * h)))
      } else {
        f <- factor(x); lv <- levels(f)
        if (length(lv) < 2L) next
        d0 <- mf; d0[[v]] <- factor(lv[1], levels = lv)
        m0 <- mu(d0)
        for (l in lv[-1]) {
          d1 <- mf; d1[[v]] <- factor(l, levels = lv)
          out[[length(out) + 1L]] <- mk(v, l, "contrast", avg(mu(d1) - m0))
        }
      }
    }
    if (!length(out)) NULL else do.call(rbind, out)
  }

  base <- eff_of(object)
  if (is.null(base)) return(NULL)
  b0 <- suppressWarnings(stats::coef(object, full = TRUE))
  V <- tryCatch(suppressWarnings(stats::vcov(object, full = TRUE)),
                error = function(e) NULL)
  se <- rep(NA_real_, nrow(base))
  if (!is.null(V) && length(b0) == nrow(V)) {
    G <- matrix(NA_real_, nrow(base), length(b0))
    for (j in seq_along(b0)) {
      hj <- eps * max(abs(b0[j]), 1)
      bp <- b0; bp[j] <- bp[j] + hj
      bm <- b0; bm[j] <- bm[j] - hj
      ep <- tryCatch(eff_of(ilm_rebuild(object, bp)), error = function(e) NULL)
      em <- tryCatch(eff_of(ilm_rebuild(object, bm)), error = function(e) NULL)
      if (is.null(ep) || is.null(em)) next
      G[, j] <- (ep$value - em$value) / (2 * hj)
    }
    ok <- apply(is.finite(G), 1L, all)
    vv <- rep(NA_real_, nrow(base))
    if (any(ok)) {
      gq <- G[ok, , drop = FALSE]
      vv[ok] <- rowSums((gq %*% V) * gq)
    }
    se <- sqrt(pmax(vv, 0))
  }
  crit <- if (isTRUE(object$exact_df)) stats::qt(0.975, object$resid_df)
          else stats::qnorm(0.975)
  out <- data.frame(term = base$term, level = base$level, kind = base$kind,
                    estimate = base$value, se = se,
                    lower = base$value - crit * se, upper = base$value + crit * se,
                    stringsAsFactors = FALSE, row.names = NULL)
  if (!is.null(cats))
    out <- cbind(out[1:2], category = base$category, out[-(1:2)],
                 stringsAsFactors = FALSE)
  out
}

## ---- the interpreter -------------------------------------------------------

#' @keywords internal
#' @noRd
ilm_scale_words <- function(fam) {
  switch(fam,
    gaussian = list(unit = "units of the response", link = "identity",
                    ratio = NULL),
    binomial = list(unit = "percentage points of probability", link = "logit",
                    ratio = "odds ratio"),
    poisson = list(unit = "counts", link = "log", ratio = "rate ratio"),
    nbinom = list(unit = "counts", link = "log", ratio = "rate ratio"),
    multinomial = list(unit = "percentage points of probability",
                       link = "logit", ratio = "odds ratio"),
    ## an ordered outcome's coefficients shift the latent scale the thresholds
    ## cut up; exp() of a logit one is the odds of a HIGHER category
    ordinal = list(unit = "percentage points of probability",
                   link = "cumulative logit", ratio = "cumulative odds ratio"),
    ordinal_probit = list(unit = "percentage points of probability",
                          link = "cumulative probit", ratio = NULL),
    ordinal_cloglog = list(unit = "percentage points of probability",
                           link = "cumulative complementary log-log",
                           ratio = NULL),
    beta = list(unit = "units of the proportion", link = "logit", ratio = NULL),
    weibull = , lognormal = , loglogistic =
      list(unit = "units of log time", link = "log", ratio = "time ratio"),
    rp = , rp_odds = , rp_normal =
      list(unit = "units of the linear predictor", link = "log cumulative hazard",
           ratio = "hazard ratio"),
    list(unit = "units of the response", link = "identity", ratio = NULL))
}

#' Interpret a fitted model in words
#'
#' Writes out what a fit says: each effect on the scale the response is measured
#' on, how strong the evidence is, what the diagnostics found, and what that
#' implies for taking the estimates at face value. Intended for reporting to
#' people who will not read a coefficient table.
#'
#' @section Causal language is licensed, not assumed:
#'
#' A regression coefficient is an association. This says "associated with"
#' unless the object carries a design that identifies an effect -- an
#' [ilm_dag_model()] with a valid adjustment set, an [ilm_did()], an
#' [ilm_rdd()] -- in which case it says so and names what licenses it. That is
#' the single most common error in reporting a model, and the one place an
#' automatic interpreter could do real damage, so it is handled explicitly.
#'
#' @section What it will not do:
#'
#' It says nothing about bias that a diagnostic did not measure. Every such
#' sentence traces to a check that was actually run, and names the remedy that
#' check names. Where no diagnostic has spoken, it is silent rather than
#' reassuring: a clean report from this function means the checks that ran
#' passed, not that the model is right.
#'
#' The prose is templated. The same fit gives the same words every time.
#'
#' @param object An [ilm_model()], [ilm_dag_model()], [ilm_did()] or
#'   [ilm_rdd()]; or the result of [ilm_power()] or [ilm_power_design()], which
#'   is written up as a power analysis, or of [ilm_contrast()], written up as
#'   the comparisons it makes.
#' @param causal Force causal or associational language. `NULL` decides from
#'   the design, which is what you want.
#' @param ame Report average marginal effects on the response scale. Costs a
#'   delta-method calculation; worth it for anything with a link function.
#' @param digits Rounding.
#' @param ... Unused.
#' @return An object of class `"ilm_interpretation"`: a list of sections, which
#'   `print()` renders as wrapped text.
#' @seealso [ilm_ame()], [ilm_appraise()].
#' @examples
#' set.seed(1); n <- 300
#' d <- data.frame(x = rnorm(n), g = factor(sample(c("a", "b"), n, TRUE)))
#' d$y <- 0.5 * d$x + 0.8 * (d$g == "b") + rnorm(n)
#' fit <- ilm_model(y ~ x + g, data = d, family = "gaussian", verbose = FALSE)
#' ilm_interpret(fit, ame = FALSE)
#' @export
ilm_interpret <- function(object, ...) UseMethod("ilm_interpret")

#' @rdname ilm_interpret
#' @export
ilm_interpret.ilm_model <- function(object, causal = NULL, ame = TRUE,
                                    digits = 3, ...) {
  fam <- if (!is.null(object$family)) object$family$name else "gaussian"
  sw <- ilm_scale_words(fam)
  is_causal <- isTRUE(causal)
  link_word <- if (is_causal) "leads to" else "is associated with"

  n <- tryCatch(stats::nobs(object), error = function(e) NA_integer_)
  tl <- attr(stats::terms(object), "term.labels")
  fixed <- setdiff(tl, grep("\\|", tl, value = TRUE))
  sec <- list()

  ## ---- what was fitted -----------------------------------------------------
  hdr <- sprintf("%s %s model of %s, fitted to %s observations.",
                 ilm_cap(ilm_article(fam)), fam, deparse(object$formula[[2]]),
                 if (is.na(n)) "an unknown number of" else format(n, big.mark = ","))
  re <- names(object$re)
  if (length(re))
    hdr <- paste(hdr, sprintf(
      "Observations are grouped by %s, and that grouping is modelled, so the estimates below already allow for it.",
      paste(re, collapse = " and ")))
  if (!identical(sw$link, "identity"))
    hdr <- paste(hdr, sprintf(
      "The model works on the %s scale, so its coefficients are not in %s; the effects below are converted.",
      sw$link, sw$unit))
  ## a likelihood the user did not choose is worth saying they did not
  if (!is.null(object$family_inferred))
    hdr <- paste(hdr, sprintf(
      "The %s family was read off the response (%s) rather than specified; if that is not how %s arises, refit with `family` set.",
      fam, object$family_inferred, deparse(object$formula[[2]])))
  sec$model <- hdr

  ## ---- the effects ---------------------------------------------------------
  ct <- ilm_coef_table(object)
  crit <- if (isTRUE(object$exact_df)) stats::qt(0.975, object$resid_df)
          else stats::qnorm(0.975)
  am <- if (isTRUE(ame)) tryCatch(ilm_ame(object), error = function(e) NULL) else NULL

  respname <- deparse(object$formula[[2]])
  ## A multinomial fit has one coefficient per predictor PER CATEGORY -- the
  ## first C categories; the last is minus their sum -- laid out category by
  ## category. This loop used to walk the design columns only, so it
  ## described the first category's coefficients and never the others, and
  ## quoted one category's marginal effect beside all of them.
  multi <- identical(fam, "multinomial")
  pX <- ncol(object$X)
  cats <- if (multi) object$ylevels[seq_len(object$C)] else NA_character_
  xn <- colnames(object$X)
  lines <- character()
  for (v in fixed) {
    k <- ilm_term_cols(object, v)
    if (!length(k)) next
    ## the label keeps its backticks; the model frame knows the plain name
    vn <- ilm_unbq(v)
    xv <- if (!is.null(object$model) && vn %in% names(object$model))
            object$model[[vn]] else NULL
    is_fac <- !is.null(xv) && !is.numeric(xv)
    for (cc in seq_along(cats)) for (i in k) {
      ii <- if (multi) (cc - 1L) * pX + i else i
      nm <- rownames(ct)[ii]
      est <- ct[ii, 1]; se <- ct[ii, 2]; p <- ct[ii, 4]
      lo <- est - crit * se; hi <- est + crit * se
      ev <- ilm_evidence(p)
      ## the level this coefficient stands for, when the term is a factor --
      ## read off the DESIGN column, since a multinomial coefficient's own
      ## name starts with its category
      lvl <- if (is_fac) {
        x1 <- xn[i]
        if (startsWith(x1, v)) substring(x1, nchar(v) + 1L) else x1
      } else NA_character_
      subj <- if (is_fac && nzchar(lvl)) sprintf("being %s rather than %s", lvl,
                    levels(factor(xv))[1])
              else sprintf("a higher %s", vn)
      s <- if (multi)
        sprintf("%s: %s that %s %s %s odds of %s, relative to the average of the categories (estimate %s, 95%% interval %s to %s, p = %s).",
                nm, ev, subj, link_word,
                if (est >= 0) "higher" else "lower", sQuote(cats[cc], FALSE),
                ilm_fmt(est, digits), ilm_fmt(lo, digits),
                ilm_fmt(hi, digits), format.pval(p, digits = 2, eps = 1e-4))
      else
        sprintf("%s: %s that %s %s a %s %s of %s (estimate %s, 95%% interval %s to %s, p = %s).",
                nm, ev, subj, link_word,
                if (est >= 0) "higher" else "lower",
                if (isTRUE(object$ordinal)) "category" else "value", respname,
                ilm_fmt(est, digits), ilm_fmt(lo, digits),
                ilm_fmt(hi, digits), format.pval(p, digits = 2, eps = 1e-4))
      ## and what it means where the response lives
      if (!is.null(sw$ratio) && fam %in% c("binomial", "poisson", "nbinom"))
        s <- paste(s, sprintf("On the %s scale that is %s.", sw$ratio,
                              ilm_fmt(exp(est), digits)))
      if (identical(fam, "ordinal"))
        s <- paste(s, sprintf(
          "As a cumulative odds ratio that is %s: the odds of being in a higher category rather than a lower one, at every cut point, are multiplied by it.",
          ilm_fmt(exp(est), digits)))
      if (!is.null(am)) {
        ## match on the term AND the level, not on a name prefix: a vectorised
        ## grepl here silently used only the first level and dropped the rest
        j <- if (is_fac) which(am$term == vn & !is.na(am$level) & am$level == lvl)
             else which(am$term == vn & is.na(am$level))
        ## and on the category: its own for a multinomial coefficient; for an
        ## ordinal one, which shifts every category at once, the highest,
        ## named as such
        catj <- NULL
        if (!is.null(am$category)) {
          catj <- if (multi) cats[cc] else utils::tail(unique(am$category), 1L)
          j <- j[am$category[j] == catj]
        }
        if (length(j) == 1L && is.finite(am$estimate[j])) {
          sc <- if (fam %in% c("binomial", "multinomial") || !is.null(catj))
            100 else 1
          s <- paste(s, sprintf(
            "%s: %s %s%s on average (%s to %s).",
            if (is.null(catj)) "In the units of the response"
            else sprintf("In the probability of %s%s", sQuote(catj, FALSE),
                         if (multi) "" else ", the highest category"),
            if (am$estimate[j] >= 0) "an increase of" else "a decrease of",
            ilm_fmt(abs(am$estimate[j]) * sc, digits),
            if (sc == 100) " percentage points" else "",
            ilm_fmt(min(am$lower[j], am$upper[j]) * sc, digits),
            ilm_fmt(max(am$lower[j], am$upper[j]) * sc, digits)))
        }
      }
      if (p >= 0.05)
        s <- paste(s, "The interval includes zero, which means the data are",
                   "consistent with no effect -- not that there is none.")
      lines <- c(lines, s)
    }
  }
  sec$effects <- lines

  ## ---- what the checks said ------------------------------------------------
  ck <- object$checks
  dl <- character()
  if (!is.null(ck) && nrow(ck)) {
    ## a check that could not reach a verdict is not one that passed
    bad <- ck[ck$status %in% c("WARN", "FAIL", "INCONCLUSIVE"), , drop = FALSE]
    bnd <- any(ck$status == "BOUNDARY")
    if (!nrow(bad) && !bnd) {
      dl <- c(dl, sprintf(
        "All %d fitting checks passed: the optimiser converged, the gradient is at zero and the information matrix is usable. These say the fit is sound, not that the model is right -- for that, run ilm_appraise().",
        nrow(ck)))
    } else if (nrow(bad)) {
      for (i in seq_len(nrow(bad))) {
        s <- sprintf("%s -- %s: %s.", bad$status[i], bad$check[i],
                     bad$detail[i])
        if (nzchar(bad$cause[i])) s <- paste0(s, " ", ilm_cap(bad$cause[i]), ".")
        if (nzchar(bad$suggestion[i]))
          s <- paste0(s, " What to do: ", bad$suggestion[i], ".")
        dl <- c(dl, s)
      }
      if (any(bad$status == "FAIL"))
        dl <- c(dl, paste(
          "At least one check FAILED. Standard errors and intervals above rest",
          "on the model being adequate, so treat them as provisional until",
          "that is resolved."))
    }
    ## a covariance at its boundary: which parts of the fit stand, in the same
    ## words summary() uses
    if (bnd)
      dl <- c(dl, paste("BOUNDARY --",
                        ilm_trust_text(object$hessian_held, ilm_boundary_at(object),
                                       avoided = identical(object$boundary, "avoid"))))
  }
  sec$diagnostics <- dl

  ## ---- how far to trust it -------------------------------------------------
  cav <- character()
  if (!is_causal)
    cav <- c(cav, paste(
      "These are associations. Nothing here rules out a common cause of a",
      "predictor and the response, so an effect could differ in size or sign",
      "from what is reported. To say more you need a design that identifies",
      "one: ilm_dag_model() with an adjustment set, ilm_did(), ilm_rdd()."))
  if (length(re) && !identical(sw$link, "identity"))
    cav <- c(cav, paste(
      "With a link function and a random effect, a coefficient is conditional",
      "on the group -- the effect for units in the same group -- while the",
      "marginal effects above are averaged over groups. The two answer",
      "different questions and will not match."))
  if (isTRUE(object$exact_df))
    cav <- c(cav, paste(
      "This model has no random or smooth terms, so its t and F tests are",
      "exact rather than large-sample approximations."))
  ## Few groups: the Wald tests' reference is a large-sample one, and for an
  ## effect that varies between groups it is the NUMBER of groups that counts.
  ## Measured on multinomial designs with 60 groups, a nominal 5% test of a
  ## between-group effect rejected a true null 5.8% and 7.8% of the time.
  gk <- which(vapply(object$re, function(e) !identical(e$kind, "basis"), TRUE))
  if (length(gk) && !isTRUE(object$exact_df)) {
    ng <- min(object$nlk[gk])
    if (is.finite(ng) && ng < 100L)
      cav <- c(cav, sprintf(
        "The tests above use a large-sample reference, and with %d groups of `%s` one for an effect that varies between groups can run somewhat liberal. For a term the conclusions rest on, ilm_pb_lrt() calibrates the p-value by simulation%s.",
        ng, names(object$re)[gk][which.min(object$nlk[gk])],
        if (identical(fam, "gaussian"))
          ", and ilm_denom_df() gives finite degrees of freedom" else ""))
  }
  sec$caveats <- cav

  structure(list(sections = sec, family = fam, causal = is_causal,
                 object_class = class(object)[1]),
            class = "ilm_interpretation")
}

#' @keywords internal
#' @noRd
ilm_cap <- function(s) {
  if (!nzchar(s)) return(s)
  paste0(toupper(substring(s, 1, 1)), substring(s, 2))
}

#' @rdname ilm_interpret
#' @export
ilm_interpret.ilm_dag_model <- function(object, causal = NULL, ame = TRUE,
                                        digits = 3, ...) {
  if (!isTRUE(object$identified)) {
    return(structure(list(sections = list(
      model = sprintf(
        "The effect of %s on %s is NOT identified by these data. No set of the measured variables closes the back-door paths in the graph, so no regression on this data set estimates it.",
        object$exposure, object$outcome),
      caveats = paste(
        "This is a conclusion, not a failure. Fitting a regression anyway would",
        "produce a number, and that number would not answer the question.",
        "Measure a confounder, or revise the graph.")),
      family = NULL, causal = FALSE, object_class = "ilm_dag_model"),
      class = "ilm_interpretation"))
  }
  ## a valid adjustment set is what licenses causal language
  lic <- is.null(causal) || isTRUE(causal)
  base <- ilm_interpret(object$fits[[1]], causal = lic, ame = ame,
                        digits = digits)
  s <- base$sections
  zz <- object$sets[[1]]
  s$model <- paste(sprintf(
    "The effect of %s on %s, identified by adjusting for %s -- a minimal sufficient set under the supplied causal graph.",
    object$exposure, object$outcome,
    if (length(zz)) paste(zz, collapse = ", ") else "nothing"), s$model)

  if (!is.null(object$dag_test)) {
    v <- attr(object$dag_test, "verdict")
    s$diagnostics <- c(sprintf(
      "The graph itself was tested against the data on %d implied conditional independence(s): %s.%s",
      nrow(object$dag_test), v,
      switch(v,
        OK = " The data are consistent with the graph, which supports but does not prove it.",
        FAIL = " The data CONTRADICT the graph. Either an arrow is missing or a variable does not measure what its name supposes, and the adjustment set follows from the graph, so the estimate inherits the problem.",
        WARN = " A claim reached significance but stayed below the effect size worth acting on.",
        "")), s$diagnostics)
  }
  if (length(object$sets) > 1L && !is.null(object$effects)) {
    e <- object$effects
    sp <- diff(range(e$estimate))
    s$caveats <- c(sprintf(
      "%d different adjustment sets identify this effect and all were fitted; the estimates span %s (%s to %s). They target the same quantity, so agreement supports the graph and disagreement is worth more attention than any single number.",
      length(object$sets), ilm_fmt(sp, digits),
      ilm_fmt(min(e$estimate), digits), ilm_fmt(max(e$estimate), digits)),
      s$caveats)
  }
  if (lic)
    s$caveats <- c(s$caveats, paste(
      "Causal language here rests entirely on the graph being right. It is an",
      "assumption you supplied, not something the data established."))
  base$sections <- s; base$causal <- lic
  base$object_class <- "ilm_dag_model"
  base
}

#' @rdname ilm_interpret
#' @export
ilm_interpret.ilm_did <- function(object, causal = NULL, ame = FALSE,
                                  digits = 3, ...) {
  a <- object$att
  sec <- list()
  sec$model <- sprintf(
    "A difference-in-differences comparison of %s across %d treated and %d control units, with treatment starting at %s = %s. Each unit's own level is modelled, so the estimate is a change relative to the control group's change rather than a level difference.",
    object$y, object$n_treated, object$n_control, object$time, object$treat_time)
  sec$effects <- sprintf(
    "The treatment %s a change of %s in %s (95%% interval %s to %s, p = %s): %s.",
    if (is.null(causal) || isTRUE(causal)) "produced" else "is associated with",
    ilm_fmt(a$estimate, digits), object$y, ilm_fmt(a$lower, digits),
    ilm_fmt(a$upper, digits), format.pval(a$p_value, digits = 2, eps = 1e-4),
    ilm_evidence(a$p_value))
  dl <- character()
  if (!is.null(object$parallel)) {
    st <- object$parallel$status
    dl <- c(dl, switch(st,
      OK = sprintf("Parallel trends holds as far as it can be checked: over %d pre-treatment periods the two groups' slopes differed by %s, which is not distinguishable from zero (p = %s). That supports the design without proving it, since the assumption concerns a period that never happened.",
                   object$parallel$n_pre, ilm_fmt(object$parallel$diff_slope, digits),
                   format.pval(object$parallel$p_value, digits = 2, eps = 1e-4)),
      FAIL = sprintf("Parallel trends FAILS: the groups' pre-treatment slopes differed by %s (p = %s). They were already diverging before the treatment, so part of the estimate above is that divergence continuing rather than an effect. The estimate should not be read as a treatment effect until this is addressed.",
                     ilm_fmt(object$parallel$diff_slope, digits),
                     format.pval(object$parallel$p_value, digits = 2, eps = 1e-4)),
      UNTESTED = sprintf("Parallel trends could NOT be checked: there %s only %d pre-treatment period%s. The assumption is doing all the work and no part of it has been verified.",
                         if (object$parallel$n_pre == 1L) "is" else "are",
                         object$parallel$n_pre,
                         if (object$parallel$n_pre == 1L) "" else "s"),
      ""))
  }
  if (!is.null(object$event)) {
    bad <- sum(object$event$rel_time < 0 &
               (object$event$lower > 0 | object$event$upper < 0))
    dl <- c(dl, sprintf(
      "The event study covers %d periods around treatment; %d pre-treatment coefficient%s exclude zero. Under parallel trends they should all sit near it, so %s",
      nrow(object$event), bad, if (bad == 1L) "" else "s",
      if (bad == 0L) "this is what the design predicts."
      else "this is a further sign the groups were not moving together."))
  }
  sec$diagnostics <- dl
  cav <- character()
  if (isTRUE(object$staggered))
    cav <- c(cav, paste(
      "Adoption is STAGGERED. When units are treated at different times and",
      "the effect varies, this pooled estimate is not an average treatment",
      "effect: already-treated units act as controls for later-treated ones.",
      "Treat the number above as indicative only."))
  cav <- c(cav, paste(
    "Difference in differences identifies an effect only if the groups would",
    "have moved together without the treatment. That cannot be verified, only",
    "made plausible."))
  if (!object$ar)
    cav <- c(cav, paste(
      "Repeated observations of the same unit are correlated, which makes",
      "standard errors too small if unmodelled. A unit random intercept is",
      "fitted; if the correlation decays with time rather than being constant,",
      "ar = TRUE adds AR(1) on top."))
  sec$caveats <- cav
  structure(list(sections = sec, family = object$family,
                 causal = is.null(causal) || isTRUE(causal),
                 object_class = "ilm_did"), class = "ilm_interpretation")
}

#' @rdname ilm_interpret
#' @export
ilm_interpret.ilm_rdd <- function(object, causal = NULL, ame = FALSE,
                                  digits = 3, ...) {
  j <- object$jump
  sec <- list()
  sec$model <- sprintf(
    "A regression discontinuity in %s at %s = %s, fitted locally on each side within a bandwidth of %s (%d observations below the cutoff, %d at or above).",
    object$y, object$running, ilm_fmt(object$cutoff, digits),
    ilm_fmt(object$h, 4), object$n_below, object$n_above)
  sec$effects <- sprintf(
    "At the cutoff, %s jumps by %s (95%% interval %s to %s, p = %s): %s. This is a local effect -- it applies to units near the cutoff and says nothing about units far from it.",
    object$y, ilm_fmt(j$estimate, digits), ilm_fmt(j$lower, digits),
    ilm_fmt(j$upper, digits), format.pval(j$p_value, digits = 2, eps = 1e-4),
    ilm_evidence(j$p_value))
  dl <- character()
  dl <- c(dl, switch(object$density$status,
    OK = "The density of the running variable is continuous at the cutoff, so there is no sign that units moved themselves across it.",
    FAIL = "The density of the running variable JUMPS at the cutoff. That is what manipulation looks like: units just above may differ from units just below in ways the design assumes they do not, and no choice of bandwidth repairs it.",
    "The density could not be checked."))
  if (!is.null(object$balance)) {
    nb <- sum(object$balance$status == "FAIL")
    dl <- c(dl, sprintf(
      "Of %d covariate%s checked for a jump at the cutoff, %d jump%s.%s",
      nrow(object$balance), if (nrow(object$balance) == 1L) "" else "s", nb,
      if (nb == 1L) "s" else "",
      if (nb) " A covariate that jumps means something other than treatment changes at the cutoff, and the estimate absorbs it." else
        " Nothing besides treatment appears to change there."))
  }
  if (!is.null(object$placebo)) {
    np <- sum(object$placebo$status == "FAIL")
    dl <- c(dl, sprintf("Of %d placebo cutoffs where nothing happens, %d show a jump.%s",
      nrow(object$placebo), np,
      if (np) " Finding jumps where there are none suggests the method is picking up noise here." else ""))
  }
  if (!is.null(object$bandwidth)) {
    bw <- object$bandwidth
    dl <- c(dl, sprintf(
      "Across bandwidths from %sx to %sx the estimate runs %s to %s, and the interval %s zero throughout.",
      ilm_fmt(min(bw$multiplier), 2), ilm_fmt(max(bw$multiplier), 2),
      ilm_fmt(min(bw$estimate), digits), ilm_fmt(max(bw$estimate), digits),
      if (all(bw$lower > 0) || all(bw$upper < 0)) "excludes" else "does not consistently exclude"))
  }
  sec$diagnostics <- dl
  cav <- character()
  if (!isTRUE(object$sharp))
    cav <- c(cav, paste(
      "Assignment is FUZZY: crossing the cutoff changes the chance of",
      "treatment without determining it. The jump above is therefore the",
      "effect of being eligible, not of being treated. Converting between the",
      "two needs an instrumental-variables estimator, which illume does not",
      "have."))
  cav <- c(cav, paste(
    "The interval is the ordinary one for a weighted local fit and does not",
    "carry the bias correction that a wide bandwidth calls for. The rdrobust",
    "package implements those intervals."))
  sec$caveats <- cav
  structure(list(sections = sec, family = object$family,
                 causal = is.null(causal) || isTRUE(causal),
                 object_class = "ilm_rdd"), class = "ilm_interpretation")
}

#' @rdname ilm_interpret
#' @param target For a power analysis, the power a study is to reach.
#' @export
ilm_interpret.ilm_power <- function(object, causal = NULL, ame = FALSE,
                                    digits = 2, target = 0.8, ...) {
  d <- as.data.frame(object); class(d) <- "data.frame"
  unit <- attr(object, "unit")
  scol <- if (!is.null(d$n_unit)) "n_unit" else "n"
  who <- if (!is.null(d$n_unit) && !is.null(unit) && !identical(unit, "row"))
    paste0(unit, "s") else "observations"
  test <- attr(object, "test"); if (is.null(test)) test <- "Wald z"
  joint <- identical(attr(object, "effect_scale"), "multiple")
  pct <- function(v) paste0(ilm_fmt(100 * v, 0), "%")
  sec <- list()
  sec$model <- paste(
    sprintf("Power for %s, estimated by simulating %d studies at each size and analysing each the way the analysis will be: %s, at a two-sided level of %s.",
            attr(object, "term"), d$sims[1L],
            if (joint) sprintf("a joint %s test of its %d coefficients", test,
                               length(attr(object, "coefs")))
            else sprintf("a %s test", test),
            format(attr(object, "alpha"))),
    if (isTRUE(attr(object, "redrawn")))
      "Each study is a fresh draw of the planned design: the allocation balanced as the protocol would balance it, and new participants every time."
    else if (isTRUE(attr(object, "grouped")))
      "Each study resamples the fitted data's clusters whole, so the number of groups moves with the size."
    else "Each study resamples the fitted data's rows.")
  lines <- character()
  pn <- if (length(unique(d$n)) >= 2L)
    tryCatch(ilm_power_n(object, target), error = function(e) NULL) else NULL
  for (e in unique(d$effect)) {
    z <- d[d$effect == e, , drop = FALSE]; z <- z[order(z[[scol]]), , drop = FALSE]
    s <- sprintf("%s, the chance of detecting it is %s.",
                 if (joint) sprintf("At %s times the assumed effect", ilm_fmt(e, digits))
                 else sprintf("At an effect of %s on the link scale", ilm_fmt(e, 3)),
                 paste(sprintf("%s with %s %s (%s to %s)", pct(z$power),
                               format(z[[scol]]), who, pct(z$mc_lower),
                               pct(z$mc_upper)), collapse = "; "))
    r <- if (is.null(pn)) NULL else pn[pn$effect == e, , drop = FALSE]
    if (!is.null(r) && nrow(r)) {
      v <- if (!is.null(r$n_unit)) c(r$n_unit, r$n_unit_lower, r$n_unit_upper)
           else c(r$n, r$n_lower, r$n_upper)
      s <- paste(s, if (!is.finite(v[1]))
        sprintf("The sizes tried do not bracket %s power, so the size that reaches it lies outside them.",
                pct(target))
      else sprintf("Reaching %s power takes about %s %s%s.", pct(target),
                   format(round(v[1])), who,
                   if (all(is.finite(v[2:3])))
                     sprintf(" -- somewhere between %s and %s, given the Monte Carlo error",
                             format(round(min(v[2:3]))), format(round(max(v[2:3]))))
                   else if (is.finite(v[2]))
                     sprintf(" -- at least %s given the Monte Carlo error, and the upper end of that range lies beyond the sizes tried",
                             format(round(v[2])))
                   else if (is.finite(v[3]))
                     sprintf(" -- at most %s given the Monte Carlo error, and the lower end lies below the sizes tried",
                             format(round(v[3])))
                   else ""))
    }
    lines <- c(lines, s)
  }
  sec$effects <- lines
  dl <- character()
  if (any(d$converged < 1, na.rm = TRUE))
    dl <- c(dl, sprintf(
      "Up to %s of the simulated studies failed to fit, and they count as missed detections: a study that will not fit has detected nothing. Dividing by the ones that did gives the more flattering power_converged column.",
      pct(max(1 - d$converged, na.rm = TRUE))))
  dft <- attr(object, "design_fails")
  if (!is.null(dft) && nrow(dft))
    for (i in seq_len(nrow(dft)))
      dl <- c(dl, sprintf(
        "FAIL -- %s: %s, in %s of the simulated studies. Every analysis of a study like this will report it.%s",
        dft$check[i], dft$detail[i], pct(dft$share[i]),
        if (nzchar(dft$suggestion[i])) paste0(" What to do: ", dft$suggestion[i], ".") else ""))
  sec$diagnostics <- dl
  cav <- c(
    "A simulated power is itself an estimate. The intervals above are its Monte Carlo error, which more replicates narrow and nothing else does.",
    if (joint)
      "The effect is a multiple of the assumed coefficients, and the power at 1 is only as good as the assumptions behind them; smaller multiples show how fast it falls if the effect is smaller than hoped."
    else
      "The power holds for the effect assumed. A smaller true effect needs a larger study; smaller values of `effect` show how fast the power falls.")
  if (isTRUE(attr(object, "grouped")) && !test %in% c("t", "F"))
    cav <- c(cav, paste(
      "The test counted uses a large-sample reference, which with few groups",
      "can run somewhat liberal for an effect that varies between them -- and",
      "flatter the power a little by the same amount."))
  sec$caveats <- cav
  structure(list(sections = sec, family = attr(object, "family"),
                 causal = FALSE, object_class = "ilm_power"),
            class = "ilm_interpretation")
}

#' @rdname ilm_interpret
#' @export
ilm_interpret.ilm_contrast <- function(object, causal = NULL, ame = FALSE,
                                       digits = 3, ...) {
  d <- as.data.frame(object); class(d) <- "data.frame"
  adj <- d$adjust[1L]; lev <- attr(object, "level"); if (is.null(lev)) lev <- 0.95
  fam <- attr(object, "family"); if (is.null(fam)) fam <- "gaussian"
  catg <- identical(fam, "multinomial") || startsWith(fam, "ordinal")
  pp <- catg && identical(attr(object, "type"), "response")
  sc <- if (pp) 100 else 1
  unit <- if (pp) " percentage points of probability"
          else if (startsWith(fam, "ordinal")) " on the latent scale"
          else if (identical(attr(object, "type"), "link") && fam != "gaussian")
            " on the link scale" else ""
  sec <- list()
  sec$model <- sprintf("%d comparison%s of marginal means%s, %s.", nrow(d),
    if (nrow(d) == 1L) "" else "s",
    if (pp) ", as differences in the probability of each category" else "",
    switch(adj,
      max_t = sprintf("with intervals and p-values adjusted so that all %d hold together at %s (the single-step studentized maximum)",
                      nrow(d), paste0(100 * lev, "%")),
      bonferroni = "adjusted by Bonferroni, which is conservative",
      "unadjusted, which suits comparisons chosen before the data were seen"))
  sec$effects <- sprintf("%s: a difference of %s%s (%s%% interval %s to %s, %sp = %s) -- %s.",
    d$contrast, ilm_fmt(sc * d$estimate, digits), unit, 100 * lev,
    ilm_fmt(sc * d$lower, digits), ilm_fmt(sc * d$upper, digits),
    if (adj == "none") "" else "adjusted ",
    format.pval(d$p_adj, digits = 2, eps = 1e-4),
    ## the graded phrase is written for the middle of a sentence
    sub(",$", "", vapply(d$p_adj, ilm_evidence, "")))
  cav <- character()
  if (identical(attr(object, "weights"), "cells"))
    cav <- c(cav, paste(
      "The means were averaged with cell weights, so each group carries its",
      "own mix of the other variables and a difference includes whatever that",
      "difference in composition contributes -- it is not an adjusted",
      "comparison."))
  if (catg && pp)
    cav <- c(cav, paste(
      "Within a group the probabilities sum to one, so a rise in one category",
      "is a fall somewhere else: the differences for a pair of groups sum to",
      "zero across the categories."))
  cav <- c(cav, paste(
    "An interval that includes zero means the data are consistent with no",
    "difference -- not that there is none."))
  sec$caveats <- cav
  structure(list(sections = sec, family = fam, causal = FALSE,
                 object_class = "ilm_contrast"), class = "ilm_interpretation")
}

#' @export
print.ilm_interpretation <- function(x, width = 76L, ...) {
  s <- x$sections
  head <- switch(x$object_class,
    ilm_dag_model = "INTERPRETATION (DAG-identified effect)",
    ilm_did = "INTERPRETATION (difference in differences)",
    ilm_rdd = "INTERPRETATION (regression discontinuity)",
    ilm_power = "INTERPRETATION (power analysis)",
    ilm_contrast = "INTERPRETATION (comparisons)",
    "INTERPRETATION")
  cat(head, "\n", strrep("=", nchar(head)), "\n\n", sep = "")
  blk <- function(title, txt) {
    if (is.null(txt) || !length(txt)) return(invisible())
    cat(title, "\n", sep = "")
    for (t in txt) cat(ilm_wrap(t, width, "  "), "\n\n", sep = "")
  }
  blk("The model", s$model)
  blk("What it says", s$effects)
  blk("What the checks found", s$diagnostics)
  blk("How far to trust it", s$caveats)
  invisible(x)
}
