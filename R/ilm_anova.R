## illume: analysis of deviance for fixed effects.
##
## WHY THIS EXISTS: car::Anova.default reads attr(model.matrix, "assign"), which
## has one entry per COLUMN of X.  A multinomial fit has p*C coefficients -- C
## per column -- so car matches a length-p assign against a length-p*C
## coefficient vector and silently tests only the first category, reporting
## df = 1 where the correct joint test has df = C.  It produces a plausible
## table that answers the wrong question, so ilm_model registers its own method.
##
## Each term spanning d columns of X is tested jointly across all C category
## dimensions: df = d * C.

#' Coefficient positions for one model term
#'
#' A term occupying `d` columns of the design matrix has `d * C` coefficients,
#' because each column has one coefficient per category dimension. This returns
#' all of them, which is what a joint test of that term requires.
#'
#' @param object A fitted `"ilm_model"` object.
#' @param t Integer term index.
#' @return Integer vector of coefficient positions.
#' @keywords internal
#' @noRd
ilm_term_idx <- function(object, t) {
  p <- ncol(object$X); C <- object$C
  cols <- which(!is.na(object$assign) & object$assign == t)
  as.vector(outer(cols, seq_len(C) - 1L, function(r, c0) c0 * p + r))
}

#' Higher-order terms containing a given term
#'
#' For the principle of marginality: `x1` is contained in `x1:grp`, so a Type II
#' test of `x1` is made after removing the interaction.
#'
#' @param mt A `terms` object.
#' @param j Integer term index.
#' @return Integer vector of term indices that contain term `j`.
#' @keywords internal
#' @noRd
ilm_relatives_of <- function(mt, j) {
  fo <- attr(mt, "factors")
  if (is.null(fo) || !ncol(fo)) return(integer(0))
  vj <- rownames(fo)[fo[, j] > 0]
  which(vapply(seq_len(ncol(fo)), function(k)
    k != j && length(vj) && all(vj %in% rownames(fo)[fo[, k] > 0]), TRUE))
}

#' Refit with selected terms removed
#'
#' Drops the relevant **columns** of the design matrix rather than rebuilding the
#' formula. This is exact and avoids re-evaluating transformations against the
#' model frame. Smooth null-space columns are always kept.
#'
#' @param object A fitted `"ilm_model"` object.
#' @param drop_terms Integer vector of term indices to remove.
#' @param restarts Integer. Optimiser restarts.
#' @return A refitted `"ilm_model"` object.
#' @keywords internal
#' @noRd
ilm_refit_drop <- function(object, drop_terms, restarts = 2L) {
  keep <- is.na(object$assign) | !(object$assign %in% drop_terms)
  ## Everything that makes this the same model -- family, weights, censoring,
  ## the dispersion model, the baseline spline -- travels through
  ## ilm_refit_like(). A reduced fit that quietly drops one of them is not a
  ## reduced model, and the test built on it is meaningless: omitting the
  ## censoring put the null rejection rate at 100%.
  f <- ilm_refit_like(object, X = object$X[, keep, drop = FALSE],
                      keep = keep, restarts = restarts)
  f$assign <- object$assign[keep]
  f$term_labels <- object$term_labels
  f
}

#' Compact summary of a reduced-model refit
#'
#' Returns only the quantities the tests need, so that fitted objects never have
#' to be sent back from parallel workers.
#'
#' @param object A fitted `"ilm_model"` object.
#' @param drop_terms Integer vector of term indices to remove.
#' @param restarts Integer. Optimiser restarts.
#' @return A list with the log-likelihood, coefficients, covariance and term map,
#'   or `NULL` if the refit failed.
#' @keywords internal
#' @noRd
ilm_drop_summary <- function(object, drop_terms, restarts = 2L, ddf = NULL) {
  f <- try(ilm_refit_drop(object, drop_terms, restarts), silent = TRUE)
  if (inherits(f, "try-error")) return(NULL)
  out <- list(ll = as.numeric(logLik(f)), conv = f$opt$convergence, pd = isTRUE(f$sdr$pdHess),
       b = setNames(f$opt$par, f$pnames)[seq_len(ncol(f$X) * f$C)],
       V = tryCatch(f$sdr$cov.fixed[seq_len(ncol(f$X) * f$C),
                                    seq_len(ncol(f$X) * f$C), drop = FALSE],
                    error = function(e) NULL),
       assign = f$assign, p = ncol(f$X), C = f$C)
  ## An F test of a term in this reduced model needs its denominator df, and
  ## the fit that gives them stays in the worker, so they are computed here,
  ## one per remaining term: the Wald F on them (and Kenward-Roger's scaled F
  ## on its adjusted covariance), keyed by term index.
  if (!is.null(ddf)) {

    tt <- sort(unique(stats::na.omit(f$assign)))
    parts <- if (identical(ddf, "kenward-roger")) ilm_kr_parts(f) else NULL
    gr <- if (!identical(ddf, "kenward-roger")) ilm_vbeta_grad(f, NULL) else NULL
    out$ftab <- lapply(stats::setNames(tt, tt), function(t) {
      idx <- which(!is.na(f$assign) & f$assign == t)
      L <- diag(length(out$b))[idx, , drop = FALSE]
      ilm_anova_f(f, out$b, L, ddf, parts = parts, gr = gr)
    })
  }
  out
}

## The Wald F test of the rows of L -- a term's block -- on finite
## denominator df: F = W / q with W = (Lb)' (L V L')^-1 (Lb), and for
## Kenward-Roger W on their adjusted covariance, times their scale.
#' @keywords internal
#' @noRd
ilm_anova_f <- function(fit, b, L, df, parts = NULL, gr = NULL) {
  d <- ilm_table_ddf(fit, L, df, parts = parts, gr = gr)
  V <- if (!is.null(d$V)) d$V else suppressWarnings(vcov(fit))
  Lb <- L %*% b; LVL <- L %*% V %*% t(L)
  W <- tryCatch(as.numeric(t(Lb) %*% solve(LVL, Lb)), error = function(e) NA_real_)
  q <- nrow(L)
  list(F = d$scale * W / q, num = q, den = d$df, method = d$method,
       reason = d$reason)
}

#' Run several reduced-model refits, optionally in parallel
#' @param object A fitted `"ilm_model"` object.
#' @param drop_sets List of integer vectors of term indices.
#' @param ncores Integer. Worker processes.
#' @param restarts Integer. Optimiser restarts.
#' @return A list of [ilm_drop_summary()] results.
#' @keywords internal
#' @noRd
ilm_refit_drops <- function(object, drop_sets, ncores = 1L, restarts = 2L,
                            ddf = NULL) {
  ## Everything that makes the reduced fit the SAME model, built by the one
  ## function that knows the list (see ilm_refit_stub()). A field missing here
  ## is a field the refit silently does without: leaving the censoring out
  ## once put the likelihood-ratio test's null rejection rate at 100%, and
  ## leaving the zero part out made every such test of a zero-inflated model
  ## compare it against a model without one.
  stub <- ilm_refit_stub(object)
  cl <- ilm_pool(ncores)
  on.exit(if (!is.null(cl)) try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
  one <- function(i) ilm_drop_summary(stub, drop_sets[[i]], restarts, ddf)
  if (!is.null(cl))
    parallel::clusterExport(cl, "stub", envir = environment())
  ilm_lapply(cl, seq_along(drop_sets), one)
}

#' Joint Wald statistic for a block of coefficients
#' @param b Numeric coefficient vector.
#' @param V Covariance matrix.
#' @param idx Integer positions of the block.
#' @return A chi-square statistic.
#' @keywords internal
#' @noRd
ilm_wald_block <- function(b, V, idx) {
  bb <- b[idx]; VV <- V[idx, idx, drop = FALSE]
  tryCatch(as.numeric(t(bb) %*% solve(VV, bb)), error = function(e) NA_real_)
}

## Type III main effects depend on how the factors are coded, and R's default
## coding is not the one the label implies. Measured on y ~ g * x + z with
## 600 rows, the x row came back Chisq 167 under contr.treatment and 1086 under
## contr.sum -- same model, same data, same term. Under treatment coding the x
## row tests the slope at the REFERENCE level of g; only under a coding whose
## columns sum to zero does it test the average slope, which is what a Type III
## main effect is supposed to mean. car::Anova() warns for the same reason.
#' @keywords internal
#' @noRd
ilm_orth_contr <- function(z) {
  if (is.character(z))
    return(z %in% c("contr.sum", "contr.helmert", "contr.poly"))
  if (is.matrix(z)) return(all(abs(colSums(z)) < 1e-8))
  FALSE
}

## Which factors in an interaction are coded in a way that stops a Type III
## main effect from meaning what it says. Separated from the warning because
## the useful response is to fix it, not only to say it.
#' @keywords internal
#' @noRd
ilm_type3_bad_coding <- function(labs, contrasts, mf = NULL) {
  inter <- grep(":", labs, fixed = TRUE, value = TRUE)
  if (!length(inter)) return(character(0))
  involved <- unique(unlist(strsplit(inter, ":", fixed = TRUE)))
  bad <- if (length(contrasts))
    names(contrasts)[names(contrasts) %in% involved &
                     !vapply(contrasts, ilm_orth_contr, TRUE)] else character(0)
  ## A numeric predictor in an interaction has exactly the same problem and no
  ## contrast to blame it on: the main effect of x in x:z is the slope where
  ## z = 0, and only when z is centred is that the average slope. The existing
  ## warning has always said so; this is what finds it.
  if (!is.null(mf)) {
    num <- involved[vapply(involved, function(v) {
      z <- mf[[v]]
      !is.null(z) && is.numeric(z) && !is.matrix(z) &&
        abs(mean(z, na.rm = TRUE)) > 1e-8 * max(1, stats::sd(z, na.rm = TRUE))
    }, TRUE)]
    bad <- unique(c(bad, num))
  }
  bad
}

#' @keywords internal
#' @noRd
ilm_warn_type3_coding <- function(labs, contrasts, mf = NULL) {
  bad <- ilm_type3_bad_coding(labs, contrasts, mf)
  if (!length(bad)) return(invisible(FALSE))
  warning("this model has interactions and ", paste(sQuote(bad), collapse = ", "),
          if (length(bad) > 1L) " use " else " uses ",
          "a coding whose columns do not sum to zero, so the main-effect rows ",
          "are not Type III tests: each one is measured at the reference level ",
          "of the factor it interacts with, not averaged over it. Refit with ",
          "contrasts = list(", bad[1], " = \"contr.sum\"), or use type = 2, ",
          "which does not depend on the coding. Numeric predictors in an ",
          "interaction have the same issue unless they are centred.",
          call. = FALSE)
  invisible(TRUE)
}

## Refit the same model with sum-to-zero coding on the named factors, so that
## a Type III main effect is averaged over the levels it interacts with rather
## than measured at one of them. afex does this and says so; saying nothing
## and reporting the other test would be worse, and so would refusing to
## report anything.
##
## Only the coding changes. The likelihood, the fitted values and every test
## that does not involve a main effect inside an interaction are identical.
#' @keywords internal
#' @noRd
ilm_recode_sum <- function(object, which) {
  ## A smooth's null-space columns and a Royston-Parmar baseline are appended
  ## to the design outside model.matrix(), so rebuilding from the terms alone
  ## would drop them. Those models keep the warning instead of a silent
  ## half-refit.
  if (length(object$smooths) || !is.null(object$rp)) return(NULL)
  ctr <- object$contrasts
  mf <- object$model
  for (v in which) {
    if (!is.null(ctr[[v]])) ctr[[v]] <- "contr.sum"
    else if (!is.null(mf[[v]]) && is.numeric(mf[[v]]) && !is.matrix(mf[[v]]))
      ## centring a numeric does for it what sum coding does for a factor:
      ## it puts zero at the average, so the other term's main effect is
      ## measured there rather than at an arbitrary origin
      mf[[v]] <- mf[[v]] - mean(mf[[v]], na.rm = TRUE)
  }
  Xn <- try(ilm_fit_cols(
    stats::model.matrix(object$terms, mf, contrasts.arg = ctr),
    object), silent = TRUE)
  if (inherits(Xn, "try-error") || ncol(Xn) != ncol(object$X)) return(NULL)
  f <- try(suppressWarnings(ilm_refit_like(object, X = Xn, restarts = 1L)),
           silent = TRUE)
  if (inherits(f, "try-error") || !isTRUE(f$ok)) return(NULL)
  ## ilm_fit() knows nothing about formulas, so the pieces every downstream
  ## accessor reads have to be carried across by hand
  f$assign <- attr(Xn, "assign")
  f$contrasts <- attr(Xn, "contrasts")
  for (nm in c("call", "formula", "fixed_formula", "terms", "xlev", "model",
               "smooths", "bars", "na.action", "n_dropped", "rows", "ylevels",
               "term_labels"))
    f[[nm]] <- object[[nm]]
  f
}

#' Analysis of deviance for fixed effects
#'
#' Tests each fixed-effect term jointly: all of its coefficients at once, and
#' for a multinomial outcome **across all category dimensions**.
#'
#' @section Why this exists rather than car::Anova:
#' `car::Anova()` reads the `assign` attribute of the model matrix, which has one
#' entry per design column. A multinomial fit has `C` coefficients per column, so
#' `car` would match a short vector against a long one and silently test only the
#' first category, reporting `df = 1` where the correct joint test has `df = C`.
#' It returns a plausible table that answers the wrong question. `illume`
#' registers its own method so that `car::Anova()` dispatches here and gets the
#' right answer.
#'
#' @section Type II versus Type III:
#' Type III tests each term with every other term in the model. Type II tests
#' each term after all terms that do **not** contain it, respecting marginality.
#' They agree when there are no interactions.
#'
#' With interactions present, Type II is usually preferred for main effects, and
#' it has the practical advantage of not depending on how factors are coded. Type
#' III does depend on the coding and is only well defined with orthogonal
#' contrasts such as [stats::contr.sum()].
#'
#' @section Wald versus likelihood ratio:
#' Wald tests are fast and need no refitting for Type III. Likelihood-ratio tests
#' refit a reduced model for each term, which is slower but avoids the
#' Hauck-Donner effect, where a Wald statistic can shrink for very strong
#' effects. If a term is important to your conclusions, confirm it with the LRT;
#' for small samples use [ilm_pb_lrt()] instead, since both rely on
#' large-sample approximations.
#'
#' @param object A fitted `"ilm_model"` object, fitted through the formula interface.
#' @param type `2` or `"II"` (the default), or `3` or `"III"`. Type II tests
#'   each term against everything not containing it, which does not depend on
#'   how the factors are coded. Type III tests each term against every other,
#'   which does -- see `recode`.
#' @param recode When `type = 3` and a term inside an interaction is coded in a
#'   way that makes its main-effect row something other than a Type III test,
#'   refit with [stats::contr.sum()] for those factors and say so. `FALSE`
#'   tests the model exactly as coded and warns instead. The fit passed in is
#'   never modified either way.
#' @param test `"Wald"` or `"LRT"`.
#' @param ncores Integer. Worker processes for the refits.
#' @param restarts Integer. Optimiser restarts in refits.
#' @param df For a gaussian mixed model's Wald tests, the denominator degrees
#'   of freedom of each F: `"auto"` (Satterthwaite's, as `lmerTest` reports),
#'   `"satterthwaite"`, `"kenward-roger"` (their df and scaled F, on their
#'   adjusted covariance; REML fits) or `"asymptotic"`. Other fits are
#'   unaffected.
#' @param statistic `"F"` (the default) or `"Chisq"`. For a gaussian mixed
#'   model's Wald tests, `"F"` gives an F test on the denominator df above and
#'   `"Chisq"` the large-sample chi-square the table gave before; with few
#'   groups the chi-square is too liberal. Every other fit keeps its own
#'   reference: exact F with nothing integrated out, chi-square otherwise.
#' @return An `"anova"` data frame with one row per fixed-effect term. The
#'   columns are `Df`, `Chisq` and `Pr(>Chisq)` in general; `Df`, `F value`
#'   and `Pr(>F)` when the model admits exact inference -- a gaussian model with
#'   no random or smooth terms, where the residual variance is estimated rather
#'   than assumed known (see [ilm_model()]); and `NumDF`, `DenDF`, `F value` and
#'   `Pr(>F)` for a gaussian mixed model's Wald tests with `statistic = "F"`. A model with no terms to test, such
#'   as an intercept-only model, returns a table with zero rows rather than an
#'   error: "there is nothing to test" is an answer, not a failure.
#' @references
#' Fox, J., & Weisberg, S. (2019). *An R Companion to Applied Regression*,
#' 3rd ed. Sage. (Chapter 5 explains Type II and Type III tests.)
#'
#' Hauck, W. W., & Donner, A. (1977). Wald's test as applied to hypotheses in
#' logit analysis. *Journal of the American Statistical Association*, 72(360),
#' 851--853.
#' @seealso [ilm_pb_lrt()], [ilm_coef_table()].
#' @export
ilm_anova <- function(object, type = 2, test = c("Wald", "LRT"),
                        ncores = 1L, restarts = 2L, recode = TRUE,
                        df = "auto", statistic = c("F", "Chisq")) {
  test <- match.arg(test); statistic <- match.arg(statistic)
  if (test == "LRT") object <- ilm_as_ml(object, "ilm_anova(test = \"LRT\")")
  type <- toupper(as.character(type)[1])
  if (!type %in% c("3", "III", "2", "II")) stop("type must be 2 / \"II\" or 3 / \"III\"")
  type3 <- type %in% c("3", "III")
  if (is.null(object$assign))
    stop("no term map: ilm_anova() needs a model fitted through the formula interface")
  ## a term held at its boundary leaves the fixed-effect tests usable
  if (test == "Wald" && !ilm_fixed_usable(object))
    warning("Hessian is not positive definite; these Wald tests are not usable. ",
            "See object$checks.", call. = FALSE)

  mt <- object$terms; labs <- object$term_labels; nt <- length(labs)
  if (type3) {
    bad <- ilm_type3_bad_coding(labs, object$contrasts, object$model)
    alt <- if (length(bad) && isTRUE(recode)) ilm_recode_sum(object, bad) else NULL
    if (!is.null(alt)) {
      fac <- bad[bad %in% names(object$contrasts)]
      num <- setdiff(bad, fac)
      message("ilm_anova: refitted for these tests",
              if (length(fac)) paste0(", with contrasts = list(",
                paste(sprintf("%s = \"contr.sum\"", fac), collapse = ", "),
                ")") else "",
              if (length(num)) paste0(
                if (length(fac)) " and with " else ", with ",
                paste(num, collapse = ", "), " centred") else "",
              ". A Type III main effect of a term that is also in an ",
              "interaction is the average over the other term only when ",
              "that term has zero at its own average -- sum-to-zero ",
              "coding for a factor, a centred value for a numeric. ",
              "Otherwise it is the effect at that factor's reference ",
              "level, or where the numeric is zero. Nothing else about ",
              "the model changes, and ",
              "the fit you passed in is untouched. Use type = 2 to avoid the ",
              "question, or recode = FALSE to test the model as coded.")
      object <- alt
      mt <- object$terms; labs <- object$term_labels
    } else if (length(bad))
      ilm_warn_type3_coding(labs, object$contrasts, object$model)
  }
  rel <- lapply(seq_len(nt), function(j) ilm_relatives_of(mt, j))
  ## a term whose every column was dropped as aliased (ilm_model(aliased =
  ## "drop")) has nothing to test: no row, and the heading says so
  gone <- which(vapply(seq_len(nt), function(j)
    !any(!is.na(object$assign) & object$assign == j), TRUE))
  if (type3) rel <- lapply(rel, function(z) integer(0))   # III conditions on all

  ## ---- assemble the refits needed, de-duplicated ---------------------------
  need <- list()
  key <- function(s) paste(sort(unique(s)), collapse = ",")
  add_need <- function(s) { k <- key(s); if (!k %in% names(need)) need[[k]] <<- s; k }
  keyM1 <- keyM0 <- character(nt)
  for (j in setdiff(seq_len(nt), gone)) {
    if (test == "LRT") {
      keyM1[j] <- if (length(rel[[j]])) add_need(rel[[j]]) else ""
      keyM0[j] <- add_need(c(rel[[j]], j))
    } else if (length(rel[[j]])) {
      keyM1[j] <- add_need(rel[[j]])
    }
  }
  ## a gaussian mixed model's Wald tests as F on finite denominator df (D-DF3);
  ## the method is resolved once, so the refits and the full fit agree
  fin <- test == "Wald" && identical(statistic, "F") && ilm_finite_df(object, df)
  dmeth <- if (!fin) NULL else if (identical(df, "auto")) "satterthwaite"
           else match.arg(as.character(df), c("satterthwaite", "kenward-roger"))
  if (fin && identical(dmeth, "kenward-roger")) {
    ok <- ilm_kr_applicable(object)
    if (!isTRUE(ok))
      stop("Kenward-Roger is not available for this model: ", ok,
           ". `df = \"satterthwaite\"` does apply here.", call. = FALSE)
  }
  fits <- if (length(need)) ilm_refit_drops(object, need, ncores, restarts,
                                            ddf = dmeth) else list()
  names(fits) <- names(need)
  kr_full <- if (fin && identical(dmeth, "kenward-roger")) ilm_kr_parts(object) else NULL
  gr_full <- if (fin && !identical(dmeth, "kenward-roger")) ilm_vbeta_grad(object, NULL)
             else NULL
  fell <- character(0)

  ## logLik(), which takes any boundary penalty back out of the objective
  ll_full <- as.numeric(logLik(object))
  b_full <- coef(object); V_full <- suppressWarnings(vcov(object))

  rows <- lapply(seq_len(nt), function(j) {
    if (j %in% gone) return(NULL)
    ncol_j <- sum(!is.na(object$assign) & object$assign == j)
    nd <- ncol_j * object$C
    ## a gaussian mixed model: F on the term's own denominator df, from the
    ## full fit, or from the reduced fit the term is tested in
    if (fin) {
      ft <- if (!nzchar(keyM1[j])) {
        L <- diag(length(b_full))[ilm_term_idx(object, j), , drop = FALSE]
        ilm_anova_f(object, b_full, L, dmeth, parts = kr_full, gr = gr_full)
      } else {
        f1 <- fits[[keyM1[j]]]
        if (is.null(f1) || is.null(f1$ftab)) NULL else f1$ftab[[as.character(j)]]
      }
      if (is.null(ft) || !is.finite(ft$F))
        return(data.frame(NumDF = nd, DenDF = NA_real_, `F value` = NA_real_,
                          `Pr(>F)` = NA_real_, row.names = labs[j], check.names = FALSE))
      if (!is.finite(ft$den)) fell <<- c(fell, labs[j])
      return(data.frame(NumDF = ft$num, DenDF = ft$den, `F value` = ft$F,
                        `Pr(>F)` = stats::pf(ft$F, ft$num, ft$den, lower.tail = FALSE),
                        row.names = labs[j], check.names = FALSE))
    }
    df <- nd
    if (test == "Wald") {
      if (!nzchar(keyM1[j])) {
        stat <- ilm_wald_block(b_full, V_full, ilm_term_idx(object, j))
      } else {
        f1 <- fits[[keyM1[j]]]
        if (is.null(f1) || is.null(f1$V)) return(data.frame(Df = df, Chisq = NA_real_,
          `Pr(>Chisq)` = NA_real_, row.names = labs[j], check.names = FALSE))
        stub <- list(X = matrix(0, 1, f1$p), assign = f1$assign, C = f1$C)
        stat <- ilm_wald_block(f1$b, f1$V, ilm_term_idx(stub, j))
      }
    } else {
      ll1 <- if (!nzchar(keyM1[j])) ll_full else {
        z <- fits[[keyM1[j]]]; if (is.null(z)) NA_real_ else z$ll }
      z0 <- fits[[keyM0[j]]]
      ll0 <- if (is.null(z0)) NA_real_ else z0$ll
      stat <- 2 * (ll1 - ll0)
      if (!is.na(stat) && stat < 0) stat <- NA_real_   # a refit failed to converge
    }
    ## With nothing integrated out the exact reference is F rather than
    ## chi-square.  The chi-square test treats the residual variance as known,
    ## which it is not, and is anti-conservative in small samples.  Dividing the
    ## Wald statistic by its degrees of freedom gives the usual F.
    if (isTRUE(object$exact_df)) {
      fv <- stat / df
      data.frame(Df = df, `F value` = fv,
                 `Pr(>F)` = stats::pf(fv, df, object$resid_df, lower.tail = FALSE),
                 row.names = labs[j], check.names = FALSE)
    } else
      data.frame(Df = df, Chisq = stat,
                 `Pr(>Chisq)` = stats::pchisq(stat, df, lower.tail = FALSE),
                 row.names = labs[j], check.names = FALSE)
  })
  ## An intercept-only model has no terms to test.  Return an empty table
  ## rather than failing: "there is nothing to test" is a valid answer.
  rows <- rows[!vapply(rows, is.null, TRUE)]
  out <- if (length(rows)) do.call(rbind, rows) else if (fin) {
    data.frame(NumDF = integer(0), DenDF = numeric(0), `F value` = numeric(0),
               `Pr(>F)` = numeric(0), check.names = FALSE)
  } else {
    e <- data.frame(Df = integer(0), stat = numeric(0), p = numeric(0),
                    check.names = FALSE)
    names(e) <- if (isTRUE(object$exact_df)) c("Df", "F value", "Pr(>F)")
                else c("Df", "Chisq", "Pr(>Chisq)")
    e
  }
  attr(out, "heading") <- c(
    sprintf("Analysis of Deviance Table (Type %s %s tests)",
            if (type3) "III" else "II",
            if (isTRUE(object$exact_df)) "F"
            else if (fin) "Wald F"
            else if (test == "Wald") "Wald chi-square"
            else "likelihood-ratio chi-square"),
    if (fin)
      sprintf("Denominator df: %s%s", if (identical(dmeth, "kenward-roger"))
                "Kenward-Roger's, with their scaled F on the adjusted covariance"
              else "Satterthwaite's (Fai and Cornelius for several columns), as lmerTest",
              if (length(fell)) sprintf("; %s tested against the chi-square, as the df could not be formed",
                                         paste(fell, collapse = ", ")) else "")
    else NULL,
    if (object$C > 1L)
      sprintf("Response: %s   (%d categories, %d contrast dimensions)",
              deparse(object$formula[[2]]), object$J, object$C)
    else sprintf("Response: %s   (family: %s)", deparse(object$formula[[2]]),
                 if (!is.null(object$family)) object$family$name else "gaussian"),
    if (object$C > 1L)
      "Each term is tested jointly across all category dimensions: Df = (columns) x C"
    else NULL,
    if (test == "LRT")
      "LRT: random structure held fixed across models, so the Laplace error largely cancels."
    else NULL,
    if (length(gone))
      sprintf("Not tested, every column dropped as aliased: %s",
              paste(labs[gone], collapse = ", "))
    else NULL)
  class(out) <- c("anova", "data.frame")
  out
}

#' car::Anova method
#'
#' Lets `car::Anova(fit, type = 3)` dispatch to [ilm_anova()], so it uses the
#' correct joint blocking across categories instead of `car`'s default, which
#' would test only one category.
#'
#' Registered with car's generic when car is loaded. It used to be exported as
#' an ordinary function, which R's method lookup does not reach, so
#' `car::Anova()` fell back to its own tests: on a multinomial fit with a
#' numeric predictor and a three-level factor, 1 and 2 degrees of freedom
#' where the joint tests have 2 and 4.
#'
#' @param mod A fitted `"ilm_model"` object.
#' @param type `"II"`, `"III"`, `2` or `3`.
#' @param test.statistic `"Chisq"` (car's default) or `"F"`: for a gaussian
#'   mixed model, the chi-square or the F test on finite denominator df; see
#'   `statistic` in [ilm_anova()].
#' @param ... Passed to [ilm_anova()].
#' @return An `"anova"` data frame.
#' @exportS3Method car::Anova
Anova.ilm_model <- function(mod, type = c("II", "III", 2, 3), test.statistic = "Chisq", ...) {
  type <- as.character(type)[1]
  ## car's own argument chooses the reference, as it does for lmer fits:
  ## "Chisq" (car's default) or "F"
  st <- if (identical(as.character(test.statistic)[1], "F")) "F" else "Chisq"
  ilm_anova(mod, type = if (type %in% c("3", "III")) "3" else "2", statistic = st, ...)
}
