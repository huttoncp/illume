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
  ## ilm_fit() defaults to the multinomial family, so an internal refit that
  ## does not pass the fitted family is silently a DIFFERENT model.  Weights
  ## must carry over for the same reason.
  f <- ilm_fit(object$X[, keep, drop = FALSE], object$y, object$J,
                ilm_re_list_of(object), object$re_struct, object$ar,
                ylevels = object$ylevels, weights = object$weights,
                family = object$family, verbose = FALSE, restarts = restarts)
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
ilm_drop_summary <- function(object, drop_terms, restarts = 2L) {
  f <- try(ilm_refit_drop(object, drop_terms, restarts), silent = TRUE)
  if (inherits(f, "try-error")) return(NULL)
  list(ll = -f$opt$objective, conv = f$opt$convergence, pd = isTRUE(f$sdr$pdHess),
       b = setNames(f$opt$par, f$pnames)[seq_len(ncol(f$X) * f$C)],
       V = tryCatch(f$sdr$cov.fixed[seq_len(ncol(f$X) * f$C),
                                    seq_len(ncol(f$X) * f$C), drop = FALSE],
                    error = function(e) NULL),
       assign = f$assign, p = ncol(f$X), C = f$C)
}

#' Run several reduced-model refits, optionally in parallel
#' @param object A fitted `"ilm_model"` object.
#' @param drop_sets List of integer vectors of term indices.
#' @param ncores Integer. Worker processes.
#' @param restarts Integer. Optimiser restarts.
#' @return A list of [ilm_drop_summary()] results.
#' @keywords internal
#' @noRd
ilm_refit_drops <- function(object, drop_sets, ncores = 1L, restarts = 2L) {
  X <- object$X; y <- object$y; J <- object$J; asg <- object$assign
  rl <- ilm_re_list_of(object); rs <- object$re_struct; arr <- object$ar
  yl <- object$ylevels; tl <- object$term_labels
  fm <- object$family; wt <- object$weights; Cc <- object$C
  ## ilm_drop_summary() calls ilm_re_list_of(), which needs $re; hand over the already
  ## normalised terms so workers do not re-derive them.
  rl_re <- object$re
  cl <- ilm_pool(ncores)
  on.exit(if (!is.null(cl)) try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
  one <- function(i) {
    ## C = J - 1L is wrong for every family but the multinomial, where J is
    ## NULL and the subtraction yields a zero-length vector.
    stub <- list(X = X, y = y, J = J, assign = asg, re = rl_re, re_struct = rs,
                 ar = arr, ylevels = yl, term_labels = tl, C = Cc,
                 family = fm, weights = wt)
    ilm_drop_summary(stub, drop_sets[[i]], restarts)
  }
  if (!is.null(cl)) parallel::clusterExport(cl, "rl_re", envir = environment())
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

#' Analysis of deviance for fixed effects
#'
#' Tests each fixed-effect term, **jointly across all category dimensions**.
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
#' @param type `3` or `"III"`, or `2` or `"II"`.
#' @param test `"Wald"` or `"LRT"`.
#' @param ncores Integer. Worker processes for the refits.
#' @param restarts Integer. Optimiser restarts in refits.
#' @return An `"anova"` data frame with one row per fixed-effect term. The
#'   columns are `Df`, `Chisq` and `Pr(>Chisq)` in general, or `Df`, `F value`
#'   and `Pr(>F)` when the model admits exact inference -- a gaussian model with
#'   no random or smooth terms, where the residual variance is estimated rather
#'   than assumed known (see [ilm_model()]). A model with no terms to test, such
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
ilm_anova <- function(object, type = 3, test = c("Wald", "LRT"),
                        ncores = 1L, restarts = 2L) {
  test <- match.arg(test)
  type <- toupper(as.character(type)[1])
  if (!type %in% c("3", "III", "2", "II")) stop("type must be 2 / \"II\" or 3 / \"III\"")
  type3 <- type %in% c("3", "III")
  if (is.null(object$assign))
    stop("no term map: ilm_anova() needs a model fitted through the formula interface")
  if (test == "Wald" && !isTRUE(object$sdr$pdHess))
    warning("Hessian is not positive definite; these Wald tests are not usable. ",
            "See object$checks.", call. = FALSE)

  mt <- object$terms; labs <- object$term_labels; nt <- length(labs)
  rel <- lapply(seq_len(nt), function(j) ilm_relatives_of(mt, j))
  if (type3) rel <- lapply(rel, function(z) integer(0))   # III conditions on all

  ## ---- assemble the refits needed, de-duplicated ---------------------------
  need <- list()
  key <- function(s) paste(sort(unique(s)), collapse = ",")
  add_need <- function(s) { k <- key(s); if (!k %in% names(need)) need[[k]] <<- s; k }
  keyM1 <- keyM0 <- character(nt)
  for (j in seq_len(nt)) {
    if (test == "LRT") {
      keyM1[j] <- if (length(rel[[j]])) add_need(rel[[j]]) else ""
      keyM0[j] <- add_need(c(rel[[j]], j))
    } else if (length(rel[[j]])) {
      keyM1[j] <- add_need(rel[[j]])
    }
  }
  fits <- if (length(need)) ilm_refit_drops(object, need, ncores, restarts) else list()
  names(fits) <- names(need)

  ll_full <- -object$opt$objective
  b_full <- coef(object); V_full <- suppressWarnings(vcov(object))

  rows <- lapply(seq_len(nt), function(j) {
    ncol_j <- sum(!is.na(object$assign) & object$assign == j)
    df <- ncol_j * object$C
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
  out <- if (length(rows)) do.call(rbind, rows) else {
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
            else if (test == "Wald") "Wald chi-square"
            else "likelihood-ratio chi-square"),
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
#' @param mod A fitted `"ilm_model"` object.
#' @param type `"II"`, `"III"`, `2` or `3`.
#' @param test.statistic Ignored; present for compatibility.
#' @param ... Passed to [ilm_anova()].
#' @return An `"anova"` data frame.
#' @export
Anova.ilm_model <- function(mod, type = c("II", "III", 2, 3), test.statistic = "Chisq", ...) {
  type <- as.character(type)[1]
  ilm_anova(mod, type = if (type %in% c("3", "III")) "3" else "2", ...)
}
