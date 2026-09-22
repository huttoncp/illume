## ---------------------------------------------------------------------------
## Factorial and repeated-measures ANOVA, specified the way it is taught.
##
## You name the participant column, the outcome, which factors vary BETWEEN
## participants and which vary WITHIN them, and any covariates. Nothing has to
## be written as a formula, and nothing has to be reshaped first.
##
## WHY THIS IS NOT A WRAPPER AROUND ilm_model(). A mixed model and a classical
## repeated-measures ANOVA answer the same question on a balanced complete
## design, but they report it differently and one of them reports things the
## other cannot. ilm_model() gives a Wald chi-square, because with random
## effects integrated out there is no exact residual degrees of freedom. A
## repeated-measures ANOVA gives an F on a denominator df, a mean squared
## error, a generalized eta squared, and -- when sphericity fails -- fractional
## degrees of freedom from a Greenhouse-Geisser or Huynh-Feldt correction.
## None of those fall out of the mixed fit.
##
## So the omnibus table is computed classically, from the multivariate
## formulation: the responses are arranged as one row per participant and one
## column per within-participant cell, and each effect is tested by projecting
## that matrix onto the contrast space for its within-participant part and
## running the between-participants model on the result. That is what
## car::Anova(idata=, idesign=) does, and what afex calls through to.
##
## ilm_model() is still fitted, for everything the classical route cannot do:
## estimated marginal means, contrasts, simple slopes, diagnostics, and an
## analysis that keeps participants with incomplete data instead of dropping
## them. `engine = "mixed"` makes it the primary rather than the companion.
##
## WHAT IT ADDS OVER afex. Every warning here names a remedy that exists in
## this package. Sphericity fails -> the unstructured mixed model, which does
## not assume sphericity and so needs no correction. Participants have missing
## cells -> the mixed model, which uses them. A covariate varies within
## participant -> the mixed model with the covariate split, because a single
## coefficient for a time-varying covariate estimates a blend of two different
## effects.
## ---------------------------------------------------------------------------

#' Check and reshape the data an ANOVA specification describes
#'
#' The specification is in terms of columns; the classical computation needs a
#' participants-by-cells matrix. This does the translation and refuses the
#' cases where the translation would quietly lose something.
#'
#' @param data A long data frame.
#' @param id,dv,between,within,covariate Column names.
#' @param fun_aggregate Function used when a participant has more than one row
#'   per cell.
#' @return A list with `Y`, `bdat`, `idata`, `dropped`, `aggregated`, `nsub`.
#' @keywords internal
#' @noRd
ilm_aov_marshal <- function(data, id, dv, between, within, covariate,
                            fun_aggregate = NULL) {
  need <- c(id, dv, between, within, covariate)
  miss <- setdiff(need, names(data))
  if (length(miss))
    stop("column(s) not in the data: ", paste(miss, collapse = ", "),
         call. = FALSE)
  d <- data[, unique(need), drop = FALSE]
  d[[id]] <- factor(d[[id]])
  for (v in c(between, within)) d[[v]] <- factor(d[[v]])
  if (!is.numeric(d[[dv]]))
    stop("`dv` (", dv, ") is a ", class(d[[dv]])[1], ". An ANOVA needs a ",
         "numeric outcome. For a categorical outcome use ilm_model() with the ",
         "binomial, multinomial or ordinal family.", call. = FALSE)
  for (v in covariate) if (!is.numeric(d[[v]]))
    stop("covariate `", v, "` is a ", class(d[[v]])[1], ". A covariate is ",
         "numeric; a categorical predictor that varies between participants ",
         "is a `between` factor.", call. = FALSE)

  ## rows with any missing value in the specification cannot be placed
  ok <- stats::complete.cases(d)
  n_na <- sum(!ok)
  d <- d[ok, , drop = FALSE]

  ## ---- a covariate must be constant within participant --------------------
  ## A time-varying covariate carries two effects at once: how participants
  ## who score higher ON AVERAGE differ, and what happens when a participant
  ## is higher THAN USUAL. They can point in opposite directions, and a single
  ## coefficient estimates a variance-weighted blend of them that answers
  ## neither question. The classical decomposition has no stratum for it at
  ## all -- it assumes one number per participant per cell.
  for (v in covariate) {
    nu <- tapply(d[[v]], d[[id]], function(z) length(unique(z)))
    if (any(nu > 1, na.rm = TRUE)) {
      k <- sum(nu > 1, na.rm = TRUE)
      stop("covariate `", v, "` varies within participant for ", k,
           " of ", length(nu), " participants, so it is a time-varying ",
           "covariate. A classical repeated-measures ANCOVA has no stratum ",
           "for one: its effect is part between-participant (people who score ",
           "higher on average) and part within-participant (a person scoring ",
           "higher than usual), and those can differ or even have opposite ",
           "signs, so one coefficient answers neither question.\n",
           "  Use engine = \"mixed\", where the covariate can be split into ",
           "its participant mean and the deviation from it, and the two ",
           "effects are estimated separately.", call. = FALSE)
    }
  }

  ## ---- one row per participant per within-cell ----------------------------
  cellkey <- if (length(within))
    interaction(d[within], drop = FALSE, sep = "\r") else
    factor(rep("", nrow(d)))
  tab <- table(d[[id]], cellkey)
  aggregated <- FALSE
  if (any(tab > 1)) {
    if (is.null(fun_aggregate)) fun_aggregate <- mean
    aggregated <- TRUE
    keep <- c(id, between, covariate)
    agg <- stats::aggregate(d[[dv]],
                            by = c(d[keep], list(.cell = cellkey)),
                            FUN = fun_aggregate)
    names(agg)[ncol(agg)] <- dv
    d <- agg
    cellkey <- factor(agg$.cell, levels = levels(cellkey))
    d$.cell <- NULL
  }

  ## ---- the participants-by-cells matrix -----------------------------------
  lv <- levels(cellkey)
  subs <- levels(droplevels(d[[id]]))
  Y <- matrix(NA_real_, length(subs), length(lv),
              dimnames = list(subs, lv))
  Y[cbind(match(as.character(d[[id]]), subs),
          match(as.character(cellkey), lv))] <- d[[dv]]

  ## Classical repeated measures needs every cell. A participant missing one
  ## occasion contributes nothing to any within-participant contrast, so the
  ## only classical option is to drop them entirely -- which is exactly the
  ## situation the mixed model exists for.
  full <- stats::complete.cases(Y)
  dropped <- sum(!full)
  Y <- Y[full, , drop = FALSE]
  subs <- subs[full]

  ## one row of between-participant information per participant
  bkeep <- unique(c(id, between, covariate))
  bdat <- d[match(subs, as.character(d[[id]])), bkeep, drop = FALSE]
  rownames(bdat) <- NULL
  for (v in between) bdat[[v]] <- droplevels(factor(bdat[[v]]))

  idata <- if (length(within)) {
    z <- as.data.frame(do.call(rbind, strsplit(lv, "\r", fixed = TRUE)),
                       stringsAsFactors = FALSE)
    names(z) <- within
    for (v in within) z[[v]] <- factor(z[[v]], levels = levels(data[[v]] <-
      factor(data[[v]])))
    z
  } else data.frame(.one = factor("a"))

  list(Y = Y, bdat = bdat, idata = idata, dropped = dropped, n_na = n_na,
       aggregated = aggregated, nsub = nrow(Y), cells = lv)
}

#' Within-participant contrast matrices, one per within-subject term
#'
#' The columns of the within-participant model matrix, grouped by the term they
#' belong to. `contr.sum` throughout, so that each term's columns span the
#' contrast space for that term alone and the projections are orthogonal to the
#' grand mean.
#'
#' @param idata One row per within-participant cell.
#' @param within Names of the within-participant factors.
#' @return A named list of contrast matrices, with `"(Intercept)"` first.
#' @keywords internal
#' @noRd
ilm_aov_within_P <- function(idata, within) {
  if (!length(within))
    return(list(`(Intercept)` = matrix(1, nrow(idata), 1)))
  ctr <- stats::setNames(rep(list("contr.sum"), length(within)), within)
  fo <- stats::as.formula(paste("~", paste(ilm_bq(within), collapse = " * ")))
  mm <- stats::model.matrix(fo, idata, contrasts.arg = ctr)
  asg <- attr(mm, "assign")
  labs <- c("(Intercept)", attr(stats::terms(fo), "term.labels"))
  out <- lapply(seq_along(labs) - 1L, function(j)
    mm[, asg == j, drop = FALSE])
  names(out) <- labs
  out
}

#' Sums of squares for one effect, in the multivariate formulation
#'
#' Projects the participants-by-cells matrix onto the within-participant
#' contrast space `P`, then tests the between-participant columns `idx` on the
#' result. With no within part `P` is the unit column and this reduces to an
#' ordinary between-participants ANOVA.
#'
#' @param Y Participants by cells.
#' @param X Between-participants model matrix.
#' @param idx Columns of `X` the effect occupies; the intercept for a pure
#'   within-participant effect.
#' @param P Within-participant contrast matrix.
#' @return A list with `SS`, `SSE`, `df_h`, `df_e`, the error SSCP `E`, and `P`.
#' @keywords internal
#' @noRd
ilm_aov_ss <- function(Y, X, idx, P) {
  YP <- Y %*% P                                   # n x q
  q <- ncol(YP)
  XtXi <- tryCatch(solve(crossprod(X)), error = function(e) NULL)
  if (is.null(XtXi))
    stop("the between-participants design is singular: two of its terms carry ",
         "the same information. Check for a factor that is constant, or two ",
         "that are duplicates of each other.", call. = FALSE)
  B <- XtXi %*% crossprod(X, YP)                  # p x q
  R <- YP - X %*% B
  E <- crossprod(R)                               # q x q error SSCP
  df_e_b <- nrow(X) - ncol(X)
  Bi <- B[idx, , drop = FALSE]
  Ai <- tryCatch(solve(XtXi[idx, idx, drop = FALSE]), error = function(e) NULL)
  if (is.null(Ai)) return(NULL)
  H <- t(Bi) %*% Ai %*% Bi                        # q x q hypothesis SSCP
  ## The univariate F is a ratio of TRACES, and a trace is only basis-free in
  ## an orthonormal basis. `contr.sum` columns are neither unit-length nor
  ## mutually orthogonal, so both traces have to be taken in the (P'P)^-1
  ## metric -- which is the same correction car applies inside its
  ## Greenhouse-Geisser calculation, for the same reason. Without it the
  ## contrasts are weighted unequally and every within-participant F is wrong,
  ## while every between-participant one stays right, because there the
  ## projection is a single column and any scaling cancels.
  Wm <- solve(crossprod(P))
  list(SS = sum(diag(H %*% Wm)), SSE = sum(diag(E %*% Wm)),
       df_h = q * length(idx), df_e = q * df_e_b,
       E = E, P = P, q = q, df_e_b = df_e_b)
}

#' Greenhouse-Geisser, Huynh-Feldt and Mauchly
#'
#' Arithmetic follows `car:::summary.Anova.mlm` exactly, so the numbers are
#' directly comparable with `car` and `afex`. Only meaningful when the effect
#' has more than one within-participant contrast; with one there is nothing for
#' sphericity to be violated about.
#'
#' @param E Error SSCP for the transformed responses.
#' @param P The within-participant contrast matrix.
#' @param df_e_b Between-participants error degrees of freedom.
#' @return A list with `gg`, `hf`, `W` and `p_mauchly`.
#' @keywords internal
#' @noRd
ilm_aov_sphericity <- function(E, P, df_e_b) {
  q <- nrow(E)
  if (q < 2) return(list(gg = NA_real_, hf = NA_real_, W = NA_real_,
                         p_mauchly = NA_real_))
  lam <- eigen(E %*% solve(crossprod(P)), only.values = TRUE)$values
  lam <- Re(lam); lam <- lam[lam > 0]
  gg <- ((sum(lam) / q)^2) / (sum(lam^2) / q)
  hf <- ((df_e_b + 1) * q * gg - 2) / (q * (df_e_b - q * gg))
  ## Mauchly, with the higher-order term car carries
  Psi <- crossprod(P)
  U <- solve(Psi, E)
  n <- df_e_b
  logW <- log(det(U)) - q * log(sum(diag(U)) / q)
  rho <- 1 - (2 * q^2 + q + 2) / (6 * q * n)
  w2 <- (q + 2) * (q - 1) * (q - 2) *
    (2 * q^3 + 6 * q^2 + 3 * nrow(P) + 2) / (288 * (n * q * rho)^2)
  z <- -n * rho * logW
  f <- q * (q + 1) / 2 - 1
  p1 <- stats::pchisq(z, f, lower.tail = FALSE)
  p2 <- stats::pchisq(z, f + 4, lower.tail = FALSE)
  list(gg = min(max(gg, 1 / q), 1), hf = min(hf, 1), W = exp(logW),
       p_mauchly = p1 + w2 * (p2 - p1))
}

#' Simple effects: what a two-way interaction is actually made of
#'
#' Comparing every cell of a 3-by-3 gives thirty-six differences, and most of
#' them move both factors at once, so they answer a question nobody asked. What
#' takes an interaction apart is the simple effect -- hold one factor fixed and
#' compare the levels of the other -- and it is worth having in both
#' directions, because "the arms differ at week 6" and "the drug arm changes
#' over time" are different findings and a study can have either without the
#' other.
#'
#' Each family is adjusted within itself. Holding `b` at one level and
#' comparing the levels of `a` is one family of comparisons; doing it at the
#' next level of `b` is another.
#'
#' @param fit The companion mixed model.
#' @param fp Exactly two factor names.
#' @param m The marshalled data, for the levels.
#' @return A list with one entry per direction, each a named list of
#'   [ilm_contrast()] results.
#' @keywords internal
#' @noRd
ilm_aov_simple <- function(fit, fp, m) {
  lev_of <- function(v) {
    z <- if (v %in% names(m$bdat)) m$bdat[[v]] else m$idata[[v]]
    levels(droplevels(factor(z)))
  }
  one_way <- function(a, b) {
    if (length(lev_of(a)) < 2L) return(NULL)
    out <- list()
    for (L in lev_of(b)) {
      r <- tryCatch({
        at <- stats::setNames(list(factor(L, levels = lev_of(b))), b)
        ilm_contrast(ilm_emmeans(fit, a, at = at))
      }, error = function(e) NULL)
      if (!is.null(r)) out[[paste0(b, " = ", L)]] <- r
    }
    if (length(out)) out else NULL
  }
  res <- list()
  d1 <- one_way(fp[1], fp[2])
  d2 <- one_way(fp[2], fp[1])
  if (!is.null(d1)) res[[paste0(fp[1], " within each ", fp[2])]] <- d1
  if (!is.null(d2)) res[[paste0(fp[2], " within each ", fp[1])]] <- d2
  if (!length(res)) return(NULL)
  structure(res, class = "ilm_simple_effects")
}

#' @export
print.ilm_simple_effects <- function(x, ...) {
  for (dir in names(x)) {
    cat("  ~ ", dir, " ~\n", sep = "")
    for (at in names(x[[dir]])) {
      cat("    [", at, "]\n", sep = "")
      print(x[[dir]][[at]])
      cat("\n")
    }
  }
  invisible(x)
}
