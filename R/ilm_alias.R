## ---- aliasing: fixed columns that are not all separable (item 4) -----------
##
## A model whose fixed-effect columns are linearly dependent has no unique
## fit. Some of it was caught: a numeric `by` variable beside its own smooth
## (ilm_model(), with its own message). Most was not: s(x, by = z) beside
## s(w, by = z), both of whose unpenalised parts contain z; poly(z, 2) beside
## s(x, by = z); two plain columns one a multiple of the other. Those fitted,
## came back with a failed Hessian and standard errors of NaN, and the checks
## pointed elsewhere. As Craig ruled, the whole fixed design is checked for
## rank before the fit, the smooths' unpenalised (null-space) columns with
## it, and a shortfall is named: which columns, which terms, and the fix.
##
## Aliasing that involves a smooth always stops: there is no sensible column
## of a smooth's null space to drop. For ordinary columns -- two collinear
## covariates, an interaction with an empty cell -- the fit stops by default,
## and ilm_model(aliased = "drop") drops the dependent columns as lm() and
## lme4 do, saying which (Craig's ruling, item 214).

## The column -> term map: a formula term's label, "(Intercept)", or the
## smooth whose null space the column is.
#' @keywords internal
#' @noRd
ilm_alias_terms <- function(X, asgn, term_labels) {
  cn <- colnames(X)
  vapply(seq_along(cn), function(j) {
    a <- asgn[j]
    if (is.na(a)) sub("\\.f[0-9]+$", "", cn[j])
    else if (a == 0L) "(Intercept)" else term_labels[a]
  }, "")
}

## The dependent columns, each with the columns it is a combination of. On
## standardised columns, so a column in dollars is not mistaken for a
## dependent one, with lm()'s tolerance. NULL when the design has full rank.
#' @keywords internal
#' @noRd
ilm_alias_find <- function(X, asgn, term_labels, tol = 1e-7) {
  if (is.null(X) || !ncol(X) || nrow(X) < 2L) return(NULL)
  Xs <- ilm_apply_scales(X, ilm_col_scales(X))
  q <- qr(Xs, tol = tol)
  if (q$rank >= ncol(X)) return(NULL)
  dep <- q$pivot[seq.int(q$rank + 1L, ncol(X))]
  ind <- q$pivot[seq_len(q$rank)]
  tm <- ilm_alias_terms(X, asgn, term_labels)
  smooth_cols <- is.na(asgn)
  qi <- qr(Xs[, ind, drop = FALSE], tol = tol)
  groups <- lapply(dep, function(d) {
    b <- tryCatch(qr.coef(qi, Xs[, d]), error = function(e) rep(NA_real_, length(ind)))
    with <- ind[!is.na(b) & abs(b) > 1e-6]
    cols <- c(d, with)
    ## a column's name as a reader knows it: a smooth's null-space column is
    ## "the unpenalised part of" that smooth, not its internal name
    say <- function(j) if (smooth_cols[j]) sprintf("the unpenalised part of %s", tm[j])
                       else if (identical(tm[j], "(Intercept)")) "the intercept"
                       else sprintf("`%s`", colnames(X)[j])
    list(column = colnames(X)[d], said = say(d), with_said = vapply(with, say, ""),
         zero = all(X[, d] == 0), terms = unique(tm[cols]),
         smooth = any(smooth_cols[cols]), all_smooth = all(smooth_cols[cols]))
  })
  list(groups = groups, smooth = any(vapply(groups, `[[`, TRUE, "smooth")),
       dependent = colnames(X)[dep])
}

## The empty cells of an interaction of two factors, the usual reason an
## interaction's columns are dependent. Character(0) when the term is not one.
#' @keywords internal
#' @noRd
ilm_alias_empty_cells <- function(mf, term) {
  vs <- ilm_unbq(strsplit(term, ":", fixed = TRUE)[[1]])
  if (length(vs) != 2L || !all(vs %in% names(mf))) return(character(0))
  a <- mf[[vs[1]]]; b <- mf[[vs[2]]]
  if (!(is.factor(a) || is.character(a)) || !(is.factor(b) || is.character(b)))
    return(character(0))
  tb <- table(factor(a), factor(b))
  z <- which(tb == 0, arr.ind = TRUE)
  if (!nrow(z)) return(character(0))
  paste0(rownames(tb)[z[, 1]], ":", colnames(tb)[z[, 2]])
}

## One clause per dependent column.
#' @keywords internal
#' @noRd
ilm_alias_lines <- function(al) {
  unique(vapply(al$groups, function(g) {
    if (g$zero) sprintf("%s is zero in every row", g$said)
    else sprintf("%s is a linear combination of %s", g$said,
                 if (length(g$with_said)) ilm_and(unique(utils::head(g$with_said, 6)))
                 else "the other columns")
  }, ""))
}

## What is said when the dependent columns are dropped, at fitting and in
## summary().
#' @keywords internal
#' @noRd
ilm_alias_drop_note <- function(al) {
  n <- length(al$dependent)
  paste0("the fixed-effect columns were not all separable, so ", n,
         if (n == 1L) " column was" else " columns were",
         " dropped as aliased, as aliased = \"drop\" asks: ",
         paste(ilm_alias_lines(al), collapse = "; "),
         ". The model is the same, written with fewer columns; ",
         if (n == 1L) "the dropped coefficient has" else "the dropped coefficients have",
         " no estimate.")
}

## The message. One sentence per dependent column, then the fix.
#' @keywords internal
#' @noRd
ilm_alias_message <- function(al, mf, aliased = "stop") {
  lines <- ilm_alias_lines(al)
  terms <- unique(unlist(lapply(al$groups, `[[`, "terms")))
  cells <- unique(unlist(lapply(terms, function(tt) ilm_alias_empty_cells(mf, tt))))
  fix <- if (al$smooth && all(vapply(al$groups, `[[`, TRUE, "all_smooth")))
    paste0("The smooths' unpenalised parts overlap: two smooths with the same ",
           "numeric `by` variable each contain that variable's main effect. ",
           "Give them different `by` variables, or combine them into one smooth")
  else if (al$smooth)
    paste0("A smooth's unpenalised part already holds its linear trend (and ",
           "a numeric `by` variable's main effect), so a parametric term in ",
           "the same direction duplicates it: drop the parametric term, or ",
           "the overlap between the smooths, and the smooth carries the effect")
  else if (length(cells))
    paste0("The interaction has combinations with no rows (",
           paste(utils::head(cells, 6), collapse = ", "),
           if (length(cells) > 6L) ", ..." else "",
           "), so some of its coefficients have nothing to estimate them: ",
           "merge levels, drop the interaction, or fit the cells as one factor")
  else
    paste0("Drop one of the terms involved (", ilm_and(sprintf("`%s`", terms)),
           "), or rewrite them so they no longer overlap")
  ## the other way out, where there is one
  alt <- if (al$smooth) {
    if (identical(aliased, "drop"))
      paste0(" aliased = \"drop\" does not apply: a smooth's unpenalised part ",
             "has no column that can sensibly be dropped.") else ""
  } else paste0(" Or set aliased = \"drop\" to drop the dependent ",
                if (length(al$dependent) == 1L) "column" else "columns",
                ", as lm() does.")
  paste0("The model's fixed-effect columns are not all separable, so it has ",
         "no unique fit (every standard error would come out NaN): ",
         paste(lines, collapse = "; "), ". ", fix, ".", alt)
}
