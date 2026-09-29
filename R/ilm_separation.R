## ---- separation: a level with no variation in the outcome --------------------
##
## When the outcome does not vary within a level of a categorical predictor --
## every count zero in one region, every trial a success in one arm -- that
## level's coefficient has no finite estimate: the likelihood keeps rising as
## the coefficient runs to minus (or plus) infinity. The optimiser stops
## somewhere out there, where the objective is flat, and every check can pass.
## Measured: three all-zero rows in one level of a Poisson fit gave a
## coefficient near minus infinity, all checks OK, and a predicted count of
## 0.00000005 for that level. The estimate, its standard error and anything
## predicted for the level mean nothing, and nothing said so.
##
## Two catches, as Craig ruled (item 5):
##   before the fit, from the data alone: a level (or a cell of two factors)
##     whose outcome takes no second value, named with what it is;
##   after the fit, cheaply: a large coefficient whose objective is flat in
##     the direction it ran, which also catches separation by a numeric
##     predictor or a combination no single level shows.
## Either gives the `separation` check FAIL, and ilm_interpret() says the
## term's effect is not estimable. The linear-programming test that finds
## every quasi-complete separation stays unbuilt unless a study shows these
## two miss binomial cases that matter.

## Which outcome would make a level separated, for this family. NULL where a
## level cannot be separated in this sense (gaussian, the survival families)
## or where a zero part already explains an all-zero level.
#' @keywords internal
#' @noRd
ilm_sep_rule <- function(fam, has_zero) {
  nm <- fam$name
  if (isTRUE(fam$ordinal)) return("ordinal")
  if (identical(nm, "binomial")) return("binomial")
  if (identical(nm, "multinomial")) return("multinomial")
  if (nm %in% c("poisson", "nbinom") && !has_zero) return("count")
  NULL
}

## The levels (and cells of two factors) whose outcome takes no second value.
## `y` is the response as the fit receives it; `w` the weights (trials for a
## binomial proportion). One row per separated level: term, level, what the
## outcome was there, and its rows.
#' @keywords internal
#' @noRd
ilm_sep_levels <- function(mf, y, w, rule, term_labels, J = NULL) {
  none <- data.frame(term = character(), level = character(), outcome = character(),
                     n = integer(), stringsAsFactors = FALSE)
  if (is.null(rule) || !length(term_labels)) return(none)
  if (is.null(w)) w <- rep(1, length(y))
  is_cat <- function(v) v %in% names(mf) && (is.factor(mf[[v]]) || is.character(mf[[v]]) ||
                                               is.logical(mf[[v]]))
  ## what the outcome was within one group of rows, or NULL when it varies
  flat <- function(rows) {
    yy <- y[rows]; ww <- w[rows]
    keep <- ww > 0 & !is.na(yy); yy <- yy[keep]; ww <- ww[keep]
    if (!length(yy)) return(NULL)
    switch(rule,
      binomial = if (all(yy <= 0)) "every trial a failure"
                 else if (all(yy >= 1)) "every trial a success" else NULL,
      count = if (all(yy == 0)) "every count zero" else NULL,
      ordinal = if (all(yy == min(y, na.rm = TRUE))) "every response in the lowest category"
                else if (all(yy == max(y, na.rm = TRUE))) "every response in the highest category"
                else NULL,
      multinomial = {
        seen <- unique(yy); miss <- setdiff(sort(unique(y)), seen)
        if (length(miss)) paste0("no response in ",
          if (length(miss) == 1L) "category " else "categories ",
          paste(miss, collapse = ", ")) else NULL
      })
  }
  out <- list()
  for (tt in term_labels) {
    vs <- ilm_unbq(strsplit(tt, ":", fixed = TRUE)[[1]])
    if (length(vs) > 2L || !all(vapply(vs, is_cat, TRUE))) next
    g <- if (length(vs) == 1L) as.character(mf[[vs]])
         else paste(as.character(mf[[vs[1]]]), as.character(mf[[vs[2]]]), sep = ":")
    for (lv in unique(g[!is.na(g)])) {
      rows <- which(g == lv)
      o <- flat(rows)
      if (!is.null(o)) out[[length(out) + 1L]] <- data.frame(term = tt, level = lv,
        outcome = o, n = length(rows), stringsAsFactors = FALSE)
    }
  }
  if (length(out)) do.call(rbind, out) else none
}

## After the fit: a coefficient far out on the link scale whose objective is
## flat in the direction it ran. Moving it five units further changes the
## objective by less than 1e-3 where it has run off to infinity; a real
## estimate that large moves it by far more. One evaluation per candidate, so
## it costs nothing for a fit with no large coefficient.
#' @keywords internal
#' @noRd
ilm_sep_flat <- function(fit, big = 8, step = 5, tol = 1e-3) {
  none <- data.frame(coef = character(), estimate = numeric(), stringsAsFactors = FALSE)
  obj <- fit$obj; pe <- fit$opt$par
  if (is.null(obj) || is.null(pe)) return(none)
  jb <- which(names(pe) == "beta")
  if (!length(jb)) return(none)
  cand <- jb[abs(pe[jb]) > big & is.finite(pe[jb])]
  if (!length(cand)) return(none)
  f0 <- tryCatch(obj$fn(pe), error = function(e) NA_real_)
  if (!is.finite(f0)) return(none)
  nms <- ilm_beta_names(fit)
  hit <- vapply(cand, function(j) {
    p1 <- pe; p1[j] <- p1[j] + step * sign(p1[j])
    f1 <- tryCatch(obj$fn(p1), error = function(e) NA_real_)
    is.finite(f1) && abs(f1 - f0) < tol
  }, TRUE)
  invisible(tryCatch(obj$fn(pe), error = function(e) NULL))   # the tape back at the optimum
  if (!any(hit)) return(none)
  k <- cand[hit] - min(jb) + 1L
  data.frame(coef = if (length(nms) >= max(k)) nms[k] else paste0("beta[", k, "]"),
             estimate = unname(pe[cand[hit]]), stringsAsFactors = FALSE)
}

## The fixed coefficients' names in the order of beta: a multinomial's are
## the design's columns for each category in turn.
#' @keywords internal
#' @noRd
ilm_beta_names <- function(fit) {
  xn <- colnames(fit$X)
  if (is.null(xn)) return(character())
  C <- fit$C %||% 1L
  if (identical(fit$family$name, "multinomial") && C > 1L) {
    ## as the fit names them: the first C categories, each against the average
    cats <- fit$ylevels[seq_len(C)]
    if (length(cats) != C || anyNA(cats)) cats <- as.character(seq_len(C))
    as.vector(outer(xn, cats, function(a, b) paste0(b, ":", a)))
  } else xn
}

## The check row and the record, added to a fit. `lv` from ilm_sep_levels(),
## `fl` from ilm_sep_flat(). Only a FAIL is added: a model with nothing
## separated gets an OK row, so the check is always there to be read.
#' @keywords internal
#' @noRd
ilm_sep_check <- function(fit, lv, fl) {
  ## a flat coefficient that is a separated level's own is said once, as the level
  if (nrow(lv) && nrow(fl)) {
    own <- paste0(ilm_unbq(lv$term), lv$level)
    fl <- fl[!vapply(fl$coef, function(cf) any(endsWith(cf, own)), TRUE), , drop = FALSE]
  }
  any_sep <- nrow(lv) > 0L || nrow(fl) > 0L
  detail <- if (!any_sep) "no level without variation in the outcome, and no coefficient on a flat likelihood"
  else paste(c(
    if (nrow(lv)) sprintf("%s = %s: %s (%d row%s)", lv$term, lv$level, lv$outcome,
                          lv$n, ifelse(lv$n == 1L, "", "s")),
    if (nrow(fl)) sprintf("coefficient %s = %s with the likelihood flat beyond it",
                          fl$coef, formatC(fl$estimate, format = "g", digits = 3))),
    collapse = "; ")
  what <- if (nrow(lv)) sprintf("%s '%s'", lv$term[1], lv$level[1]) else sprintf("coefficient %s", fl$coef[1])
  fit$checks <- ilm_add_check(fit$checks, "separation", if (any_sep) "FAIL" else "OK", detail,
    if (any_sep) paste0("the outcome does not vary within ", what,
                        if (nrow(lv) + nrow(fl) > 1L) " (and the others listed)" else "",
                        ", so the coefficient has no finite estimate: it runs off towards infinity,",
                        " and its estimate, its standard error and anything predicted for it mean nothing")
    else "",
    if (any_sep) paste0("merge ", if (nrow(lv)) sprintf("'%s' with a neighbouring level of %s", lv$level[1], lv$term[1])
                        else "the sparse level with a neighbouring one",
                        ", drop its rows, or remove the term; the other coefficients are then estimable")
    else "")
  fit$separation <- list(levels = lv, flat = fl)
  fit$ok <- !any(fit$checks$status == "FAIL")
  fit
}

## The levels a factor declared that no row used. The model frame drops them
## (drop.unused.levels), so they have no coefficient and cannot be predicted
## for; that used to happen without a word.
#' @keywords internal
#' @noRd
ilm_empty_levels <- function(data, mf, vars) {
  if (!is.data.frame(data)) return(list())
  out <- list()
  for (v in intersect(vars, names(data))) {
    if (!is.factor(data[[v]]) || !v %in% names(mf) || !is.factor(mf[[v]])) next
    e <- setdiff(levels(data[[v]]), levels(mf[[v]]))
    if (length(e)) out[[v]] <- e
  }
  out
}

#' @keywords internal
#' @noRd
ilm_sep_flat_none <- function()
  data.frame(coef = character(), estimate = numeric(), stringsAsFactors = FALSE)

## One sentence for the message at the fit: what is separated.
#' @keywords internal
#' @noRd
ilm_sep_words <- function(sep) {
  lv <- sep$levels; fl <- sep$flat
  bits <- c(if (nrow(lv)) sprintf("%s '%s' (%s)", lv$term, lv$level, lv$outcome),
            if (nrow(fl)) sprintf("coefficient %s, on a flat likelihood", fl$coef))
  paste0("the outcome does not vary within ", ilm_and(utils::head(bits, 4)),
         if (length(bits) > 4L) sprintf(" and %d more", length(bits) - 4L) else "",
         ", so ", if (length(bits) == 1L) "that coefficient has" else "those coefficients have",
         " no finite estimate.")
}
