## ---- remedies the code can act on -------------------------------------------
##
## Every check names its remedy in words, in the `suggestion` column of
## fit$checks. This file turns those words into changes a refit can make:
## each remedy is a set of ilm_model() arguments, or none when only a person
## can make the change -- pooling levels, merging categories. The rules read
## the fit itself, not the words, so the words can be edited without breaking
## a remedy; test-remedies.R fails when a check has no rule.
##
## Each remedy has a TIER, which says what applying it would change:
##   numerical   the same model, fitted harder (more optimiser restarts)
##   structural  a different random-effect or variance structure, or a
##               different way of estimating the variances, whose fixed
##               effects mean what they meant before: a term at a variance of
##               zero removed, a covariance of lower rank, a dispersion or
##               zero part, a penalty that keeps a covariance off its edge,
##               REML in place of ML
##   estimand    a change to what the fixed effects estimate, or to what their
##               standard errors account for: a random effect whose variance
##               is not zero removed, categories merged, levels pooled
## Remedies can also come from another package's diagnostics, through
## ilm_remedy_table(); they are applied the same way, and the tiers mean the
## same thing there.
## The tier is for whoever decides which remedies may be tried without being
## asked: the person running the analysis. ilm_apply_remedy() applies
## whichever remedy it is given.
##
## Each remedy also has a KEY, `kind/target`, which names the change itself
## and not the check that asked for it: one remedy often answers several
## checks, and an id is renumbered whenever tables are combined or a
## standalone check is included. A script that applies "drop/g" applies the
## same remedy under the next version or another combination, or stops. The
## target is the term as the fit names it -- as its checks name it -- since a
## nested (1 | a/b) is two terms from one written bar. A change that can be
## read off its arguments alone takes its key from them (ilm_rem_args_key()),
## so the same change has the same key whichever package lists it.
##
## The generics, c(), print() and the table's assembly are illumex's, shared
## with its data checks; illume supplies the methods for a fitted model.

ilm_rem_tiers <- c("numerical", "structural", "estimand")

## what print() says under a model's table, and in place of a change for a
## remedy made by hand
ilm_rem_tier_note <- paste0(
  "Refit with one by ilm_apply_remedy(fit, <this list>, id or key). ",
  "Numerical is the same model fitted harder; structural changes the ",
  "random-effect or variance structure, or how the variances are ",
  "estimated, and not what the fixed effects mean; ",
  "estimand changes what they estimate or what their standard errors ",
  "account for, so apply one of those only by choice.")
ilm_rem_by_hand <- "by hand: not something a refit can do"

## the checks that have a rule below, by the name before any `[term]`
ilm_rem_known <- c("category_counts", "weights_type", "re_levels",
                   "obs_per_level", "latent_budget", "obs_per_ar_latent",
                   "optimizer", "gradient", "hessian", "variance_boundary",
                   "smooth_shrinkage", "sigma_rank", "sigma_within",
                   "rho_boundary", "dispersion_limit", "parameter_aliasing",
                   "separation")

## one remedy: its tier, what it does in words, the ilm_model() arguments
## that make it -- NULL when it has to be done by hand -- and its key, which
## a change read off its arguments need not be given
#' @keywords internal
#' @noRd
ilm_rem <- function(tier, remedy, args = NULL, key = NULL) {
  if (!tier %in% ilm_rem_tiers) stop("internal: unknown remedy tier ", tier)
  if (is.null(key)) key <- ilm_rem_args_key(args)
  if (is.null(key)) stop("internal: a remedy without a key: ", remedy)
  list(tier = tier, remedy = remedy, args = args, key = key)
}

## The key of a change that says what it is by its arguments alone, or NULL.
## illume's own rules and another package's table both key these changes
## here, so c() lists a change both name once, as it did when changes were
## matched by their text.
#' @keywords internal
#' @noRd
ilm_rem_args_key <- function(args) {
  if (is.null(args) || !length(args)) return(NULL)
  nm <- sort(names(args))
  one <- function(v) is.character(v) && length(v) == 1L && !is.na(v)
  rhs <- function(f) paste(deparse(f[[length(f)]], width.cutoff = 500L), collapse = " ")
  if (identical(nm, "restarts")) return(illumex::ilm_remedy_key("restarts", "optimiser"))
  if (identical(nm, "boundary") && one(args$boundary))
    return(illumex::ilm_remedy_key("boundary", args$boundary))
  if (identical(nm, "family") && one(args$family))
    return(illumex::ilm_remedy_key("family", args$family))
  if (identical(nm, "reml") && is.logical(args$reml))
    return(illumex::ilm_remedy_key("reml", if (isTRUE(args$reml)) "on" else "off"))
  if (identical(nm, "dispformula") && inherits(args$dispformula, "formula"))
    return(illumex::ilm_remedy_key("disp", rhs(args$dispformula)))
  if (all(nm %in% c("ziformula", "zi_type")) && inherits(args$ziformula, "formula")) {
    ty <- if (is.null(args$zi_type)) "inflated" else args$zi_type
    r <- rhs(args$ziformula)
    return(illumex::ilm_remedy_key("zi", ty, if (r != "1") r))
  }
  if (identical(nm, "ar") && is.null(args$ar))
    return(illumex::ilm_remedy_key("drop", "ar"))
  NULL
}

## ---- reading the fit ----------------------------------------------------------

#' @keywords internal
#' @noRd
ilm_rem_groups <- function(fit)
  names(fit$re)[vapply(fit$re, function(e) !identical(e$kind, "basis"), TRUE)]

## The bar each grouping term came from. ilm_model() names a bar's term after
## its grouping factor as it builds them, in the order the bars were written,
## and after every smooth -- so the i-th grouping term is the i-th bar.
#' @keywords internal
#' @noRd
ilm_rem_bars <- function(fit) {
  grp <- ilm_rem_groups(fit)
  if (!length(fit$bars) || length(fit$bars) != length(grp)) return(list())
  stats::setNames(fit$bars, grp)
}

## The term of the formula AS WRITTEN that a fitted term came from, or NULL.
## A smooth's label -- s(x) -- is not how s(x, k = 5) was written, so terms
## are matched by what they parse to rather than by their text. For a bar the
## answer also carries the OTHER bars the same written term expands to: a
## nested (1 | a/b) is the two terms (1 | b:a) and (1 | a) to the fit, and
## one of them can only be removed by writing the other out.
#' @keywords internal
#' @noRd
ilm_rem_term <- function(fit, term) {
  fo <- fit$formula
  if (!inherits(fo, "formula")) return(NULL)
  tl <- attr(stats::terms(fo), "term.labels")
  bar <- ilm_rem_bars(fit)[[term]]
  if (!is.null(bar)) {
    for (t in tl[grepl("|", tl, fixed = TRUE)]) {
      ex <- tryCatch(ilm_findbars(stats::as.formula(paste0("~ (", t, ")"))),
                     error = function(e) NULL)
      hit <- vapply(ex, function(b) identical(b, bar), TRUE)
      if (sum(hit) == 1L)
        return(structure(t, rest = vapply(ex[!hit], function(b)
          paste0("(", ilm_term_text(b), ")"), "")))
    }
    return(NULL)
  }
  ## a smooth: evaluate each smooth term with mgcv's constructors and compare
  ## the label it builds with the one the fit carries
  env <- new.env(parent = environment(fo) %||% globalenv())
  for (f in c("s", "te", "ti", "t2")) assign(f, getExportedValue("mgcv", f), envir = env)
  for (t in tl) {
    e <- tryCatch(str2lang(t), error = function(err) NULL)
    if (!is.call(e) || !sub("^mgcv::", "", deparse1(e[[1L]])) %in% c("s", "te", "ti", "t2"))
      next
    lab <- tryCatch(ilm_smooth_label(eval(e, env)), error = function(err) NULL)
    if (identical(lab, term)) return(t)
  }
  NULL
}

## The formula with one term taken out and, optionally, text put in its
## place, or NULL when the term cannot be found as written. update() keeps an
## offset and the formula's environment, both of which rebuilding the formula
## from its term labels would lose.
#' @keywords internal
#' @noRd
ilm_rem_formula <- function(fit, term, add = NULL) {
  t <- ilm_rem_term(fit, term)
  if (is.null(t)) return(NULL)
  fo <- fit$formula
  add <- c(attr(t, "rest"), add)
  t <- as.character(t)
  chg <- paste0(". ~ . - ", if (grepl("|", t, fixed = TRUE)) paste0("(", t, ")") else t,
                if (length(add)) paste0(" + ", paste(add, collapse = " + ")) else "")
  out <- tryCatch(stats::update(fo, stats::as.formula(chg)), error = function(e) NULL)
  if (is.null(out) || t %in% attr(stats::terms(out), "term.labels")) return(NULL)
  environment(out) <- environment(fo)
  out
}

## re_struct as a person would write it: every term that differs from the
## default, and the one being changed. A structure the fit was given for
## another term is carried along rather than lost.
#' @keywords internal
#' @noRd
ilm_rem_struct <- function(fit, term, spec) {
  rs <- lapply(fit$re_struct, function(s) {
    out <- list(type = s$type)
    if (identical(s$type, "rr")) out$rank <- as.integer(s$rank)
    if (isFALSE(s$d_cor)) out$d_cor <- FALSE
    out
  })
  new <- rs[[term]]
  for (k in names(spec)) new[[k]] <- spec[[k]]
  if (!identical(new$type, "rr")) new$rank <- NULL
  rs[[term]] <- new
  keep <- vapply(names(rs), function(nm)
    nm == term || !identical(rs[[nm]], list(type = "us")), TRUE)
  rs[keep]
}

## latent values per term, as the latent_budget check counts them
#' @keywords internal
#' @noRd
ilm_rem_latent <- function(fit) {
  lat <- vapply(names(fit$re), function(nm) {
    e <- fit$re[[nm]]
    as.numeric(e$nl * e$d * ilm_str_width(fit$re_struct[[nm]], fit$C))
  }, 0)
  if (!is.null(fit$ar)) lat["ar"] <- ilm_ar_nlat(fit$ar) * fit$C
  lat
}

## whether every standard deviation of a grouping term is at zero, by the
## threshold the variance_boundary check uses
#' @keywords internal
#' @noRd
ilm_rem_at_zero <- function(fit, term) {
  S <- fit$Sigma[[term]]
  !is.null(S) && all(sqrt(pmax(diag(as.matrix(S)), 0)) < 1e-3)
}

## ---- remedies shared by several checks ---------------------------------------

#' @keywords internal
#' @noRd
ilm_rem_restarts <- function(fit) {
  cur <- tryCatch(eval(fit$call$restarts, environment(fit$formula)),
                  error = function(e) NULL)
  if (!is.numeric(cur) || length(cur) != 1L || !is.finite(cur)) cur <- 3L
  nr <- as.integer(max(10L, 2L * cur))
  ilm_rem("numerical",
    sprintf("refit with %d optimiser restarts from the solution reached, not %d: the same model, fitted harder",
            nr, as.integer(cur)),
    list(restarts = nr))
}

## boundary = "avoid", named with what it costs, as the checks name it; not
## offered to a fit that already used it, or that has no covariance to keep
## off its boundary
#' @keywords internal
#' @noRd
ilm_rem_avoid <- function(fit) {
  if (identical(fit$boundary, "avoid") || !length(ilm_rem_groups(fit))) return(list())
  list(ilm_rem("structural",
    paste0("refit with boundary = \"avoid\", a small penalty that keeps every ",
           "random-effect covariance inside its range. Each variance is then ",
           "assumed nonzero rather than estimated at zero, so do not test ",
           "whether it is; it comes out larger, and for a binary or ",
           "categorical outcome the fixed effects come out a little further ",
           "from zero -- markedly so when a category is rare"),
    list(boundary = "avoid")))
}

## Remove a grouping term's bar. At a variance of zero the fit is the same
## without it, which the boundary tests check; anywhere else it changes what
## the standard errors account for, and for a non-gaussian response what the
## fixed effects mean. Which of the two is read off the fit, whichever check
## asked: a term can fail its level count and sit at zero as well.
#' @keywords internal
#' @noRd
ilm_rem_drop <- function(fit, term) {
  at_zero <- ilm_rem_at_zero(fit, term)
  fo <- ilm_rem_formula(fit, term)
  why <- if (at_zero)
    sprintf("drop '%s': its variance is estimated at zero, so the fixed effects are the same without it and their standard errors barely move. Keep it instead if the design calls for it -- repeated measures, say -- since a zero estimate is not evidence of no clustering", term)
  else sprintf("drop '%s' from the model. Its variance is not zero, so the standard errors stop accounting for the grouping, and for a non-gaussian response the fixed effects change meaning", term)
  ilm_rem(if (at_zero) "structural" else "estimand", why,
          if (is.null(fo)) NULL else list(formula = fo),
          illumex::ilm_remedy_key("drop", term))
}

## keep a grouping term's intercept and remove its random slopes, which is no
## change when every slope variance is at zero, by sigma_within's threshold
#' @keywords internal
#' @noRd
ilm_rem_noslope <- function(fit, term) {
  bar <- ilm_rem_bars(fit)[[term]]
  if (is.null(bar) || fit$re[[term]]$d < 2L) return(list())
  Sd <- fit$Sigma_d[[term]]
  sv <- if (is.null(Sd)) numeric(0) else sqrt(pmax(diag(as.matrix(Sd)), 0))
  at_zero <- length(sv) > 1L && all(sv[-1L] / sv[1L] < 1e-3)
  lhs <- stats::terms(stats::as.formula(call("~", bar[[2L]])))
  fo <- if (attr(lhs, "intercept") == 1L)
    ilm_rem_formula(fit, term, paste0("(1 | ", ilm_term_text(bar[[3L]]), ")"))
  why <- if (at_zero)
    sprintf("remove the random slope of '%s' and keep its random intercept: the slope variance is estimated at zero, so the fit is the same without it", term)
  else sprintf("remove the random slope of '%s' and keep its random intercept. The slope variance is not zero, so the standard error of the fixed slope then understates its uncertainty", term)
  list(ilm_rem(if (at_zero) "structural" else "estimand", why,
               if (is.null(fo)) NULL else list(formula = fo),
               illumex::ilm_remedy_key("noslope", term)))
}

#' @keywords internal
#' @noRd
ilm_rem_dcor <- function(fit, term) {
  s <- fit$re_struct[[term]]
  if (fit$re[[term]]$d < 2L || isFALSE(s$d_cor)) return(list())
  list(ilm_rem("structural",
    sprintf("drop the correlation between the intercept and the slope of '%s' (d_cor = FALSE), keeping both variances", term),
    list(re_struct = ilm_rem_struct(fit, term, list(d_cor = FALSE))),
    illumex::ilm_remedy_key("nocor", term)))
}

#' @keywords internal
#' @noRd
ilm_rem_pool <- function(term)
  ilm_rem("estimand", sprintf(
    "pool levels of '%s' that belong together, so that each level carries more observations; the grouping then means something different, so choose the pooling by what the levels are",
    term), key = illumex::ilm_remedy_key("pool", term))

## ---- the rules -----------------------------------------------------------------

## The remedies for one row of fit$checks. Only a row that is not OK has any;
## INCONCLUSIVE says a check could not run, which is not a problem to remedy.
#' @keywords internal
#' @noRd
ilm_rem_rules <- function(fit, check, status, detail = "", suggestion = "") {
  if (!status %in% c("WARN", "FAIL", "BOUNDARY")) return(list())
  base <- sub("\\[.*$", "", check)
  term <- if (grepl("\\[", check)) sub("^[^[]*\\[(.*)\\]$", "\\1", check) else NA_character_
  C <- fit$C
  if (!base %in% ilm_rem_known)
    ## a check this file does not know yet: its own words, done by hand, and
    ## never tried without being asked
    return(list(ilm_rem("estimand", suggestion, key = illumex::ilm_remedy_key("advice", check))))
  switch(base,

  category_counts = {
    w <- if (is.null(fit$weights)) rep(1, length(fit$y)) else fit$weights
    tabc <- vapply(seq_len(fit$J), function(j) sum(w[as.integer(fit$y) == j]), 0)
    k <- which.min(tabc)
    list(ilm_rem("estimand", sprintf(
      "merge the rarest category, '%s' (%g of %g), with another before refitting. The model then describes a different outcome, so choose the merge by what the categories mean",
      fit$ylevels[k], tabc[k], sum(tabc)),
      key = illumex::ilm_remedy_key("merge_category", fit$ylevels[k])))
  },

  weights_type = list(ilm_rem("estimand", paste0(
    "if the weights count replicate observations, aggregate them so that they ",
    "are whole numbers of rows; if they are sampling weights, pass design = ",
    "ilm_design(weights = , ids = , strata = ) and read the standard errors ",
    "from ilm_svy_coef() -- which needs the clusters and strata, which only ",
    "you know",
    if (length(ilm_rem_groups(fit))) ", and refuses a model with random effects" else ""),
    key = illumex::ilm_remedy_key("weights", "declare"))),

  re_levels = {
    e <- fit$re[[term]]; s <- fit$re_struct[[term]]; out <- list()
    if (C >= 2L && s$type %in% c("us", "diag")) {
      r <- max(1L, ilm_rec_rank(C, e$nl))
      out <- list(ilm_rem("structural", sprintf(
        "give '%s' a reduced-rank category covariance, rr(%d): %d parameters instead of %d",
        term, r, ilm_nrr(C, r), ilm_str_npar(s, C)),
        list(re_struct = ilm_rem_struct(fit, term, list(type = "rr", rank = r))),
        illumex::ilm_remedy_key("rank", term, r)))
    } else if (C >= 2L && isTRUE(s$rank > 1L)) {
      r <- max(1L, min(s$rank - 1L, ilm_rec_rank(C, e$nl)))
      out <- list(ilm_rem("structural", sprintf(
        "lower the rank of the category covariance of '%s' from %d to %d",
        term, s$rank, r),
        list(re_struct = ilm_rem_struct(fit, term, list(type = "rr", rank = r))),
        illumex::ilm_remedy_key("rank", term, r)))
    }
    c(out, ilm_rem_dcor(fit, term), ilm_rem_noslope(fit, term),
      list(ilm_rem_pool(term), ilm_rem_drop(fit, term)))
  },

  obs_per_level = list(ilm_rem_pool(term), ilm_rem_drop(fit, term)),

  latent_budget = {
    lat <- ilm_rem_latent(fit); big <- names(lat)[which.max(lat)]
    if (identical(big, "ar")) list(
      ilm_rem("structural", "coarsen the time grid given to the AR term, so that several observations share each latent value",
              key = illumex::ilm_remedy_key("coarsen", "ar")),
      ilm_rem("estimand", "drop the AR term: the correlation over time is then not modelled, and the standard errors stop accounting for it",
              list(ar = NULL)))
    else {
      e <- fit$re[[big]]; s <- fit$re_struct[[big]]; w <- ilm_str_width(s, C)
      N <- if (is.null(fit$weights)) length(fit$y) else sum(fit$weights)
      keep <- max(1L, as.integer(floor((lat[[big]] - (sum(lat) - N / 5)) / (e$nl * e$d))))
      out <- list()
      if (keep < w)
        out <- list(ilm_rem("structural", sprintf(
          "give '%s', the term with the most latent values, a category covariance of rank %d (rr(%d)), which cuts its latent values from %d to %d",
          big, keep, keep, as.integer(lat[[big]]), as.integer(e$nl * e$d * keep)),
          list(re_struct = ilm_rem_struct(fit, big, list(type = "rr", rank = keep))),
          illumex::ilm_remedy_key("rank", big, keep)))
      if (identical(e$kind, "basis"))
        c(out, list(ilm_rem("structural", sprintf(
          "lower the basis dimension k of the smooth '%s', which has the most latent values", big),
          key = illumex::ilm_remedy_key("basis", big))))
      else c(out, ilm_rem_noslope(fit, big),
             list(ilm_rem_pool(big), ilm_rem_drop(fit, big)))
    }
  },

  obs_per_ar_latent = list(
    ilm_rem("structural", "coarsen the time grid given to the AR term, so that several observations share each latent value",
            key = illumex::ilm_remedy_key("coarsen", "ar")),
    ## s(time) describes the past and does not forecast: beyond the data a
    ## smooth runs on its last slope, where a correlated process reverts,
    ## so it is offered as a description only
    ilm_rem("estimand", "to DESCRIBE the past only, not to forecast: replace the AR term with a smooth of time, s(time), plus a random slope on time for the group -- a trend and each group's departure from it in place of a correlated process. Beyond the data a smooth runs on at its last slope, so do not use it for forecasts",
            key = illumex::ilm_remedy_key("smooth_time", "ar"))),

  ## where every start broke down, more restarts from the same place is
  ## advice the fit has already taken
  optimizer = if (status == "FAIL")
    c(ilm_rem_rescale(fit),
      if (!isTRUE(fit$opt$all_starts)) list(ilm_rem_restarts(fit))) else list(),

  gradient = c(ilm_rem_rescale(fit),
    if (!isTRUE(fit$opt$all_starts)) list(ilm_rem_restarts(fit)),
    list(ilm_rem("structural", "simplify the random-effect structure; the other checks that are not OK name the term",
                 key = illumex::ilm_remedy_key("simplify", "random_effects")))),

  hessian = if (status == "BOUNDARY") {
    held <- intersect(fit$hessian_held, ilm_rem_groups(fit))
    ## a dispersion held alone is dispersion_limit's to remedy, and
    ## boundary = "avoid" does nothing for it
    if (!length(setdiff(fit$hessian_held, "dispersion"))) list()
    ## a term held at a variance of zero comes out; one held at a correlation
    ## of +/-1 is sigma_rank's or sigma_within's to simplify
    else c(lapply(held[vapply(held, function(h) ilm_rem_at_zero(fit, h), TRUE)],
                  function(h) ilm_rem_drop(fit, h)),
           ilm_rem_avoid(fit))
  } else c(ilm_rem_rescale(fit), list(ilm_rem_restarts(fit)), ilm_rem_avoid(fit),
           list(ilm_rem("structural", "simplify the random-effect structure; the other checks that are not OK name the term",
                 key = illumex::ilm_remedy_key("simplify", "random_effects")))),

  variance_boundary = {
    zero <- Filter(function(h) ilm_rem_at_zero(fit, h), ilm_rem_groups(fit))
    c(lapply(zero, function(h) ilm_rem_drop(fit, h)), ilm_rem_avoid(fit))
  },

  smooth_shrinkage = {
    sob <- fit$smooths[[term]]; t <- ilm_rem_term(fit, term)
    nx <- if (is.null(sob) || is.null(sob$Xf)) NA_integer_ else ncol(sob$Xf)
    v <- if (is.null(sob)) character(0) else sob$sm$term
    by <- if (is.null(sob) || is.null(sob$sm$by)) "NA" else sob$sm$by
    ## mgcv writes the variable as code -- `my x`, log(x) -- so it goes into
    ## the formula as it is, and is looked up in the model frame without its
    ## backticks. The unpenalised column is checked to BE that straight line
    ## rather than assumed to be: a higher-order penalty leaves a curve.
    lin <- FALSE
    if (!is.null(t) && identical(nx, 1L) && length(v) == 1L && identical(by, "NA")) {
      x <- fit$model[[ilm_unbq(v)]]
      lin <- is.numeric(x) && length(x) == nrow(sob$Xf) &&
        isTRUE(abs(stats::cor(sob$Xf[, 1L], x)) > 1 - 1e-6)
    }
    if (lin) {
      fo <- ilm_rem_formula(fit, term, v)
      list(ilm_rem("structural", sprintf(
        "replace '%s' with %s: the smooth was shrunk to a straight line in %s, so this is the same fit with one parameter fewer",
        term, v, v), if (is.null(fo)) NULL else list(formula = fo),
        illumex::ilm_remedy_key("linear", term)))
    } else if (!is.null(t) && identical(nx, 0L)) {
      fo <- ilm_rem_formula(fit, term)
      list(ilm_rem("structural", sprintf(
        "drop '%s': it was shrunk to nothing and has no unpenalised part, so this is the same fit without it",
        term), if (is.null(fo)) NULL else list(formula = fo),
        illumex::ilm_remedy_key("drop", term)))
    } else list(ilm_rem("structural", sprintf(
      "replace the smooth that '%s' belongs to with its unpenalised part, the only part of it these data support",
      term), key = illumex::ilm_remedy_key("unpenalised", term)))
  },

  sigma_rank = {
    s <- fit$re_struct[[term]]
    ev <- sort(eigen(as.matrix(fit$Sigma[[term]]), symmetric = TRUE,
                     only.values = TRUE)$values, decreasing = TRUE)
    kp <- if (identical(s$type, "rr")) seq_len(s$rank) else seq_along(ev)
    eff <- max(1L, sum(ev[kp] / max(ev[kp]) >= 1e-3))
    cur <- if (identical(s$type, "rr")) s$rank else C
    c(if (eff < cur) list(ilm_rem("structural", sprintf(
        "give '%s' a category covariance of rank %d (rr(%d)) in place of %s: the %s dropped %s no variance",
        term, eff, eff, ilm_str_label(s),
        if (cur - eff > 1L) "directions" else "direction",
        if (cur - eff > 1L) "have" else "has"),
        list(re_struct = ilm_rem_struct(fit, term, list(type = "rr", rank = eff))),
        illumex::ilm_remedy_key("rank", term, eff))),
      ilm_rem_avoid(fit))
  },

  sigma_within = {
    Sd <- as.matrix(fit$Sigma_d[[term]]); sv <- sqrt(pmax(diag(Sd), 0))
    rel <- sv[-1L] / sv[1L]
    if (any(rel < 1e-3))
      c(ilm_rem_noslope(fit, term), ilm_rem_avoid(fit))
    else c(ilm_rem_dcor(fit, term), ilm_rem_noslope(fit, term),
           ilm_rem_avoid(fit))
  },

  ## a negative binomial at k = infinity is the Poisson; a beta at phi =
  ## infinity has no simpler family, only a term to look at
  dispersion_limit = if (identical(fit$family$name, "nbinom"))
    list(ilm_rem("structural", paste0(
      "refit with family = \"poisson\": the data show no overdispersion ",
      "beyond the model's other terms, so this is the same fit with one ",
      "parameter fewer"), list(family = "poisson")))
  else list(ilm_rem("structural", paste0(
      "look at the term absorbing the variation: coarsen a correlation over ",
      "time to a grid several observations share, or drop a random effect ",
      "with one observation per level"), key = illumex::ilm_remedy_key("inspect", "dispersion"))),

  rho_boundary = list(
    ilm_rem("structural", "coarsen the time grid given to the AR term, so that neighbouring observations share a latent value",
            key = illumex::ilm_remedy_key("coarsen", "ar")),
    ilm_rem("estimand", "drop the AR term: the correlation over time is then not modelled, and the standard errors stop accounting for it",
            list(ar = NULL))),

  separation = {
    sp <- fit$separation
    lv <- if (!is.null(sp) && nrow(sp$levels)) sp$levels[1, ] else NULL
    list(
      ilm_rem("estimand", if (!is.null(lv))
        sprintf("merge '%s' with a neighbouring level of %s, or drop its rows: the level whose outcome does not vary then has a coefficient that can be estimated, or none",
                lv$level, lv$term)
        else "merge the level whose outcome does not vary with a neighbouring one, or drop its rows",
        key = if (!is.null(lv)) illumex::ilm_remedy_key("merge_level", lv$term, lv$level)
              else illumex::ilm_remedy_key("merge_level", "separated")),
      ilm_rem("structural", if (!is.null(lv))
        sprintf("remove %s from the model, if its effect is not the question", lv$term)
        else "remove the term whose coefficient runs off, if its effect is not the question",
        key = illumex::ilm_remedy_key("remove", if (!is.null(lv)) lv$term else "separated")))
  },

  parameter_aliasing = {
    ## from the FIRST parenthesis: a parameter name can hold its own,
    ## as beta.(Intercept) does
    pair <- if (grepl("\\(.*<->.*\\)$", detail)) sub("^[^(]*\\((.*)\\)$", "\\1", detail) else ""
    list(ilm_rem("structural", paste0(
      "remove one of the two terms whose parameters these data cannot separate",
      if (nzchar(pair)) paste0(" (", pair, ")") else ""),
      key = illumex::ilm_remedy_key("unalias", if (nzchar(pair)) pair else "terms")))
  })
}

## The remedies for a standalone check's result. These checks simulate, so
## they are run by choice and passed in rather than run here.
#' @keywords internal
#' @noRd
ilm_rem_standalone <- function(fit, which, res) {
  if (!res$status %in% c("WARN", "FAIL")) return(list())
  fam <- fit$family$name
  has_disp <- isTRUE(fit$family$n_disp > 0)
  disp_add <- function(term) {
    old <- attr(fit$Zd, "formula")
    f <- if (is.null(old)) stats::as.formula(paste("~", term))
         else stats::update(old, stats::as.formula(paste("~ . +", term)))
    environment(f) <- environment(fit$formula)
    f
  }
  switch(which,

  dispersion = c(
    if (identical(fam, "poisson")) list(ilm_rem("structural",
      "refit as a negative binomial, family = \"nbinom\", whose variance can exceed its mean: the coefficients keep their meaning, and the standard errors take the extra variation into account",
      list(family = "nbinom"))),
    if (identical(fam, "binomial")) list(ilm_rem("structural",
      "with several trials per row, add a random effect for each row -- a row identifier as a grouping factor, (1 | row) -- to take up variation between rows",
      key = illumex::ilm_remedy_key("row_effect", "rows"))),
    list(ilm_rem("estimand",
      "look for a predictor the model leaves out, which makes a response look more variable than its family allows; ilm_check_omitted() tests candidates",
      key = illumex::ilm_remedy_key("omitted", "predictor")))),

  zeros = if (is.null(fit$Zzi)) list(
    ilm_rem("structural",
      "add a zero part with a constant excess-zero probability, ziformula = ~ 1, as a mixture in which some units were never at risk (zi_type = \"inflated\"); the count coefficients then describe the units that were -- see ilm_zi_coef()",
      list(ziformula = stats::as.formula("~ 1", env = environment(fit$formula)))),
    ilm_rem("structural",
      "the same zero part as a hurdle (zi_type = \"hurdle\"): every zero from the zero process, and a count that cannot be zero for the rest. Choose between the two by what the zeros mean, not by which fits better",
      list(ziformula = stats::as.formula("~ 1", env = environment(fit$formula)),
           zi_type = "hurdle")),
    ilm_rem("structural",
      "let the excess-zero probability depend on a predictor thought to drive it: ziformula = ~ x",
      key = illumex::ilm_remedy_key("zi", "predictor")))
  else list(ilm_rem("structural",
    "let the zero part depend on a predictor thought to drive the excess: a ziformula with that predictor in it",
    key = illumex::ilm_remedy_key("zi", "predictor"))),

  variance = {
    out <- list()
    if (isTRUE(res$p_trend < 0.05)) out <- c(out,
      if (has_disp) list(ilm_rem("structural",
        "let the spread change with the fitted mean, dispformula = ~ mu, which makes the dispersion a power of it",
        list(dispformula = disp_add("mu"))))
      else if (identical(fam, "poisson")) list(ilm_rem("structural",
        "refit as a negative binomial, family = \"nbinom\", whose variance grows faster than its mean",
        list(family = "nbinom")))
      else list(ilm_rem("structural",
        "use a family whose variance grows with its mean",
        key = illumex::ilm_remedy_key("family", "variance"))))
    if (isTRUE(res$p_ratio < 0.05)) {
      col <- if (is.null(res$by)) NA_character_ else res$by
      out <- c(out, if (has_disp && !is.na(col)) list(ilm_rem("structural",
        sprintf("give each level of %s its own spread, dispformula = ~ %s", col, ilm_bq(col)),
        list(dispformula = disp_add(ilm_bq(col)))))
      else list(ilm_rem("structural",
        "give each group its own spread, with a dispersion formula in the variable that defines the groups",
        key = illumex::ilm_remedy_key("disp", "groups"))))
    }
    out
  })
}

## ---- the functions ------------------------------------------------------------

## what makes a list of remedies belong to one fit
#' @keywords internal
#' @noRd
ilm_rem_id <- function(fit)
  paste(format(fit$opt$objective, digits = 15), nrow(fit$X),
        paste(deparse(fit$call), collapse = " "))

#' @keywords internal
#' @noRd
ilm_rem_change <- function(args) {
  if (is.null(args)) return("")
  paste(vapply(names(args), function(a)
    paste(a, "=", paste(deparse(args[[a]], width.cutoff = 500L), collapse = " ")), ""),
    collapse = ", ")
}

#' Remedies for what a model's checks found
#'
#' Lists a remedy for every check of a fitted model that is not `"OK"`, each
#' one written out as the change to [ilm_model()] that makes it, so that it
#' can be made with [ilm_apply_remedy()] rather than retyped. The checks
#' made at every fit are read from `object$checks`. The ones that simulate --
#' [ilm_check_dispersion()], [ilm_check_zeros()] and [ilm_check_variance()]
#' -- are passed in, since running them is a choice.
#'
#' Each remedy has a tier, which says what applying it would change:
#' \describe{
#'   \item{`numerical`}{The same model, fitted harder: more optimiser
#'     restarts.}
#'   \item{`structural`}{A different random-effect or variance structure, or
#'     a different way of estimating the variances, whose fixed effects mean
#'     what they meant before: a term whose variance is estimated at zero
#'     removed (the fit is the same without it), a covariance of lower rank,
#'     a dispersion model, a zero part, a negative binomial in place of a
#'     poisson, the boundary-avoiding penalty with its measured costs, or
#'     REML in place of maximum likelihood.}
#'   \item{`estimand`}{A change to what the fixed effects estimate, or to
#'     what their standard errors account for: a random effect whose variance
#'     is not zero removed, categories merged, levels pooled. Apply one only
#'     because the question calls for it, not because a check did.}
#' }
#' Remedies are listed numerical first, then structural, then estimand. Some
#' can only be made by hand -- which categories to merge is a question about
#' what they mean -- and those have no `change`. A remedy is a candidate, not
#' a cure: refit, and read the checks of the new fit.
#'
#' @section Keys:
#' Each remedy also has a `key`, `kind/target`, naming the change itself:
#' `"drop/g"` removes the grouping term `g`, `"family/nbinom"` refits as a
#' negative binomial, `"rank/g:2"` gives `g` a category covariance of rank 2.
#' The target is the term as the fit's checks name it. The `id` numbers the
#' rows of one table and changes when tables are combined with `c()` or a
#' standalone check is included; the key does not, so a script that gives
#' [ilm_apply_remedy()] a key makes the same remedy in another session, or
#' stops if the fit no longer calls for it.
#'
#' `ilm_remedies()` is a generic, shared with 'illumex', whose data checks
#' list remedies to data in the same table. Another package whose
#' diagnostics name remedies that are refits of an illume model gives its
#' own results a method, built with [ilm_remedy_table()], and
#' [ilm_apply_remedy()] applies them as it applies these.
#'
#' @param object A fitted `"ilm_model"` object.
#' @param dispersion,zeros,variance Optional results of
#'   [ilm_check_dispersion()], [ilm_check_zeros()] and
#'   [ilm_check_variance()] run on `object`, whose remedies are then listed
#'   too.
#' @param ... Unused.
#' @return A data frame of class `"ilm_remedies"`, one row per remedy, with
#'   `id`, `key`, the `check` (or checks) it answers and its `status`, the
#'   `tier`, the `remedy` in words, and the `change`: the [ilm_model()]
#'   arguments that make it, as code, or `""` when it is made by hand. An
#'   empty data frame when every check is OK.
#' @seealso [ilm_apply_remedy()] to refit with one, [summary.ilm_model()] for
#'   the checks themselves, [ilm_remedy_table()] for remedies from another
#'   package's diagnostics.
#' @examples
#' ## x and y vary within groups only, so the groups differ by nothing and
#' ## the random intercept's variance is estimated at zero
#' set.seed(1)
#' d <- data.frame(g = factor(rep(1:20, each = 10)), x = rnorm(200), e = rnorm(200))
#' d$x <- d$x - ave(d$x, d$g)
#' d$y <- 1 + 0.5 * d$x + d$e - ave(d$e, d$g)
#' f <- ilm_model(y ~ x + (1 | g), data = d, verbose = FALSE)
#' ilm_remedies(f)
#' @exportS3Method illumex::ilm_remedies
ilm_remedies.ilm_model <- function(object, dispersion = NULL, zeros = NULL,
                                   variance = NULL, ...) {
  rows <- list()
  add <- function(check, status, rs)
    for (r in rs) rows[[length(rows) + 1L]] <<- c(list(check = check, status = status), r)
  ck <- object$checks
  for (i in seq_len(nrow(ck)))
    add(ck$check[i], ck$status[i],
        ilm_rem_rules(object, ck$check[i], ck$status[i], ck$detail[i], ck$suggestion[i]))
  ## a standalone result is recognised by what only it carries
  std <- list(dispersion = dispersion, zeros = zeros, variance = variance)
  mark <- list(dispersion = "ratio", zeros = "expected", variance = "p_trend")
  for (nm in names(std)) {
    res <- std[[nm]]
    if (is.null(res)) next
    if (!is.list(res) || is.null(res$status) ||
        (!identical(res$status, "INCONCLUSIVE") && is.null(res[[mark[[nm]]]])))
      stop("`", nm, "` must be the result of ilm_check_", nm, "() on this model.",
           call. = FALSE)
    add(paste0("ilm_check_", nm), res$status, ilm_rem_standalone(object, nm, res))
  }
  ilm_rem_assemble(object, rows)
}

## The table from its rows -- each a list of check, status, tier, remedy,
## key and args -- for ilm_remedies() and ilm_remedy_table() alike, through
## illumex's assembly, so a remedy from another package is listed, merged,
## ordered and tied to its fit exactly as illume's own are.
#' @keywords internal
#' @noRd
ilm_rem_assemble <- function(object, rows) {
  rows <- lapply(rows, function(r)
    list(check = r$check, status = r$status, tier = r$tier, remedy = r$remedy,
         key = r$key, change = ilm_rem_change(r$args), payload = r$args))
  illumex::ilm_remedy_assemble(rows, target_id = ilm_rem_id(object),
                               tiers = ilm_rem_tiers,
                               tier_note = ilm_rem_tier_note,
                               by_hand = ilm_rem_by_hand)
}

## what a check can say when it names a remedy: an OK check names none
ilm_rem_statuses <- c("WARN", "FAIL", "BOUNDARY", "INCONCLUSIVE")

#' A table of remedies from another package's diagnostics
#'
#' For a package whose diagnostics find something wrong with an illume model
#' and name a remedy that is a refit of it -- a negative binomial for counts
#' more variable than a poisson allows, a dispersion or zero part -- so that
#' its users apply those remedies as they apply illume's own, with
#' [ilm_apply_remedy()]. Its `ilm_remedies()` method builds the table here:
#'
#' ```
#' ilm_remedies.my_calibration <- function(object, ...)
#'   ilm_remedy_table(object$fit, check = "pit_shape", status = "FAIL",
#'                    tier = "structural", remedy = "...",
#'                    args = list(list(family = "nbinom")))
#' ```
#'
#' The table is tied to `object` as illume's own are: listed, merged when
#' several checks name the same change, and ordered by tier the same way, and
#' refused by [ilm_apply_remedy()] for any other fit. The `change` column is
#' written from `args`, never passed in, so the change a person reads is the
#' refit that is made.
#'
#' A remedy that is not a refit of the model -- more simulations, another
#' kind of fold, a recalibration -- is not something [ilm_apply_remedy()] can
#' make. List it by hand, with `args = NULL` for its row, or leave it to the
#' package's own output.
#'
#' `ilm_remedy_table()` is a generic, shared with 'illumex'; this is its
#' method for a fitted model.
#'
#' @param object The fitted `"ilm_model"` the remedies would refit.
#' @param check Character: the name of the check each remedy answers. Not one
#'   of the fit's own checks, nor an `ilm_check_` name -- those are
#'   [ilm_remedies()]'s, and the report after a refit reads them from the fit.
#' @param status Character: what the check found, `"WARN"`, `"FAIL"`,
#'   `"BOUNDARY"` or `"INCONCLUSIVE"`.
#' @param tier Character: `"numerical"`, `"structural"` or `"estimand"`, as
#'   in [ilm_remedies()].
#' @param remedy Character: the remedy in words, a sentence a person can act
#'   on.
#' @param args A list with one element per row: a named list of
#'   [ilm_model()] arguments that make the remedy, or `NULL` for one made by
#'   hand. `NULL` makes every row by hand. The data and `verbose` are not
#'   arguments a remedy sets. A formula is evaluated where the model's own
#'   formula was written, so names in it mean what they mean there.
#' @param key Optional character, one per row: each remedy's key, `kind/target`
#'   as [illumex::ilm_remedy_key()] builds it. By default a change that says
#'   what it is by its arguments -- a family, a zero part, a dispersion
#'   formula -- takes the key illume gives the same change, so `c()` lists
#'   it once; any other refit takes `"change/"` and its change, and a remedy
#'   made by hand takes `"by_hand/"` and its check, numbered when one check
#'   has several.
#' @param ... Unused.
#' @return A data frame of class `"ilm_remedies"`, as from [ilm_remedies()].
#' @seealso [ilm_remedies()], [ilm_apply_remedy()].
#' @examples
#' set.seed(2)
#' d <- data.frame(x = rnorm(300))
#' d$y <- rnbinom(300, mu = exp(1 + 0.5 * d$x), size = 1.5)
#' f <- ilm_model(y ~ x, data = d, family = "poisson", verbose = FALSE)
#' ## what another package's diagnostic of this fit would hand back
#' rem <- ilm_remedy_table(f, check = "pit_shape", status = "FAIL",
#'   tier = "structural",
#'   remedy = "the intervals are too narrow: refit as a negative binomial",
#'   args = list(list(family = "nbinom")))
#' rem
#' f2 <- ilm_apply_remedy(f, rem, "family/nbinom")
#'
#' ## one list, illume's remedies and the other package's together: tables
#' ## for the same fit combine with c(), numbered afresh, a change several
#' ## lists name listed once
#' c(ilm_remedies(f), rem)
#' @exportS3Method illumex::ilm_remedy_table
ilm_remedy_table.ilm_model <- function(object, check, status, tier, remedy,
                                       args = NULL, key = NULL, ...) {
  txt <- list(check = check, status = status, tier = tier, remedy = remedy)
  for (nm in names(txt)) {
    v <- txt[[nm]]
    if (!is.character(v) || anyNA(v) || !all(nzchar(trimws(v))))
      stop("`", nm, "` must be character, with no missing or empty entries.",
           call. = FALSE)
  }
  n <- length(check)
  if (any(lengths(txt) != n))
    stop("`check`, `status`, `tier` and `remedy` must have one entry per ",
         "remedy, and have lengths ", paste(lengths(txt), collapse = ", "), ".",
         call. = FALSE)
  bad <- setdiff(status, ilm_rem_statuses)
  if (length(bad))
    stop("`status` must be one of ", paste(ilm_rem_statuses, collapse = ", "),
         " -- a check that is OK names no remedy -- not ",
         paste(unique(bad), collapse = ", "), ".", call. = FALSE)
  bad <- setdiff(tier, ilm_rem_tiers)
  if (length(bad))
    stop("`tier` must be one of ", paste(ilm_rem_tiers, collapse = ", "),
         ", not ", paste(unique(bad), collapse = ", "), ".", call. = FALSE)
  ## a merged row lists its checks with commas between them
  if (any(grepl(",", check, fixed = TRUE)))
    stop("`check` names must not contain a comma.", call. = FALSE)
  ## the report after a refit reads the fit's own checks from the new fit;
  ## another check of the same name would be reported as that one
  own <- unique(c(object$checks$check, ilm_rem_known))
  bad <- check[check %in% own | startsWith(check, "ilm_check_")]
  if (length(bad))
    stop("`check` names a check of illume's own (", paste(unique(bad), collapse = ", "),
         "): its remedies come from ilm_remedies(). Name the check that ",
         "found the problem.", call. = FALSE)
  if (!is.null(key) && (!is.character(key) || length(key) != n || anyNA(key) ||
                        !all(grepl("^[^/]+/.+", key))))
    stop("`key` must be character, one per remedy, each `kind/target` as ",
         "ilm_remedy_key() builds it -- or NULL to have them made.", call. = FALSE)

  if (is.null(args)) args <- vector("list", n)
  if (!is.list(args) || length(args) != n)
    stop("`args` must be a list with one element per remedy -- a named list ",
         "of ilm_model() arguments, or NULL for one made by hand.", call. = FALSE)
  ok <- setdiff(names(formals(ilm_model_formula)), c("data", "verbose"))
  env <- environment(object$formula)
  for (i in seq_len(n)) {
    a <- args[[i]]
    if (is.null(a)) next
    if (!is.list(a) || !length(a) || is.null(names(a)) || !all(nzchar(names(a))))
      stop("`args[[", i, "]]` must be a named list of ilm_model() arguments, ",
           "or NULL for a remedy made by hand: `args` holds one such list per ",
           "remedy, so a single remedy is args = list(list(family = ",
           "\"nbinom\")).", call. = FALSE)
    no <- setdiff(names(a), ok)
    if (length(no))
      stop("`args[[", i, "]]` sets ", paste(no, collapse = ", "), ", which ",
           if (length(no) == 1L) "is not an argument" else "are not arguments",
           " of ilm_model() that a remedy can set. A setting of another ",
           "package's own belongs in its output, or in a remedy made by hand ",
           "(args = NULL).", call. = FALSE)
    for (k in names(a))
      if (inherits(a[[k]], "formula") && !is.null(env)) environment(a[[k]]) <- env
    args[[i]] <- a
  }
  if (is.null(key)) {
    hand <- vapply(args, is.null, TRUE)
    nth <- ave(seq_len(n), check, hand, FUN = seq_along)
    many <- ave(seq_len(n), check, hand, FUN = length) > 1L
    key <- vapply(seq_len(n), function(i) {
      if (hand[i])
        return(illumex::ilm_remedy_key("by_hand", check[i], if (many[i]) nth[i]))
      ilm_rem_args_key(args[[i]]) %||%
        illumex::ilm_remedy_key("change", ilm_rem_change(args[[i]]))
    }, "")
  }
  rows <- lapply(seq_len(n), function(i)
    list(check = check[i], status = status[i], tier = tier[i],
         remedy = remedy[i], args = args[[i]], key = key[i]))
  ilm_rem_assemble(object, rows)
}

#' Refit a model with one of its remedies
#'
#' Makes one remedy from [ilm_remedies()]: refits the model through its own
#' call with only that change, so everything else about it -- the family, the
#' weights, a zero part, a dispersion model, REML -- comes along, and the new
#' fit's call is the one a person would have written, ready for [update()] or
#' a further remedy. Says what the checks that called for the remedy say now,
#' and which checks are not OK after the refit.
#'
#' The data are found where the model's formula was created, as `update()`
#' would; pass `data` when they are not there, e.g. when the model was fitted
#' inside a function from a formula written outside it.
#'
#' `ilm_apply_remedy()` is a generic, shared with 'illumex', which applies
#' remedies to data; this is its method for a fitted model.
#'
#' @param object A fitted `"ilm_model"` object, from the formula interface.
#' @param remedies The list from [ilm_remedies()] for `object` that the
#'   remedy was chosen from -- whether illume's own, or another package's
#'   built with [ilm_remedy_table()]. Required rather than recomputed, so the
#'   remedy made is always the one that was read.
#' @param which The remedy to make: its `id` in `remedies`, or its `key`.
#'   An id numbers the rows of this one list, and a list that includes the
#'   standalone checks numbers its remedies differently from one that does
#'   not; a key names the change itself, so a script that gives one makes
#'   the same remedy, or stops if the list no longer has it. Spaces in a key
#'   are ignored.
#' @param data The data the model was fitted to, when it cannot be found.
#' @param verbose Logical. Print the new fit's checks as it is fitted.
#' @param reason Optional: why this remedy is being made, in a sentence --
#'   "the zero check failed and the excess zeros are structural", "a
#'   reviewer asked for the simpler model". It is kept in the new fit's
#'   `remedy_log` beside the remedy, so the record of how the model was
#'   reached says why each change was made as well as what it was.
#' @param ... Unused.
#' @return The refitted model. Its `remedy_log` holds every remedy applied
#'   to reach it, in order: the key, check, status, tier, remedy and change,
#'   and the `reason` given for it (`NA` where none was).
#' @seealso [ilm_remedies()].
#' @examples
#' set.seed(1)
#' d <- data.frame(g = factor(rep(1:20, each = 10)), x = rnorm(200), e = rnorm(200))
#' d$x <- d$x - ave(d$x, d$g)
#' d$y <- 1 + 0.5 * d$x + d$e - ave(d$e, d$g)   # no group effect at all
#' f <- ilm_model(y ~ x + (1 | g), data = d, verbose = FALSE)
#' rem <- ilm_remedies(f)
#' rem
#' f2 <- ilm_apply_remedy(f, rem, "drop/g",
#'                        reason = "the groups were not expected to differ")
#' f2$remedy_log
#' @exportS3Method illumex::ilm_apply_remedy
ilm_apply_remedy.ilm_model <- function(object, remedies, which, data = NULL,
                                       verbose = FALSE, reason = NULL, ...) {
  if (!is.null(reason)) {
    if (!is.character(reason) || length(reason) != 1L || is.na(reason) ||
        !nzchar(trimws(reason)))
      stop("`reason` must be a single sentence saying why the remedy is made, ",
           "or NULL.", call. = FALSE)
    reason <- trimws(reason)
  }
  cl <- object$call
  if (is.null(cl) || !inherits(object$formula, "formula"))
    stop("a remedy is made by refitting through the formula interface, and ",
         "this model did not come from it. Fit it with ilm_model() and a ",
         "formula first.", call. = FALSE)
  if (!inherits(remedies, "ilm_remedies"))
    stop("`remedies` must come from ilm_remedies() or ilm_remedy_table().",
         call. = FALSE)
  if (is.null(attr(remedies, "target_id")) || is.null(attr(remedies, "payload")))
    stop("`remedies` has lost what ties it to its fit, as a subset of it ",
         "does. Pass the whole list from ilm_remedies() and choose with ",
         "`which`.", call. = FALSE)
  if (!identical(attr(remedies, "target_id"), ilm_rem_id(object)))
    stop("`remedies` was listed for a different fit. List them for this one ",
         "with ilm_remedies().", call. = FALSE)
  if (!nrow(remedies))
    stop("there is nothing to remedy: every check of this fit is OK.",
         call. = FALSE)
  i <- illumex::ilm_remedy_find(remedies, which)
  id <- remedies$id[i]
  args <- attr(remedies, "payload")[[as.character(id)]]
  if (is.null(args))
    stop("remedy ", id, " (", remedies$key[i], ") is made by hand, not by a ",
         "refit: ", remedies$remedy[i], ".", call. = FALSE)

  env0 <- environment(object$formula)
  if (is.null(env0)) env0 <- parent.frame()
  if (is.null(data)) {
    data <- tryCatch(eval(cl$data, env0), error = function(e) NULL)
    if (!is.data.frame(data))
      stop("the data this model was fitted to (`", deparse1(cl$data), "`) ",
           "cannot be found from where its formula was created. Pass them ",
           "as `data`.", call. = FALSE)
  } else if (!is.data.frame(data))
    stop("`data` must be a data frame.", call. = FALSE)

  new <- cl
  new[[1L]] <- as.name("ilm_model_formula")
  new$formula <- object$formula
  for (a in names(args)) new[a] <- list(args[[a]])
  new$data <- quote(.ilm_remedy_data)
  new$verbose <- verbose
  ## the formula's own environment, for every other name in the call, with
  ## the data and the (internal) fitting function added on top -- as
  ## ilm_model() evaluates, and for the same reason
  env <- new.env(parent = env0)
  assign(".ilm_remedy_data", data, envir = env)
  assign("ilm_model_formula", ilm_model_formula, envir = env)
  fit <- tryCatch(eval(new, env), error = function(e)
    stop("the refit with ", remedies$change[i], " failed: ",
         conditionMessage(e), call. = FALSE))
  ## the call as it would have been written, so that update() and the next
  ## remedy find the data where this one did
  fit$call$data <- cl$data
  fit$call$verbose <- cl$verbose
  ## a log from before reasons, or keys, were kept has no column for them
  old <- object$remedy_log
  if (!is.null(old) && is.null(old$reason)) old$reason <- NA_character_
  if (!is.null(old) && is.null(old$key)) old$key <- NA_character_
  if (!is.null(old)) old <- old[c("key", setdiff(names(old), "key"))]
  fit$remedy_log <- rbind(old, data.frame(
    key = remedies$key[i],
    check = remedies$check[i], status = remedies$status[i],
    tier = remedies$tier[i], remedy = remedies$remedy[i],
    change = remedies$change[i],
    reason = if (is.null(reason)) NA_character_ else reason,
    stringsAsFactors = FALSE))
  if (nrow(fit$X) != nrow(object$X))
    warning("the refit used ", nrow(fit$X), " rows against the original ",
            nrow(object$X), ": a variable the remedy brought in has missing ",
            "values, or `data` is not the data the model was fitted to. The ",
            "two fits cannot be compared.", call. = FALSE)

  ## what the checks that called for the remedy say now
  trig <- trimws(strsplit(remedies$check[i], ",", fixed = TRUE)[[1L]])
  st0 <- trimws(strsplit(remedies$status[i], ",", fixed = TRUE)[[1L]])
  ck <- fit$checks
  now <- vapply(seq_along(trig), function(k) {
    ## a check the fit did not make itself -- a standalone ilm_check_*(), or
    ## another package's -- is not in the new fit either; it has to be run
    if (startsWith(trig[k], "ilm_check_") || !trig[k] %in% object$checks$check)
      return(sprintf("  %s: %s before; run it again on the new fit", trig[k], st0[k]))
    a <- ck$status[ck$check == trig[k]]
    sprintf("  %s: %s -> %s", trig[k], st0[k],
            if (length(a)) a[1L] else "no longer checked, as the term is gone")
  }, "")
  bad <- ck[ck$status != "OK", , drop = FALSE]
  message("ilm_apply_remedy(): refitted with ", remedies$change[i], ".\n",
          if (!is.null(reason)) paste0("  reason given: ", reason, "\n"),
          paste(unique(now), collapse = "\n"), "\n",
          if (nrow(bad)) paste0("  not OK after the refit: ",
                                paste(sprintf("%s (%s)", bad$check, bad$status), collapse = ", "))
          else "  every check is OK after the refit")
  fit
}

## Rescaling, named first when a fit fails. The fit rescales its fixed-effect
## columns itself (ilm_rescale.R), but not a random slope's covariate, whose
## scale sits in the random-effect covariance: one far from unit scale (an SD
## outside 1e-3 to 1e3) is the likeliest cause of a failed optimisation, and
## the fix is by hand. None when every slope covariate is on a sane scale.
#' @keywords internal
#' @noRd
ilm_rem_rescale <- function(fit) {
  bars <- fit$bars
  if (is.null(bars) || !length(bars) || is.null(fit$model)) return(list())
  vs <- unique(unlist(lapply(bars, function(b)
    tryCatch(all.vars(b[[2]]), error = function(e) character()))))
  out <- list()
  for (v in intersect(vs, names(fit$model))) {
    x <- fit$model[[v]]
    if (!is.numeric(x)) next
    sdv <- stats::sd(x, na.rm = TRUE)
    if (!is.finite(sdv) || sdv <= 0 || (sdv >= 1e-3 && sdv <= 1e3)) next
    k <- 10^round(log10(sdv))
    out[[length(out) + 1L]] <- ilm_rem("numerical", sprintf(
      "rescale %s by hand -- its SD is %s; %s it by %s so it is near 1. The fit rescales its fixed-effect columns itself, but not the covariate of a random slope, whose scale sits in the random-effect covariance",
      v, formatC(sdv, format = "g", digits = 3), if (k > 1) "divide" else "multiply",
      formatC(if (k > 1) k else 1 / k, format = "g", digits = 3)),
      key = illumex::ilm_remedy_key("rescale", v))
  }
  out
}
