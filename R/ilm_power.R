## ---------------------------------------------------------------------------
## A-priori power, by simulation.
##
## Take a model as the account of how the data arise, generate datasets from
## it at a given size and effect, refit each one, and count how often the term
## of interest is detected. That works for any family this package fits --
## zero-inflated counts, ordinal and multinomial outcomes, mixed models --
## because nothing in it depends on a closed-form variance existing.
##
## Two things about a simulated power estimate are routinely dropped, and both
## change what the number means.
##
## IT IS ITSELF AN ESTIMATE. A power of 0.80 from 200 replicates has a standard
## error of 2.8%, so it is 0.74 to 0.86. Reading a required sample size off it
## to three figures is false precision, and reporting the interval is the
## difference between "you need 340 people" and "somewhere between 300 and
## 400, and here is how to narrow that".
##
## SOME FITS DO NOT CONVERGE. Dropping those and dividing by the rest gives
## power CONDITIONAL ON CONVERGENCE, which is not power: a study that fails to
## fit has not detected anything. Both numbers are reported, along with the
## failure rate, and the default counts a failure as a non-detection because
## that is what it is.
##
## And three things about the simulation itself, each of which was once wrong.
##
## THE STUDY. A fitted model's rows are resampled, whole clusters at a time,
## because the observed covariates are the best account of who will be
## recruited. A scaffold's design is DRAWN AGAIN for every replicate: allocation
## balanced the way the protocol balances it, and a covariate given as a
## function drawn fresh, because a new study recruits new people. Resampling a
## scaffold's grid instead left cells of a small factorial empty -- 18
## participants in a 3 x 3 design converged 28.5% of the time, where the
## planned design has two per cell and is always estimable.
##
## THE TEST. What is counted is the test the analysis will report: t on the
## residual degrees of freedom when nothing is integrated out, as summary()
## does, and the Wald z otherwise; a term with several coefficients -- a factor
## with three levels, or any term of a multinomial model -- jointly, as F or
## chi-square, the way ilm_anova() tests it. A z test in place of the t put the
## power of a 20-person two-arm trial 3.7 points too high.
##
## THE MODEL. Each replicate is refitted from the fitted model's own design --
## the rows the study drew and the response it produced -- through the refit
## machinery every other test in the package uses. Rebuilding
## ilm_model(formula, data, family) by hand dropped the contrasts, so a model
## fitted with sum-to-zero coding came back with differently named
## coefficients and a power of exactly zero; it also dropped the zero part, the
## dispersion model and the censoring, and failed outright on a transformed
## response or predictor, which the model frame holds only transformed.
##
## References:
##   Arnold, B. F., Hogan, D. R., Colford, J. M. and Hubbard, A. E. (2011).
##     Simulation methods to estimate design power. BMC Medical Research
##     Methodology 11, 94.
##   Green, P. and MacLeod, C. J. (2016). SIMR: an R package for power analysis
##     of generalised linear mixed models by simulation. Methods in Ecology and
##     Evolution 7, 493-498.
## ---------------------------------------------------------------------------

#' Resample a design to a given size, keeping clusters whole
#'
#' Rows are drawn with replacement. When the model has a grouping factor the
#' CLUSTERS are drawn instead, because a mixed design's power depends far more
#' on how many groups there are than on how many rows sit inside each, and
#' resampling rows would quietly hold the number of groups fixed. The rows
#' drawn are kept as the `"rows"` attribute.
#'
#' @keywords internal
#' @noRd
ilm_power_resample <- function(mf, n, group = NULL) {
  if (is.null(group) || !group %in% names(mf)) {
    i <- sample.int(nrow(mf), n, replace = TRUE)
    out <- mf[i, , drop = FALSE]
    rownames(out) <- NULL
    attr(out, "rows") <- i
    return(out)
  }
  g <- factor(mf[[group]])
  r <- ilm_power_rows(n, as.integer(g), seq_len(nrow(mf)))
  out <- mf[r$rows, , drop = FALSE]
  ## a cluster drawn twice has to become two clusters, or the design has
  ## fewer independent groups than it appears to
  out[[group]] <- factor(r$copy)
  rownames(out) <- NULL
  attr(out, "rows") <- r$rows
  out
}

## Which rows a simulated study of about `n` rows draws, from the positions in
## `base`, and which drawn copy of a cluster each belongs to. With no unit the
## rows are drawn singly and each is its own copy.
#' @keywords internal
#' @noRd
ilm_power_rows <- function(n, unit, base) {
  if (is.null(unit)) {
    i <- sample.int(length(base), n, replace = TRUE)
    return(list(rows = base[i], copy = seq_len(n)))
  }
  idx <- split(seq_along(base), unit[base])
  per <- length(base) / length(idx)
  ng <- max(2L, round(n / per))
  pick <- sample.int(length(idx), ng, replace = TRUE)
  i <- unlist(idx[pick], use.names = FALSE)
  list(rows = base[i], copy = rep(seq_len(ng), lengths(idx[pick])))
}

## The unit a simulated study recruits: the first grouping factor that is a
## population rather than a smooth, else the units an AR(1) chain runs within,
## else the rows themselves (NULL).
##
## The first grouping factor, and not simply the first random term: a smooth's
## penalised part sits in the same list, ahead of the bars, and taking it as
## the unit resampled a mixed model with a smooth row by row.
#' @keywords internal
#' @noRd
ilm_power_unit <- function(object) {
  gk <- which(vapply(object$re, function(e) !identical(e$kind, "basis"), TRUE))
  if (length(gk))
    return(list(codes = object$re[[gk[1L]]]$group, name = names(object$re)[gk[1L]]))
  if (!is.null(object$ar))
    return(list(codes = ilm_ar_group(object$ar), name = "the AR(1) units"))
  NULL
}

## Each row's unit in an AR(1) or CAR(1) specification.
#' @keywords internal
#' @noRd
ilm_ar_group <- function(ar) {
  if (identical(ar$type, "car1")) findInterval(ar$idx, ar$first)
  else (ar$idx - 1L) %/% ar$Tt + 1L
}

## Is every level of `g` inside a single level of `unit`? A grouping nested in
## the resampled unit -- classes within a school -- has to be relabelled with
## the copy, or two copies of one school would share their classes' random
## effects. One that crosses it -- items seen by every participant -- keeps its
## labels, because the items do not change when participants are redrawn.
#' @keywords internal
#' @noRd
ilm_nested_in <- function(g, unit) {
  if (is.null(unit)) return(FALSE)
  all(vapply(split(unit, g), function(v) length(unique(v)) == 1L, TRUE))
}

## A censoring specification restricted to some rows, attributes and all.
#' @keywords internal
#' @noRd
ilm_censor_rows <- function(spec, rows) {
  ct <- attr(spec, "ctime")
  structure(unclass(spec)[rows], lower = attr(spec, "lower"),
            upper = attr(spec, "upper"),
            ctime = if (is.null(ct)) NULL else ct[rows], class = class(spec))
}

## An AR(1) or CAR(1) specification for the rows a simulated study drew, with
## its units relabelled. The times come back from the specification itself --
## an AR(1) slot, or a CAR(1) time relative to the unit's first, rebuilt from
## the gaps it kept; only differences in time enter the correlation.
#' @keywords internal
#' @noRd
ilm_ar_rows <- function(ar, rows, group) {
  if (identical(ar$type, "car1")) {
    tc <- numeric(ar$n_cell)
    for (k in seq_along(ar$rest)) tc[ar$rest[k]] <- tc[ar$prev[k]] + ar$gap[k]
    ilm_car1(tc[ar$idx][rows], group, verbose = FALSE)
  } else {
    ilm_ar1(((ar$idx - 1L) %% ar$Tt + 1L)[rows], group, verbose = FALSE)
  }
}

## The fitted model's design restricted to the rows a study drew. Everything
## row-indexed moves together -- the fixed design, the dispersion and zero
## designs, the censoring, a smooth's basis, each random term's slopes and
## groups -- so the replicate is the same model on different rows, and no
## formula is evaluated again.
#' @keywords internal
#' @noRd
ilm_power_stub <- function(st, rows, copy, nested, ar_nested) {
  s <- st
  keep_attr <- function(M) {
    out <- M[rows, , drop = FALSE]
    attr(out, "formula") <- attr(M, "formula")
    out
  }
  s$X <- st$X[rows, , drop = FALSE]
  s$y <- st$y[rows]
  if (!is.null(st$weights)) s$weights <- st$weights[rows]
  if (!is.null(st$Zd))  s$Zd  <- keep_attr(st$Zd)
  if (!is.null(st$Zzi)) s$Zzi <- keep_attr(st$Zzi)
  if (!is.null(st$censor)) s$censor <- ilm_censor_rows(st$censor, rows)
  if (!is.null(st$rp)) s$rp$D <- st$rp$D[rows, , drop = FALSE]
  relabel <- function(g, nest) {
    gi <- g[rows]
    as.integer(factor(if (nest) paste(gi, copy) else gi))
  }
  for (k in seq_along(st$re)) {
    e <- st$re[[k]]
    if (identical(e$kind, "basis")) {
      e$basis <- e$basis[rows, , drop = FALSE]
    } else {
      e$group <- relabel(e$group, nested[k])
      e$Z <- e$Z[rows, , drop = FALSE]
      e$nl <- max(e$group)
    }
    s$re[[k]] <- e
  }
  if (!is.null(st$ar))
    s$ar <- ilm_ar_rows(st$ar, rows, relabel(ilm_ar_group(st$ar), ar_nested))
  s
}

## The fitted model's specification on a newly drawn design -- a scaffold's
## fresh draw of its plan. The fixed design is built from the draw with the
## fit's own terms, levels and contrasts, a smooth's basis is evaluated on
## it, and each random term's groups and slopes are read off it; everything
## else is the fit's. Refitting that is the same analysis as refitting through
## the formula, without the formula front end's cost on every replicate.
#' @keywords internal
#' @noRd
ilm_power_stub_grid <- function(object, st, d) {
  nd <- ilm_newX(object, d)
  s <- st
  s$X <- nd$X
  s$y <- numeric(nrow(d))
  s$weights <- NULL
  for (k in seq_along(st$re)) {
    e <- st$re[[k]]
    if (identical(e$kind, "basis")) {
      ## the same lookup ilm_eta() makes for a smooth's penalised block
      lab <- sub("\\.[0-9]+$", "", names(object$re)[k])
      sdl <- nd$smooths[[lab]]
      if (is.null(sdl)) sdl <- nd$smooths[[names(object$re)[k]]]
      i <- suppressWarnings(as.integer(sub(".*\\.", "", names(object$re)[k])))
      e$basis <- if (length(sdl$Xr) == 1L) sdl$Xr[[1L]] else
        sdl$Xr[[if (is.na(i)) 1L else i]]
    } else {
      nm <- names(object$re)[k]
      e$group <- as.integer(factor(d[[nm]]))
      e$Z <- ilm_re_design(object, k, d, e$d)
      e$nl <- max(e$group)
    }
    s$re[[k]] <- e
  }
  s
}

## Fresh random effects for the groups of a simulated study, and each row's
## share of its group's. One matrix-normal draw per group, U = A Z B, so a
## random slope keeps its own variance and a multinomial term its covariance
## across categories. The factorisation is shared with predict(marginal =
## TRUE) through ilm_re_factors(), so the two cannot drift; the difference is
## only WHOSE effect is drawn -- here one per simulated group, there one
## population draw applied to every row and averaged over.
#' @keywords internal
#' @noRd
ilm_power_re <- function(object, k, g, Z) {
  f <- ilm_re_factors(object, k)
  d <- nrow(f$A); C <- ncol(f$B); ng <- max(g)
  if (is.null(Z) || ncol(Z) != d)
    stop("the random-effect design for '", names(object$re)[k], "' could not ",
         "be rebuilt for the simulated study.", call. = FALSE)
  U <- ilm_re_draws(f$A, f$B, ng, nrow(f$B))
  out <- matrix(0, length(g), C)
  for (cc in seq_len(C)) {
    Uc <- matrix(vapply(U, function(u) u[, cc], numeric(d)), ncol = d,
                 byrow = TRUE)
    out[, cc] <- rowSums(Z * Uc[g, , drop = FALSE])
  }
  out
}

## The linear predictor of a resampled study: its fixed part at the effect
## being simulated, the fitted smooths -- mean structure, not a population to
## redraw -- and fresh random effects and serial correlation.
#' @keywords internal
#' @noRd
ilm_power_eta_stub <- function(object, s, B) {
  keep <- seq_len(nrow(B))
  ## a flexible baseline's columns are functions of the response; its part
  ## of the model enters through the simulator, not the linear predictor
  if (!is.null(object$rp)) keep <- keep[-object$rp$cols]
  eta <- s$X[, keep, drop = FALSE] %*% B[keep, , drop = FALSE]
  for (k in seq_along(object$re)) {
    e <- s$re[[k]]
    if (identical(e$kind, "basis")) {
      ctb <- e$basis %*% ilm_Bhat_term(object, k)
      if (identical(object$re_struct[[k]]$type, "rr"))
        ctb <- ctb %*% t(object$Lambda[[k]])
      eta <- eta + ctb
    } else {
      eta <- eta + ilm_power_re(object, k, e$group, e$Z)
    }
  }
  if (!is.null(s$ar))
    eta <- eta + ilm_ar_draw(s$ar, object$rho, object$Sigma[["ar"]], ncol(B))
  eta
}

## The linear predictor of a study laid out as a data frame -- a scaffold's
## freshly drawn design -- built through the formula the model was fitted
## with.
#' @keywords internal
#' @noRd
ilm_power_eta_grid <- function(object, d, B) {
  eta <- ilm_eta(object, ilm_newX(object, d), beta = B)
  for (k in seq_along(object$re)) {
    if (identical(object$re[[k]]$kind, "basis")) next
    nm <- names(object$re)[k]
    if (!nm %in% names(d))
      stop("the grouping factor '", nm, "' is not in the simulated design.",
           call. = FALSE)
    dk <- if (is.null(object$dk)) 1L else object$dk[k]
    eta <- eta + ilm_power_re(object, k, as.integer(factor(d[[nm]])),
                              ilm_re_design(object, k, d, dk))
  }
  eta
}

## A response drawn from a linear predictor, the way the model says it arises.
## Categories come back as integer codes, 1 to J.
#' @keywords internal
#' @noRd
ilm_power_response <- function(object, eta, logsig = NULL, Zzi = NULL,
                               w = NULL) {
  fam <- object$family; N <- nrow(eta)
  if (is.null(fam) || identical(fam$name, "multinomial")) {
    J <- object$J
    P <- exp(eta %*% t(stats::contr.sum(J))); P <- P / rowSums(P)
    cp <- P %*% (upper.tri(diag(J), diag = TRUE) * 1)   # cumulative, by row
    ## pmin: rounding can leave the last cumulative a hair under 1
    return(pmin(as.integer(rowSums(stats::runif(N) > cp)) + 1L, J))
  }
  if (isTRUE(object$ordinal)) {
    ## the latent variable, and the band it falls in: the model's own account
    ## of where a category comes from
    z <- as.numeric(eta[, 1L]) + fam$qfun(stats::runif(N))
    return(as.integer(rowSums(outer(z, as.numeric(object$zeta), `>`))) + 1L)
  }
  if (!is.null(object$rp)) return(ilm_rp_sim(object, off = as.vector(eta[, 1L])))
  dsp <- if (!is.null(object$dispersion)) log(unname(object$dispersion)) else
    numeric(0)
  y <- as.numeric(fam$sim(eta, if (is.null(w)) rep(1, N) else w, dsp,
                          logsig = logsig))
  if (!is.null(object$Zzi) && !is.null(Zzi)) {
    mu <- as.numeric(fam$linkinv(eta[, 1L]))
    dv <- if (!is.null(logsig)) exp(logsig) else
      if (length(dsp)) rep(exp(dsp[1L]), N) else rep(1, N)
    y <- ilm_zi_rng(object, mu, dv, ilm_zi_p(object, Zzi), y)
  }
  y
}

## Integer category codes as the factor a formula refit expects.
#' @keywords internal
#' @noRd
ilm_power_as_response <- function(object, y) {
  fam <- object$family
  if (is.null(fam) || identical(fam$name, "multinomial") ||
      isTRUE(object$ordinal))
    return(factor(object$ylevels[y], levels = object$ylevels,
                  ordered = isTRUE(object$ordinal)))
  y
}

## A response for the rows of `d`, at coefficients `beta`: one draw of the
## study `d` lays out. Used by ilm_scaffold() to give its grid a response
## that really is a sample from the parameters it claims.
#' @keywords internal
#' @noRd
ilm_power_draw <- function(object, d, beta) {
  B <- matrix(beta, ncol(object$X), object$C)
  ilm_power_as_response(object,
                        ilm_power_response(object, ilm_power_eta_grid(object, d, B)))
}

## What is being tested: a single coefficient, or every coefficient of a term.
##
## A variable or term name resolves through the design's `assign` map, not by
## matching name prefixes -- an exposure `x` would otherwise collect a
## covariate called `xray` -- and takes every column the term produced, in
## every category of a multinomial model.
#' @keywords internal
#' @noRd
ilm_power_target <- function(object, term) {
  b <- stats::coef(object); nm <- names(b)
  p <- ncol(object$X); C <- if (is.null(object$C)) 1L else object$C
  asg <- object$assign; tl <- object$term_labels
  cols <- function(j) if (is.null(asg)) integer(0) else
    which(!is.na(asg) & asg == j)
  pos <- function(k) as.vector(outer(k, (seq_len(C) - 1L) * p, "+"))
  real <- if (is.null(asg)) which(nm != "(Intercept)") else
    which(rep(!is.na(asg) & asg > 0L, C))
  if (!length(real))
    stop("the model has no term to find: it is an intercept only.",
         call. = FALSE)
  if (is.null(term)) {
    ## one coefficient per term in a single-equation model, as it always was;
    ## a multinomial model's first coefficient is a category's intercept, so
    ## there it is the first term, across every category
    if (C == 1L || is.null(asg)) return(list(label = nm[real[1L]], idx = real[1L]))
    j <- asg[real[1L]]
    return(list(label = tl[j], idx = pos(cols(j))))
  }
  if (length(term) != 1L)
    stop("`term` names one coefficient or one term.", call. = FALSE)
  if (term %in% nm) return(list(label = term, idx = match(term, nm)))
  j <- if (is.null(tl)) NA_integer_ else match(ilm_as_label(term, tl), tl)
  k <- if (is.na(j)) integer(0) else cols(j)
  if (!length(k))
    stop("`term` (", term, ") is not a coefficient or a term in this model. ",
         if (length(tl)) paste0("Terms: ", paste(tl, collapse = ", "), ". ") else "",
         "Coefficients: ", paste(nm[real], collapse = ", "), call. = FALSE)
  idx <- pos(k)
  list(label = if (length(idx) == 1L) nm[idx] else tl[j], idx = idx)
}

## Did a replicate's fit work? Judged on what fitting can get wrong -- the
## optimiser, the gradient, the Hessian -- and on the fixed effects being
## usable. A pre-fit check is a verdict on the DESIGN, the same for every
## study drawn from it: counting it against each replicate put the power of a
## three-category multinomial with four visits per participant at exactly
## zero, because two observations per latent value fails that check whatever
## the data. Those verdicts are reported with the result instead.
#' @keywords internal
#' @noRd
ilm_power_usable <- function(fit) {
  ck <- fit$checks
  np <- if (is.null(fit$n_precheck)) 0L else fit$n_precheck
  post <- ck$status[seq_len(nrow(ck)) > np]
  !any(post == "FAIL") && ilm_fixed_usable(fit)
}

## The pre-fit checks a replicate failed: their names, details and remedies.
#' @keywords internal
#' @noRd
ilm_power_design_fails <- function(fit) {
  ck <- fit$checks
  np <- if (is.null(fit$n_precheck)) 0L else fit$n_precheck
  ck[seq_len(nrow(ck)) <= np & ck$status == "FAIL",
     c("check", "detail", "suggestion"), drop = FALSE]
}

## The p-value the analysis would report for these coefficients: t or z for
## one, F or chi-square for several, exactly as summary() and ilm_anova() do.
#' @keywords internal
#' @noRd
ilm_power_p <- function(fit, nm) {
  cf <- stats::coef(fit); k <- match(nm, names(cf))
  if (anyNA(k)) return(NA_real_)
  V <- tryCatch(suppressWarnings(stats::vcov(fit)), error = function(e) NULL)
  if (is.null(V)) return(NA_real_)
  exact <- isTRUE(fit$exact_df) && isTRUE(is.finite(fit$resid_df))
  if (length(k) == 1L) {
    se <- sqrt(V[k, k])
    if (!is.finite(se) || se <= 0) return(NA_real_)
    z <- unname(cf[k] / se)
    return(if (exact) 2 * stats::pt(-abs(z), fit$resid_df) else
      2 * stats::pnorm(-abs(z)))
  }
  W <- ilm_wald_block(cf, V, k); q <- length(k)
  if (!is.finite(W) || W < 0) return(NA_real_)
  if (exact) stats::pf(W / q, q, fit$resid_df, lower.tail = FALSE) else
    stats::pchisq(W, q, lower.tail = FALSE)
}

#' Power for a term, by simulating from a fitted model
#'
#' Treats a fitted model as the truth, generates datasets of the requested
#' sizes, refits each, and counts how often the named term is detected. It
#' works for any family [ilm_model()] fits, because it never needs a
#' closed-form variance.
#'
#' @section The studies that are simulated:
#'
#' From a model fitted to data, each replicate resamples the fitted rows --
#' whole clusters at a time when the model has a grouping factor, so the number
#' of groups moves with `n` -- and draws a new response from the model: fresh
#' random effects for the new groups, fresh serial correlation, the fitted
#' smooths, the zero part, the dispersion model and the censoring. It is then
#' refitted as the same model, contrasts and all.
#'
#' From an [ilm_scaffold()], which has no data behind it, each replicate is a
#' fresh draw of the planned design instead: the allocation balanced as the
#' protocol would balance it, and any covariate given as a function drawn
#' again, because a new study recruits new people.
#'
#' @section The test that is counted:
#'
#' The one the analysis will report. A single coefficient gets `summary()`'s
#' test: t on the residual degrees of freedom when nothing is integrated out,
#' the Wald z otherwise. A term with several coefficients -- a factor with
#' three levels, or any term of a multinomial model, which has one coefficient
#' per category -- is tested jointly, F or chi-square, as [ilm_anova()] tests
#' it. Name one coefficient (`"b:x"`) to follow a single category instead.
#'
#' @section The estimate has a standard error:
#'
#' Power from `sims` replicates is a proportion, so it carries a standard error
#' of `sqrt(p (1 - p) / sims)` -- 2.8% at 0.80 from 200 replicates. The
#' returned interval is a Wilson interval on that, and the sample size implied
#' by a target power is a range rather than a number. Raising `sims` narrows
#' it, and nothing else does.
#'
#' @section Fits that do not converge:
#'
#' A replicate that fails to fit has not detected anything, so `power` counts
#' it as a non-detection. `power_converged` divides by the replicates that
#' worked, which is the number most software reports and is **power
#' conditional on convergence** -- a different and more flattering quantity.
#' When `converged` is below one the two differ, and the gap is the size of the
#' problem rather than something to smooth over.
#'
#' A check the analysis makes BEFORE fitting -- too few observations per random
#' effect, say -- is a verdict on the design, the same for every study drawn
#' from it, and is not counted against each replicate. The result reports it
#' instead, with its remedy, because every analysis of such a study will print
#' it as a failure and the plan should hear it first.
#'
#' @param object A fitted [ilm_model()], or an [ilm_scaffold()], to treat as
#'   the truth.
#' @param n Sample sizes to try. Defaults to a spread around the fitted size.
#'   For a model with a grouping factor this is the number of ROWS, and the
#'   number of clusters moves with it.
#' @param term What to test: one coefficient, named as the coefficient
#'   (`"armtreatment"`, or `"b:x"` in a multinomial model), or a whole term,
#'   named as the variable or term that produced it (`"arm"`, `"x"`). A term
#'   that produced a single coefficient is that coefficient; one that produced
#'   several is tested jointly. Defaults to the first coefficient that is not
#'   an intercept, and in a multinomial model to the first term.
#' @param effect For a single coefficient, values for it on the LINK scale;
#'   defaults to the fitted value. For a term tested jointly, MULTIPLES of its
#'   fitted coefficients -- `c(0.5, 1)` asks what happens if the effect is half
#'   what was assumed; defaults to 1. Several values trace power across effect
#'   sizes.
#' @param sims Replicates per cell.
#' @param alpha Two-sided level.
#' @param seed Random seed.
#' @param progress Show a progress bar; see [illumex::ilm_progress_arg].
#' @return An object of class `"ilm_power"`: one row per `n` by `effect` cell
#'   with `power`, its Monte Carlo interval, `power_converged` and `converged`.
#' @seealso [ilm_power_n()] to read off the size for a target power,
#'   [plot.ilm_power()] for the curve, [ilm_power_design()] for a study with no
#'   data yet, [ilm_simulate()] for the generator.
#' @references Arnold, B. F., Hogan, D. R., Colford, J. M. and Hubbard, A. E.
#'   (2011). Simulation methods to estimate design power. *BMC Medical Research
#'   Methodology* 11, 94.
#' @examples
#' set.seed(1); n <- 200
#' d <- data.frame(x = rnorm(n))
#' d$y <- rbinom(n, 1, plogis(-0.5 + 0.5 * d$x))
#' f <- ilm_model(y ~ x, data = d, family = "binomial", verbose = FALSE)
#' ilm_power(f, n = c(200, 400), sims = 50)
#' @export
ilm_power <- function(object, n = NULL, term = NULL, effect = NULL,
                      sims = 200L, alpha = 0.05, seed = 1L, progress = NULL) {
  if (!inherits(object, "ilm_model"))
    stop("`object` must be a fitted ilm_model, not ", class(object)[1],
         call. = FALSE)
  if (!is.null(object$design))
    stop("this model was fitted to a complex sample, and ilm_power() draws a ",
         "study as a simple random sample of rows or of clusters -- not the ",
         "stratified, weighted design the data came from, so the power it ",
         "reported would be for a different survey. Size a survey from its ",
         "design effect: the ratio of the ilm_svy_coef() variances to the ",
         "model-based ones, applied to ilm_power() on the unweighted model.",
         call. = FALSE)
  scaffold <- inherits(object, "ilm_scaffold") && !is.null(object$scaffold)
  tgt <- ilm_power_target(object, term)
  b <- stats::coef(object)
  nm_t <- names(b)[tgt$idx]
  joint <- length(tgt$idx) > 1L
  n0 <- nrow(object$X)
  if (is.null(n)) n <- unique(round(n0 * c(0.5, 1, 2)))
  n <- sort(unique(as.integer(n)))
  if (any(n < 10L)) stop("`n` below 10 is not a study", call. = FALSE)
  if (is.null(effect)) effect <- if (joint) 1 else unname(b[tgt$idx])
  grid <- expand.grid(n = n, effect = effect, KEEP.OUT.ATTRS = FALSE)
  fam <- if (is.null(object$family)) "multinomial" else object$family$name
  p <- ncol(object$X); C <- if (is.null(object$C)) 1L else object$C

  if (scaffold) {
    s <- object$scaffold
    rpu <- max(1L, as.integer(s$rows_per_unit))
    st <- ilm_refit_stub(object)
    one <- function(nn, B) {
      d <- ilm_scaffold_grid(s$design, max(4L, as.integer(round(nn / rpu))),
                             s$within, s$group)
      y <- ilm_power_response(object, ilm_power_eta_grid(object, d, B))
      ilm_refit_like(ilm_power_stub_grid(object, st, d), y = y, restarts = 1L)
    }
    grouped <- !is.null(s$group)
  } else {
    st <- ilm_refit_stub(object)
    unit <- ilm_power_unit(object)
    ucodes <- if (is.null(unit)) NULL else unit$codes
    nested <- vapply(st$re, function(e) !identical(e$kind, "basis") &&
                       ilm_nested_in(e$group, ucodes), TRUE)
    ar_nested <- !is.null(st$ar) && ilm_nested_in(ilm_ar_group(st$ar), ucodes)
    ## Frequency weights say a row stands for that many observations, and a
    ## new study draws each of them separately; simulating one outcome and
    ## weighting it would claim they all came out the same. Binomial weights
    ## on a proportion are numbers of trials instead, and travel with the row.
    w <- st$weights
    trials <- !is.null(w) && identical(fam, "binomial") &&
      any(!(st$y %in% c(0, 1)))
    base <- seq_len(n0)
    if (!is.null(w) && !trials) {
      base <- rep(base, pmax(round(w), 0L))
      st$weights <- NULL
    }
    dispv <- if (!is.null(object$Zd) || isTRUE(object$disp_mu))
      log(ilm_disp_vec(object)) else NULL
    one <- function(nn, B) {
      r <- ilm_power_rows(nn, ucodes, base)
      sr <- ilm_power_stub(st, r$rows, r$copy, nested, ar_nested)
      y <- ilm_power_response(object, ilm_power_eta_stub(object, sr, B),
                              logsig = if (is.null(dispv)) NULL else dispv[r$rows],
                              Zzi = sr$Zzi, w = sr$weights)
      if (!is.null(sr$censor)) y <- ilm_censor_apply(sr$censor, y)
      ilm_refit_like(sr, y = y, restarts = 1L)
    }
    grouped <- !is.null(unit)
  }

  set.seed(seed)
  pb <- ilm_progress(nrow(grid) * sims, progress, "simulating studies")
  tick <- 0L
  res <- vector("list", nrow(grid))
  dfail <- list()                       # design verdicts: check -> detail, n
  nfit <- 0L
  for (g in seq_len(nrow(grid))) {
    bb <- b
    bb[tgt$idx] <- if (joint) b[tgt$idx] * grid$effect[g] else grid$effect[g]
    B <- matrix(bb, p, C)
    rej <- ok <- logical(sims)
    for (s_i in seq_len(sims)) {
      fit <- try(suppressWarnings(suppressMessages(one(grid$n[g], B))),
                 silent = TRUE)
      tick <- tick + 1L; pb$tick(tick)
      if (inherits(fit, "try-error")) next
      nfit <- nfit + 1L
      fl <- ilm_power_design_fails(fit)
      for (f_i in seq_len(nrow(fl))) {
        k <- fl$check[f_i]
        dfail[[k]] <- if (is.null(dfail[[k]]))
          list(detail = fl$detail[f_i], suggestion = fl$suggestion[f_i], n = 1L)
        else utils::modifyList(dfail[[k]], list(n = dfail[[k]]$n + 1L))
      }
      if (!ilm_power_usable(fit)) next
      pv <- ilm_power_p(fit, nm_t)
      if (!is.finite(pv)) next
      ok[s_i] <- TRUE
      rej[s_i] <- pv < alpha
    }
    nok <- sum(ok)
    ci <- ilm_wilson(sum(rej), sims, alpha = 0.05)
    res[[g]] <- data.frame(
      n = grid$n[g], effect = grid$effect[g],
      power = mean(rej),                 # failures count as non-detections
      mc_lower = ci[1L], mc_upper = ci[2L],
      power_converged = if (nok) sum(rej) / nok else NA_real_,
      converged = nok / sims, sims = sims, row.names = NULL)
  }
  pb$done()
  exact <- isTRUE(object$exact_df)
  test <- if (joint) (if (exact) "F" else "Wald chi-square") else
    (if (exact) "t" else "Wald z")
  out <- do.call(rbind, res)
  df_tab <- if (length(dfail))
    data.frame(check = names(dfail),
               detail = vapply(dfail, `[[`, "", "detail"),
               suggestion = vapply(dfail, `[[`, "", "suggestion"),
               share = vapply(dfail, `[[`, 0L, "n") / max(nfit, 1L),
               row.names = NULL, stringsAsFactors = FALSE) else NULL
  structure(out, class = c("ilm_power", "data.frame"), term = tgt$label,
            coefs = nm_t, test = test, df = length(nm_t),
            effect_scale = if (joint) "multiple" else "link",
            alpha = alpha, family = fam, formula = object$formula,
            grouped = grouped, redrawn = scaffold, fitted_n = n0,
            design_fails = df_tab)
}

#' Wilson interval for a proportion
#'
#' Rather than the Wald one, which at a power of 0.98 from 200 replicates puts
#' the upper end above 1 and is worst exactly where power estimates live.
#'
#' @keywords internal
#' @noRd
ilm_wilson <- function(x, n, alpha = 0.05) {
  z <- stats::qnorm(1 - alpha / 2)
  p <- x / n
  d <- 1 + z^2 / n
  ctr <- (p + z^2 / (2 * n)) / d
  hw <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / d
  c(max(ctr - hw, 0), min(ctr + hw, 1))
}

#' The sample size a target power implies
#'
#' Interpolates the simulated curve, and reports the range the Monte Carlo
#' error allows rather than a single number the simulation cannot support.
#'
#' @param object An [ilm_power()] result.
#' @param target Power to reach.
#' @return A data frame with one row per effect size: the interpolated `n`, and
#'   the `n_lower`/`n_upper` implied by the Monte Carlo interval. For a result
#'   from [ilm_power_design()] it also carries `n_unit`, the same figure in
#'   participants rather than rows.
#' @seealso [ilm_power()].
#' @export
ilm_power_n <- function(object, target = 0.8) {
  if (!inherits(object, "ilm_power"))
    stop("`object` must be an ilm_power() result, not ", class(object)[1],
         call. = FALSE)
  d <- as.data.frame(object); class(d) <- "data.frame"
  if (length(unique(d$n)) < 2L)
    stop("at least two sample sizes are needed to interpolate one; rerun ",
         "ilm_power() with several values of `n`.", call. = FALSE)
  cross <- function(x, y) {
    o <- order(x); x <- x[o]; y <- y[o]
    if (all(y < target) || all(y > target)) return(NA_real_)
    stats::approx(y, x, xout = target, ties = "ordered")$y
  }
  out <- do.call(rbind, lapply(split(d, d$effect), function(z)
    data.frame(effect = z$effect[1L], n = cross(z$n, z$power),
               n_lower = cross(z$n, z$mc_upper),
               n_upper = cross(z$n, z$mc_lower),
               row.names = NULL)))
  ## An ilm_power_design() curve is reported in PARTICIPANTS and this is in
  ## rows, so the two numbers differ by the rows per participant. Returning
  ## only the rows invites "we need 330 people" for a study of 165.
  rpu <- attr(object, "rows_per_unit")
  if (!is.null(rpu) && rpu > 1L) {
    out$n_unit <- out$n / rpu
    out$n_unit_lower <- out$n_lower / rpu
    out$n_unit_upper <- out$n_upper / rpu
  }
  out
}

#' @export
print.ilm_power <- function(x, ...) {
  test <- attr(x, "test"); if (is.null(test)) test <- "Wald z"
  cf <- attr(x, "coefs")
  joint <- identical(attr(x, "effect_scale"), "multiple")
  cat(sprintf("Simulated power for %s (%s family, alpha = %.3g)\n",
              attr(x, "term"), attr(x, "family"), attr(x, "alpha")))
  cat(if (joint)
        sprintf("  counting the joint %s test of its %d coefficients, as the analysis reports it\n",
                test, length(cf))
      else sprintf("  counting the %s test, as the analysis reports it\n", test))
  cat(sprintf("  %d replicates per cell; %s had %d rows\n", x$sims[1L],
              if (isTRUE(attr(x, "redrawn"))) "the scaffold's own grid"
              else "the fitted study", attr(x, "fitted_n")))
  if (isTRUE(attr(x, "redrawn")))
    cat("  each replicate is a fresh draw of the planned design: balanced\n",
        "  allocation, new participants\n", sep = "")
  else if (isTRUE(attr(x, "grouped")))
    cat("  clusters are resampled whole, so the number of GROUPS moves with n\n")
  if (joint)
    cat("  `effect` is a MULTIPLE of the assumed coefficients: 1 is the model\n",
        "  as given, 0.5 half of every one of them\n", sep = "")
  d <- as.data.frame(x); class(d) <- "data.frame"
  d$sims <- NULL
  for (j in c("power", "mc_lower", "mc_upper", "power_converged", "converged"))
    d[[j]] <- round(d[[j]], 3)
  d$effect <- signif(d$effect, 4)
  print(d, row.names = FALSE)
  cat("\n  mc_lower/mc_upper is a Wilson interval on the power ESTIMATE: at\n",
      "  0.80 from ", x$sims[1L], " replicates the standard error is ",
      sprintf("%.3f", sqrt(0.8 * 0.2 / x$sims[1L])),
      ", so a\n  sample size read off this curve is a range. More replicates ",
      "narrow it;\n  nothing else does.\n", sep = "")
  if (any(x$converged < 1, na.rm = TRUE))
    cat("\n  Some replicates did not converge. `power` counts those as ",
        "non-detections,\n  because a study that will not fit has not detected ",
        "anything.\n  `power_converged` divides by the ones that worked, which ",
        "is power\n  CONDITIONAL ON CONVERGENCE and is the more flattering ",
        "number.\n", sep = "")
  ## a verdict on the design itself: every analysis of such a study will
  ## print it, so the plan should hear it first
  dft <- attr(x, "design_fails")
  if (!is.null(dft) && nrow(dft)) {
    cat("\n  The DESIGN fails a check the analysis makes before fitting, and\n",
        "  every analysis of such a study will report it as FAILED:\n", sep = "")
    for (i in seq_len(nrow(dft))) {
      cat(sprintf("    %s: %s  (%.0f%% of the simulated studies)\n",
                  dft$check[i], dft$detail[i], 100 * dft$share[i]))
      if (nzchar(dft$suggestion[i]))
        cat("      try: ", dft$suggestion[i], "\n", sep = "")
    }
    cat("  The power above counts the fits that converged regardless, so\n",
        "  read it with that verdict attached.\n", sep = "")
  }
  invisible(x)
}

#' Power curve
#'
#' @param x An [ilm_power()] result.
#' @param target Draw a line at this power.
#' @param ... Passed to the plotting function.
#' @return `x`, invisibly.
#' @seealso [ilm_power()].
#' @export
plot.ilm_power <- function(x, target = 0.8, ...) {
  d <- as.data.frame(x); class(d) <- "data.frame"
  eff <- unique(d$effect)
  op <- graphics::par(no.readonly = TRUE); on.exit(graphics::par(op))
  graphics::plot(range(d$n), c(0, 1), type = "n", xlab = "sample size",
                 ylab = "power",
                 main = sprintf("power for %s", attr(x, "term")), ...)
  graphics::abline(h = target, lty = 2, col = "grey50")
  for (i in seq_along(eff)) {
    z <- d[d$effect == eff[i], ]
    z <- z[order(z$n), ]
    ## the Monte Carlo band, because the curve is itself estimated
    graphics::polygon(c(z$n, rev(z$n)), c(z$mc_lower, rev(z$mc_upper)),
                      col = grDevices::adjustcolor(i + 1, alpha.f = 0.18),
                      border = NA)
    graphics::lines(z$n, z$power, col = i + 1, lwd = 2)
    graphics::points(z$n, z$power, col = i + 1, pch = 19)
  }
  if (length(eff) > 1L)
    graphics::legend("bottomright",
                     legend = if (identical(attr(x, "effect_scale"), "multiple"))
                       sprintf("effect x %.3g", eff) else sprintf("effect %.3g", eff),
                     col = seq_along(eff) + 1, lwd = 2, bty = "n")
  invisible(x)
}
