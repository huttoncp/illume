## ---------------------------------------------------------------------------
## The parts of a fitted model, for code built on it.
##
## Each piece is laid out BY NAME -- which block of the parameter vector, which
## level of which grouping factor, which cell of a correlation over time -- so
## code reading a fit never repeats illume's index arithmetic, which is how two
## copies of it come apart. The positions returned are positions in the fit's
## full parameter vector, object$obj$env$par: the order the joint precision
## holds, fixed and random interleaved as the template declares them.
## ---------------------------------------------------------------------------

#' Expectation of a function of a normal variable
#'
#' `E[h(eta + sd * Z)]` for `Z` standard normal, by the Gauss-Hermite
#' quadrature `predict(groups = "population")` uses for a population average
#' -- the same nodes and the same rule for how many -- so code built on a fit
#' averages over a latent spread exactly as the fit's own predictions do.
#'
#' @details
#' The number of nodes grows with the spread: 40 times the smallest whole
#' number at least `sd^2`, and no more than 800, with one rule for the whole
#' call, set by the largest `sd`. The rules are computed once and cached.
#' Measured against [stats::integrate()], through an inverse link that is
#' better than 1e-12 with a logit up to `sd = 6`, and better than 1e-9 with a
#' complementary log-log up to `sd = 5`. For any other integrand, check the
#' accuracy it needs; `n` sets the number of nodes directly.
#'
#' @param h A vectorised function of one argument.
#' @param eta,sd The mean and the standard deviation of the normal, recycled
#'   to a common length.
#' @param n Optional number of nodes, in place of the rule.
#' @return The expectation for each element: a vector shaped as `h(eta)`.
#' @seealso [predict.ilm_model()], whose `groups = "population"` uses the
#'   same rule.
#' @examples
#' ## a population-averaged probability through a logit
#' ilm_normal_expect(stats::plogis, eta = 0.5, sd = 1.2)
#' ## the probability of a count of 3 in a new group, under a log link
#' ilm_normal_expect(function(e) stats::dpois(3, exp(e)), eta = 1, sd = 0.6)
#' @export
ilm_normal_expect <- function(h, eta, sd, n = NULL) {
  if (!is.function(h))
    stop("`h` must be a function", call. = FALSE)
  m <- max(length(eta), length(sd))
  if (!m) return(h(numeric(0)))
  eta <- rep_len(as.numeric(eta), m); sd <- rep_len(as.numeric(sd), m)
  if (anyNA(sd) || any(sd < 0))
    stop("`sd` must be non-negative and not missing", call. = FALSE)
  if (is.null(n)) n <- ilm_gh_n(max(sd))
  if (!is.numeric(n) || length(n) != 1L || n < 1)
    stop("`n` must be a single positive number of nodes", call. = FALSE)
  gh <- ilm_gh(as.integer(n))
  out <- 0
  for (q in seq_along(gh$z)) out <- out + gh$w[q] * h(eta + sd * gh$z[q])
  out
}

#' Random effects of a fitted model
#'
#' The predicted random effects -- the conditional modes -- labelled by the
#' groups they belong to: one row per level of each grouping factor and each
#' of its coefficients, one per basis function of a smooth, and one per cell
#' of a correlation over time.
#'
#' @details
#' `sd` is each mode's conditional standard deviation with the variance
#' components at their estimates, from the inner Hessian of the Laplace
#' approximation: the quantity `lme4` calls `condVar`. With `reml = TRUE` the
#' fixed effects are integrated out alongside the random ones, so there it
#' also carries their uncertainty. A fit read back from disk no longer has the
#' compiled objective, and gives `NA`.
#'
#' `row` is each value's position in the fit's full parameter vector,
#' `object$obj$env$par`, which is also the order of its joint precision. Under
#' REML the fixed effects sit in that vector's random block too; they are not
#' random effects, and are not listed.
#'
#' For a reduced-rank term the rows are its latent factor scores; multiplied by
#' the loadings, `object$Lambda`, they give the effect on each category. A
#' random walk's first cell in each group is held at zero rather than
#' estimated, so it has no `row`, and a mode and `sd` of zero.
#'
#' With `nlme` or `lme4` attached, `ranef(fit)` gives the same.
#'
#' @param object A fitted `"ilm_model"`.
#' @param ... Unused.
#' @return A data frame of class `"ilm_ranef"`, with columns
#'   \describe{
#'     \item{`type`}{`"re"`, `"smooth"`, or the correlation over time:
#'       `"ar1"`, `"car1"`, `"rw1"`.}
#'     \item{`term`}{the term's name in the fit, as in `object$Sigma`.}
#'     \item{`factor`}{the grouping variable.}
#'     \item{`level`}{the group, a factor with the fitted levels; a smooth's
#'       basis functions are `b1`, `b2`, ...}
#'     \item{`dim`}{the coefficient -- `"(Intercept)"`, a slope's variable --
#'       prefixed by the category for a multinomial outcome.}
#'     \item{`time`, `cell`}{for a correlation over time, the cell's time and
#'       its row of [ilm_cells()]; `NA` otherwise.}
#'     \item{`row`}{the position in `object$obj$env$par`.}
#'     \item{`mode`, `sd`}{the conditional mode and its standard deviation.}
#'   }
#' @seealso [ilm_varcorr()] for the variance components, [ilm_cells()].
#' @examples
#' set.seed(1)
#' d <- data.frame(id = factor(rep(sprintf("s%02d", 1:10), each = 6)),
#'                 x = rnorm(60))
#' d$y <- 1 + 0.5 * d$x + rnorm(10)[d$id] + rnorm(60)
#' fit <- ilm_model(y ~ x + (1 | id), data = d, family = "gaussian",
#'                  verbose = FALSE)
#' head(ilm_ranef(fit))
#' @export
ilm_ranef <- function(object, ...) {
  if (!inherits(object, "ilm_model"))
    stop("`object` must be a fitted ilm_model, not ", class(object)[1],
         call. = FALSE)
  ilm_ranef_table(object, sd = TRUE)
}

## The table itself; `sd = FALSE` leaves out the conditional SDs, which cost a
## sparse inverse, for callers that want only the layout.
#' @keywords internal
#' @noRd
ilm_ranef_table <- function(object, sd = TRUE) {
  env <- object$obj$env
  pn <- names(env$par)
  csd <- if (sd) ilm_cond_sd(object) else rep(NA_real_, length(pn))
  pr <- object$sdr$par.random
  C <- object$C
  cats <- if (is.null(object$ylevels)) paste0("cat", seq_len(C))
          else object$ylevels[seq_len(C)]
  parts <- list()

  ## the random terms: bvec in term order, each an (levels x dims) by width
  ## block stored column-major, dimension-major within a column
  bv_rows <- which(pn == "bvec")
  bvec <- unname(pr[names(pr) == "bvec"])
  for (k in seq_along(object$re)) {
    e <- object$re[[k]]; nm <- names(object$re)[k]
    idx <- object$b_idx[[k]]
    nl <- object$nlk[k]; d <- object$dk[k]; w <- object$wk[k]
    basis <- identical(e$kind, "basis")
    lev <- if (basis) paste0("b", seq_len(nl))
           else if (length(e$levels) == nl) e$levels
           else as.character(seq_len(nl))
    dn <- if (basis) "(penalised)"
          else if (!is.null(colnames(e$Z))) colnames(e$Z)
          else if (d == 1L) "(Intercept)" else paste0("z", seq_len(d))
    rr <- identical(object$re_struct[[k]]$type, "rr")
    cn <- if (rr) paste0("factor", seq_len(w)) else cats
    li <- rep(rep(seq_len(nl), d), w)
    di <- rep(rep(seq_len(d), each = nl), w)
    ci <- rep(seq_len(w), each = nl * d)
    parts[[length(parts) + 1L]] <- list(
      type = if (basis) "smooth" else "re", term = nm,
      factor = if (basis) nm else if (!is.null(e$factor)) e$factor else nm,
      level = lev[li], lev_all = lev,
      dim = if (C > 1L || rr) paste0(cn[ci], ":", dn[di]) else dn[di],
      time = NULL, cell = NA_integer_,
      row = bv_rows[idx], mode = bvec[idx], sd = csd[bv_rows[idx]])
  }

  ## the cells of a correlation over time, anchors included
  ar <- object$ar
  if (!is.null(ar)) {
    cl <- ilm_cells(object)
    nlat <- attr(cl, "n_latent")
    ba_rows <- which(pn == "B_ar")
    Ba <- ilm_Bar_hat(object)
    nc <- nrow(cl)
    ci <- rep(seq_len(C), each = nc)
    ri <- rep(seq_len(nc), C)
    pos <- ifelse(is.na(cl$index[ri]), NA_integer_,
                  ba_rows[cl$index[ri] + (ci - 1L) * nlat])
    grp <- if (!is.null(ar$vars)) ar$vars[["group"]] else "group"
    parts[[length(parts) + 1L]] <- list(
      type = ar$type, term = "ar", factor = grp,
      level = as.character(cl$group[ri]), lev_all = levels(cl$group),
      dim = if (C > 1L) paste0(cats[ci], ":(Intercept)")
            else rep("(Intercept)", length(ri)),
      time = cl$time[ri], cell = ri,
      row = pos,
      mode = if (is.null(Ba)) rep(NA_real_, length(ri)) else Ba[cbind(ri, ci)],
      sd = if (sd) ifelse(is.na(pos), 0, csd[pos]) else rep(NA_real_, length(pos)))
  }

  n <- vapply(parts, function(p) length(p$mode), 1L)
  levs <- unique(unlist(lapply(parts, `[[`, "lev_all")))
  tm <- NULL
  if (!is.null(ar)) {
    ## a time column on the scale the times were given on, NA off the cells
    tm <- ilm_cor_as_time(rep(NA_real_, sum(n)), ar)
    at <- sum(n[-length(n)])
    tm[at + seq_len(n[length(n)])] <- parts[[length(parts)]]$time
  }
  out <- data.frame(
    type = rep(vapply(parts, `[[`, "", "type"), n),
    term = rep(vapply(parts, `[[`, "", "term"), n),
    factor = rep(vapply(parts, `[[`, "", "factor"), n),
    level = factor(unlist(lapply(parts, `[[`, "level")), levels = levs),
    dim = unlist(lapply(parts, `[[`, "dim")),
    stringsAsFactors = FALSE)
  out$time <- if (is.null(tm)) rep(NA, sum(n)) else tm
  out$cell <- unlist(lapply(seq_along(parts), function(i)
    rep_len(parts[[i]]$cell, n[i])))
  out$row <- unlist(lapply(parts, `[[`, "row"))
  out$mode <- unlist(lapply(parts, `[[`, "mode"))
  out$sd <- unlist(lapply(parts, `[[`, "sd"))
  rownames(out) <- NULL
  class(out) <- c("ilm_ranef", "data.frame")
  out
}

#' @rdname ilm_ranef
#' @exportS3Method nlme::ranef
ranef.ilm_model <- function(object, ...) ilm_ranef(object, ...)

#' @export
print.ilm_ranef <- function(x, digits = 4, ...) {
  cat("<ilm_ranef>", nrow(x), "predicted random effects (conditional modes)\n")
  tab <- table(factor(paste(x$type, x$term, sep = ": "),
                      levels = unique(paste(x$type, x$term, sep = ": "))))
  cat(paste0("  ", names(tab), ": ", as.integer(tab), collapse = "\n"), "\n\n")
  d <- as.data.frame(unclass(x), stringsAsFactors = FALSE)
  if (all(is.na(d$time))) d$time <- NULL
  if (all(is.na(d$cell))) d$cell <- NULL
  num <- vapply(d, is.numeric, TRUE) & !names(d) %in% c("row", "cell")
  d[num] <- lapply(d[num], signif, digits)
  print(utils::head(d, 10L), row.names = FALSE)
  if (nrow(d) > 10L) cat("  ...", nrow(d) - 10L, "more rows\n")
  invisible(x)
}

#' @param x An `"ilm_ranef"` object.
#' @param row.names,optional Unused; for the generic.
#' @param lme4 Logical. Give `lme4`'s long-format names -- `grpvar`, `term`,
#'   `grp`, `condval`, `condsd` -- in place of illume's own.
#' @rdname ilm_ranef
#' @export
as.data.frame.ilm_ranef <- function(x, row.names = NULL, optional = FALSE,
                                    lme4 = FALSE, ...) {
  d <- as.data.frame(unclass(x), stringsAsFactors = FALSE)
  if (!lme4) return(d)
  data.frame(grpvar = d$factor, term = d$dim, grp = d$level,
             condval = d$mode, condsd = d$sd, stringsAsFactors = FALSE)
}

## Conditional standard deviations of every random-block entry, from the inner
## Hessian at the optimum, laid out along object$obj$env$par (NA for fixed
## parameters). The Laplace approximation's own curvature, so it is lme4's
## condVar for a gaussian model; sdr$diag.cov.random is not, because it adds
## the uncertainty propagated from the fixed parameters.
#' @keywords internal
#' @noRd
ilm_cond_sd <- function(object) {
  env <- object$obj$env
  out <- rep(NA_real_, length(env$par))
  ri <- env$random
  if (!length(ri)) return(out)
  ## the full vector at the optimum, block by block BY NAME: under REML the
  ## stored opt$par has the fixed effects prepended, so it is not the vector
  ## TMB optimised, and filling by position put them where the variance
  ## parameters go
  full <- ilm_full_par(object)
  H <- tryCatch(env$spHess(full, random = TRUE), error = function(e) NULL)
  if (is.null(H)) return(out)
  v <- ilm_diag_inv(H)
  out[ri] <- sqrt(pmax(v, 0))
  out
}

## The fit's full parameter vector at the optimum, in object$obj$env$par order,
## each block taken by name: a random block from par.random, a fixed one from
## opt$par. The joint precision and ilm_ranef()'s rows are laid out along it.
#' @keywords internal
#' @noRd
ilm_full_par <- function(object) {
  full <- object$obj$env$par
  rn <- names(full)
  pf <- object$opt$par; pr <- object$sdr$par.random
  for (nm in unique(rn)) {
    src <- if (nm %in% names(pr)) pr[names(pr) == nm] else pf[names(pf) == nm]
    if (length(src) != sum(rn == nm))
      stop("internal: the fit's '", nm, "' block has ", length(src),
           " values for ", sum(rn == nm), " places", call. = FALSE)
    full[rn == nm] <- src
  }
  full
}

## The diagonal of a sparse symmetric matrix's inverse, a block of columns at
## a time through its sparse Cholesky factor, so the dense inverse of a large
## random block is never formed.
#' @keywords internal
#' @noRd
ilm_diag_inv <- function(H, chunk = 512L) {
  if (!inherits(H, "sparseMatrix")) H <- Matrix::Matrix(H, sparse = TRUE)
  H <- Matrix::forceSymmetric(H)
  n <- nrow(H)
  ch <- tryCatch(Matrix::Cholesky(H, LDL = FALSE, perm = TRUE),
                 error = function(e) NULL)
  if (is.null(ch)) return(rep(NA_real_, n))
  out <- numeric(n)
  for (s in seq(1L, n, by = chunk)) {
    j <- s:min(n, s + chunk - 1L)
    E <- Matrix::sparseMatrix(i = j, j = seq_along(j), x = 1,
                              dims = c(n, length(j)))
    X <- Matrix::solve(ch, E, system = "A")
    out[j] <- as.numeric(X[cbind(j, seq_along(j))])
  }
  out
}

#' Variance components of a fitted model
#'
#' Each random term's covariance, a correlation over time's parameters and
#' the family's dispersion, all on their natural scale: what `lme4::VarCorr()`
#' reports, for every structure illume fits.
#'
#' @details
#' **Random terms.** For one linear predictor, a term's covariance is that of
#' its coefficients: an intercept's variance, or an intercept and slopes'
#' variances and covariances. For a multinomial outcome it is the covariance
#' across categories, with the within-group structure beside it as the
#' `"Sigma_d"` attribute, whose `[1, 1]` the fit holds at one because the
#' category covariance carries the scale.
#'
#' **Over time.** An AR(1) term's `Sigma` is its stationary covariance and
#' `rho` the correlation one grid step apart; a CAR(1) term's `rho` is per unit
#' of time, and `range` its reciprocal decay; a random walk has
#' `var_per_time`, the variance of a step one time unit long. The words in
#' `meaning` are those of [ilm_cells()].
#'
#' **Dispersion.** On the natural scale and named for what it is. For a
#' gaussian model it is `sigma`, the residual standard DEVIATION, not the
#' variance. With a dispersion model it varies by row, and this is its value
#' at the median row.
#'
#' These are the variance parameters transformed exactly as a draw of them
#' would be -- through the one transform, `ilm_rebuild()` -- so that anything
#' summing over draws and this function cannot disagree at the estimate.
#'
#' With `nlme` or `lme4` attached, `VarCorr(fit)` gives the same.
#'
#' @param object,x A fitted `"ilm_model"`.
#' @param sigma Unused; for the generic.
#' @param ... Unused.
#' @return An object of class `"ilm_VarCorr"`: a list with `re` (one
#'   covariance matrix per random term, with `"stddev"` and `"correlation"`
#'   attributes as `lme4` gives them), `latent` (the correlation over time, or
#'   `NULL`), `dispersion` (`value`, `meaning`, `modelled`, or `NULL`) and
#'   `family`. `as.data.frame()` gives `lme4`'s columns: `grp`, `var1`,
#'   `var2`, `vcov`, `sdcor`.
#' @seealso [ilm_ranef()], [ilm_cells()].
#' @examples
#' set.seed(1)
#' d <- data.frame(id = factor(rep(1:12, each = 6)), t = rep(0:5, 12))
#' d$y <- 1 + rnorm(12)[d$id] + rnorm(12, 0, 0.3)[d$id] * d$t + rnorm(72)
#' fit <- ilm_model(y ~ t + (1 + t | id), data = d, family = "gaussian",
#'                  verbose = FALSE)
#' ilm_varcorr(fit)
#' as.data.frame(ilm_varcorr(fit))
#' @export
ilm_varcorr <- function(object) {
  if (!inherits(object, "ilm_model"))
    stop("`object` must be a fitted ilm_model, not ", class(object)[1],
         call. = FALSE)
  structure(ilm_natural(object), class = "ilm_VarCorr")
}

#' @rdname ilm_varcorr
#' @exportS3Method nlme::VarCorr
VarCorr.ilm_model <- function(x, sigma = 1, ...) ilm_varcorr(x)

## The variance parameters on their natural scale, from a fit or from a
## parameter vector in opt$par order. ONE transform, ilm_rebuild(), so the
## variance components at the estimate and those of any draw agree.
#' @keywords internal
#' @noRd
ilm_natural <- function(object, par = NULL) {
  o <- if (is.null(par)) object else ilm_rebuild(object, par)
  C <- o$C
  cats <- if (is.null(o$ylevels)) paste0("cat", seq_len(C))
          else o$ylevels[seq_len(C)]
  re <- lapply(seq_along(o$re), function(k) {
    e <- o$re[[k]]; nm <- names(o$re)[k]
    S <- as.matrix(o$Sigma[[nm]])
    if (C == 1L) {
      d <- o$dk[k]
      V <- S[1, 1] * (if (d > 1L) as.matrix(o$Sigma_d[[nm]]) else matrix(1, 1L, 1L))
      lab <- if (identical(e$kind, "basis")) "(penalised)"
             else if (!is.null(colnames(e$Z))) colnames(e$Z)
             else if (d == 1L) "(Intercept)" else paste0("z", seq_len(d))
    } else {
      V <- S; lab <- cats
      if (!is.null(o$Sigma_d[[nm]])) {
        Sd <- as.matrix(o$Sigma_d[[nm]])
        zl <- colnames(e$Z)
        if (!is.null(zl)) dimnames(Sd) <- list(zl, zl)
        attr(V, "Sigma_d") <- Sd
      }
    }
    dimnames(V) <- list(lab, lab)
    ilm_vc_attrs(V, factor = if (identical(e$kind, "basis")) nm
                             else if (!is.null(e$factor)) e$factor else nm,
                 kind = if (identical(e$kind, "basis")) "smooth" else "re")
  })
  names(re) <- names(o$re)
  latent <- NULL
  if (!is.null(o$ar)) {
    Sa <- as.matrix(o$Sigma[["ar"]])
    lab <- if (C > 1L) cats else "(Intercept)"
    dimnames(Sa) <- list(lab, lab)
    Sa <- ilm_vc_attrs(Sa, factor = if (!is.null(o$ar$vars)) o$ar$vars[["group"]]
                                    else "group", kind = o$ar$type)
    latent <- switch(o$ar$type,
      rw1 = list(type = "rw1", var_per_time = Sa),
      car1 = list(type = "car1", Sigma = Sa, rho = unname(o$rho),
                  range = unname(o$ar_range)),
      list(type = "ar1", Sigma = Sa, rho = unname(o$rho), step = o$ar$step))
    latent$meaning <- ilm_ar_words(o$ar)
  }
  disp <- NULL
  if (!is.null(o$dispersion)) {
    v <- unname(o$dispersion)
    names(v) <- sub("^log_", "", names(o$dispersion))
    disp <- list(value = v, meaning = ilm_disp_words(o$family$name),
                 modelled = !is.null(o$disp_formula) || isTRUE(o$disp_mu))
  }
  list(re = re, latent = latent, dispersion = disp,
       family = if (is.null(o$family)) "multinomial" else o$family$name)
}

## lme4's attributes on a covariance matrix: the standard deviations and the
## correlations, with a zero variance giving a correlation of zero rather than
## NaN.
#' @keywords internal
#' @noRd
ilm_vc_attrs <- function(V, factor, kind) {
  s <- sqrt(pmax(diag(V), 0))
  R <- V / outer(s, s)
  R[!is.finite(R)] <- 0
  diag(R) <- 1
  attr(V, "stddev") <- stats::setNames(s, rownames(V))
  attr(V, "correlation") <- R
  attr(V, "factor") <- factor
  attr(V, "kind") <- kind
  V
}

## What a family's dispersion parameter is, in words.
#' @keywords internal
#' @noRd
ilm_disp_words <- function(fam) {
  switch(fam,
    gaussian = "sigma, the residual standard deviation (not the variance)",
    nbinom = "k, the negative binomial size: Var(y) = mu + mu^2 / k",
    beta = "phi, the beta precision: Var(y) = mu (1 - mu) / (1 + phi)",
    weibull = , lognormal = , loglogistic =
      "the scale of the error on the log-time scale",
    "the family's dispersion parameter, on its natural scale")
}

#' @export
print.ilm_VarCorr <- function(x, digits = 4, ...) {
  d <- as.data.frame(x)
  fmt <- function(v) formatC(signif(v, digits), format = "fg", digits = digits)
  cat("Variance components (natural scale)\n")
  vv <- d[is.na(d$var2), , drop = FALSE]
  out <- data.frame(Groups = ifelse(duplicated(vv$grp), "", vv$grp),
                    Name = ifelse(is.na(vv$var1), "", vv$var1),
                    Variance = trimws(fmt(vv$vcov)),
                    Std.Dev. = trimws(fmt(vv$sdcor)), stringsAsFactors = FALSE)
  print(out, row.names = FALSE, right = FALSE)
  cc <- d[!is.na(d$var2), , drop = FALSE]
  if (nrow(cc)) {
    cat("\nCorrelations\n")
    print(data.frame(Groups = cc$grp, Between = paste(cc$var1, "and", cc$var2),
                     Corr = round(cc$sdcor, 3), stringsAsFactors = FALSE),
          row.names = FALSE, right = FALSE)
  }
  if (!is.null(x$latent)) {
    lt <- x$latent
    cat("\nOver time:", switch(lt$type,
      rw1 = "a random walk; var_per_time is the variance of a one-unit step",
      car1 = sprintf("CAR(1), rho = %s per unit of time (range %s)",
                     fmt(lt$rho), fmt(lt$range)),
      sprintf("AR(1), rho = %s per grid step of %s", fmt(lt$rho), fmt(lt$step))),
      "\n")
  }
  if (!is.null(x$dispersion)) {
    cat("\nDispersion:", paste(names(x$dispersion$value),
                               trimws(fmt(x$dispersion$value)), sep = " = ",
                               collapse = ", "),
        "--", x$dispersion$meaning,
        if (isTRUE(x$dispersion$modelled)) "(modelled; the value at the median row)" else "",
        "\n")
  }
  invisible(x)
}

#' @rdname ilm_varcorr
#' @param row.names,optional Unused; for the generic.
#' @export
as.data.frame.ilm_VarCorr <- function(x, row.names = NULL, optional = FALSE,
                                      ...) {
  rows <- list()
  add <- function(grp, V) {
    lab <- rownames(V); s <- attr(V, "stddev"); R <- attr(V, "correlation")
    for (i in seq_along(lab))
      rows[[length(rows) + 1L]] <<- data.frame(grp = grp, var1 = lab[i],
        var2 = NA_character_, vcov = V[i, i], sdcor = unname(s[i]),
        stringsAsFactors = FALSE)
    if (length(lab) > 1L)
      for (j in seq_along(lab)[-length(lab)]) for (i in (j + 1L):length(lab))
        rows[[length(rows) + 1L]] <<- data.frame(grp = grp, var1 = lab[j],
          var2 = lab[i], vcov = V[i, j], sdcor = R[i, j],
          stringsAsFactors = FALSE)
  }
  for (nm in names(x$re)) add(nm, x$re[[nm]])
  if (!is.null(x$latent))
    add("ar", if (identical(x$latent$type, "rw1")) x$latent$var_per_time
              else x$latent$Sigma)
  if (identical(x$family, "gaussian") && !is.null(x$dispersion) &&
      !isTRUE(x$dispersion$modelled)) {
    s <- unname(x$dispersion$value[1])
    rows[[length(rows) + 1L]] <- data.frame(grp = "Residual", var1 = NA_character_,
      var2 = NA_character_, vcov = s^2, sdcor = s, stringsAsFactors = FALSE)
  }
  if (!length(rows))
    return(data.frame(grp = character(0), var1 = character(0),
                      var2 = character(0), vcov = numeric(0),
                      sdcor = numeric(0)))
  do.call(rbind, rows)
}
