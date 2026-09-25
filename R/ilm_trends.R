## ---------------------------------------------------------------------------
## Simple slopes: taking an interaction apart when one side of it is continuous.
##
## An omnibus test says a treatment-by-time interaction exists. It does not say
## what is driving it, and the two questions people actually have are different
## from each other:
##
##   1. Is the slope on time different BETWEEN groups?   -- the interaction
##      restated one pair at a time
##   2. Is the slope on time different from zero WITHIN each group? -- whether
##      anything is happening in that arm at all
##
## They can disagree in both directions. Two arms can have slopes that differ
## significantly while neither differs from zero, and both can be strongly
## non-zero while not differing from each other. Reporting only one of them is
## how "the treatment group improved and the control group did not" gets
## written when what was tested never compared the two.
##
## A marginal MEAN cannot answer either. Averaging the fitted values at the
## mean of the covariate collapses exactly the thing being asked about, which
## is why emmeans needs emtrends alongside it.
##
## THE CONSTRUCTION. A marginal mean is L %*% beta where the rows of L are
## averaged rows of the grid's model matrix. A marginal SLOPE is the same thing
## with a different L: the derivative of those rows with respect to the
## covariate. Everything downstream -- the covariance L V L', the simultaneous
## adjustment, the contrasts between rows -- is then identical, which is why
## ilm_trends() returns the same class ilm_emmeans() does and ilm_contrast()
## works on it without knowing the difference.
##
## The derivative is taken numerically, by a central difference, for the same
## reason emtrends does: it is exact for any term linear in the covariate, and
## it keeps working for polynomials, splines and interactions without anyone
## having to write down their derivatives.
## ---------------------------------------------------------------------------

#' Estimated marginal slopes
#'
#' The slope of a continuous predictor within each level of a factor, with the
#' other predictors averaged over. This is what takes a significant
#' interaction between a factor and a covariate apart.
#'
#' @section The two questions, which are not the same question:
#' A treatment-by-time interaction says the slope on time depends on the arm.
#' Two different follow-ups are then available, and they answer different
#' things:
#'
#' * **Do the slopes differ between arms?** Pass the result to
#'   [ilm_contrast()], which differences the rows. This is the interaction
#'   restated one pair at a time.
#' * **Is the slope different from zero within an arm?** Read the
#'   `statistic` and `p.value` columns of this table, which test each row
#'   against zero.
#'
#' They can disagree in both directions: two arms can have slopes that differ
#' significantly while neither is distinguishable from zero, and both can be
#' strongly non-zero while not differing from one another. Reporting one and
#' describing it as the other is the common way a within-arm result gets
#' written up as a between-arm claim.
#'
#' @section Why a marginal mean cannot do this:
#' [ilm_emmeans()] averages predictions with the covariate held at its mean,
#' which collapses the very thing the interaction is about. The slope has to be
#' estimated as a slope. The two share all their machinery -- the same
#' reference grid, the same weighting, the same covariance -- and differ only
#' in what fills the contrast matrix.
#'
#' @param object A fitted `"ilm_model"` object.
#' @param specs Factor(s) to compute the slope within, as a character vector.
#' @param var The continuous predictor whose slope is wanted.
#' @param at Named list fixing other predictors, as in [ilm_emmeans()]. The
#'   slope of `var` is constant across these unless `var` enters non-linearly
#'   or interacts with them, in which case where you evaluate it matters and
#'   the result says so.
#' @param weights How to average over factors not in `specs`; see
#'   [ilm_emmeans()].
#' @param delta Step for the numerical derivative. The default is 1/1000 of
#'   the observed range of `var`, which is exact for any linear term and small
#'   enough elsewhere.
#' @param level Confidence level.
#' @param df Denominator degrees of freedom for the within-level tests: one of
#'   `"auto"`, `"satterthwaite"`, `"kenward-roger"`, `"residual"`,
#'   `"asymptotic"`, or a single number. See [ilm_denom_df()]. `"auto"` gives
#'   an exact t where nothing was integrated out, Satterthwaite for a gaussian
#'   mixed model, and a z test otherwise.
#' @return A data frame of class `"ilm_emm"`, one row per level of `specs`,
#'   with `estimate`, `se`, `df`, `statistic`, `p.value`, `lower` and `upper`.
#'   Being an `"ilm_emm"` it can be passed straight to [ilm_contrast()].
#' @seealso [ilm_emmeans()] for means, [ilm_contrast()] for differences
#'   between the slopes, [ilm_anova()] for the omnibus interaction this takes
#'   apart.
#' @examples
#' set.seed(1)
#' d <- data.frame(arm = factor(rep(c("ctl", "trt"), each = 60)),
#'                 week = rep(0:5, 20))
#' d$y <- 2 + 0.1 * d$week + 0.4 * d$week * (d$arm == "trt") + rnorm(120)
#' f <- ilm_model(y ~ arm * week, data = d, family = "gaussian",
#'                verbose = FALSE)
#' tr <- ilm_trends(f, "arm", var = "week")
#' tr                      # is each arm's slope different from zero?
#' ilm_contrast(tr)        # do the two arms' slopes differ from each other?
#' @export
ilm_trends <- function(object, specs, var, at = NULL,
                       weights = c("equal", "proportional", "cells"),
                       delta = NULL, level = 0.95, df = "auto") {
  weights <- match.arg(weights)
  if (!inherits(object, "ilm_model"))
    stop("`object` must be a fitted ilm_model, not ", class(object)[1],
         call. = FALSE)
  mf <- ilm_data(object)
  if (is.null(mf)) stop("the fit did not keep its model frame", call. = FALSE)
  var <- as.character(var)
  if (length(var) != 1L)
    stop("`var` must name a single continuous predictor, not ", length(var),
         call. = FALSE)
  if (!var %in% names(mf))
    stop("`var` (", var, ") is not a predictor in the model. The model has: ",
         paste(setdiff(names(mf), names(mf)[1]), collapse = ", "),
         call. = FALSE)
  if (!is.numeric(mf[[var]]))
    stop("`var` (", var, ") is a ", class(mf[[var]])[1], ", and a slope is ",
         "only defined for a continuous predictor. For a factor the ",
         "comparison you want is between its levels: use ilm_emmeans() and ",
         "ilm_contrast().", call. = FALSE)
  specs <- as.character(specs)
  miss <- setdiff(specs, names(mf))
  if (length(miss))
    stop("variable(s) not in the model: ", paste(miss, collapse = ", "),
         call. = FALSE)
  if (var %in% specs)
    stop("`var` (", var, ") is also in `specs`. The slope is computed WITHIN ",
         "the levels of `specs`, so the covariate cannot be one of them.",
         call. = FALSE)
  if (names(mf)[1L] %in% c(specs, var))
    stop("`", names(mf)[1L], "` is the response; slopes are taken along a ",
         "predictor, within the levels of others.", call. = FALSE)

  ## the grid, and the same grid nudged along `var` in both directions
  g <- ilm_ref_grid(object, at)
  rng <- range(mf[[var]], na.rm = TRUE)
  h <- if (!is.null(delta)) as.numeric(delta)[1] else {
    sp <- diff(rng)
    if (!is.finite(sp) || sp <= 0) 1e-4 else sp / 1000
  }
  gp <- g; gp[[var]] <- g[[var]] + h / 2
  gm <- g; gm[[var]] <- g[[var]] - h / 2

  tt <- stats::delete.response(stats::terms(object))
  mmof <- function(z) ilm_drop_intercept(
    stats::model.matrix(tt, data = z, contrasts.arg = object$contrasts), object)
  ## the derivative of the model matrix with respect to `var`; exact for any
  ## term linear in it, and a good approximation for anything else
  mmd <- (mmof(gp) - mmof(gm)) / h

  b <- stats::coef(object)
  mn <- is.null(object$family) || identical(object$family$name, "multinomial")
  if (ncol(mmd) * (if (mn) object$C else 1L) != length(b))
    stop("the reference grid does not match the fitted coefficients; a term ",
         "here is not a plain variable (a smooth or a matrix column), and a ",
         "marginal slope is not defined for it", call. = FALSE)
  V <- suppressWarnings(as.matrix(stats::vcov(object)))

  w <- ilm_emm_cellw(g, mf, weights)
  av <- ilm_emm_avg(mmd, g, specs, w)
  if (mn) return(ilm_trends_multinom(object, av, specs, V, b, level, weights,
                                     var, h))
  L <- av$L; lv <- av$lv
  est <- as.numeric(L %*% b)
  Vem <- L %*% V %*% t(L)
  se <- sqrt(pmax(diag(Vem), 0))

  ## Degrees of freedom for the within-level tests. A slope being different
  ## from zero is the question here, unlike a marginal mean, so the reference
  ## has to be right rather than merely available.
  ddf <- ilm_trend_df(object, L, df)
  ## Kenward-Roger's t is the estimate over its ADJUSTED standard error, on
  ## its df; the df alone, beside the unadjusted one, is neither method
  if (!is.null(attr(ddf, "V"))) {
    Vem <- L %*% attr(ddf, "V") %*% t(L)
    se <- sqrt(pmax(diag(Vem), 0))
  }
  stat <- est / se
  pval <- if (is.finite(ddf[1])) 2 * stats::pt(-abs(stat), ddf) else
    2 * stats::pnorm(-abs(stat))
  crit <- if (is.finite(ddf[1])) stats::qt(1 - (1 - level) / 2, ddf) else
    stats::qnorm(1 - (1 - level) / 2)

  out <- if (length(specs))
    as.data.frame(do.call(rbind, strsplit(lv, "\r", fixed = TRUE)),
                  stringsAsFactors = FALSE)
  else data.frame(.all = "", stringsAsFactors = FALSE)
  names(out) <- if (length(specs)) specs else ".all"
  out$estimate <- est; out$se <- se
  out$df <- as.numeric(ddf)          # drop the method attribute; it lives on the object
  out$statistic <- stat; out$p.value <- pval
  out$lower <- est - crit * se; out$upper <- est + crit * se
  rownames(out) <- NULL

  fam <- if (!is.null(object$family)) object$family$name else "gaussian"
  structure(out, class = c("ilm_trends", "ilm_emm", "data.frame"),
            L = L, V = Vem, specs = specs, weights = weights, type = "link",
            level = level, family = fam, object = object, var = var,
            delta = h, df_method = attr(ddf, "method"))
}

## A multinomial fit's slopes: one per category in every level of `specs`, the
## slope of the category's centred log-odds -- exact and linear, like its
## marginal means, so ilm_contrast() compares them within a category. A slope
## in PROBABILITY is not linear in the coefficients and depends on where it is
## taken; ilm_ame() averages it over the sample, which is the usual question.
#' @keywords internal
#' @noRd
ilm_trends_multinom <- function(object, av, specs, V, b, level, weights, var,
                                h) {
  J <- object$J; C <- object$C; p <- ncol(av$L)
  Tc <- stats::contr.sum(J); cats <- object$ylevels
  nl <- length(av$lv)
  G <- matrix(0, nl * J, p * C)
  for (i in seq_len(nl)) for (j in seq_len(J))
    G[(i - 1L) * J + j, ] <- kronecker(Tc[j, , drop = FALSE],
                                       av$L[i, , drop = FALSE])
  est <- as.numeric(G %*% b)
  Vem <- G %*% V %*% t(G)
  se <- sqrt(pmax(diag(Vem), 0))
  crit <- stats::qnorm(1 - (1 - level) / 2)
  sp <- if (length(specs))
    as.data.frame(do.call(rbind, strsplit(av$lv, "\r", fixed = TRUE)),
                  stringsAsFactors = FALSE)
  else data.frame(.all = "", stringsAsFactors = FALSE)
  names(sp) <- if (length(specs)) specs else ".all"
  out <- sp[rep(seq_len(nl), each = J), , drop = FALSE]
  out$category <- factor(rep(cats, nl), levels = cats)
  out$estimate <- est; out$se <- se; out$df <- Inf
  out$statistic <- ifelse(se > 0, est / se, NA_real_)
  out$p.value <- 2 * stats::pnorm(-abs(out$statistic))
  out$lower <- est - crit * se; out$upper <- est + crit * se
  rownames(out) <- NULL
  structure(out, class = c("ilm_trends", "ilm_emm", "data.frame"),
            L = G, V = Vem, est_c = est, category = as.character(out$category),
            specs = specs, weights = weights, type = "link", level = level,
            family = "multinomial", object = object, var = var, delta = h,
            df_method = "asymptotic")
}

#' Degrees of freedom for each row of a trends table
#'
#' Each row is its own one-degree contrast, so each gets its own df. Under
#' `"auto"` this is the exact residual df where nothing was integrated out,
#' Satterthwaite for a gaussian mixed model, and infinite otherwise -- which
#' recovers the z test the package used before finite df existed.
#'
#' @param object A fitted `"ilm_model"` object.
#' @param L The contrast matrix, one row per level.
#' @param df What the caller asked for: a method name or a number.
#' @return A numeric vector of df, carrying a `"method"` attribute.
#' @keywords internal
#' @noRd
ilm_trend_df <- function(object, L, df) {
  if (is.numeric(df))
    return(structure(rep(as.numeric(df)[1], nrow(L)), method = "supplied"))
  meth <- match.arg(as.character(df)[1],
                    c("auto", "satterthwaite", "kenward-roger", "residual",
                      "asymptotic"))
  ## Kenward-Roger's adjusted covariance is computed once and serves every
  ## row, and it comes back with the df because the standard errors need it.
  ## Asked for where it is not available, it says why instead of quietly
  ## becoming a z test.
  if (identical(meth, "kenward-roger") &&
      identical(object$family$name, "gaussian") && !isTRUE(object$exact_df)) {
    ok <- ilm_kr_applicable(object)
    if (!isTRUE(ok))
      stop("Kenward-Roger is not available for this model: ", ok,
           ". `df = \"satterthwaite\"` does apply here.", call. = FALSE)
    parts <- ilm_kr_parts(object)
    out <- vapply(seq_len(nrow(L)), function(i)
      ilm_df_kr(object, L[i, , drop = FALSE], parts = parts)$df, numeric(1))
    return(structure(out, method = "kenward-roger", V = parts$PhiA))
  }
  out <- vapply(seq_len(nrow(L)), function(i) {
    r <- tryCatch(ilm_denom_df(object, L[i, , drop = FALSE], method = meth),
                  error = function(e) list(df = Inf, method = "asymptotic"))
    r$df
  }, numeric(1))
  used <- tryCatch(ilm_denom_df(object, L[1, , drop = FALSE],
                                method = meth)$method,
                   error = function(e) "asymptotic")
  structure(out, method = used)
}

#' @export
print.ilm_trends <- function(x, digits = 4, ...) {
  v <- attr(x, "var"); sp <- attr(x, "specs")
  cat("Estimated marginal slopes of `", v, "`", sep = "")
  if (length(sp)) cat(" within ", paste(sp, collapse = " x "), sep = "")
  cat("\n")
  dm <- attr(x, "df_method")
  cat("  averaged with ", attr(x, "weights"), " weights; ",
      if (identical(dm, "asymptotic")) "z tests" else
        paste0("t tests on ", dm, " df"), "\n", sep = "")
  d <- as.data.frame(x)
  num <- vapply(d, is.numeric, TRUE)
  d[num] <- lapply(d[num], function(z) round(z, digits))
  print(d, row.names = FALSE)
  if (identical(attr(x, "family"), "multinomial"))
    cat("\n  Slopes of each category's CENTRED log-odds -- its log-probability\n",
        "  less the average over the categories. For the change in each\n",
        "  category's probability, use ilm_ame().\n", sep = "")
  cat("\n  p-values above test each slope against ZERO. To test whether the\n",
      "  slopes DIFFER from one another, pass this to ilm_contrast().\n",
      sep = "")
  invisible(x)
}
