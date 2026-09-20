## ---------------------------------------------------------------------------
## Model-side plots: what makes this more than an exploratory set.
##
## Same principle as the data-side plots: say what the picture means, and say
## when it cannot be trusted.
##
##   coef      fixed effects with confidence intervals, zero marked
##   effect    partial dependence for one predictor, holding others fixed
##   random    caterpillar plot of the estimated random effects
##   residual  delegates to ilm_appraise(), which already owns this
##
## Both an ilm_plot_model() function and a plot() method are provided: the
## method for interactive use, the function so it can be called explicitly and
## found by prefix completion.
## ---------------------------------------------------------------------------

ILM_MODEL_PLOTS <- c("coef", "effect", "random", "residual", "acf")

#' Diagnostic and summary plots for a fitted model
#'
#' @param model A fitted `"ilm_model"` object.
#' @param what `"coef"` for fixed effects with intervals, `"effect"` for the
#'   partial effect of one predictor, `"random"` for a caterpillar plot of the
#'   estimated group effects, `"residual"` to hand off to [ilm_appraise()], or
#'   `"acf"` to hand off to [ilm_plot_acf()], which also needs `time` and
#'   `group`.
#' @param term Which predictor or grouping factor to show, for `"effect"` and
#'   `"random"`.
#' @param conf Confidence level for intervals.
#' @param colour,fill,alpha,size Appearance.
#' @param main Plot title.
#' @param ... Passed to the underlying plot.
#' @return Invisibly, the data behind the plot.
#' @seealso [ilm_appraise()] for the full residual panel.
#' @examples
#' set.seed(1)
#' d <- ilm_sim()
#' f <- ilm_model(score ~ income + grp + (1 | id), data = d,
#'                family = "gaussian", verbose = FALSE)
#' ilm_plot_model(f, "coef")
#' ilm_plot_model(f, "random")
#' @export
ilm_plot_model <- function(model, what = "coef", term = NULL, conf = 0.95,
                           colour = "black", fill = "grey70", alpha = NULL,
                           size = 1, main = NULL, ...) {
  if (!inherits(model, "ilm_model"))
    stop("`model` must be a fitted ilm_model object, not ", class(model)[1],
         call. = FALSE)
  if (length(what) != 1L || !what %in% ILM_MODEL_PLOTS)
    stop("unknown `what`: ", paste(sQuote(what), collapse = ", "),
         ". Options are ", paste(sQuote(ILM_MODEL_PLOTS), collapse = ", "), ".",
         call. = FALSE)
  switch(what,
    coef     = ilm_plot_coef(model, conf, colour, size, main, ...),
    effect   = ilm_plot_effect(model, term, conf, colour, fill, alpha, main, ...),
    random   = ilm_plot_random(model, term, conf, colour, size, main, ...),
    residual = ilm_appraise(model, ...),
    acf      = ilm_plot_acf(model, colour = colour, fill = fill, alpha = alpha,
                            size = size, main = main, ...))
}

#' Plot a fitted model
#'
#' A shorthand for [ilm_plot_model()], which is the fuller interface.
#'
#' @param x A fitted `"ilm_model"` object.
#' @param what Which plot to draw; see [ilm_plot_model()].
#' @param ... Passed to [ilm_plot_model()].
#' @return Invisibly, the data behind the plot.
#' @examples
#' set.seed(1)
#' d <- ilm_sim()
#' f <- ilm_model(score ~ income + (1 | id), data = d, family = "gaussian",
#'                verbose = FALSE)
#' plot(f)
#' @exportS3Method base::plot
plot.ilm_model <- function(x, what = "coef", ...) ilm_plot_model(x, what = what, ...)

## ---- coefficients ----------------------------------------------------------

#' @keywords internal
#' @noRd
ilm_plot_coef <- function(model, conf = 0.95, colour = "black", size = 1,
                          main = NULL, intercept = FALSE, ...) {
  b <- coef(model)
  s <- suppressWarnings(sqrt(diag(vcov(model))))
  if (!length(b) || all(!is.finite(s)))
    stop("this model has no usable fixed-effect standard errors; ",
         "see model$checks", call. = FALSE)
  keep <- if (intercept) rep(TRUE, length(b)) else !grepl("(Intercept)", names(b), fixed = TRUE)
  if (!any(keep))
    stop("the model has only an intercept; use `intercept = TRUE` to show it",
         call. = FALSE)
  b <- b[keep]; s <- s[keep]
  ## a gaussian model with nothing integrated out has an exact t reference
  crit <- if (isTRUE(model$exact_df)) stats::qt(1 - (1 - conf) / 2, model$resid_df)
          else stats::qnorm(1 - (1 - conf) / 2)
  lo <- b - crit * s; hi <- b + crit * s
  i <- seq_along(b)
  op <- graphics::par(mar = c(4, max(6, max(nchar(names(b))) * 0.65), 3, 2))
  on.exit(graphics::par(op), add = TRUE)
  graphics::plot(b, i, xlim = range(c(lo, hi, 0), na.rm = TRUE), yaxt = "n",
                 ylab = "", xlab = "estimate", pch = 19, cex = size, col = colour,
                 main = main %||% sprintf("Fixed effects (%.0f%% %s)", conf * 100,
                   if (isTRUE(model$exact_df)) "t interval" else "Wald interval"),
                 ylim = c(0.5, length(b) + 0.5), ...)
  graphics::axis(2, i, names(b), las = 1, cex.axis = 0.85)
  graphics::abline(v = 0, lty = 2, col = "grey50")
  graphics::segments(lo, i, hi, i, col = colour, lwd = 1.6)
  ## intervals excluding zero are the ones a reader will act on
  sig <- is.finite(lo) & (lo > 0 | hi < 0)
  if (any(sig)) graphics::points(b[sig], i[sig], pch = 19, cex = size * 1.25, col = colour)
  invisible(data.frame(term = names(b), estimate = b, lower = lo, upper = hi,
                       excludes_zero = sig, row.names = NULL))
}

`%||%` <- function(a, b) if (is.null(a)) b else a

## ---- partial dependence ----------------------------------------------------

## One predictor varied across its range with the others held at typical
## values. This is a conditional effect, not a marginal one: it answers "what
## does the model say about x, all else equal", which is the question a
## coefficient answers too but on the response scale.
#' @keywords internal
#' @noRd
ilm_plot_effect <- function(model, term = NULL, conf = 0.95, colour = "black",
                            fill = "grey70", alpha = NULL, main = NULL,
                            n_grid = 60L, data = NULL, ...) {
  ## A partial-dependence grid needs the observed predictors, to hold the ones
  ## not being varied at typical values. Three sources, in order of reliability:
  ## the frame the fit carried, the data named in the original call, or a frame
  ## the caller supplies.
  ## ilm_model() stores the model frame as $model, the same name lm() uses
  mf <- model$model
  if (!is.null(mf)) mf <- as.data.frame(mf)
  if (is.null(mf) && !is.null(data)) mf <- as.data.frame(data)
  if (is.null(mf)) {
    cl <- tryCatch(stats::getCall(model), error = function(e) NULL)
    if (!is.null(cl) && !is.null(cl$data))
      mf <- tryCatch(as.data.frame(eval(cl$data, parent.frame())),
                     error = function(e) NULL)
  }
  if (is.null(mf)) {
    tl <- model$term_labels
    stop("an effect plot needs the observed predictors, and this model does ",
         "not carry them. Pass the data used to fit it, e.g. ",
         "ilm_plot_model(model, \"effect\", data = your_data)",
         if (!is.null(tl)) paste0(". Terms: ", paste(tl, collapse = ", ")), ".",
         call. = FALSE)
  }
  resp <- if (!is.null(model$formula)) all.vars(model$formula[[2]]) else character(0)
  tl <- model$term_labels
  vars <- intersect(names(mf), if (!is.null(tl)) tl else setdiff(names(mf), resp))
  if (!length(vars)) vars <- setdiff(names(mf), resp)
  if (is.null(term)) {
    term <- vars[1]
    message("ilm_plot_model: no `term` given, showing ", sQuote(term),
            "; available: ", paste(vars, collapse = ", "))
  }
  if (!term %in% vars)
    stop("`term` (", term, ") is not a predictor in this model. Available: ",
         paste(vars, collapse = ", "), call. = FALSE)

  v <- mf[[term]]
  grid <- if (is.numeric(v) && !is.logical(v))
    seq(min(v, na.rm = TRUE), max(v, na.rm = TRUE), length.out = n_grid)
    else sort(unique(v[!is.na(v)]))
  nd <- mf[rep(1L, length(grid)), , drop = FALSE]
  for (cn in vars) if (!identical(cn, term)) {
    cv <- mf[[cn]]
    nd[[cn]] <- if (is.numeric(cv) && !is.logical(cv)) stats::median(cv, na.rm = TRUE)
                else { tb <- sort(table(cv), decreasing = TRUE); as.character(names(tb)[1]) }
    if (is.factor(mf[[cn]])) nd[[cn]] <- factor(nd[[cn]], levels = levels(mf[[cn]]))
  }
  nd[[term]] <- grid
  pr <- tryCatch(stats::predict(model, newdata = nd, se.fit = TRUE),
                 error = function(e)
                   stop("predict() failed for this model: ", conditionMessage(e),
                        call. = FALSE))
  fit <- if (is.list(pr)) pr$fit else pr
  se  <- if (is.list(pr) && !is.null(pr$se.fit)) pr$se.fit else NULL
  fit <- as.matrix(fit)[, 1]
  crit <- stats::qnorm(1 - (1 - conf) / 2)

  ttl <- main %||% paste("Partial effect of", term)
  if (is.numeric(grid)) {
    ylim <- if (!is.null(se)) range(c(fit - crit * as.matrix(se)[, 1],
                                      fit + crit * as.matrix(se)[, 1])) else range(fit)
    graphics::plot(grid, fit, type = "n", xlab = term, ylab = "predicted",
                   main = ttl, ylim = ylim, ...)
    if (!is.null(se)) {
      s1 <- as.matrix(se)[, 1]
      graphics::polygon(c(grid, rev(grid)),
                        c(fit - crit * s1, rev(fit + crit * s1)),
                        col = grDevices::adjustcolor(fill, alpha %||% 0.5),
                        border = NA)
    }
    graphics::lines(grid, fit, col = colour, lwd = 2)
  } else {
    graphics::plot(seq_along(grid), fit, xaxt = "n", xlab = term,
                   ylab = "predicted", main = ttl, pch = 19, col = colour,
                   ylim = if (!is.null(se)) range(c(fit - crit * as.matrix(se)[,1],
                                                    fit + crit * as.matrix(se)[,1]))
                          else range(fit), ...)
    graphics::axis(1, seq_along(grid), as.character(grid))
    if (!is.null(se)) { s1 <- as.matrix(se)[, 1]
      graphics::segments(seq_along(grid), fit - crit * s1,
                         seq_along(grid), fit + crit * s1, col = colour) }
  }
  graphics::mtext(sprintf("other predictors held at median / most common level"),
                  side = 3, line = 0.2, cex = 0.7, col = "grey35")
  invisible(data.frame(term = term, value = as.character(grid), fit = fit,
                       row.names = NULL))
}

## ---- random effects --------------------------------------------------------

## A caterpillar plot: estimated group effects in rank order. A flat line means
## the grouping factor is doing nothing; a few points far from zero are the
## groups driving the variance component.
#' @keywords internal
#' @noRd
ilm_plot_random <- function(model, term = NULL, conf = 0.95, colour = "black",
                            size = 0.7, main = NULL, ...) {
  if (!length(model$re))
    stop("this model has no random effects to plot", call. = FALSE)
  nms <- names(model$re)
  grp <- vapply(model$re, function(e) identical(e$kind, "group"), TRUE)
  if (!any(grp))
    stop("this model has only smooth (basis) terms, which are not caterpillar ",
         "plots; use what = \"effect\"", call. = FALSE)
  nms <- nms[grp]
  if (is.null(term)) {
    term <- nms[1]
    if (length(nms) > 1L)
      message("ilm_plot_model: no `term` given, showing ", sQuote(term),
              "; available: ", paste(nms, collapse = ", "))
  }
  if (!term %in% nms)
    stop("`term` (", term, ") is not a grouping factor in this model. ",
         "Available: ", paste(nms, collapse = ", "), call. = FALSE)
  k <- match(term, names(model$re))
  B <- ilm_Bhat_term(model, k)
  b <- as.matrix(B)[, 1]
  o <- order(b); b <- b[o]
  graphics::plot(b, seq_along(b), pch = 19, cex = size, col = colour,
                 xlab = "estimated effect", ylab = paste(term, "(ranked)"),
                 main = main %||% paste("Random effects:", term), ...)
  graphics::abline(v = 0, lty = 2, col = "grey50")
  invisible(data.frame(level = o, effect = b, row.names = NULL))
}
