## ---------------------------------------------------------------------------
## Seeing a moderation, rather than only testing it.
##
## NOT the partial dependence plot. ilm_plot_effect() varies one predictor with
## the others HELD AT TYPICAL VALUES, which is the right picture for "what does
## the model say about x, all else equal" and exactly the wrong one here:
## pinning the moderator at its mean or modal level collapses the thing being
## shown. So this is a separate function rather than another `what =`.
##
## It is built on ilm_emmeans() and ilm_trends() rather than on new estimation
## code, so the picture and the p-value come from the same fit through the same
## estimator. Two shapes, chosen by what `x` is:
##
##   x categorical  ->  cell means of x within each level (or at values) of the
##                      moderator. The classic interaction plot, which is what
##                      the audience for this already reads.
##   x continuous   ->  the SLOPE of x within each level of the moderator, from
##                      ilm_trends(), which is the quantity that is actually
##                      being said to differ.
## ---------------------------------------------------------------------------

#' See a moderation
#'
#' Draws the effect of `x` at each level, or several values, of a moderator --
#' the picture behind an [ilm_moderation()] row.
#'
#' For a categorical `x` this is the cell means, the interaction plot. For a
#' continuous `x` it is the **slope** of `x` within each level of the
#' moderator, from [ilm_trends()], because a slope is the quantity being said
#' to differ. Both carry intervals.
#'
#' This is deliberately not [ilm_plot_model()]'s effect plot, which holds the
#' other predictors at typical values -- pinning the moderator is exactly what
#' would hide the moderation.
#'
#' @param x An [ilm_moderation()] result, or a fitted [ilm_model()] that
#'   already contains the interaction.
#' @param moderator Which moderator to draw. Defaults to the strongest row.
#' @param exposure The treatment or exposure, when `x` is a plain model.
#' @param at For a continuous moderator, the values to evaluate at. Defaults to
#'   its quartiles, so the picture spans the data rather than a nominal range.
#' @param level Confidence level.
#' @param data Needed only if a refit is required and the moderator is not in
#'   the model frame.
#' @param main,colour,pch Passed through; `pch` takes a name.
#' @param ... Passed to [tinyplot::tinyplot()].
#' @return Invisibly, the data frame that was plotted.
#' @seealso [ilm_moderation()], [ilm_trends()], [ilm_emmeans()].
#' @examples
#' \donttest{
#' set.seed(1)
#' d <- data.frame(tx = rbinom(300, 1, 0.5), age = rnorm(300),
#'                 site = factor(sample(c("a", "b", "c"), 300, TRUE)))
#' d$y <- 0.4 * d$tx + 0.8 * d$tx * (d$site == "c") + rnorm(300)
#' fit <- ilm_model(y ~ tx + age + site, data = d, verbose = FALSE)
#' m <- ilm_moderation(fit, x = "tx", progress = FALSE)
#' ilm_plot_moderation(m)
#' }
#' @export
ilm_plot_moderation <- function(x, moderator = NULL, exposure = NULL,
                                at = NULL, level = 0.95, data = NULL,
                                main = NULL, colour = NULL, pch = NULL, ...) {
  if (!requireNamespace("tinyplot", quietly = TRUE))
    stop("the tinyplot package is needed to draw this.", call. = FALSE)
  z <- ilm_plot_mod_resolve(x, moderator, exposure, data)
  fit <- z$fit; ex <- z$x; mod <- z$moderator
  mf <- fit$model
  if (!mod %in% names(mf))
    stop("the moderator '", mod, "' is not in the refitted model's frame.",
         call. = FALSE)
  mv <- mf[[mod]]; xv <- mf[[ex]]
  x_cat <- is.factor(xv) || is.character(xv) || is.logical(xv) ||
    length(unique(stats::na.omit(xv))) <= 2L
  m_cat <- is.factor(mv) || is.character(mv) || is.logical(mv)
  atl <- list()
  if (!m_cat)
    atl[[mod]] <- if (is.null(at))
      stats::quantile(mv, c(.25, .5, .75), na.rm = TRUE) else at
  ## A 0/1 treatment stored as a NUMBER is categorical in meaning and numeric
  ## in type, so ilm_emmeans() would hold it at its mean -- giving one curve
  ## where two are wanted, and silently. Its two values are named explicitly.
  if (x_cat && is.numeric(xv))
    atl[[ex]] <- sort(unique(stats::na.omit(xv)))
  if (!length(atl)) atl <- NULL

  if (x_cat) {
    ## cell means: the interaction plot
    e <- as.data.frame(ilm_emmeans(fit, c(ex, mod), at = atl, level = level))
    e[[ex]] <- factor(e[[ex]]); e[[mod]] <- factor(round_lab(e[[mod]]))
    ttl <- main %||% paste0("Effect of ", ex, " by ", mod)
    ylb <- "estimated marginal mean"
    d <- data.frame(xv = e[[mod]], y = e$estimate, lo = e$lower, hi = e$upper,
                    by = e[[ex]], stringsAsFactors = FALSE)
    xlb <- mod
  } else {
    ## the SLOPE of x within each level of the moderator
    tr <- as.data.frame(ilm_trends(fit, mod, ex, at = atl, level = level))
    lab <- tr[[mod]]
    d <- data.frame(xv = factor(round_lab(lab)), y = tr$estimate,
                    lo = tr$lower, hi = tr$upper, by = factor("slope"),
                    stringsAsFactors = FALSE)
    ttl <- main %||% paste0("Slope of ", ex, " by ", mod)
    ylb <- paste0("slope of ", ex); xlb <- mod
  }
  ## a reference line at no effect, because "do the intervals cross zero" is
  ## the question a reader brings to this picture
  args <- list(x = d$xv, y = d$y, ymin = d$lo, ymax = d$hi,
               type = "pointrange", xlab = xlb, ylab = ylb, main = ttl)
  if (nlevels(d$by) > 1L) {
    args$by <- d$by
    ## The legend title must be given explicitly. Building the call with
    ## do.call() puts the EVALUATED `by` vector into it, and tinyplot deparses
    ## its arguments to title a legend -- so the title becomes the deparsed
    ## factor and the width computation throws "invalid graphics state". It
    ## depends on how long the level NAMES are, which is why a test using
    ## "a"/"b"/"c" passes while "north"/"central"/"south" does not.
    args$legend <- list(title = ex)
  }
  if (!is.null(colour)) args$col <- colour
  if (!is.null(pch)) args$pch <- ilm_pch(pch)
  do.call(tinyplot::tinyplot, c(args, list(...)))
  if (!x_cat) graphics::abline(h = 0, lty = 2, col = "grey50")
  invisible(d)
}

#' @keywords internal
#' @noRd
round_lab <- function(v) if (is.numeric(v)) signif(v, 3) else as.character(v)

## Work out which fit, exposure and moderator are being drawn, refitting only
## if the moderation object did not retain the one that was asked for.
#' @keywords internal
#' @noRd
ilm_plot_mod_resolve <- function(x, moderator, exposure, data) {
  if (inherits(x, "ilm_moderation")) {
    d <- as.data.frame(x); class(d) <- "data.frame"
    ex <- attr(x, "x")
    if (is.null(moderator)) moderator <- d$moderator[1]
    if (!moderator %in% d$moderator)
      stop("`moderator` (", moderator, ") was not among the candidates: ",
           paste(d$moderator, collapse = ", "), call. = FALSE)
    fits <- attr(x, "fits")
    i <- match(moderator, d$moderator)
    fit <- if (i <= length(fits) && !is.null(fits[[i]])) fits[[i]] else NULL
    if (is.null(fit))
      stop("the refit for '", moderator, "' was not retained -- ",
           "ilm_moderation() keeps the strongest `n_keep` of them. Re-run ",
           "with a larger `n_keep`, or pass the model itself.", call. = FALSE)
    return(list(fit = fit, x = ex, moderator = moderator))
  }
  if (!inherits(x, "ilm_model"))
    stop("`x` must be an ilm_moderation() result or a fitted ilm_model(); ",
         "it is ", class(x)[1], call. = FALSE)
  tl <- x$term_labels
  it <- grep(":", tl, value = TRUE)
  if (!length(it))
    stop("this model has no interaction to draw. Fit one, or pass an ",
         "ilm_moderation() result.", call. = FALSE)
  if (is.null(exposure) || is.null(moderator)) {
    ## plain names, because they are looked up in the model frame next
    p <- ilm_unbq(strsplit(it[1], ":", fixed = TRUE)[[1]])
    if (is.null(exposure)) exposure <- p[1]
    if (is.null(moderator)) moderator <- setdiff(p, exposure)[1]
    message("ilm_plot_moderation(): drawing `", exposure, "` by `", moderator,
            "`; name `exposure` and `moderator` to choose another.")
  }
  list(fit = x, x = exposure, moderator = moderator)
}
