## ---------------------------------------------------------------------------
## Adaptive plotting.
##
## Three ideas distinguish this from a generic EDA plotting set:
##
##   1. The plot carries the DIAGNOSTIC VERDICT. A histogram of a skewed
##      variable says so in its subtitle, using the same machinery that
##      produces `gauss` and `gauss_note` in ilm_describe(). The number tells
##      you something is wrong; the plot shows you what kind.
##
##   2. Rendering is N-AWARE, and says so. A scatter of a million points is a
##      black rectangle, so above a threshold it becomes a binned density and
##      announces the substitution rather than silently misleading.
##
##   3. The geom CHOICE is a separate, inspectable function. `ilm_pick_geom()`
##      returns what would be drawn and why, without drawing anything, so the
##      dispatch can be tested and explained.
##
## Backend is tinyplot, which is an extension of base graphics with no
## dependencies of its own, so these compose with illume's existing base
## graphics diagnostics on the same device.
## ---------------------------------------------------------------------------

ILM_GEOMS <- c("auto", "histogram", "density", "bar", "box", "violin",
               "point", "bin2d", "line", "spine")

## ---- geom choice -----------------------------------------------------------

## Returns the geom that would be drawn and the reason, without drawing. Kept
## separate so the dispatch is testable and so a user can ask "what will this
## do?" before committing to a plot of ten million rows.
#' Which geom would be drawn, and why
#'
#' Returns the choice without drawing anything, so the dispatch can be checked
#' before committing to a plot of many million rows.
#'
#' @param x A vector.
#' @param y An optional second vector.
#' @param n_max Above this many points a scatter becomes a binned density.
#' @return A list with `geom`, `reason` and `n`.
#' @examples
#' d <- ilm_sim()
#' ilm_pick_geom(d$score)
#' ilm_pick_geom(d$score, d$income)
#' @export
ilm_pick_geom <- function(x, y = NULL, n_max = 5000L) {
  cx <- ilm_class_of(x)
  cy <- if (is.null(y)) NA_character_ else ilm_class_of(y)
  n <- if (is.null(y)) sum(!is.na(x)) else sum(!is.na(x) & !is.na(y))

  if (is.null(y)) {
    g <- switch(cx,
      numeric     = "histogram",
      categorical = "bar",
      logical     = "bar",
      time        = "line")
    return(list(geom = g, reason = sprintf("one %s variable", cx), n = n))
  }
  key <- paste(cx, cy, sep = "/")
  g <- switch(key,
    "numeric/numeric"         = if (n > n_max) "bin2d" else "point",
    "numeric/categorical"     = "box",
    "categorical/numeric"     = "box",
    "numeric/logical"         = "box",
    "logical/numeric"         = "box",
    "categorical/categorical" = "spine",
    "time/numeric"            = "line",
    "numeric/time"            = "line",
    "time/categorical"        = "line",
    NULL)
  if (is.null(g))
    return(list(geom = "point", reason = sprintf("%s vs %s (no better default)", cx, cy), n = n))
  reason <- if (identical(g, "bin2d"))
    sprintf("%d points exceeds n_max = %d, so density is drawn instead of points", n, n_max)
    else sprintf("%s vs %s", cx, cy)
  list(geom = g, reason = reason, n = n)
}

## ---- plotting characters by name -------------------------------------------
##
## Base R's `pch` is 26 integers nobody remembers. Nothing in the argument says
## that 16 is a filled circle, and the difference between 16, 19, 20 and 21 is
## not guessable from the numbers -- so a plot gets whichever code the author
## happened to recall, and a reader comparing two plots cannot tell whether a
## difference in the markers was meant.
##
## Names are also checkable, which numbers are not: a pch outside 0:25 is
## accepted by graphics and quietly draws nothing, while a name that is not in
## the table stops here and lists the ones that are.
##
## A SINGLE character is left alone, because base R draws it literally --
## pch = "x" means the letter x, and translating it would silently turn the
## plot into crosses.

#' Plotting characters, by name
#'
#' The lookup behind the `pch` argument of the plotting functions. Numbers pass
#' through after a range check, a single character passes through (base R draws
#' it literally), and a name becomes the number base R wants.
#'
#' @param x Numeric `pch` codes, single characters, or names such as
#'   `"filled circle"`. Case, spaces, underscores, hyphens and dots are all
#'   ignored, so `"filled_circle"` and `"Filled Circle"` are the same thing.
#' @return An integer `pch` vector, or the input unchanged where it was already
#'   numeric or a single character.
#' @keywords internal
#' @noRd
ilm_pch <- function(x) {
  if (is.null(x) || !length(x)) return(x)
  if (is.numeric(x)) {
    bad <- x[is.finite(x) & (x < 0 | x > 25 | x != as.integer(x))]
    if (length(bad))
      stop("`pch` codes run from 0 to 25; got ",
           paste(unique(bad), collapse = ", "),
           ". Names work too, such as \"filled circle\".", call. = FALSE)
    return(x)
  }
  if (!is.character(x)) return(x)
  tab <- ilm_pch_table()
  key <- tolower(gsub("[ _.-]+", "", trimws(x)))
  out <- vector("list", length(x))
  for (i in seq_along(x)) {
    if (is.na(x[i])) { out[[i]] <- NA_integer_; next }
    if (nchar(x[i]) == 1L) { out[[i]] <- x[i]; next }   # a literal glyph
    j <- match(key[i], names(tab))
    if (is.na(j))
      stop("'", x[i], "' is not a plotting character illume knows. ",
           "The names are: ", paste(ilm_pch_names(), collapse = ", "), ".",
           call. = FALSE)
    out[[i]] <- tab[[j]]
  }
  ## a mix of glyphs and codes has to stay character, since that is the only
  ## vector that can carry both
  if (any(vapply(out, is.character, TRUE))) as.character(unlist(out))
  else as.integer(unlist(out))
}

## The names, in code order, primary name first. Aliases follow it, because a
## person reaching for "solid circle" should not have to find out that illume
## calls it something else.
#' @keywords internal
#' @noRd
ilm_pch_spec <- function() {
  list(
    c("0",  "open square", "square", "hollow square", "empty square"),
    c("1",  "open circle", "circle", "hollow circle", "empty circle"),
    c("2",  "open triangle", "triangle", "triangle up", "hollow triangle"),
    c("3",  "plus"),
    c("4",  "cross", "times"),
    c("5",  "open diamond", "diamond", "hollow diamond"),
    c("6",  "open triangle down", "triangle down", "down triangle"),
    c("7",  "square cross", "crossed square"),
    c("8",  "star", "asterisk"),
    c("9",  "diamond plus"),
    c("10", "circle plus"),
    c("11", "star of david", "double triangle"),
    c("12", "square plus"),
    c("13", "circle cross"),
    c("14", "square triangle"),
    c("15", "filled square", "solid square"),
    c("16", "filled circle", "solid circle", "point", "dot"),
    c("17", "filled triangle", "solid triangle"),
    c("18", "filled diamond", "solid diamond"),
    c("19", "large filled circle", "bold circle"),
    c("20", "small filled circle", "bullet", "small dot"),
    c("21", "circle fill", "filled circle outline", "bg circle"),
    c("22", "square fill", "filled square outline", "bg square"),
    c("23", "diamond fill", "filled diamond outline", "bg diamond"),
    c("24", "triangle fill", "filled triangle outline", "bg triangle"),
    c("25", "triangle down fill", "filled triangle down"))
}

#' @keywords internal
#' @noRd
ilm_pch_table <- function() {
  out <- list()
  for (s in ilm_pch_spec()) {
    code <- as.integer(s[1])
    for (nm in s[-1]) out[[gsub("[ _.-]+", "", nm)]] <- code
  }
  out
}

## The primary name of each code, spelled the way a person would write it --
## the lookup keys have had their spaces stripped and would read badly here.
#' @keywords internal
#' @noRd
ilm_pch_names <- function() {
  vapply(ilm_pch_spec(), function(s) sprintf("%s (%s)", s[2], s[1]), "")
}

## NOT a do.call() funnel, though that was tried first and is the obvious
## design. `do.call()` puts the EVALUATED arguments into the call, and tinyplot
## deparses its arguments to title the legend -- so a numeric `by` column came
## back as a deparsed 32-element vector, and the legend width computation threw
## "invalid graphics state". It passed every test and failed R CMD check, on an
## example that had been in the package for months.
##
## So `pch` is a real argument on the functions that draw a marker, placed
## after `...` where it cannot be matched positionally, and translated on the
## way through. Everything else reaches tinyplot exactly as it always did.

## ---- the verdict shown on the plot -----------------------------------------

#' @keywords internal
#' @noRd
ilm_plot_verdict <- function(v) {
  if (is.numeric(v) && !is.logical(v)) {
    g <- ilm_gauss_check(v[is.finite(v)])
    if (is.na(g$gauss)) return(g$gauss_note)
    return(sprintf("gauss %.2f%s", g$gauss,
                   if (nzchar(g$gauss_note)) paste0(" \u2014 ", g$gauss_note) else ""))
  }
  if (is.logical(v)) {
    p <- fmean(v)
    return(if (is.finite(p) && (p < 0.01 || p > 0.99))
      sprintf("p(TRUE) = %.3f \u2014 near-constant (separation risk)", p) else "")
  }
  d <- ilm_describe_cat(v)
  d$note
}

## many levels make an unreadable bar chart, so the tail is pooled and the
## pooling is stated rather than done quietly
#' @keywords internal
#' @noRd
ilm_lump <- function(v, max_levels = 20L) {
  tb <- sort(table(as.character(v)[!is.na(v)]), decreasing = TRUE)
  if (length(tb) <= max_levels) return(list(v = v, note = ""))
  keep <- names(tb)[seq_len(max_levels - 1L)]
  w <- as.character(v); w[!is.na(w) & !(w %in% keep)] <- "(other)"
  list(v = factor(w, levels = c(keep, "(other)")),
       note = sprintf("%d levels pooled into (other)", length(tb) - max_levels + 1L))
}

## ---- the main entry point --------------------------------------------------
## ---- what each geom needs --------------------------------------------------
##
## A single source of truth for the requirements, used both to validate a call
## and to answer "what does this geom need?" at the console. Documentation that
## lives only in a help page goes stale; this cannot, because the checks read
## the same table.
##
##   geom       x            y            size controls
##   ---------  -----------  -----------  -------------
##   histogram  numeric      not used     border width
##   density    numeric      not used     line width
##   bar        categorical  not used     border width
##   box        cat or num   required     line width
##   violin     cat or num   required     line width
##   point      numeric      required     point size
##   bin2d      numeric      required     (not used)
##   line       time or num  required     line width
##   spine      categorical  required     border width

ILM_GEOM_SPEC <- list(
  histogram = list(needs_y = FALSE, x = "numeric",                  size = "lwd"),
  density   = list(needs_y = FALSE, x = "numeric",                  size = "lwd"),
  bar       = list(needs_y = FALSE, x = c("categorical", "logical"), size = "lwd"),
  box       = list(needs_y = TRUE,  x = c("categorical", "logical", "numeric"),
                   y = c("numeric", "categorical", "logical"),      size = "lwd"),
  violin    = list(needs_y = TRUE,  x = c("categorical", "logical", "numeric"),
                   y = c("numeric", "categorical", "logical"),      size = "lwd"),
  point     = list(needs_y = TRUE,  x = "numeric", y = "numeric",   size = "cex"),
  bin2d     = list(needs_y = TRUE,  x = "numeric", y = "numeric",   size = NA),
  line      = list(needs_y = TRUE,  x = c("time", "numeric"),
                   y = c("numeric", "categorical"),                 size = "lwd"),
  spine     = list(needs_y = TRUE,  x = "categorical", y = "categorical", size = "lwd")
)

## A queryable version of the table above, so the requirements can be read
## without opening the help page.
#' What each geom requires
#'
#' The table the argument checks read, so the documented requirements and the
#' enforced ones cannot drift apart.
#'
#' @param geom Optional geom name or names; all of them if `NULL`.
#' @return A data frame with `geom`, `needs_y`, `x_class`, `y_class` and
#'   `size_controls`.
#' @examples
#' ilm_geom_spec()
#' ilm_geom_spec("point")
#' @export
ilm_geom_spec <- function(geom = NULL) {
  g <- if (is.null(geom)) names(ILM_GEOM_SPEC) else geom
  bad <- setdiff(g, names(ILM_GEOM_SPEC))
  if (length(bad))
    stop("unknown `geom`: ", paste(sQuote(bad), collapse = ", "),
         ". Options are ", paste(sQuote(names(ILM_GEOM_SPEC)), collapse = ", "),
         ".", call. = FALSE)
  do.call(rbind, lapply(g, function(k) {
    s <- ILM_GEOM_SPEC[[k]]
    data.frame(geom = k, needs_y = s$needs_y,
               x_class = paste(s$x, collapse = " / "),
               y_class = if (is.null(s$y)) "not used" else paste(s$y, collapse = " / "),
               size_controls = if (is.na(s$size)) "-" else
                 c(lwd = "line/border width", cex = "point size")[[s$size]],
               stringsAsFactors = FALSE)
  }))
}

## ---- the main entry point --------------------------------------------------

#' Adaptive plot of one or two variables
#'
#' Chooses a plot appropriate to the classes of what you give it, annotates it
#' with the same diagnostic verdict [ilm_describe()] reports, and switches to a
#' binned density when there are too many points to show individually.
#'
#' @section Which arguments each geom needs:
#' `geom = "auto"` picks from the data. When you name a geom instead, some
#' require both `x` and `y` and some use `x` alone. [ilm_geom_spec()] returns
#' the full table, and it is the same table the argument checks read, so it
#' cannot drift from the behaviour:
#'
#' \tabular{llll}{
#'   **geom** \tab **x** \tab **y** \tab **`size` controls** \cr
#'   histogram \tab numeric \tab not used \tab border width \cr
#'   density \tab numeric \tab not used \tab line width \cr
#'   bar \tab categorical or logical \tab not used \tab border width \cr
#'   box \tab categorical or numeric \tab required \tab line width \cr
#'   violin \tab categorical or numeric \tab required \tab line width \cr
#'   point \tab numeric \tab required numeric \tab point size \cr
#'   bin2d \tab numeric \tab required numeric \tab not used \cr
#'   line \tab time or numeric \tab required \tab line width \cr
#'   spine \tab categorical \tab required categorical \tab border width
#' }
#'
#' @param data A data frame.
#' @param x Name of the variable on the horizontal axis.
#' @param y Optional second variable; required by some geoms.
#' @param by Optional grouping variable, drawn as colour with a legend.
#' @param geom `"auto"` or one of the geoms in the table above.
#' @param colour,color Colour of lines, points and borders. Synonyms; give one.
#' @param fill Fill colour for bars, boxes, violins and densities.
#' @param alpha Opacity between 0 (transparent) and 1 (opaque).
#' @param size Point size for `"point"`, line or border width elsewhere.
#' @param palette Colours for groups: a vector of colour names, or the name of
#'   a palette such as `"Dark 2"`.
#' @param theme A tinyplot theme name, such as `"clean"`. With tinyplot 0.7.0
#'   or later the available names are listed by `tinyplot::tinytheme_list()`
#'   and an unknown one is reported here; on earlier versions tinyplot reports
#'   it instead.
#' @param n_max Above this many points a scatter becomes a binned density.
#' @param max_levels Categorical levels beyond this are pooled into `(other)`.
#' @param verdict Annotate the plot with the diagnostic verdict.
#' @param main Plot title.
#' @param ... Passed to [tinyplot::tinyplot()].
#' @param pch Plotting character. Takes a NAME as well as a number:
#'   `"filled circle"` is 16, and every code from 0 to 25 has one. Case,
#'   spaces, underscores and hyphens are ignored. A single character is
#'   drawn literally, so `pch = "x"` is still the letter x.
#' @return Invisibly, a list with the geom used, the reason, and the note shown.
#' @seealso [ilm_geom_spec()], [ilm_pick_geom()], [ilm_plot_all()],
#'   [ilm_plot_model()].
#' @examples
#' d <- ilm_sim()
#' ilm_plot(d, "score")
#' ilm_plot(d, "grp", "score", geom = "violin", fill = "lightblue")
#' ilm_plot(d, "score", "income", by = "grp", alpha = 0.6)
#' @export
ilm_plot <- function(data, x, y = NULL, by = NULL, geom = "auto",
                     colour = NULL, color = NULL, fill = NULL, alpha = NULL,
                     size = NULL, palette = NULL, theme = NULL,
                     n_max = 5000L, max_levels = 20L, verdict = TRUE,
                     main = NULL, ..., pch = NULL) {
  if (!is.data.frame(data)) stop("`data` must be a data frame", call. = FALSE)
  if (length(geom) != 1L || !geom %in% ILM_GEOMS) {
    auto <- if (is.character(x) && length(x) == 1L && x %in% names(data))
      ilm_pick_geom(data[[x]], if (!is.null(y) && y %in% names(data)) data[[y]],
                    n_max)$geom else NA
    stop("unknown `geom`: ", paste(sQuote(geom), collapse = ", "),
         ". Options are ", paste(sQuote(ILM_GEOMS), collapse = ", "),
         if (!is.na(auto)) paste0(". \"auto\" would use ", sQuote(auto)), ".",
         call. = FALSE)
  }
  miss <- setdiff(c(x, y, by), names(data))
  if (length(miss))
    stop("column(s) not found in the data: ", paste(miss, collapse = ", "),
         ". Available: ", paste(utils::head(names(data), 12), collapse = ", "),
         if (ncol(data) > 12) ", ..." else "", call. = FALSE)

  ## `color` is accepted as a synonym so the call works either way
  if (!is.null(color)) {
    if (!is.null(colour) && !identical(colour, color))
      stop("`colour` and `color` were both given with different values; ",
           "use one", call. = FALSE)
    colour <- color
  }
  if (!is.null(theme)) {
    if (length(theme) != 1L || !is.character(theme))
      stop("`theme` must be a single theme name, not ", length(theme), " ",
           class(theme)[1], " value(s)", call. = FALSE)
    ## tinytheme_list() arrived in tinyplot 0.7.0. Older tinyplot still accepts
    ## `theme` and still has tinytheme(); it just cannot enumerate the names.
    ## Validate against the list where there is one and let tinyplot object
    ## where there is not, rather than refusing a theme that would have worked.
    ok <- tryCatch(unlist(getExportedValue("tinyplot", "tinytheme_list")(),
                          use.names = FALSE), error = function(e) NULL)
    if (!is.null(ok) && !theme %in% ok)
      stop("unknown `theme`: ", sQuote(theme), ". Options are ",
           paste(sQuote(ok), collapse = ", "), ".", call. = FALSE)
  }
  if (!is.null(alpha) && (!is.numeric(alpha) || length(alpha) != 1L ||
                          alpha < 0 || alpha > 1))
    stop("`alpha` must be a single number between 0 (transparent) and 1 (opaque)",
         call. = FALSE)
  if (!is.null(size) && (!is.numeric(size) || length(size) != 1L || size <= 0))
    stop("`size` must be a single positive number", call. = FALSE)

  xv <- data[[x]]; yv <- if (is.null(y)) NULL else data[[y]]
  pick <- ilm_pick_geom(xv, yv, n_max)
  g <- if (identical(geom, "auto")) pick$geom else geom
  spec <- ILM_GEOM_SPEC[[g]]

  ## requirements, checked against the same table the documentation shows
  if (spec$needs_y && is.null(y)) {
    need <- names(ILM_GEOM_SPEC)[vapply(ILM_GEOM_SPEC, `[[`, TRUE, "needs_y")]
    stop("geom ", sQuote(g), " needs both `x` and `y`. Geoms requiring `y`: ",
         paste(sQuote(need), collapse = ", "), ". Geoms using `x` alone: ",
         paste(sQuote(setdiff(names(ILM_GEOM_SPEC), need)), collapse = ", "),
         ".", call. = FALSE)
  }
  if (!spec$needs_y && !is.null(y))
    warning("geom ", sQuote(g), " uses `x` only; `y` (", y, ") is ignored",
            call. = FALSE)
  cx <- ilm_class_of(xv)
  if (!cx %in% spec$x) {
    alt <- names(ILM_GEOM_SPEC)[vapply(ILM_GEOM_SPEC,
             function(s) cx %in% s$x, TRUE)]
    stop("geom ", sQuote(g), " needs `x` to be ", paste(spec$x, collapse = " or "),
         ", but ", sQuote(x), " is ", cx,
         if (length(alt)) paste0(". Geoms that accept a ", cx, " `x`: ",
                                 paste(sQuote(alt), collapse = ", ")), ".",
         call. = FALSE)
  }
  if (spec$needs_y && !is.null(spec$y)) {
    cy <- ilm_class_of(yv)
    if (!cy %in% spec$y)
      stop("geom ", sQuote(g), " needs `y` to be ", paste(spec$y, collapse = " or "),
           ", but ", sQuote(y), " is ", cy, ".", call. = FALSE)
  }

  notes <- character(0)
  if (identical(geom, "auto") && identical(g, "bin2d")) notes <- c(notes, pick$reason)
  if (g %in% c("bar", "box", "violin", "spine")) {
    cv <- if (ilm_class_of(xv) == "categorical") "x" else
          if (!is.null(yv) && ilm_class_of(yv) == "categorical") "y" else NA
    if (!is.na(cv)) {
      lp <- ilm_lump(if (cv == "x") xv else yv, max_levels)
      if (nzchar(lp$note)) {
        notes <- c(notes, lp$note)
        if (cv == "x") xv <- lp$v else yv <- lp$v
      }
    }
  }
  sub <- if (verdict) ilm_plot_verdict(if (is.null(yv)) xv else yv) else ""
  sub <- paste(c(sub, notes)[nzchar(c(sub, notes))], collapse = "  |  ")
  if (nzchar(sub)) message("ilm_plot: ", sub)

  ttl <- if (!is.null(main)) main else if (is.null(y)) x else paste(y, "by", x)

  if (identical(g, "bin2d")) return(invisible(
    ilm_bin2d(xv, yv, xlab = x, ylab = y, main = ttl, sub = sub,
              palette = palette, ...)))

  tt <- switch(g, histogram = "histogram", density = "density", bar = "barplot",
               box = "boxplot", violin = "violin", point = "p", line = "l",
               spine = "spineplot")
  ## The formula interface is used rather than passing vectors. Passing vectors
  ## through do.call() silently produces an EMPTY plot -- no error, no warning,
  ## just a blank device -- while the formula path is unaffected. Lumped or
  ## coerced columns are written back into a working copy so the formula can
  ## refer to them by name.
  bt <- function(nm) paste0("`", nm, "`")
  dwork <- data
  dwork[[x]] <- xv
  if (!is.null(y)) dwork[[y]] <- yv
  rhs <- paste0(bt(x), if (!is.null(by)) paste0(" | ", bt(by)))
  fml <- stats::as.formula(paste(if (is.null(y)) "" else bt(y), "~", rhs))
  args <- list(fml, data = dwork, type = tt, main = ttl, sub = sub, xlab = x)
  if (!is.null(y)) args$ylab <- y
  if (!is.null(by)) {
    ## tinyplot reserves legend space by altering the device layout, which
    ## silently blanks the figure when one is already in force (par(mfrow),
    ## layout()). Inside a multi-panel layout the legend is dropped instead.
    mf <- graphics::par("mfrow")
    dots <- list(...)
    if (!identical(mf, c(1L, 1L)) && is.null(dots$legend)) {
      args$legend <- FALSE
      message("ilm_plot: legend suppressed inside a multi-panel layout ",
              "(it would reset the layout); plot on its own device for a legend")
    }
  }
  ## ggplot2-style aesthetic names mapped onto tinyplot's arguments
  if (!is.null(colour))  args$col <- colour
  if (!is.null(fill))    args$fill <- fill
  if (!is.null(alpha))   args$alpha <- alpha
  if (!is.null(palette)) args$palette <- palette
  if (!is.null(theme))   args$theme <- theme
  ## `size` means point size for a scatter and line width elsewhere, per the
  ## size_controls column of ilm_geom_spec()
  if (!is.null(size) && !is.na(spec$size)) args[[spec$size]] <- size
  if (!is.null(pch)) args$pch <- ilm_pch(pch)
  do.call(tinyplot::tinyplot, c(args, list(...)))
  invisible(list(geom = g, reason = pick$reason, note = sub))
}

## 2D binned counts, drawn with base image(). tinyplot has no hex or bin2d
## type, and a dependency for one fallback is not worth it.
#' @keywords internal
#' @noRd
ilm_bin2d <- function(x, y, bins = 60L, xlab = "", ylab = "", main = "", sub = "",
                      palette = NULL, ...) {
  ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]
  bx <- seq(min(x), max(x), length.out = bins + 1L)
  by_ <- seq(min(y), max(y), length.out = bins + 1L)
  ix <- findInterval(x, bx, rightmost.closed = TRUE, all.inside = TRUE)
  iy <- findInterval(y, by_, rightmost.closed = TRUE, all.inside = TRUE)
  m <- matrix(0L, bins, bins)
  tb <- table(ix, iy)
  m[cbind(as.integer(rownames(tb))[row(tb)], as.integer(colnames(tb))[col(tb)])] <- as.vector(tb)
  m[m == 0L] <- NA_integer_
  graphics::image(x = (bx[-1] + bx[-length(bx)]) / 2,
                  y = (by_[-1] + by_[-length(by_)]) / 2, z = m,
                  xlab = xlab, ylab = ylab, main = main,
                  col = if (is.null(palette)) grDevices::hcl.colors(32, "YlGnBu", rev = TRUE)
                        else if (length(palette) > 1L) grDevices::colorRampPalette(palette)(32)
                        else grDevices::hcl.colors(32, palette, rev = TRUE))
  if (nzchar(sub)) graphics::mtext(sub, side = 1, line = 4, cex = 0.7)
  invisible(list(geom = "bin2d", bins = bins))
}

## ---- small multiples -------------------------------------------------------

#' Small multiples of every variable
#'
#' Legends are suppressed in this layout: tinyplot reserves legend space by
#' altering the device layout, which blanks a multi-panel figure. Plot a single
#' variable with [ilm_plot()] when you need one.
#'
#' @inheritParams ilm_plot
#' @param class `"all"`, or one or more of `"numeric"`, `"categorical"`,
#'   `"logical"`, `"time"`.
#' @param max_panels Stop after this many variables.
#' @return Invisibly, a list of the per-panel results.
#' @examples
#' ilm_plot_all(ilm_sim(), class = "numeric", max_panels = 4)
#' @export
ilm_plot_all <- function(data, by = NULL, class = "all", max_panels = 12L,
                         n_max = 5000L, verdict = FALSE, ...) {
  bad <- setdiff(class, c("all", ILM_CLASSES))
  if (length(bad))
    stop("unknown `class`: ", paste(sQuote(bad), collapse = ", "),
         ". Options are ", paste(sQuote(c("all", ILM_CLASSES)), collapse = ", "),
         ".", call. = FALSE)
  cand <- setdiff(names(data), by)
  cl <- vapply(data[cand], ilm_class_of, "")
  if (!identical(class, "all")) cand <- cand[cl %in% class]
  if (!length(cand)) stop("no columns of class ",
                          paste(sQuote(class), collapse = ", "), " in the data",
                          call. = FALSE)
  if (length(cand) > max_panels) {
    message("ilm_plot_all: showing the first ", max_panels, " of ",
            length(cand), " variables; raise `max_panels` to see more")
    cand <- cand[seq_len(max_panels)]
  }
  op <- graphics::par(mfrow = grDevices::n2mfrow(length(cand)))
  on.exit(graphics::par(op), add = TRUE)
  if (!is.null(by))
    message("ilm_plot_all: legends are suppressed in small multiples; ",
            "use ilm_plot() for a single variable with a legend")
  invisible(lapply(cand, function(cn)
    suppressMessages(ilm_plot(data, cn, by = by, n_max = n_max,
                              verdict = verdict, ...))))
}

## ---- missingness -----------------------------------------------------------

#' Proportion missing, by variable
#'
#' @param data A data frame.
#' @param ... Passed to [graphics::barplot()].
#' @return Invisibly, the named vector of proportions.
#' @seealso [ilm_describe_na_all()] for the same information as a table.
#' @examples
#' ilm_plot_missing(ilm_sim())
#' @export
ilm_plot_missing <- function(data, ...) {
  if (!is.data.frame(data)) stop("`data` must be a data frame", call. = FALSE)
  p <- vapply(data, function(v) mean(is.na(v)), 1)
  if (all(p == 0)) {
    message("ilm_plot_missing: no missing values"); return(invisible(p))
  }
  ## drawn by ilm_plot_na_all() so there is one implementation and one
  ## backend: this was the last base-graphics plot in the package
  ilm_plot_na_all(data, ...)
  invisible(sort(p, decreasing = TRUE))
}
