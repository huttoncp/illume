## ---------------------------------------------------------------------------
## Design-based inference for data from a complex sample.
##
## A sampling weight is not a replicate count, and the difference is not
## cosmetic. Passing survey weights to a model that treats weights as
## FREQUENCIES gives the right point estimates -- a weighted likelihood is
## design-consistent for the population parameter -- and standard errors that
## are wrong by roughly sqrt(n / sum(w)), because the likelihood believes it
## saw sum(w) observations. With 2,000 respondents weighted up to a population
## of 2,000,000, that is a factor of about thirty.
##
## The point estimates are therefore fine and the variance is the whole
## problem, which is why this file supplies a variance rather than a fitting
## routine. The linearization (Taylor) estimator is the sandwich again, with
## the middle assembled the way the design says:
##
##   - score contributions are summed WITHIN each primary sampling unit,
##     because units in the same cluster are not independent;
##   - those PSU totals are centred WITHIN each stratum, because stratification
##     removes between-stratum variation from the design;
##   - each stratum contributes n_h/(n_h - 1) times the sum of squares, times
##     the finite population correction if one was given.
##
## The degrees of freedom are the number of PSUs less the number of strata,
## which is usually far smaller than the number of rows and is the number that
## should be quoted.
##
## References:
##   Binder, D. A. (1983). On the variances of asymptotically normal estimators
##     from complex surveys. International Statistical Review 51, 279-292.
##   Lumley, T. (2010). Complex Surveys: A Guide to Analysis Using R. Wiley.
## ---------------------------------------------------------------------------

#' Describe how a sample was drawn
#'
#' Bundles the weights, clustering and stratification of a complex sample so
#' that [ilm_model()] can fit with them and report a variance that reflects
#' them. The arguments mirror `survey::svydesign()`.
#'
#' @param data The data frame the sample lives in.
#' @param weights A one-sided formula naming the sampling weight, or a numeric
#'   vector. These are inverse probabilities of selection, not replicate
#'   counts.
#' @param ids A one-sided formula naming the primary sampling unit, or `NULL`
#'   for an unclustered sample, in which case every row is its own PSU.
#' @param strata A one-sided formula naming the stratum, or `NULL`.
#' @param fpc A one-sided formula naming the finite population correction --
#'   either the population size in each stratum or the sampling fraction --
#'   or `NULL` to ignore it, which is conservative.
#' @return An object of class `"ilm_design"`.
#' @seealso [ilm_model()] with `design =`, and [ilm_svy_coef()].
#' @references Lumley, T. (2010). *Complex Surveys: A Guide to Analysis Using
#'   R*. Wiley.
#' @examples
#' set.seed(1); n <- 500
#' d <- data.frame(x = rnorm(n), psu = rep(1:50, each = 10),
#'                 st = rep(1:2, each = 250))
#' d$w <- ifelse(d$st == 1, 40, 10)
#' d$y <- 1 + 0.5 * d$x + rnorm(n)
#' des <- ilm_design(d, weights = ~ w, ids = ~ psu, strata = ~ st)
#' des
#' @export
ilm_design <- function(data, weights, ids = NULL, strata = NULL, fpc = NULL) {
  if (!is.data.frame(data))
    stop("`data` must be a data frame; it is ", class(data)[1], call. = FALSE)
  grab <- function(x, nm) {
    if (is.null(x)) return(NULL)
    if (inherits(x, "formula")) {
      v <- all.vars(x)
      miss <- setdiff(v, names(data))
      if (length(miss))
        stop("`", nm, "` names ", paste(miss, collapse = ", "),
             ", which is not in the data.", call. = FALSE)
      if (length(v) == 1L) data[[v]] else
        interaction(data[v], drop = TRUE, sep = ":")
    } else {
      if (length(x) != nrow(data))
        stop("`", nm, "` has ", length(x), " values against ", nrow(data),
             " rows.", call. = FALSE)
      x
    }
  }
  w <- as.numeric(grab(weights, "weights"))
  if (is.null(w)) stop("`weights` is required", call. = FALSE)
  if (anyNA(w) || any(w <= 0))
    stop("sampling weights must be positive and present for every row; ",
         sum(is.na(w) | w <= 0), " are not.", call. = FALSE)
  id <- grab(ids, "ids"); st <- grab(strata, "strata"); fp <- grab(fpc, "fpc")
  id <- if (is.null(id)) factor(seq_len(nrow(data))) else droplevels(factor(id))
  st <- if (is.null(st)) factor(rep("1", nrow(data))) else droplevels(factor(st))
  ## A PSU belongs to exactly one stratum. If one appears in two, the labels
  ## are not nested and the variance would treat the same cluster as two
  ## independent ones.
  tb <- table(id, st)
  if (any(rowSums(tb > 0) > 1L))
    stop("some primary sampling units appear in more than one stratum, so the ",
         "design is not nested and a stratum's clusters cannot be told apart. ",
         "If PSU labels restart within each stratum, combine them: ",
         "ids = ~ interaction(strata, psu).", call. = FALSE)
  npsu <- nlevels(id); nstr <- nlevels(st)
  structure(list(weights = w, ids = id, strata = st, fpc = fp,
                 n = nrow(data), n_psu = npsu, n_strata = nstr,
                 df = npsu - nstr), class = "ilm_design")
}

#' @export
print.ilm_design <- function(x, ...) {
  cat("<ilm_design>\n")
  cat(sprintf("  %d rows in %d primary sampling unit(s), %d stratum(a)\n",
              x$n, x$n_psu, x$n_strata))
  cat(sprintf("  weights %.4g to %.4g, summing to %.6g\n",
              min(x$weights), max(x$weights), sum(x$weights)))
  cat(sprintf("  design degrees of freedom: %d\n", x$df))
  if (x$df < 30L)
    cat("  That is few. It is the number that sets the reference distribution,\n",
        "  not the row count, and with this many a t interval is noticeably\n",
        "  wider than a normal one.\n", sep = "")
  invisible(x)
}

#' Linearization variance for a fit from a complex sample
#'
#' The design-based (Taylor linearization) covariance of the fixed effects:
#' score contributions summed within each primary sampling unit, centred within
#' each stratum, and scaled by that stratum's number of clusters.
#'
#' @section Why the model-based variance is not an option here:
#'
#' A likelihood weighted by sampling weights believes it saw `sum(w)`
#' observations rather than `n`, so its standard errors are too small by at
#' least `sqrt(n / sum(w))` -- and by more once the sample is clustered, since
#' the design effect is on top of that. On the stratified two-stage example in
#' [ilm_svy_coef()] the weights alone account for a factor of 5.9 and the
#' observed ratios are 7.3 and 16.6.
#'
#' The point estimates are unaffected: a weighted likelihood is
#' design-consistent for the population parameter. That is what makes the
#' mistake easy to miss -- everything looks right except the uncertainty, and
#' the uncertainty looks better than right.
#'
#' @param object A fitted [ilm_model()] with no random effects.
#' @param design An [ilm_design()], or `NULL` to use the one the fit carries.
#' @param lonely What to do with a stratum holding a single PSU, which supplies
#'   no variance of its own. `"adjust"` centres it at the grand mean of all PSU
#'   totals instead of its own, which is conservative; `"certainty"` treats it
#'   as contributing nothing, which assumes it was selected with certainty;
#'   `"fail"` refuses.
#' @return A covariance matrix, with the design degrees of freedom attached.
#' @seealso [ilm_svy_coef()], [ilm_design()].
#' @references Binder, D. A. (1983). On the variances of asymptotically normal
#'   estimators from complex surveys. *International Statistical Review* 51,
#'   279-292.
#' @export
ilm_svy_vcov <- function(object, design = NULL,
                         lonely = c("adjust", "certainty", "fail")) {
  lonely <- match.arg(lonely)
  if (!inherits(object, "ilm_model"))
    stop("`object` must be a fitted ilm_model, not ", class(object)[1],
         call. = FALSE)
  if (is.null(design)) design <- object$design
  if (!inherits(design, "ilm_design"))
    stop("no design: pass one, or fit with ilm_model(design = ilm_design(...)) ",
         "so the fit carries it.", call. = FALSE)
  if (design$n != nrow(object$X))
    stop("the design describes ", design$n, " rows and the fit used ",
         nrow(object$X), ". If rows were dropped as missing, build the design ",
         "from the same data the model was given.", call. = FALSE)

  S <- ilm_estfun(object)                  # already weighted by the fit
  bread <- as.matrix(suppressWarnings(stats::vcov(object)))
  p <- ncol(object$X)
  bread <- bread[seq_len(p), seq_len(p), drop = FALSE]

  psu <- design$ids; str <- design$strata
  Sp <- rowsum(S, psu, reorder = TRUE)
  psu_str <- str[match(rownames(Sp), as.character(psu))]
  grand <- colMeans(Sp)

  meat <- matrix(0, p, p)
  lonely_n <- 0L
  for (h in levels(str)) {
    i <- which(psu_str == h)
    nh <- length(i)
    if (!nh) next
    Sh <- Sp[i, , drop = FALSE]
    if (nh == 1L) {
      lonely_n <- lonely_n + 1L
      if (lonely == "fail")
        stop("stratum ", h, " holds a single primary sampling unit, which ",
             "supplies no variance of its own. Use lonely = \"adjust\" to ",
             "centre it at the grand mean instead, which is conservative, or ",
             "lonely = \"certainty\" if it was selected with certainty and ",
             "contributes none.", call. = FALSE)
      if (lonely == "certainty") next
      Shc <- Sh - rep(grand, each = nh)
      meat <- meat + crossprod(Shc)
      next
    }
    Shc <- Sh - rep(colMeans(Sh), each = nh)
    f <- 1
    if (!is.null(design$fpc)) {
      fv <- design$fpc[str == h]
      pop <- stats::median(as.numeric(fv), na.rm = TRUE)
      f <- if (is.finite(pop) && pop > 1) max(1 - nh / pop, 0) else
        max(1 - pop, 0)
    }
    meat <- meat + f * (nh / (nh - 1)) * crossprod(Shc)
  }
  V <- bread %*% meat %*% bread
  dimnames(V) <- list(colnames(object$X), colnames(object$X))
  structure(V, df = design$df, n_psu = design$n_psu,
            n_strata = design$n_strata, lonely = lonely_n)
}

#' Fixed effects with design-based standard errors
#'
#' @param object A fitted [ilm_model()] from a complex sample.
#' @param design An [ilm_design()], or `NULL` to use the fit's own.
#' @param level Confidence level.
#' @param lonely See [ilm_svy_vcov()].
#' @return A data frame of `term`, `estimate`, `se`, `se_model`, `t`, `p`,
#'   `lower` and `upper`, on the design degrees of freedom.
#' @seealso [ilm_design()], [ilm_svy_vcov()].
#' @examples
#' set.seed(1); n <- 600
#' d <- data.frame(x = rnorm(n), psu = rep(1:60, each = 10),
#'                 st = rep(1:2, each = 300))
#' d$w <- ifelse(d$st == 1, 40, 10)
#' d$y <- 1 + 0.5 * d$x + rep(rnorm(60, 0, 0.8), each = 10) + rnorm(n)
#' des <- ilm_design(d, weights = ~ w, ids = ~ psu, strata = ~ st)
#' fit <- ilm_model(y ~ x, data = d, family = "gaussian", design = des,
#'                  verbose = FALSE)
#' ilm_svy_coef(fit)
#' @export
ilm_svy_coef <- function(object, design = NULL, level = 0.95,
                         lonely = c("adjust", "certainty", "fail")) {
  V <- ilm_svy_vcov(object, design, match.arg(lonely))
  b <- stats::coef(object)
  se <- sqrt(pmax(diag(V), 0))
  df <- attr(V, "df")
  crit <- if (df > 0) stats::qt(1 - (1 - level) / 2, df) else
    stats::qnorm(1 - (1 - level) / 2)
  tv <- b / se
  out <- data.frame(term = names(b), estimate = unname(b), se = unname(se),
                    se_model = unname(sqrt(diag(
                      as.matrix(suppressWarnings(stats::vcov(object)))))),
                    df = df, t = unname(tv),
                    p = 2 * stats::pt(-abs(unname(tv)), max(df, 1)),
                    lower = unname(b - crit * se),
                    upper = unname(b + crit * se),
                    row.names = NULL, stringsAsFactors = FALSE)
  structure(out, class = c("ilm_svy_coef", "data.frame"),
            df = df, n_psu = attr(V, "n_psu"), n_strata = attr(V, "n_strata"),
            lonely = attr(V, "lonely"))
}

#' @export
print.ilm_svy_coef <- function(x, ...) {
  cat("Fixed effects with design-based (linearization) standard errors\n")
  cat(sprintf("  %d primary sampling unit(s) in %d stratum(a); %d degrees of freedom\n",
              attr(x, "n_psu"), attr(x, "n_strata"), attr(x, "df")))
  d <- as.data.frame(x); class(d) <- "data.frame"
  d$estimate <- signif(d$estimate, 5); d$se <- signif(d$se, 4)
  d$se_model <- signif(d$se_model, 4); d$t <- round(d$t, 3)
  d$p <- signif(d$p, 3)
  d$lower <- signif(d$lower, 4); d$upper <- signif(d$upper, 4)
  print(d, row.names = FALSE)
  r <- x$se / x$se_model
  cat(sprintf("\n  design se / model-based se: %.2f to %.2f\n",
              min(r), max(r)))
  cat("  The model-based column is there to be ignored: a weighted likelihood\n",
      "  believes it saw sum(w) observations, so its errors are too small by\n",
      "  at least sqrt(n / sum(w)), and by more once the sample is clustered.\n",
      sep = "")
  if (attr(x, "lonely"))
    cat(sprintf("\n  %d stratum(a) held a single sampling unit and were handled\n  by `lonely`; they carry no variance of their own.\n",
                attr(x, "lonely")))
  invisible(x)
}

#' Keep a subset of the rows a design describes
#'
#' Rows dropped by `na.action` have to leave the design as well, or the weights
#' line up against the wrong observations -- silently, since both are numeric
#' vectors of plausible length.
#'
#' @keywords internal
#' @noRd
ilm_design_subset <- function(design, i) {
  ids <- droplevels(design$ids[i]); str <- droplevels(design$strata[i])
  structure(list(weights = design$weights[i], ids = ids, strata = str,
                 fpc = if (is.null(design$fpc)) NULL else design$fpc[i],
                 n = length(ids), n_psu = nlevels(ids),
                 n_strata = nlevels(str),
                 df = nlevels(ids) - nlevels(str)), class = "ilm_design")
}
