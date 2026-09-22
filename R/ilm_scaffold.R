## illume: a model you can reason about before you have any data.
##
## WHY THIS EXISTS: every simulation-based tool in this package -- ilm_power(),
## ilm_simulate(), ilm_emmeans(), ilm_trends() -- takes a FITTED model, and a
## study being planned has nothing to fit. The workaround is to invent a small
## data set, fit it, and overwrite the parameters with the assumptions you
## actually want. That works, and it is exactly what a separate project built
## against illume ended up doing, and it is where set_coef()'s dispersion bug
## was found: the workaround is easy to get subtly wrong and nothing tells you.
##
## So this does it once, properly. ilm_scaffold() returns an ordinary
## "ilm_model" whose parameters are the assumptions rather than estimates, and
## whose data is a design grid rather than observations. Everything downstream
## then works unchanged, which is the point: a user can ask what their
## assumptions IMPLY -- the marginal means, the slopes, a simulated data set --
## before committing to a power calculation built on them.
##
## The object is marked, and its print method says plainly that nothing here
## was estimated, because an "ilm_model" that was never fitted to real data is
## exactly the sort of thing that gets mistaken for a result.

#' A model built from assumptions instead of data
#'
#' Builds a design grid from a study specification, attaches the parameters you
#' assume, and returns a fitted-shaped `"ilm_model"`. Use it to plan a study:
#' pass it to [ilm_power()] for a power curve, to [ilm_simulate()] for data a
#' pipeline can be rehearsed on, or to [ilm_emmeans()] and [ilm_trends()] to
#' see what your assumptions actually imply before you rely on them.
#'
#' @section Saying what you assume:
#'
#' Give **either** `coefs` or `cells`, not both.
#'
#' `coefs` is a named vector of regression coefficients, on the link scale,
#' named as the model matrix names them (`"armtreatment"`, not `"arm"`). It is
#' exact and works for every model, including continuous predictors and
#' interactions with them.
#'
#' `cells` is the expected mean in each cell of the design, on the **response**
#' scale, which is how most people actually hold an assumption: "controls
#' average 12, the treated group averages 14.5". illume solves back to
#' coefficients. If the cell means you give cannot be produced by the formula
#' you gave -- crossed means under an additive formula, say -- it says so and
#' names the term that is missing, rather than quietly fitting the closest
#' thing it can.
#'
#' @section Saying how big the study is:
#'
#' `n_unit` counts **independent units**: participants when the formula has a
#' grouping bar, rows when it does not. A variable named in `within` is crossed
#' inside each unit; everything else is allocated across units as evenly as the
#' numbers allow. A variable that appears in a bar's left-hand side -- the
#' `time` of `(1 + time | id)` -- is within by construction and does not need
#' naming.
#'
#' @section What is NOT here:
#'
#' The standard errors on a scaffold come from one realisation of the design at
#' its own size. They are not a property of your assumptions, and reading power
#' off them would be reading one coin flip. That is what [ilm_power()] is for:
#' it refits many simulated studies and counts.
#'
#' @param formula A model formula, with random-effect bars if the design has
#'   repeated measures.
#' @param design Named list describing each variable in the design. A character
#'   vector becomes a factor with those levels; a numeric vector of length two
#'   or more becomes numeric values to cross; a function is called with the
#'   number of values needed and must return that many, which is the escape
#'   hatch for any distribution you like.
#' @param n_unit Integer. Independent units: participants if the formula has a
#'   grouping bar, otherwise rows.
#' @param family Family name or an [ilm_family()] object.
#' @param coefs Named numeric vector of coefficients on the link scale.
#' @param cells Expected cell means on the response scale: either a named
#'   numeric vector whose names are the factor levels joined by `"."`, or a
#'   data frame with one column per design factor and a column `mean`.
#' @param sd **Residual** standard deviation, for families that have one --
#'   not the standard deviation of the outcome. The random effects sit on top
#'   of it, so with `icc` the outcome's own spread is `sd / sqrt(1 - icc)`:
#'   `sd = 6, icc = 0.5` is data whose standard deviation is 8.49. The print
#'   method reports both. Superpower and faux take the **total** standard
#'   deviation, and so does a paper you might read one off, so converting one
#'   of those is `sd = sd_total * sqrt(1 - icc)`.
#' @param re_sd Random-effect standard deviations, named by grouping factor.
#'   For a term with a random slope, give a vector of standard deviations in
#'   the order the bar lists them; see `re_cor`.
#' @param re_cor Correlation between the random effects within a term, named by
#'   grouping factor. A single number for a two-column term, or a correlation
#'   matrix. Defaults to zero.
#' @param icc Intraclass correlation, as an alternative to `re_sd` for a random
#'   intercept: the random-effect variance becomes `icc / (1 - icc)` times the
#'   residual variance.
#' @param within Character vector of design variables that vary within a unit.
#' @param contrasts Passed to [ilm_model()]; also used when solving `cells`,
#'   so the two always agree.
#' @param seed Integer seed for the grid's own realisation.
#' @param verbose Logical. Report what was built.
#' @return An `"ilm_model"` that also carries class `"ilm_scaffold"`.
#' @seealso [ilm_power_design()] for the common case in one call;
#'   [ilm_power()], [ilm_simulate()], [ilm_emmeans()], [ilm_trends()].
#' @examples
#' \donttest{
#' ## a parallel-arm trial, stated as cell means
#' s <- ilm_scaffold(y ~ arm, design = list(arm = c("control", "treatment")),
#'                   n_unit = 120, cells = c(control = 12, treatment = 14.5),
#'                   sd = 4, verbose = FALSE)
#' coef(s)
#'
#' ## a treatment-by-time design with repeated measures
#' s2 <- ilm_scaffold(y ~ arm * time + (1 | id),
#'                    design = list(arm = c("control", "treatment"),
#'                                  time = c("pre", "post")),
#'                    within = "time", n_unit = 60,
#'                    cells = c(control.pre = 12, control.post = 12.2,
#'                              treatment.pre = 12, treatment.post = 14.5),
#'                    sd = 4, icc = 0.5, verbose = FALSE)
#' ilm_emmeans(s2, c("arm", "time"))
#' }
#' @export
ilm_scaffold <- function(formula, design, n_unit, family = "gaussian",
                         coefs = NULL, cells = NULL, sd = NULL,
                         re_sd = NULL, re_cor = NULL, icc = NULL,
                         within = NULL, contrasts = NULL, seed = 1L,
                         verbose = TRUE) {
  if (!inherits(formula, "formula"))
    stop("`formula` must be a formula, not ", class(formula)[1], call. = FALSE)
  if (!is.list(design) || !length(design) || is.null(names(design)) ||
      any(!nzchar(names(design))))
    stop("`design` must be a named list, one entry per variable in the ",
         "design; it is ", class(design)[1], call. = FALSE)
  if (is.null(coefs) == is.null(cells))
    stop("give exactly one of `coefs` (coefficients on the link scale) and ",
         "`cells` (expected means on the response scale).", call. = FALSE)
  n_unit <- as.integer(n_unit)[1]
  if (is.na(n_unit) || n_unit < 4L)
    stop("`n_unit` below 4 is not a study; it is ", n_unit, call. = FALSE)
  fam <- if (inherits(family, "ilm_family")) family else ilm_family(family)
  if (!requireNamespace("lme4", quietly = TRUE))
    stop("the formula interface needs lme4 to read random-effect bars; ",
         "install.packages(\"lme4\")", call. = FALSE)

  bars  <- lme4::findbars(formula)
  fform <- lme4::nobars(formula)
  resp  <- all.vars(formula)[1L]
  group <- if (length(bars)) ilm_mf_name(bars[[1L]][[3L]]) else NULL
  ## a variable a bar varies over is within-unit whether or not it was named
  slope_vars <- unlist(lapply(bars, function(b) all.vars(b[[2L]])))
  within <- union(as.character(within), intersect(slope_vars, names(design)))
  ## The grouping factor is BUILT from n_unit -- it is the thing being counted,
  ## not a variable with levels to allocate. Listing it in `design` would have
  ## been silently overwritten, which is a slow way to find out that n_unit is
  ## what decides how many participants there are.
  if (length(bars) && group %in% names(design))
    stop("'", group, "' is the grouping factor, so it is built from `n_unit` ",
         "rather than described in `design`. Drop it from `design` and set ",
         "`n_unit` to the number of ", group, "s you are planning for.",
         call. = FALSE)

  set.seed(seed)
  grid <- ilm_scaffold_grid(design, n_unit, within, group)
  N <- nrow(grid)

  ## ---- the coefficients ----------------------------------------------------
  ## the response does not exist yet -- that is the whole point -- so the
  ## design is built from the right-hand side alone
  mt  <- stats::terms(fform, data = grid)
  mtr <- stats::delete.response(mt)
  Xf <- stats::model.matrix(mtr, stats::model.frame(mtr, grid),
                            contrasts.arg = contrasts)
  bet <- ilm_scaffold_beta(coefs, cells, Xf, design, grid, fam, contrasts, mt)

  ## ---- the variance parameters --------------------------------------------
  sdv <- ilm_scaffold_sd(sd, fam)
  rev <- ilm_scaffold_resd(re_sd, re_cor, icc, sdv, bars, group)

  ## ---- a first fit, only to get an object of the right shape ---------------
  ## It is fitted to data drawn crudely from the assumptions rather than to
  ## noise, because a mixed model fitted to noise lands on the boundary and a
  ## boundary fit is a bad thing to start rebuilding from.
  grid[[resp]] <- ilm_scaffold_y0(Xf, bet, grid, group, rev, sdv, fam)
  fit <- suppressWarnings(suppressMessages(
    ilm_model(formula, data = grid, family = fam$name, contrasts = contrasts,
              verbose = FALSE, restarts = 1L)))
  fit <- ilm_scaffold_impose(fit, bet, rev, sdv, fam)

  ## ---- and now draw the data properly, through the package's own simulator,
  ## so the grid the scaffold carries really is a sample from the parameters it
  ## claims. Refit and impose again: the refit is only there to give every
  ## downstream field a consistent object, and the second impose is what makes
  ## the parameters exact.
  grid[[resp]] <- ilm_power_draw(fit, grid, bet)
  fit2 <- try(suppressWarnings(suppressMessages(
    ilm_model(formula, data = grid, family = fam$name, contrasts = contrasts,
              verbose = FALSE, restarts = 1L))), silent = TRUE)
  if (inherits(fit2, "try-error"))
    ## the PARAMETERS are still exactly what was asked for, because the impose
    ## above set them; it is the grid's response column that is the cruder
    ## draw. Worth saying, because a model this design cannot fit once is a
    ## model ilm_power() is about to try to fit several hundred times.
    warning("the scaffold's model would not refit on its own simulated data: ",
            sub("^Error[^:]*: ", "", conditionMessage(attr(fit2, "condition"))),
            " The assumed parameters are unaffected, but a power curve on ",
            "this design will lose replicates -- check `converged` in the ",
            "result.", call. = FALSE)
  else fit <- ilm_scaffold_impose(fit2, bet, rev, sdv, fam)

  fit$scaffold <- list(design = design, within = within, n_unit = n_unit,
                       rows_per_unit = if (is.null(group)) 1L else N %/% n_unit,
                       group = group, cells = cells, coefs = bet,
                       re_sd = rev, sd = sdv, seed = seed)
  class(fit) <- c("ilm_scaffold", class(fit))
  if (verbose) print(fit)
  fit
}

## Build the design grid.
##
## The unit of allocation is the point. A between-unit variable is assigned to
## UNITS and then repeated down that unit's rows, because assigning it to rows
## would put a participant in both arms of the trial.
#' @keywords internal
#' @noRd
ilm_scaffold_grid <- function(design, n_unit, within, group) {
  design <- design[!vapply(design, is.null, TRUE)]
  wv <- intersect(within, names(design))
  bv <- setdiff(names(design), wv)
  ## a factor's levels are its own, in its own order; anything else keeps the
  ## order it was written in, so a cell name the user types lines up with the
  ## level the model matrix builds
  lev <- function(x) if (is.factor(x)) levels(x) else unique(x)

  ## the within-unit part is crossed in full inside every unit
  wgrid <- if (length(wv)) {
    w <- lapply(design[wv], function(x)
      if (is.function(x)) stop("a `within` variable cannot be a function: it ",
                               "has to take the same values in every unit, ",
                               "or the units are not comparable.",
                               call. = FALSE) else lev(x))
    expand.grid(w, stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)
  } else NULL
  per <- if (is.null(wgrid)) 1L else nrow(wgrid)

  ## between-unit variables are allocated across units, balanced and shuffled
  ub <- data.frame(.unit = seq_len(n_unit))
  for (v in bv) {
    x <- design[[v]]
    ub[[v]] <- if (is.function(x)) {
      val <- x(n_unit)
      if (length(val) != n_unit)
        stop("`design$", v, "` returned ", length(val), " value(s) for ",
             n_unit, " unit(s); a design function must return as many values ",
             "as it is asked for.", call. = FALSE)
      val
    } else {
      l <- lev(x)
      sample(rep_len(l, n_unit))          # balanced to within one unit
    }
  }

  out <- if (is.null(wgrid)) ub else
    cbind(ub[rep(seq_len(n_unit), each = per), , drop = FALSE],
          wgrid[rep(seq_len(per), times = n_unit), , drop = FALSE])
  rownames(out) <- NULL
  if (!is.null(group)) out[[group]] <- factor(out$.unit)
  out$.unit <- NULL
  ## character columns become factors with the levels in the order given, so a
  ## cell name the user writes matches the level the model matrix builds
  for (v in names(design)) {
    x <- design[[v]]
    if (!is.function(x) && (is.character(x) || is.factor(x)))
      out[[v]] <- factor(as.character(out[[v]]), levels = lev(x))
  }
  out
}

## Resolve `coefs` or `cells` into a coefficient vector aligned to Xf.
#' @keywords internal
#' @noRd
ilm_scaffold_beta <- function(coefs, cells, Xf, design, grid, fam, contrasts,
                              mt) {
  nmX <- colnames(Xf)
  if (!is.null(coefs)) {
    if (is.null(names(coefs)))
      stop("`coefs` must be named. This model's coefficients are: ",
           paste(nmX, collapse = ", "), call. = FALSE)
    miss <- setdiff(nmX, names(coefs))
    extra <- setdiff(names(coefs), nmX)
    if (length(extra))
      stop("`coefs` names this model does not have: ",
           paste(extra, collapse = ", "), ". It has: ",
           paste(nmX, collapse = ", "), call. = FALSE)
    b <- stats::setNames(rep(0, length(nmX)), nmX)
    b[names(coefs)] <- as.numeric(coefs)
    if (length(miss) && !identical(miss, "(Intercept)"))
      message("coefficients not given are taken as zero: ",
              paste(setdiff(miss, "(Intercept)"), collapse = ", "))
    return(b)
  }

  ## ---- cell means ----------------------------------------------------------
  fv <- names(design)[vapply(design, function(x)
    !is.function(x) && (is.character(x) || is.factor(x)), TRUE)]
  fv <- intersect(fv, all.vars(mt))
  if (!length(fv))
    stop("`cells` needs at least one factor in the design to have cells; ",
         "give `coefs` instead.", call. = FALSE)
  cg <- expand.grid(lapply(grid[fv], function(x) levels(factor(x))),
                    stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)
  names(cg) <- fv
  ## anything else in the model sits at its average, so the cells describe the
  ## factors alone rather than the factors at one arbitrary covariate value
  for (v in setdiff(all.vars(mt), c(fv, all.vars(mt)[1L])))
    if (v %in% names(grid))
      cg[[v]] <- if (is.numeric(grid[[v]])) mean(grid[[v]]) else
        factor(levels(factor(grid[[v]]))[1L], levels = levels(factor(grid[[v]])))
  for (v in fv) cg[[v]] <- factor(cg[[v]], levels = levels(factor(grid[[v]])))

  key <- do.call(paste, c(cg[fv], list(sep = ".")))
  mu <- ilm_scaffold_cellvec(cells, cg, fv, key)
  lf <- ilm_scaffold_link(fam)
  bad <- which(!is.finite(suppressWarnings(lf(mu))))
  if (length(bad))
    stop("cell mean(s) '", paste(key[bad], collapse = "', '"),
         "' cannot be put on the ", fam$link, " scale: ",
         paste(mu[bad], collapse = ", "), ". ",
         if (identical(fam$link, "logit"))
           "`cells` is a proportion here, strictly between 0 and 1."
         else "`cells` is a mean on the response scale, and this link needs ",
         if (identical(fam$link, "logit")) "" else "it to be positive.",
         call. = FALSE)
  Xc <- stats::model.matrix(stats::delete.response(mt),
                            stats::model.frame(stats::delete.response(mt), cg),
                            contrasts.arg = contrasts)
  eta <- lf(mu)
  b <- tryCatch(qr.solve(Xc, eta), error = function(e)
    stop("the cell means could not be solved for coefficients: ",
         conditionMessage(e), call. = FALSE))
  ## A least-squares solve ALWAYS returns something. If the formula cannot
  ## reproduce the means that were asked for, the closest thing it can produce
  ## is a different study from the one being planned, so this stops.
  gap <- max(abs(as.numeric(Xc %*% b) - eta))
  if (gap > 1e-6 * max(1, max(abs(eta)))) {
    worst <- key[which.max(abs(as.numeric(Xc %*% b) - eta))]
    stop("this formula cannot produce the cell means given: cell '", worst,
         "' is off by ", format(gap, digits = 3), " on the link scale. ",
         "Cell means that differ by more than the terms in the formula allow ",
         "need the interaction between them -- ", paste(fv, collapse = " * "),
         " rather than ", paste(fv, collapse = " + "), ".", call. = FALSE)
  }
  stats::setNames(as.numeric(b), colnames(Xc))
}

## The forward link. The family objects carry `linkinv` because that is what
## fitting and prediction need; going the other way is only needed here, where
## the user states a mean and illume has to find the coefficient.
#' @keywords internal
#' @noRd
ilm_scaffold_link <- function(fam) {
  switch(fam$link,
         identity = function(m) m,
         log      = function(m) log(m),
         logit    = function(m) log(m / (1 - m)),
         probit   = function(m) stats::qnorm(m),
         cloglog  = function(m) log(-log(1 - m)),
         inverse  = function(m) 1 / m,
         sqrt     = function(m) sqrt(m),
         stop("`cells` is not supported for the ", fam$link, " link, because ",
              "illume has no way back from a mean to a coefficient here. ",
              "Give `coefs` on the link scale instead.", call. = FALSE))
}

## Read `cells` in either of the two forms into one vector aligned to cg.
#' @keywords internal
#' @noRd
ilm_scaffold_cellvec <- function(cells, cg, fv, key) {
  if (is.data.frame(cells)) {
    if (!"mean" %in% names(cells))
      stop("a `cells` data frame needs a column called `mean`; it has: ",
           paste(names(cells), collapse = ", "), call. = FALSE)
    if (!all(fv %in% names(cells)))
      stop("a `cells` data frame needs one column per design factor. ",
           "Missing: ", paste(setdiff(fv, names(cells)), collapse = ", "),
           call. = FALSE)
    k <- do.call(paste, c(lapply(cells[fv], as.character), list(sep = ".")))
    i <- match(key, k)
    if (anyNA(i))
      stop("`cells` does not cover every cell of the design. Missing: ",
           paste(key[is.na(i)], collapse = ", "), call. = FALSE)
    return(as.numeric(cells$mean[i]))
  }
  if (is.null(names(cells)))
    stop("`cells` must be named, one per cell of the design: ",
         paste(key, collapse = ", "), call. = FALSE)
  i <- match(key, names(cells))
  if (anyNA(i))
    stop("`cells` is missing ", sum(is.na(i)), " cell(s) of the design. ",
         "Expected names: ", paste(key, collapse = ", "),
         ". Got: ", paste(names(cells), collapse = ", "), call. = FALSE)
  as.numeric(cells)[i]
}

## The residual scale, checked against what the family actually has.
#' @keywords internal
#' @noRd
ilm_scaffold_sd <- function(sd, fam) {
  if (fam$n_disp == 0L) {
    if (!is.null(sd))
      stop("the ", fam$name, " family has no separate dispersion, so `sd` ",
           "has nothing to set: its spread follows from the mean.",
           call. = FALSE)
    return(NULL)
  }
  if (is.null(sd))
    stop("`sd` is needed: for a ", fam$name, " model the residual spread is ",
         "what decides power, and there is no data here to take it from.",
         call. = FALSE)
  sd <- as.numeric(sd)[1]
  if (!is.finite(sd) || sd <= 0)
    stop("`sd` must be a positive number; it is ", sd, call. = FALSE)
  sd
}

## Random-effect covariances, one per bar, as plain covariance matrices.
#' @keywords internal
#' @noRd
ilm_scaffold_resd <- function(re_sd, re_cor, icc, sdv, bars, group) {
  if (!length(bars)) {
    if (!is.null(re_sd) || !is.null(icc))
      stop("`re_sd`/`icc` describe a grouping factor, and this formula has ",
           "no random-effect bar. Add one, as in `+ (1 | id)`.", call. = FALSE)
    return(NULL)
  }
  if (!is.null(icc)) {
    if (!is.null(re_sd))
      stop("give `icc` or `re_sd`, not both: they set the same thing.",
           call. = FALSE)
    if (is.null(sdv))
      stop("`icc` is a ratio to the residual variance, so it needs `sd` too.",
           call. = FALSE)
    if (icc <= 0 || icc >= 1)
      stop("`icc` must be strictly between 0 and 1; it is ", icc, call. = FALSE)
    re_sd <- stats::setNames(list(sdv * sqrt(icc / (1 - icc))), group)
  }
  if (is.null(re_sd))
    stop("`re_sd` or `icc` is needed: this formula has a random-effect bar, ",
         "and how much units differ from each other is what decides power in ",
         "a design with repeated measures.", call. = FALSE)
  if (!is.list(re_sd)) re_sd <- as.list(re_sd)
  gv <- vapply(bars, function(b) ilm_mf_name(b[[3L]]), "")
  out <- vector("list", length(bars)); names(out) <- gv
  for (j in seq_along(bars)) {
    nm <- gv[j]
    d <- length(all.vars(bars[[j]][[2L]])) +
      as.integer(attr(stats::terms(ilm_one_sided(bars[[j]][[2L]])),
                      "intercept"))
    s <- re_sd[[nm]]
    if (is.null(s))
      stop("`re_sd` has nothing for the grouping factor '", nm, "'. It ",
           "names: ", paste(names(re_sd), collapse = ", "), call. = FALSE)
    s <- as.numeric(s)
    if (length(s) == 1L && d > 1L)
      stop("'", nm, "' has ", d, " random effects (", deparse(bars[[j]][[2L]]),
           ") but `re_sd` gives one number. Give ", d,
           ", in the order the bar lists them.", call. = FALSE)
    if (length(s) != d)
      stop("`re_sd$", nm, "` has ", length(s), " value(s) for ", d,
           " random effect(s).", call. = FALSE)
    if (any(!is.finite(s) | s <= 0))
      stop("`re_sd$", nm, "` must be positive; it is ",
           paste(s, collapse = ", "), call. = FALSE)
    R <- diag(d)
    if (!is.null(re_cor) && !is.null(re_cor[[nm]])) {
      rc <- re_cor[[nm]]
      R <- if (is.matrix(rc)) rc else {
        m <- diag(d); m[lower.tri(m)] <- rc; m[upper.tri(m)] <- t(m)[upper.tri(m)]; m
      }
      if (!identical(dim(R), c(d, d)))
        stop("`re_cor$", nm, "` is not ", d, " by ", d, call. = FALSE)
    }
    out[[j]] <- diag(s, d) %*% R %*% diag(s, d)
  }
  out
}

## Crude draw, used only so the first fit starts somewhere sensible.
#' @keywords internal
#' @noRd
ilm_scaffold_y0 <- function(Xf, bet, grid, group, rev, sdv, fam) {
  N <- nrow(Xf)
  eta <- as.numeric(Xf %*% bet)
  if (!is.null(group) && !is.null(rev) && group %in% names(grid)) {
    g <- factor(grid[[group]])
    V <- rev[[1L]]
    u <- matrix(stats::rnorm(nlevels(g) * nrow(V)), nlevels(g)) %*% ilm_msqrt(V)
    eta <- eta + u[as.integer(g), 1L]
  }
  mu <- fam$linkinv(eta)
  s <- if (is.null(sdv)) 1 else sdv
  switch(fam$name,
         gaussian = stats::rnorm(N, mu, s),
         poisson  = ,
         nbinom1  = ,
         nbinom2  = stats::rpois(N, pmax(mu, 1e-8)),
         binomial = stats::rbinom(N, 1L, pmin(pmax(mu, 1e-6), 1 - 1e-6)),
         Gamma    = stats::rgamma(N, shape = 2, scale = pmax(mu, 1e-8) / 2),
         beta     = stats::rbeta(N, pmax(mu, 1e-3) * 4,
                                 (1 - pmin(mu, 1 - 1e-3)) * 4),
         stats::rnorm(N, mu, s))
}

## Put the assumed parameters into a fitted object, exactly.
##
## Goes through the same rebuild set_coef() uses, so a scaffold and a
## marginaleffects nudge cannot drift apart: whatever one rebuilds, so does the
## other.
#' @keywords internal
#' @noRd
ilm_scaffold_impose <- function(fit, bet, rev, sdv, fam) {
  tl <- names(fit$opt$par)
  par <- as.numeric(fit$opt$par)
  nmX <- colnames(fit$X)
  ib <- which(tl == "beta")
  if (length(ib) != length(nmX))
    stop("the scaffold fit has ", length(ib), " fixed effects and the design ",
         "has ", length(nmX), "; this is a bug in ilm_scaffold().",
         call. = FALSE)
  par[ib] <- bet[match(nmX, names(bet))]
  if (!is.null(sdv) && any(tl == "logdisp")) {
    ## `sd` is the spread the DATA are to be generated with, and that is what
    ## $dispersion has to end up holding, because that is the field the
    ## simulator reads. An exact fit reports sigma on the unbiased (n - p)
    ## scale, so the parameter goes in a touch below the target for the
    ## reported value to land on it -- otherwise asking for sd = 4 on 120 rows
    ## quietly simulates at 4.03.
    ld <- log(sdv)
    if (isTRUE(fit$exact_df) && isTRUE(is.finite(fit$resid_df)))
      ld <- ld - 0.5 * log(nrow(fit$X) / fit$resid_df)
    par[tl == "logdisp"] <- ld
  }

  if (!is.null(rev) && any(tl == "theta")) {
    th <- numeric(0)
    for (k in seq_along(fit$re)) {
      if (identical(fit$re[[k]]$kind, "basis")) next
      V <- rev[[names(fit$re)[k]]]
      if (is.null(V)) V <- rev[[k]]
      dk <- if (is.null(fit$dk)) 1L else fit$dk[k]
      ## the scale of the whole term lives in Sigma and the shape in Sigma_d,
      ## whose leading entry is fixed at 1 for identifiability -- so the first
      ## random effect's variance IS the scale
      s2 <- V[1L, 1L]
      th <- c(th, log(sqrt(s2)))               # npc = 1 when C = 1
      if (dk > 1L) {
        Ld <- t(chol(V / s2))
        v <- numeric(0)
        for (j in 1:dk) for (i in j:dk) {
          if (i == 1L && j == 1L) next
          v <- c(v, if (i == j) log(Ld[i, j]) else Ld[i, j])
        }
        th <- c(th, v)
      }
    }
    if (length(th) == sum(tl == "theta")) par[tl == "theta"] <- th
    else warning("the scaffold could not match ", length(th), " variance ",
                 "parameter(s) to the ", sum(tl == "theta"), " this fit has; ",
                 "the random-effect assumptions were not imposed.",
                 call. = FALSE)
  }
  ilm_rebuild(fit, par)
}

#' @export
print.ilm_scaffold <- function(x, ...) {
  s <- x$scaffold
  cat("<ilm_scaffold>  parameters ASSUMED, not estimated\n")
  cat("  formula    : ", deparse(x$formula), "\n", sep = "")
  cat("  family     : ", x$family$name, "\n", sep = "")
  if (is.null(s$group))
    cat("  size       : ", s$n_unit, " rows\n", sep = "")
  else
    cat("  size       : ", s$n_unit, " ", s$group, "s x ", s$rows_per_unit,
        " = ", nrow(x$model), " rows\n", sep = "")
  bw <- setdiff(names(s$design), s$within)
  if (length(bw))
    cat("  between    : ", paste(bw, collapse = ", "), "\n", sep = "")
  if (length(s$within))
    cat("  within     : ", paste(s$within, collapse = ", "), "\n", sep = "")
  cat("\n  assumed coefficients\n")
  print(round(stats::coef(x), 4))
  lab <- function(t) sprintf("  %-12s: ", t)
  if (!is.null(s$sd))
    cat("\n", lab("residual sd"), format(s$sd, digits = 4), "\n", sep = "")
  if (!is.null(s$re_sd)) for (nm in names(s$re_sd))
    cat(lab(paste(nm, "sd")),
        paste(format(sqrt(diag(as.matrix(s$re_sd[[nm]]))), digits = 4),
              collapse = ", "), "\n", sep = "")
  ## `sd` is the RESIDUAL scale and the random effects sit on top of it, so
  ## the outcome's own spread is larger -- sd = 6 with icc = 0.5 is data whose
  ## standard deviation is 8.49. Every other package that takes a design
  ## specification (Superpower, faux) takes the TOTAL standard deviation, and
  ## a published SD read off a paper is the total one too. Printing both is
  ## the difference between a sample size and a wrong sample size.
  if (!is.null(s$sd) && !is.null(s$re_sd)) {
    dks <- vapply(s$re_sd, nrow, 0L)
    if (all(dks == 1L)) {
      tot <- sqrt(s$sd^2 + sum(vapply(s$re_sd, function(v) v[1L, 1L], 0)))
      cat(lab("total sd"), format(tot, digits = 4),
          "   <- the outcome's own spread; `sd` above is the residual\n",
          sep = "")
    } else {
      cat(lab("total sd"), "varies by row -- a random slope makes the\n",
          "                spread depend on where you are in the design\n",
          sep = "")
    }
  }
  cat("\n  The standard errors on this object come from ONE realisation of\n",
      "  the design at this size. For power, use ilm_power() or\n",
      "  ilm_power_design(), which refit many simulated studies.\n", sep = "")
  invisible(x)
}

#' Power for a study that has not been run
#'
#' Builds a scaffold from your assumptions and runs [ilm_power()] on it, in one
#' call. `n_unit` counts participants when the formula has a grouping bar and
#' rows when it does not, so the numbers you give are the numbers you would
#' write in a protocol.
#'
#' For anything beyond a power curve -- seeing what the assumptions imply,
#' simulating a data set, checking marginal means -- build the scaffold with
#' [ilm_scaffold()] and use it directly.
#'
#' @inheritParams ilm_scaffold
#' @param n_unit Integer vector of unit counts to trace power across.
#' @param term The term to test, named as the variable or as the coefficient;
#'   see [ilm_power()].
#' @param effect Values for that coefficient on the link scale. Defaults to
#'   whatever the assumptions imply for it.
#' @param sims,alpha,seed,progress Passed to [ilm_power()].
#' @return An [ilm_power()] result, with an extra `n_unit` column.
#' @seealso [ilm_scaffold()], [ilm_power()], [ilm_power_n()].
#' @examples
#' \donttest{
#' ilm_power_design(y ~ arm, design = list(arm = c("control", "treatment")),
#'                  n_unit = c(60, 120, 200),
#'                  cells = c(control = 12, treatment = 14.5), sd = 4,
#'                  term = "arm", sims = 50)
#' }
#' @export
ilm_power_design <- function(formula, design, n_unit, family = "gaussian",
                             coefs = NULL, cells = NULL, sd = NULL,
                             re_sd = NULL, re_cor = NULL, icc = NULL,
                             within = NULL, contrasts = NULL, term = NULL,
                             effect = NULL, sims = 200L, alpha = 0.05,
                             seed = 1L, progress = NULL, verbose = TRUE) {
  n_unit <- sort(unique(as.integer(n_unit)))
  if (any(is.na(n_unit)) || any(n_unit < 4L))
    stop("`n_unit` below 4 is not a study", call. = FALSE)
  ## the scaffold is built at the LARGEST size asked for, so every smaller one
  ## is a subsample of a grid that really contains every cell
  s <- ilm_scaffold(formula, design, n_unit = max(n_unit), family = family,
                    coefs = coefs, cells = cells, sd = sd, re_sd = re_sd,
                    re_cor = re_cor, icc = icc, within = within,
                    contrasts = contrasts, seed = seed, verbose = FALSE)
  rpu <- s$scaffold$rows_per_unit
  ## shown BEFORE the curve is simulated, not after: the simulation is the
  ## slow part, and what it assumes is what the user needs to see while it runs
  if (verbose) print(s)
  out <- ilm_power(s, n = n_unit * rpu, term = term, effect = effect,
                   sims = sims, alpha = alpha, seed = seed, progress = progress)
  out$n_unit <- out$n %/% rpu
  ## `[.data.frame` keeps names, row names and class and drops everything
  ## else, including the attributes print.ilm_power() and plot.ilm_power()
  ## read -- so reordering the columns the obvious way makes printing the
  ## result fail on its own output. They are put back by hand.
  keep <- attributes(out)
  keep <- keep[setdiff(names(keep), c("names", "row.names", "class"))]
  out <- out[, c("n_unit", setdiff(names(out), "n_unit")), drop = FALSE]
  for (nm in names(keep)) attr(out, nm) <- keep[[nm]]
  attr(out, "unit") <- if (is.null(s$scaffold$group)) "row" else
    s$scaffold$group
  attr(out, "rows_per_unit") <- rpu
  attr(out, "scaffold") <- s
  class(out) <- c("ilm_power", "data.frame")
  out
}
