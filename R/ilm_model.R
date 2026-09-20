## illume: lme4-style formula interface.
##
##   ilm_model(y ~ x1 + x2 + (1 | subj) + (1 + time | site) + s(x), data = dd)
##
## Random-effect bars are parsed with lme4::findbars; smooths with
## mgcv::interpret.gam, then reparameterised by ilm_smooth() into fixed
## null-space columns (appended to X) plus iid random basis blocks.
##
## This stores call / terms / xlev / contrasts / model frame / smooth objects.
## That is not bookkeeping for its own sake: emmeans::recover_data, any
## predict(newdata=) and car::Anova all rebuild a reference grid from exactly
## those components, and retrofitting them later is painful.

#' Fit a multinomial linear mixed model
#'
#' The main entry point. Takes an lme4-style formula and a data frame and fits a
#' mixed-effects model for a nominal categorical outcome with three or more
#' unordered categories.
#'
#' @details
#' Random-effect terms use lme4 syntax and smooths use mgcv syntax, so a model
#' can mix them freely:
#'
#' \preformatted{
#'   y ~ 1                                      intercept only
#'   y ~ x1 + grp                               ordinary regression
#'   y ~ x1 + grp + (1 | subj)                  random intercept
#'   y ~ x1 + time + (1 + time | subj)          random intercept and slope
#'   y ~ x1 + s(xs, k = 10) + (1 | subj)        penalised smooth of xs
#'   y ~ x1 + t2(lon, lat) + (1 | site)         2-D surface (spatial)
#' }
#'
#' Two rules for smooth terms. Use `t2()` rather than `te()` for tensor
#' products, because `te()` cannot be converted to the mixed-model form this
#' package relies on. And write smooths **unqualified** -- `s(x)`, not
#' `mgcv::s(x)` -- because mgcv identifies them by name, so a namespaced call
#' would be mistaken for an ordinary predictor. This matches `mgcv::gam()`.
#'
#' The outcome may be a factor or a character vector. Its levels set the category
#' labels used throughout the output.
#'
#' @section Simple models get exact inference:
#' A gaussian model with no random or smooth terms is an ordinary linear model.
#' In that case there is nothing to integrate out, so illume reports **exact**
#' t tests on `n - p` degrees of freedom and F tests in [ilm_anova()], rather
#' than the large-sample normal and chi-square approximations it must use when
#' random effects are present. Coefficients, standard errors, the residual
#' standard deviation and the F tests then agree with [stats::lm()] and
#' `car::Anova()` to numerical precision.
#'
#' This applies only where an exact reference genuinely exists. A Poisson or
#' binomial model without random effects is still a GLM, with no exact
#' small-sample analogue, so it keeps z and chi-square.
#'
#' @section Reading the coefficients:
#' Categories are coded **sum-to-zero**, so a coefficient is that category's
#' deviation from the average across categories, *not* a contrast against a
#' baseline. `summary()` prints a reminder, because this is easy to misread if
#' you are used to [nnet::multinom()].
#'
#' @section Always read the checks:
#' This model class fails quietly: a fit can return sensible-looking
#' coefficients while its covariance matrix is unusable, making the standard
#' errors meaningless. `summary()` prints the check verdicts for that reason,
#' and `fit$checks` holds the full table with a reason and a suggested remedy
#' for anything that is not `"OK"`.
#'
#' @param formula A formula with random-effect bars and optional smooth terms.
#' @param ... Arguments passed to the formula interface, listed below.
#' @param family Response distribution: one of "gaussian", "binomial",
#'   "poisson", "nbinom" or "multinomial". See [ilm_family()].
#' @param data A data frame.
#' @param re_struct Optional named list of category covariance structures, named
#'   by grouping variable. See [ilm_fit()].
#' @param ar Optional correlation over time, from [ilm_ar1()] or [ilm_car1()].
#' @param censor Optional censoring specification from [ilm_censor()], for a
#'   response with a floor, a ceiling or a detection limit.
#' @param rp_df Degrees of freedom for a flexible parametric baseline, used by
#'   `family = "rp"`, `"rp_odds"` and `"rp_normal"`. `1` is a straight line in
#'   log time, and so the corresponding parametric model; `3` is the usual
#'   default and allows two interior knots.
#' @param rp_knots Knot positions on the log-time scale, given directly in
#'   place of `rp_df`.
#' @param dispformula Optional one-sided formula for the dispersion, modelling
#'   its logarithm: `~ group` for a separate spread per level, `~ x` for one
#'   that changes with a covariate, `~ mu` for a power of the fitted mean.
#'   `mu` is a reserved name. This is the remedy for what
#'   [ilm_check_variance()] diagnoses.
#' @param weights Optional **frequency** weights: the number of replicate
#'   observations each row stands for. Evaluated inside `data`. See [ilm_fit()]
#'   for when this is valid, and why survey weights are not.
#' @param contrasts Optional contrasts for factor predictors, passed to
#'   [stats::model.matrix()]. Type III tests require orthogonal contrasts such
#'   as [stats::contr.sum()].
#' @param verbose Logical. Print checks while fitting.
#' @param restarts Integer. Optimiser restarts.
#' @param joint Logical or `NULL`. Compute the joint precision over fixed and
#'   random parameters. `NULL` (the default) switches it on when the model
#'   contains smooths, which is when [predict.ilm_model()] needs it.
#' @param na.action How to handle missing values; default [stats::na.omit()].
#'
#' @return An object of class `"ilm_model"`. Beyond the elements listed in
#'   [ilm_fit()], a formula fit also stores `call`, `terms`, `xlev`,
#'   `contrasts`, the model frame and the smooth objects -- everything needed to
#'   rebuild a reference grid for [predict.ilm_model()] and for `emmeans` or
#'   `marginaleffects`.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' n <- 600
#' dd <- data.frame(
#'   subj = factor(sample(40, n, TRUE)),
#'   x1   = rnorm(n),
#'   grp  = factor(sample(c("a", "b", "c"), n, TRUE))
#' )
#' dd$y <- factor(sample(c("low", "mid", "high"), n, TRUE),
#'                levels = c("low", "mid", "high"))
#'
#' fit <- ilm_model(y ~ x1 + grp + (1 | subj), data = dd)
#' summary(fit)
#' ilm_anova(fit, type = 3)
#' head(predict(fit))
#' }
#'
#' @seealso [summary.ilm_model()], [ilm_anova()], [predict.ilm_model()],
#'   [ilm_appraise()], [ilm_pb_lrt()].
#' @rdname ilm_model
#' @export
ilm_model <- function(formula, ...) {
  cl <- match.call()
  use_formula <- inherits(formula, "formula")
  fn_name <- if (use_formula) "ilm_model_formula" else "ilm_fit"
  fn <- if (use_formula) ilm_model_formula else ilm_fit
  cl[[1L]] <- as.name(fn_name)
  ## Evaluate in a CHILD of the caller's frame that also carries the target
  ## function.  Evaluating in the caller's frame alone would fail once the
  ## package is installed, because ilm_model_formula() is internal and the caller
  ## cannot see it -- a failure that does not appear under devtools::load_all(),
  ## where everything happens to sit in the global environment.
  env <- new.env(parent = parent.frame())
  assign(fn_name, fn, envir = env)
  eval(cl, env)
}

#' Formula front end
#'
#' Parses the formula, builds the design matrices and calls [ilm_fit()]. Called
#' by [ilm_model()]; documented separately only because [ilm_model()] dispatches to it.
#'
#' Random-effect bars are extracted with `lme4::findbars()` and smooths with
#' `mgcv::interpret.gam()`, then reparameterised by `ilm_smooth()` so that the
#' unpenalised part of each smooth joins the fixed effects and the penalised part
#' becomes a random term.
#'
#' A smooth is centred on the observed data, which changes what the INTERCEPT
#' means: it becomes the mean response at the sample average of the smooth,
#' not at a fixed point. That average moves from sample to sample, and the
#' standard error on the intercept does not carry that movement, so an interval
#' for it is narrower than its sampling spread. Measured on 2000 replicates of
#' a 600-row design, the intercept covered a fixed population value 88% of the
#' time and the sample-specific value 94.3%; the slopes were unaffected, which
#' is what the centring is for. Read the intercept of a model with a smooth as
#' a property of the sample, and take conclusions from the slopes and from
#' [ilm_plot_model()]'s effect curves.
#'
#' @return An object of class `"ilm_model"`.
#' @rdname ilm_model
#' @export
ilm_model_formula <- function(formula, data, family = "gaussian",
                         re_struct = NULL, ar = NULL,
                         weights = NULL, contrasts = NULL, verbose = TRUE,
                         restarts = 3L, joint = NULL, na.action = stats::na.omit,
                         censor = NULL, dispformula = NULL,
                         rp_df = 3L, rp_knots = NULL) {
  fam <- if (is.list(family)) family else ilm_family(family)
  cl <- match.call()
  if (!requireNamespace("lme4", quietly = TRUE)) stop("lme4 is required for the formula interface")
  if (!requireNamespace("mgcv", quietly = TRUE)) stop("mgcv is required for the formula interface")

  ## Where the terms of the formula get evaluated. lme4::nobars(),
  ## mgcv::interpret.gam() and reformulate() all hand back a formula carrying
  ## an environment of their own, so without this the model frame is built
  ## somewhere the caller's local variables do not exist -- and a term such as
  ## ilm_fourier(t, 12, K), or ns(x, df = d), written inside a function fails
  ## with "object 'K' not found".
  fenv <- environment(formula)
  if (is.null(fenv)) fenv <- parent.frame()

  bars  <- lme4::findbars(formula)
  fform <- lme4::nobars(formula)
  environment(fform) <- fenv
  gp    <- mgcv::interpret.gam(fform)
  environment(gp$pf) <- fenv
  if (!is.null(gp$fake.formula)) environment(gp$fake.formula) <- fenv
  smsp  <- gp$smooth.spec
  respn <- deparse(formula[[2]])

  ## one model frame for every part, so NA handling is consistent across
  ## fixed effects, smooths and grouping factors
  ## mgcv::interpret.gam() recognises smooths by the bare symbols s, te, t2, ti.
  ## A namespace-qualified call such as mgcv::s(x) is NOT recognised: it falls
  ## through as an ordinary parametric term and fails later with an obscure
  ## "invalid type (list)" error from model.frame().  Catch it here instead.
  ptl <- attr(stats::terms(gp$pf), "term.labels")
  qualified <- grep("^[A-Za-z0-9.]+::(s|te|t2|ti)[(]", ptl, value = TRUE)
  if (length(qualified))
    stop("smooth terms must be written unqualified, e.g. s(x) not ",
         qualified[1], ". mgcv detects smooths by name, so a namespaced call ",
         "is treated as an ordinary predictor.", call. = FALSE)

  bar_terms <- unlist(lapply(bars, function(b) c(deparse(b[[2]]), deparse(b[[3]]))))
  sm_terms  <- unlist(lapply(smsp, `[[`, "term"))
  rhs <- unique(c(attr(stats::terms(gp$pf), "term.labels"), sm_terms, bar_terms))
  rhs <- setdiff(rhs, c("1", "0", "-1"))
  ## Carry weights through the MODEL FRAME rather than evaluating them
  ## separately: that keeps them aligned when na.action drops rows, which a
  ## separate eval() would silently get wrong.  form_all is used only to build
  ## the frame -- X comes from gp$pf -- so the extra variable is harmless.
  wexpr <- cl$weights
  wnm <- if (is.null(wexpr)) NULL else deparse(wexpr)
  if (!is.null(wnm)) rhs <- c(rhs, wnm)
  ## The dispersion model's variables have to be in the same model frame as
  ## everything else, or na.action will drop different rows from each and the
  ## two designs will not line up.
  dvars <- character(0); disp_mu <- FALSE
  if (!is.null(dispformula)) {
    if (!inherits(dispformula, "formula") || length(dispformula) != 2L)
      stop("`dispformula` must be a one-sided formula, such as ~ group or ",
           "~ mu.", call. = FALSE)
    dtl <- attr(stats::terms(dispformula), "term.labels")
    disp_mu <- "mu" %in% dtl
    dvars <- setdiff(all.vars(dispformula), "mu")
    miss <- setdiff(dvars, names(data))
    if (length(miss))
      stop("`dispformula` refers to ", paste(miss, collapse = ", "),
           ", which ", if (length(miss) > 1L) "are" else "is",
           " not in the data. `mu` is the one reserved name, meaning the ",
           "fitted mean.", call. = FALSE)
    rhs <- unique(c(rhs, dvars))
  }
  if (!length(rhs)) rhs <- "1"
  form_all <- stats::reformulate(rhs, response = formula[[2]], env = fenv)
  mf <- stats::model.frame(form_all, data, na.action = na.action,
                           drop.unused.levels = TRUE)

  yraw <- stats::model.response(mf); N <- nrow(mf)
  if (fam$name == "multinomial") {
    yf <- factor(yraw)
    J <- nlevels(yf); ylevels <- levels(yf); yi <- as.integer(yf)
    if (J < 3L)
      stop("the multinomial family needs at least 3 outcome categories; ",
           "use family = \"binomial\" for 2", call. = FALSE)
  } else {
    ## A two-level factor, character or logical is a perfectly ordinary binomial
    ## response.  Characters are converted first, so that the same column behaves
    ## the same way whether it is stored as character or as a factor -- the
    ## multinomial path already converts, and falling through to as.numeric()
    ## here would silently produce NAs.
    if (is.character(yraw)) yraw <- factor(yraw)
    if (is.logical(yraw))   yraw <- factor(yraw, levels = c(FALSE, TRUE))
    if (is.factor(yraw)) {
      if (fam$name != "binomial")
        stop("family \"", fam$name, "\" needs a numeric response, but this one is ",
             "a factor with ", nlevels(yraw), " levels; use family = ",
             if (nlevels(yraw) == 2L) "\"binomial\"" else "\"multinomial\"",
             call. = FALSE)
      if (nlevels(yraw) != 2L)
        stop("the binomial family needs a response with 2 categories, but this ",
             "one has ", nlevels(yraw), "; use family = \"multinomial\"",
             call. = FALSE)
      ## second level is the modelled outcome, as in glm()
      ylevels <- levels(yraw)
      yi <- as.numeric(yraw) - 1
    } else { yi <- as.numeric(yraw); ylevels <- NULL }
    J <- NULL
  }

  ## ---- fixed effects ------------------------------------------------------
  mt <- stats::terms(gp$pf, data = mf)
  ## model.frame() is what runs makepredictcall(), and it is the only thing
  ## that does -- terms(formula, data = ) does not. Without carrying the result
  ## across, any term that stores state fitted from the data (ns, bs, poly,
  ## scale) gets recomputed by predict() from whatever rows it was handed, on a
  ## different basis. ns() errors outright on a short newdata; poly() returns
  ## numbers, quietly wrong ones. Copy the prediction calls model.frame()
  ## worked out onto the terms the model keeps.
  pvt <- attr(mf, "terms")
  pv <- attr(pvt, "predvars")
  if (!is.null(pv)) {
    src <- vapply(as.list(attr(pvt, "variables"))[-1], deparse1, "")
    want <- vapply(as.list(attr(mt, "variables"))[-1], deparse1, "")
    idx <- match(want, src)
    if (!anyNA(idx))
      attr(mt, "predvars") <- as.call(c(quote(list), as.list(pv)[-1][idx]))
  }
  X  <- stats::model.matrix(mt, mf, contrasts.arg = contrasts)
  xlev <- stats::.getXlevels(mt, mf)
  ctr  <- attr(X, "contrasts")
  ## cbind() below drops attributes, so carry the column -> term map by hand.
  ## Smooth null-space columns get NA: they belong to no parametric term.
  asgn <- attr(X, "assign")

  ## ---- smooths: null space -> X, penalised blocks -> basis terms ----------
  re_list <- list(); sm_store <- list()
  for (sp in smsp) {
    lab <- sp$label
    sob <- ilm_smooth(sp, mf)
    sm_store[[lab]] <- sob
    if (!is.null(sob$Xf) && ncol(sob$Xf)) {
      cn <- paste0(lab, ".f", seq_len(ncol(sob$Xf)))
      Xf <- sob$Xf; colnames(Xf) <- cn
      X <- cbind(X, Xf); asgn <- c(asgn, rep(NA_integer_, ncol(Xf)))
    }
    rb <- sob$rand
    for (i in seq_along(rb)) {
      nm <- if (length(rb) == 1L) lab else paste0(lab, ".", i)
      re_list[[nm]] <- list(basis = rb[[i]])
    }
  }

  ## ---- random-effect bars -------------------------------------------------
  for (b in bars) {
    gvar <- deparse(b[[3]])
    g <- mf[[gvar]]
    if (is.null(g)) stop("grouping factor '", gvar, "' not found in the data")
    Z <- stats::model.matrix(stats::as.formula(paste("~", deparse(b[[2]]))), mf)
    nm <- gvar; k <- 1L
    while (nm %in% names(re_list)) { k <- k + 1L; nm <- paste0(gvar, ".", k) }
    re_list[[nm]] <- list(group = g, Z = Z)
  }

  ## ---- fit ----------------------------------------------------------------
  ## joint defaults to TRUE when the model contains smooths.  Measured: holding
  ## the penalised coefficients fixed breaks their correlation with the smooth's
  ## null-space term and inflates interval width nearly threefold (0.358 vs
  ## 0.123), pushing intervals outside [0, 1).  Joint draws are the correct
  ## default wherever a smooth is present.
  if (is.null(joint)) joint <- length(smsp) > 0L
  w <- if (is.null(wnm)) NULL else as.numeric(mf[[wnm]])
  Zd <- NULL
  if (!is.null(dispformula)) {
    dfm <- stats::update(dispformula, ~ . )
    keep <- attr(stats::terms(dfm), "term.labels")
    keep <- setdiff(keep, "mu")
    f2 <- if (length(keep))
      stats::reformulate(keep, env = fenv) else stats::as.formula("~ 1", fenv)
    Zd <- stats::model.matrix(f2, mf)
    attr(Zd, "formula") <- dispformula
  }
  ## A flexible parametric baseline is a spline in log time, and log time is
  ## the response, so the basis is DATA: its columns simply join the model
  ## matrix. The spline coefficients are then part of beta, and the standard
  ## errors, the anova and the bootstrap all work on them unchanged.
  rp <- NULL
  if (isTRUE(fam$rp)) {
    tt <- as.numeric(yi)
    if (any(!is.na(tt) & tt <= 0))
      stop("a flexible parametric baseline is a spline in log(time), so every ",
           "time must be strictly positive; ", sum(tt <= 0, na.rm = TRUE),
           " are not.", call. = FALSE)
    ev <- if (is.null(censor)) rep(1L, length(tt))
          else as.integer(ilm_censor_for(censor, tt) == 0L)
    lt <- log(tt)
    kn <- if (!is.null(rp_knots)) sort(unique(as.numeric(rp_knots)))
          else ilm_rp_knots(lt, ev, as.integer(rp_df)[1])
    if (length(kn) < 2L)
      stop("`rp_knots` needs at least two knots, the boundaries", call. = FALSE)
    Bs <- ilm_rcs(lt, kn)
    Ds <- ilm_rcs(lt, kn, deriv = TRUE)
    ## the derivative design is zero in every column that is not the spline,
    ## so D %*% beta picks out exactly d(eta)/d(log t)
    Dfull <- cbind(Ds, matrix(0, nrow(X), ncol(X)))
    X <- cbind(Bs, X)
    colnames(Dfull) <- colnames(X)
    rp <- list(knots = kn, df = length(kn) - 1L, D = Dfull,
               cols = seq_len(ncol(Bs)), lo = min(lt) - 3, hi = max(lt) + 3)
    ## the spline columns belong to no formula term, the same as a smooth's
    ## null space, so they carry NA and no per-term test picks them up
    asgn <- c(rep(NA_integer_, ncol(Bs)), asgn)
  }
  fit <- ilm_fit(X, yi, J, re_list, re_struct = re_struct, ar = ar, censor = censor,
                 Zd = Zd, disp_mu = disp_mu, rp = rp,
                  ylevels = ylevels, weights = w, family = fam, verbose = verbose,
                  restarts = restarts, joint = joint)

  ## ---- everything the ecosystem layer reconstructs a reference grid from ---
  fit$call      <- cl
  fit$formula   <- formula
  fit$fixed_formula <- gp$pf
  fit$terms     <- mt
  fit$xlev      <- xlev
  fit$contrasts <- ctr
  fit$model     <- mf
  fit$smooths   <- sm_store
  fit$bars      <- bars
  fit$na.action <- attr(mf, "na.action")
  fit$ylevels   <- ylevels
  ## map each fixed-effect column back to its formula term, so per-term
  ## hypotheses can be blocked across the category dimension later
  fit$assign      <- asgn
  fit$term_labels <- attr(mt, "term.labels")
  fit
}

#' Standard accessors for a fitted model
#'
#' Supply the pieces that other packages expect from a model object, so that
#' `car::Anova()`, `emmeans` and `marginaleffects` can work with a fit.
#'
#' @param x,object,formula A fitted `"ilm_model"` object.
#' @param ... Unused.
#' @return The model terms, formula, model frame, model matrix or original call.
#' @rdname ilm_model-accessors
#' @export
terms.ilm_model        <- function(x, ...) x$terms

#' Standard accessors for a fitted model
#'
#' Supply the pieces that other packages expect from a model object, so that
#' `car::Anova()`, `emmeans` and `marginaleffects` can work with a fit.
#'
#' @param x,object,formula A fitted `"ilm_model"` object.
#' @param ... Unused.
#' @return The model terms, formula, model frame, model matrix or original call.
#' @rdname ilm_model-accessors
#' @export
formula.ilm_model      <- function(x, ...) x$formula

#' Standard accessors for a fitted model
#'
#' Supply the pieces that other packages expect from a model object, so that
#' `car::Anova()`, `emmeans` and `marginaleffects` can work with a fit.
#'
#' @param x,object,formula A fitted `"ilm_model"` object.
#' @param ... Unused.
#' @return The model terms, formula, model frame, model matrix or original call.
#' @rdname ilm_model-accessors
#' @export
model.frame.ilm_model  <- function(formula, ...) formula$model

#' Standard accessors for a fitted model
#'
#' Supply the pieces that other packages expect from a model object, so that
#' `car::Anova()`, `emmeans` and `marginaleffects` can work with a fit.
#'
#' @param x,object,formula A fitted `"ilm_model"` object.
#' @param ... Unused.
#' @return The model terms, formula, model frame, model matrix or original call.
#' @rdname ilm_model-accessors
#' @export
model.matrix.ilm_model <- function(object, ...) {
  X <- object$X
  attr(X, "assign") <- object$assign
  attr(X, "contrasts") <- object$contrasts
  X
}

#' Standard accessors for a fitted model
#'
#' Supply the pieces that other packages expect from a model object, so that
#' `car::Anova()`, `emmeans` and `marginaleffects` can work with a fit.
#'
#' @param x,object,formula A fitted `"ilm_model"` object.
#' @param ... Unused.
#' @return The model terms, formula, model frame, model matrix or original call.
#' @rdname ilm_model-accessors
#' @export
getCall.ilm_model <- function(x, ...) x$call
