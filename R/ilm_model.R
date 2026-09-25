## illume: lme4-style formula interface.
##
##   ilm_model(y ~ x1 + x2 + (1 | subj) + (1 + time | site) + s(x), data = dd)
##
## Random-effect bars are parsed with reformulas::findbars (lme4's own parser,
## which moved there); smooths with mgcv::interpret.gam, then reparameterised
## by ilm_smooth() into fixed null-space columns (appended to X) plus iid
## random basis blocks.
##
## This stores call / terms / xlev / contrasts / model frame / smooth objects.
## That is not bookkeeping for its own sake: emmeans::recover_data, any
## predict(newdata=) and car::Anova all rebuild a reference grid from exactly
## those components, and retrofitting them later is painful.

## The random-effect bars of a formula, and the formula without them, by
## lme4's own parser, which lme4 2.0 moved to reformulas. It is imported: every
## formula passes through it, bars or none, and with it only suggested a new
## install could not fit even y ~ x, which stopped with "the formula interface
## needs reformulas (or lme4)".
#' @keywords internal
#' @noRd
ilm_findbars <- function(f) reformulas::findbars(f)

#' @keywords internal
#' @noRd
ilm_nobars <- function(f) reformulas::nobars(f)

#' Fit generalized linear and additive mixed models
#'
#' The main entry point. Takes an lme4-style formula and a data frame and fits
#' a regression model for a gaussian, binomial, Poisson, negative binomial,
#' beta, multinomial or ordinal response, or a survival time (accelerated
#' failure time or Royston-Parmar) -- with or without random effects,
#' penalised smooths, AR(1) or CAR(1) correlation, a model for the dispersion,
#' zero inflation or a hurdle, and censoring at a floor or a ceiling. See
#' `family` for the full list.
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
#' A categorical outcome may be a factor or a character vector; its levels set
#' the category labels used throughout the output.
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
#' Every family but the multinomial has one coefficient per predictor, read as
#' in [stats::glm()]. A multinomial outcome has one per predictor per category,
#' and its categories are coded **sum-to-zero**, so a coefficient is that
#' category's deviation from the average across categories, *not* a contrast
#' against a baseline. `summary()` prints a reminder, because this is easy to
#' misread if you are used to [nnet::multinom()].
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
#'   "poisson", "nbinom", "beta", "multinomial", the ordinal families
#'   ("ordinal", "ordinal_probit", "ordinal_cloglog"), the accelerated failure
#'   time families ("weibull", "lognormal", "loglogistic") or the
#'   Royston-Parmar ones ("rp", "rp_odds", "rp_normal"). See [ilm_family()]. The default, `"auto"` (or `NULL`), reads
#'   the family off the response and says which it chose and why: a factor
#'   with 2 levels, a logical or a 0/1 variable is binomial; an unordered
#'   factor with 3 or more levels is multinomial and an ordered one ordinal;
#'   whole numbers that reach down to 0 or 1 are a count, and poisson; whole
#'   numbers that never come near zero, such as a blood pressure, are read
#'   as a measurement, and gaussian; values strictly between 0 and 1 are
#'   beta; anything else numeric is gaussian, as is a response censored at a
#'   floor or a ceiling. Where the response cannot settle it -- a proportion
#'   that touches 0 or 1 or comes with weights, a numeric variable with two
#'   values other than 0 and 1, a survival time, a date -- the fit stops and
#'   asks rather than guessing. The choice is a starting point, not a
#'   verdict: counts are often overdispersed (`"nbinom"`), and ratings on a
#'   short scale are often better read as ordinal. The chosen family is
#'   written into `fit$call`, so a refit uses it rather than guessing again.
#' @param data A data frame.
#' @param re_struct Optional named list setting the covariance structure of
#'   random terms, one element per term, named by its grouping variable as the
#'   formula names it; a term it leaves out keeps the default. Each element is
#'   itself a list:
#'   \describe{
#'     \item{`type`, and `rank` with `"rr"`}{the covariance across the
#'       outcome's categories, so it matters only for a multinomial model:
#'       `"us"` (the default, every variance and correlation free), `"diag"`
#'       (categories uncorrelated) or `"rr"` (reduced rank, with `rank` below
#'       the number of categories less one). Every other family has one linear
#'       predictor, where a rank has nothing to reduce.}
#'     \item{`d_cor = FALSE`}{drops the correlation between a random intercept
#'       and its slopes, in any family.}
#'   }
#'   For example `list(site = list(type = "rr", rank = 1))`, or
#'   `list(subj = list(d_cor = FALSE))`; see the examples, and [ilm_fit()] for
#'   what each costs. The pre-fit checks say when a structure is too rich for
#'   the data and name a rank to try.
#' @param ar Optional correlation over time, from [ilm_ar1()], [ilm_car1()] or
#'   [ilm_rw1()]. Written by name, as `ilm_car1(~ time | group)`, it takes the
#'   two columns from `data` after rows with missing values are dropped, so it
#'   cannot come out of step with the response, and the fit remembers them for
#'   [ilm_cells()] and for predictions on new rows.
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
#' @param ziformula Optional one-sided formula for the zero part of a count
#'   model, on the logit scale: `~ 1` for a constant excess-zero probability,
#'   `~ x` for one that depends on a predictor. This is the remedy for what
#'   [ilm_check_zeros()] diagnoses. Fixed effects only -- a random effect in
#'   the zero part is not supported, and a `|` here is an error rather than
#'   something quietly dropped.
#' @param zi_type `"inflated"` or `"hurdle"`, and they are different models.
#'   `"inflated"` is a mixture: some rows are structural zeros and the rest
#'   come from a count that may itself be zero, so a zero in the data could
#'   have come from either. `"hurdle"` is two processes: whether the response
#'   clears zero, and how far past it goes, the latter fitted to a count that
#'   cannot be zero. Choose by what the zeros mean, not by fit -- a structural
#'   zero is a unit that was never at risk. See [ilm_zi_coef()].
#' @param design An [ilm_design()] describing how a complex sample was drawn.
#'   Supplying one fits with the sampling weights and attaches the design, so
#'   [ilm_svy_coef()] can report a variance that reflects the clustering and
#'   stratification. Do not also pass `weights`: a sampling weight and a
#'   replicate count are different things and the design already carries one.
#' @param weights Optional **frequency** weights: the number of replicate
#'   observations each row stands for. Evaluated inside `data`. See [ilm_fit()]
#'   for when this is valid, and why survey weights are not.
#' @param contrasts How to code factor predictors. `NULL` (the default) uses
#'   R's own setting, which is treatment coding for unordered factors and
#'   polynomial for ordered ones -- the same as [stats::lm()], so a coefficient
#'   means what it means everywhere else in R. A single string (`"treatment"`,
#'   `"sum"`, `"helmert"` or `"poly"`, with or without the `contr.` prefix)
#'   applies that coding to every unordered factor at once. A named list sets
#'   them one factor at a time, exactly as [stats::lm()] takes it.
#'
#'   Type III tests of a main effect that is also in an interaction need
#'   orthogonal coding such as [stats::contr.sum()]; [ilm_anova()] checks for
#'   that and says so rather than reporting a test that is not the one it
#'   claims. Note that this is a separate matter from the sum-to-zero coding a
#'   multinomial fit uses across its outcome CATEGORIES, which is not a
#'   predictor contrast and is not affected by this argument.
#' @param verbose Logical. Print checks while fitting.
#' @param restarts Integer. Optimiser restarts.
#' @param joint Logical or `NULL`. Compute the joint precision over fixed and
#'   random parameters. `NULL` (the default) switches it on when the model
#'   contains smooths, which is when [predict.ilm_model()] needs it.
#' @param na.action How to handle missing values; default [stats::na.omit()].
#' @param reml Logical. Estimate the variance components by RESTRICTED maximum
#'   likelihood instead of maximum likelihood. Defaults to `FALSE`, and the
#'   reason is the order you work in: maximum likelihood is what lets you
#'   compare fixed-effect structures, because a restricted likelihood belongs
#'   to contrasts orthogonal to the design matrix and changing that matrix
#'   changes which data it is the likelihood of. Settle the fixed effects under
#'   the default, then refit with `reml = TRUE` for the estimates you report --
#'   maximum likelihood biases the variance components downward, and with few
#'   clusters that carries through to standard errors and to the degrees of
#'   freedom from [ilm_denom_df()]. With few groups the difference is
#'   measurable: in a study of 8 groups of 6, nominal 90% prediction intervals
#'   for new groups, built from the fitted variance components, covered 0.903
#'   with REML against 0.891 with maximum likelihood for a gaussian response,
#'   and 0.892 against 0.872 for a Poisson one. Once set, any likelihood-ratio
#'   test refuses
#'   rather than quietly comparing things that are not comparable, and so does
#'   [ilm_robust()], whose sandwich needs per-observation scores that a
#'   restricted likelihood does not have. Available for every family, but exact
#'   only for a gaussian response: for any other family it is an approximately
#'   restricted likelihood, which reduces the downward bias without REML's
#'   exact properties, and `fit$reml_exact` says which a fit has (see the
#'   *Regression models* vignette). [ilm_dag_model()] defaults to `TRUE` for a
#'   gaussian response, because there the graph fixed the adjustment set
#'   before any data were seen.
#' @param boundary What to do about a random-effect covariance at the edge of
#'   its range -- a variance of zero, or a correlation of +/-1 -- where the
#'   likelihood is flat and cannot say where in that direction the truth is.
#'   `"hold"`, the default, is maximum likelihood: an estimate that lands there
#'   is held at it in the direction that reached the edge, the rest of the
#'   fit's uncertainty is computed around it, and the fixed effects remain
#'   usable (see the BOUNDARY verdict in [summary.ilm_model()]). Their
#'   standard errors are then those of the reduced model the boundary
#'   implies -- a covariance of lower rank, or the term dropped -- on every
#'   platform. `"avoid"` adds the boundary-avoiding penalty of
#'   Chung et al. (2013, 2015) -- half the log-determinant of each grouping
#'   term's covariance -- which keeps every estimate strictly inside its
#'   range; in one dimension it is a gamma(2) prior on the standard deviation.
#'   The penalty is small against the likelihood, so it matters only where
#'   the data cannot resolve the covariance, and `logLik()` reports the
#'   likelihood of the data at the penalised estimate.
#'
#'   Measured against `"hold"` on a three-category outcome with 60 groups of
#'   8, 400 datasets per condition: with a true between-group SD of 0.05,
#'   every `"avoid"` fit was usable against 395 of 400 under `"hold"`, and the
#'   fixed effects' intervals covered at 0.949 against 0.946. The cost is in
#'   the estimates.
#'   A variance is pulled away from zero rather than estimated at it -- that
#'   SD of 0.05 came out at a median of 0.21, against 0.10 -- so a test of
#'   whether it IS zero no longer applies. And a larger between-group variance
#'   means larger within-group effects on a logit scale, so the fixed effects
#'   moved further from zero with it: by about 1% in four conditions of six,
#'   but by 11% when one outcome category was rare, where coverage fell from
#'   0.930 to 0.916, and by 20% when the groups were unbalanced and
#'   heavy-tailed as well. `"hold"` stays the default for that reason; the fit
#'   says when a boundary was reached under it, and names `"avoid"` as the
#'   alternative.
#'
#' @return An object of class `"ilm_model"`. Beyond the elements listed in
#'   [ilm_fit()], a formula fit also stores `call`, `terms`, `xlev`,
#'   `contrasts`, the model frame and the smooth objects -- everything needed to
#'   rebuild a reference grid for [predict.ilm_model()] and for `emmeans` or
#'   `marginaleffects`.
#'
#' @references
#' Chung, Y., Rabe-Hesketh, S., Dorie, V., Gelman, A., & Liu, J. (2013). A
#' nondegenerate penalized likelihood estimator for variance parameters in
#' multilevel models. *Psychometrika*, 78(4), 685--709.
#'
#' Chung, Y., Gelman, A., Rabe-Hesketh, S., Liu, J., & Dorie, V. (2015).
#' Weakly informative prior for point estimation of covariance matrices in
#' hierarchical models. *Journal of Educational and Behavioral Statistics*,
#' 40(2), 136--157.
#'
#' @examples
#' set.seed(1)
#' n <- 300
#' dd <- data.frame(subj = factor(sample(30, n, TRUE)), x1 = rnorm(n),
#'                  grp = factor(sample(c("a", "b", "c"), n, TRUE)))
#' dd$y <- 1 + 0.5 * dd$x1 + rnorm(30)[dd$subj] + rnorm(n)
#'
#' ## a linear mixed model: a random intercept for each subject
#' fit <- ilm_model(y ~ x1 + grp + (1 | subj), data = dd, family = "gaussian",
#'                  verbose = FALSE)
#' summary(fit)
#'
#' ## a random slope for x1 as well, without its correlation with the
#' ## intercept -- re_struct names the term by its grouping variable
#' fit2 <- ilm_model(y ~ x1 + (1 + x1 | subj), data = dd, family = "gaussian",
#'                   re_struct = list(subj = list(d_cor = FALSE)),
#'                   verbose = FALSE)
#'
#' \donttest{
#' ## a count, and a nominal outcome with four categories
#' dd$n_events <- rpois(n, exp(0.3 + 0.2 * dd$x1))
#' fit3 <- ilm_model(n_events ~ x1 + (1 | subj), data = dd, family = "poisson",
#'                   verbose = FALSE)
#' ## subjects that differ along one direction across the categories: more
#' ## "x" and "y" and less "w" and "z", in fixed proportions
#' u <- rnorm(30)[dd$subj]
#' eta <- cbind(-u, u + 0.4 * dd$x1, 0.5 * u, -0.5 * u)
#' dd$k <- factor(apply(exp(eta) / rowSums(exp(eta)), 1, function(p)
#'   sample(c("w", "x", "y", "z"), 1, prob = p)))
#' ## the subjects' random intercepts vary across the three category
#' ## dimensions; a reduced rank of 1 describes that covariance with one
#' ## dimension instead of three -- fewer parameters and fewer latent values
#' fit4 <- ilm_model(k ~ x1 + (1 | subj), data = dd, family = "multinomial",
#'                   re_struct = list(subj = list(type = "rr", rank = 1)),
#'                   verbose = FALSE)
#' ilm_anova(fit4, type = 3)
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
#' Random-effect bars are extracted with `findbars()` from reformulas (lme4's
#' parser; lme4 itself is used when reformulas is absent) and smooths with
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
#' @section How smooth terms are named:
#' A smooth takes the name mgcv gives the smooth it builds, and three things
#' are named from it:
#' * **the smooth itself**, as `fit$smooths` and the random terms list it:
#'   `"s(x)"` for `s(x)`, `"t2(x,z)"` for `t2(x, z)`, and `"s(x):z"` for a
#'   smooth with a numeric `by`, `s(x, by = z)`;
#' * **its unpenalised columns** among the fixed effects: the name with `.f1`,
#'   `.f2`, ... after it, so `"s(x).f1"`, or `"s(x):z.f1"` and `"s(x):z.f2"`;
#' * **its penalised part**, the random term whose variance `summary()`,
#'   [ilm_varcorr()] and `fit$Sigma` report: the name itself, or `"t2(x,z).1"`,
#'   `"t2(x,z).2"`, ... for a tensor product, which has a penalty per margin.
#'
#' Two smooths may not share a name. `s(x) + s(x, k = 5)` stops, because the
#' second would overwrite the first.
#'
#' **What changed.** A smooth with a numeric `by` used to take the name of the
#' same smooth without it: `"s(x)"`, with columns `"s(x).f1"` and
#' `"s(x).f2"`. Beside a plain `s(x)` the two names collided: one smooth's
#' penalised part was lost from the fit, and `predict()` on new rows stopped.
#' Its columns are now `"s(x):z.f1"` and `"s(x):z.f2"`, and its variance
#' `"s(x):z"`. Code that picked them out by the old names needs the new ones.
#' A smooth without a `by` is named as before.
#'
#' @return An object of class `"ilm_model"`.
#' @rdname ilm_model
#' @export
ilm_model_formula <- function(formula, data, family = "auto",
                         re_struct = NULL, ar = NULL,
                         weights = NULL, contrasts = NULL, verbose = TRUE,
                         restarts = 3L, joint = NULL, na.action = stats::na.omit,
                         censor = NULL, dispformula = NULL,
                         rp_df = 3L, rp_knots = NULL, ziformula = NULL,
                         zi_type = c("inflated", "hurdle"), design = NULL,
                         reml = FALSE, boundary = c("hold", "avoid")) {
  zi_type <- match.arg(zi_type)
  boundary <- match.arg(boundary)
  ## A survey design supplies the weights, so taking them from both places
  ## would silently apply one and ignore the other.
  if (!is.null(design)) {
    if (!inherits(design, "ilm_design"))
      stop("`design` must come from ilm_design(), not ", class(design)[1],
           call. = FALSE)
    if (!is.null(substitute(weights)))
      stop("give the weights to ilm_design() or to `weights`, not both: a ",
           "design already carries them, and they mean different things.",
           call. = FALSE)
  }
  ## "auto" is resolved once the response is in hand, below: the family is
  ## read off the response as the model will see it, after na.action
  auto <- is.null(family) || identical(family, "auto")
  fam <- if (auto) NULL else if (is.list(family)) family else ilm_family(family)
  cl <- match.call()

  ## Where the terms of the formula get evaluated. nobars(),
  ## mgcv::interpret.gam() and reformulate() all hand back a formula carrying
  ## an environment of their own, so without this the model frame is built
  ## somewhere the caller's local variables do not exist -- and a term such as
  ## ilm_fourier(t, 12, K), or ns(x, df = d), written inside a function fails
  ## with "object 'K' not found".
  fenv <- environment(formula)
  if (is.null(fenv)) fenv <- parent.frame()

  bars  <- ilm_findbars(formula)
  fform <- ilm_nobars(formula)
  environment(fform) <- fenv
  ## Not mgcv::interpret.gam() directly: it rebuilds formulas from text and
  ## fails outright on a column name that needs backticks. See ilm_names.R.
  gp    <- ilm_interpret_gam(fform)
  environment(gp$pf) <- fenv
  if (!is.null(gp$fake.formula)) environment(gp$fake.formula) <- fenv
  smsp  <- gp$smooth.spec

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
  ## A smooth with a numeric `by` is not centred, so its unpenalised part
  ## already spans the by variable itself: for s(x, by = z), z and z * x,
  ## whatever the basis. With z also a term of the model the two are one
  ## column, the fit has no unique answer, and it came back with a failed
  ## Hessian and every standard error NaN, saying nothing of why.
  if (is.data.frame(data)) for (sp in smsp) {
    bv <- sp$by
    if (is.null(bv) || identical(bv, "NA") || !bv %in% names(data) ||
        !is.numeric(data[[bv]]) || !bv %in% ilm_unbq(ptl)) next
    lab <- paste0(sub("\\)$", "", sp$label), ", by = ", bv, ")")
    stop("`", bv, "` is a term of the model and also the `by` variable of ",
         lab, ". A smooth with a numeric `by` is not centred, so its ",
         "unpenalised part already contains ", bv, "'s main effect: the two ",
         "are the same column, and the model as written has no unique fit ",
         "(every standard error would come out NaN). Drop `", bv, "` from ",
         "the formula; ", lab, " carries its effect. See Wood (2017), ",
         "Generalized Additive Models, 2nd ed., p. 326.", call. = FALSE)
  }

  ## all.vars() on the grouping side, not deparse(). A NESTED bar is expanded
  ## by findbars() into a grouping EXPRESSION rather than a name --
  ## `(1 | continent/country)` becomes `(1 | country:continent)` plus
  ## `(1 | continent)` -- and deparsing that yields "country:continent", which
  ## is not a column, so neither constituent ever reached the model frame.
  ##
  ## Every name in `rhs` becomes formula TEXT for reformulate(), so each one
  ## is quoted on the way in: all.vars() hands names back bare, and a column
  ## called `site id` would otherwise be parsed as two symbols.
  bar_terms <- unlist(lapply(bars, function(b)
    c(ilm_term_text(b[[2]]), ilm_bq(all.vars(b[[3]])))))
  ## Checked HERE, before model.frame() gets it. Those variables now go into
  ## the frame, so a missing one would otherwise surface as R's bare "object
  ## 'nope' not found" from inside eval(predvars) -- which does not say it was
  ## a grouping factor, does not name the bar, and does not say what to do.
  if (is.data.frame(data)) {
    gvars <- unique(unlist(lapply(bars, function(b) all.vars(b[[3]]))))
    miss <- setdiff(gvars, names(data))
    if (length(miss))
      stop("the grouping factor(s) ", paste(sQuote(miss), collapse = ", "),
           " could not be built from the data: ",
           if (length(miss) == 1L) "it is" else "they are",
           " not a column of `data`. The bars name ",
           paste(sQuote(gvars), collapse = ", "), ".", call. = FALSE)
  }
  ## A smooth's `by` variable is data the basis needs as much as its `term`
  ## is. Only the term used to be carried into the model frame, so a numeric
  ## `by` that appeared nowhere else in the formula never reached it and
  ## mgcv could not find it.
  sm_terms  <- ilm_bq(unlist(lapply(smsp, function(s)
    c(s$term, if (!is.null(s$by) && !identical(s$by, "NA")) s$by))))
  rhs <- unique(c(attr(stats::terms(gp$pf), "term.labels"), sm_terms, bar_terms))
  rhs <- setdiff(rhs, c("1", "0", "-1"))
  ## Carry weights through the MODEL FRAME rather than evaluating them
  ## separately: that keeps them aligned when na.action drops rows, which a
  ## separate eval() would silently get wrong.  form_all is used only to build
  ## the frame -- X comes from gp$pf -- so the extra variable is harmless.
  ## Two spellings of the one expression: quoted to go INTO the formula, and
  ## as model.frame() names the column, to be found in it afterwards.
  wexpr <- cl$weights
  wnm <- if (is.null(wexpr)) NULL else ilm_mf_name(wexpr)
  if (!is.null(wnm)) rhs <- c(rhs, ilm_term_text(wexpr))
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
    rhs <- unique(c(rhs, ilm_bq(dvars)))
  }
  zvars <- character(0)
  if (!is.null(ziformula)) {
    if (!inherits(ziformula, "formula") || length(ziformula) != 2L)
      stop("`ziformula` must be a one-sided formula, such as ~ 1 or ~ group.",
           call. = FALSE)
    if (any(grepl("|", as.character(ziformula), fixed = TRUE)))
      stop("`ziformula` takes fixed effects only. A random effect in the zero ",
           "part would need its own latent vector integrated out alongside the ",
           "count part's, which this fit does not do -- rather than drop the ",
           "bar and fit something you did not ask for, it stops here.",
           call. = FALSE)
    zvars <- all.vars(ziformula)
    miss <- setdiff(zvars, names(data))
    if (length(miss))
      stop("`ziformula` refers to ", paste(miss, collapse = ", "),
           ", which ", if (length(miss) > 1L) "are" else "is",
           " not in the data.", call. = FALSE)
    ## the zero part shares the model frame, so its rows are dropped by the
    ## same na.action and the two designs cannot come out of step
    rhs <- unique(c(rhs, ilm_bq(zvars)))
  }
  ## A correlation over time given by name, `~ time | group`, reads both
  ## columns from the model frame, so they go into it: na.action then drops a
  ## row from the structure and the response together, where vectors taken
  ## from the whole data frame beforehand fall out of step with it.
  if (inherits(ar, "ilm_cor_named")) {
    if (is.data.frame(data)) {
      miss <- setdiff(ar$vars, names(data))
      if (length(miss))
        stop("the correlation over time names ",
             paste(sQuote(miss), collapse = " and "), ", which ",
             if (length(miss) > 1L) "are not columns" else "is not a column",
             " of `data`.", call. = FALSE)
    }
    rhs <- unique(c(rhs, ilm_bq(ar$vars)))
  }
  if (!length(rhs)) rhs <- "1"
  form_all <- stats::reformulate(rhs, response = formula[[2]], env = fenv)
  mf <- stats::model.frame(form_all, data, na.action = na.action,
                           drop.unused.levels = TRUE)
  ## The variables the terms are built from, for the rows the frame kept. The
  ## frame holds a transformed term as its own column -- "log(x)", a Fourier
  ## basis -- and not the variable underneath, so everything that builds new
  ## rows from the fit's own (marginal means, average effects, scenarios)
  ## could not rebuild the term: "object 'x' not found", or, for a column
  ## called t or time, R's own t() and time() in its place. Only the columns
  ## the frame lacks are kept, matched to its rows by name.
  data_extra <- NULL
  if (is.data.frame(data)) {
    rv <- unique(c(all.vars(formula), zvars, dvars,
                   if (inherits(ar, "ilm_cor_named")) ar$vars))
    extra <- setdiff(intersect(rv, names(data)), names(mf))
    if (length(extra)) {
      idx <- match(rownames(mf), rownames(data))
      if (!anyNA(idx)) {
        data_extra <- as.data.frame(data)[idx, extra, drop = FALSE]
        rownames(data_extra) <- rownames(mf)
      }
    }
  }
  ## Dropping incomplete rows is the default and usually the right thing, but
  ## doing it SILENTLY is not: a model fitted to 61% of the data with no note
  ## of it invites conclusions the data cannot carry. Say how many went, and
  ## which columns took them, so the choice is visible at the point it is made.
  n_drop <- if (is.data.frame(data)) nrow(data) - nrow(mf) else 0L
  if (verbose && n_drop > 0L) {
    culprits <- names(mf)[vapply(names(mf), function(v)
      anyNA(data[[v]]), TRUE, USE.NAMES = FALSE)]
    culprits <- culprits[culprits %in% names(data)]
    message(sprintf(
      "%d of %d row(s) dropped for missing values (%.1f%%), leaving %d%s",
      n_drop, nrow(data), 100 * n_drop / nrow(data), nrow(mf),
      if (length(culprits))
        paste0("; missing in: ", paste(utils::head(culprits, 6), collapse = ", "))
      else ""))
    if (n_drop / nrow(data) > 0.1)
      message("  that is more than a tenth of the data. Complete cases stay ",
              "unbiased when missingness is unrelated to the OUTCOME given ",
              "the predictors, and are biased otherwise -- ilm_check_missing() ",
              "tells the two apart, and ilm_impute() is the remedy for the ",
              "second.")
  }

  if (inherits(ar, "ilm_cor_named")) ar <- ilm_cor_build(ar, mf)

  yraw <- stats::model.response(mf); N <- nrow(mf)
  ## family = "auto": chosen from the response, and SAID, because a likelihood
  ## chosen silently is exactly the kind of decision this package does not
  ## make behind the user's back. The choice is written into the stored call,
  ## so everything that refits through the call -- moderation searches, power
  ## by simulation -- refits the same family rather than guessing afresh on
  ## simulated data, where a low-rate count can come out all 0s and 1s.
  fam_why <- NULL
  if (auto) {
    g <- ilm_guess_family(yraw, ilm_mf_name(formula[[2]]),
                          zero_part = !is.null(ziformula), censor = censor,
                          weighted = !is.null(wnm))
    fam <- ilm_family(g$family)
    fam_why <- g$why
    message(sprintf("ilm_model(): family = \"%s\", inferred from `%s`: %s.",
                    g$family, ilm_mf_name(formula[[2]]), g$why),
            if (nzchar(g$hint)) paste0(" ", g$hint) else "",
            " Pass `family` to choose another.")
    cl$family <- g$family
  }
  if (isTRUE(fam$ordinal)) {
    ## An ordered factor already carries the ordering. A plain factor is taken
    ## in its level order, which is alphabetical unless someone set it, and
    ## that is worth saying out loud -- "agree", "disagree", "neutral" is a
    ## perfectly ordinary alphabetical ordering and a nonsensical one here.
    if (is.character(yraw)) yraw <- factor(yraw)
    if (is.factor(yraw)) {
      if (!is.ordered(yraw))
        message("ilm_model: the response is an unordered factor, so its ",
                "levels are taken in the order they are stored: ",
                paste(levels(yraw), collapse = " < "),
                ". If that is not the order you mean, make it an ordered ",
                "factor first.")
      ylevels <- levels(yraw); yi <- as.integer(yraw)
    } else {
      yi <- as.integer(round(as.numeric(yraw)))
      ylevels <- as.character(sort(unique(yi)))
      yi <- match(yi, sort(unique(yi)))
    }
    J <- length(ylevels)
    if (J < 3L)
      stop("an ordinal family needs at least 3 ordered categories; with 2 ",
           "there is a single cut and the model is a binomial one -- use ",
           "family = \"binomial\".", call. = FALSE)
  } else if (fam$name == "multinomial") {
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
  contrasts <- ilm_contrasts_arg(contrasts, mf, mt)
  X  <- stats::model.matrix(mt, mf, contrasts.arg = contrasts)
  xlev <- stats::.getXlevels(mt, mf)
  ctr  <- attr(X, "contrasts")
  ## The thresholds ARE the intercepts of a cumulative link model, so the
  ## design must not carry one as well: a column of ones is confounded with
  ## every threshold at once, and the optimiser wanders along that ridge
  ## rather than failing.
  ##
  ## The column is dropped HERE, after model.matrix has coded the factors,
  ## and not by putting - 1 in the formula. Those are not the same thing: a
  ## formula without an intercept makes model.matrix expand the first factor
  ## to ALL its levels instead of contrasting them, which adds a parameter the
  ## data cannot identify. The fit still reaches the same likelihood -- it is
  ## the same model written down twice -- but the Hessian is singular and
  ## every standard error comes back NaN.
  if (isTRUE(fam$ordinal)) {
    ic <- match("(Intercept)", colnames(X), nomatch = 0L)
    if (ic > 0L) {
      X <- X[, -ic, drop = FALSE]
      attr(X, "assign") <- attr(stats::model.matrix(mt, mf,
                                 contrasts.arg = contrasts), "assign")[-ic]
      attr(X, "contrasts") <- ctr
    }
    if (!ncol(X))
      stop("an ordinal model needs at least one predictor: with none, the ",
           "thresholds are the whole model and there is nothing to estimate ",
           "beyond the observed category proportions.", call. = FALSE)
  }
  ## cbind() below drops attributes, so carry the column -> term map by hand.
  ## Smooth null-space columns get NA: they belong to no parametric term.
  asgn <- attr(X, "assign")

  ## ---- smooths: null space -> X, penalised blocks -> basis terms ----------
  re_list <- list(); sm_store <- list()
  for (sp in smsp) {
    sob <- ilm_smooth(sp, mf)
    ## mgcv's own name for the smooth it built: "s(x)", or "s(x):z" for a
    ## numeric `by`. The name as written in the formula, "s(x)" for both,
    ## used to key everything, so s(x) + s(x, by = z) gave the second the
    ## first's name: it overwrote the first's penalised part, the two shared
    ## column names, and predict() on new rows stopped on a column count.
    lab <- ilm_smooth_label(sp)
    if (lab %in% names(sm_store))
      stop("two smooths are both called '", lab, "': the second would ",
           "overwrite the first. Write each smooth once -- they differ only ",
           "in their settings -- or give the second a copy of the ",
           "variable under another name.", call. = FALSE)
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
    ## the grouping side as model.frame() names a column: a lone name bare,
    ## so `site id` is found under "site id"
    gvar <- ilm_mf_name(b[[3]])
    ## EVALUATE the grouping side rather than looking it up by name. A nested
    ## bar's group is an interaction that exists only as an expression --
    ## `country:continent` from `(1 | continent/country)` -- so a lookup finds
    ## nothing and the model refuses a formula lme4 accepts. Evaluating builds
    ## the crossed factor, which is what nesting means.
    g <- if (!is.null(mf[[gvar]])) mf[[gvar]] else
      tryCatch(eval(b[[3]], mf, environment(formula)), error = function(e) NULL)
    if (is.null(g))
      stop("the grouping factor '", gvar, "' could not be built from the ",
           "data. Its variable(s) -- ", paste(all.vars(b[[3]]), collapse = ", "),
           " -- must all be columns of `data`.", call. = FALSE)
    ## an interaction of two factors comes back with every combination as a
    ## level, including the ones that never occur; unused levels would each
    ## claim a random effect that no row informs
    g <- droplevels(as.factor(g))
    ## Built from the expression, not from pasted text: a lone backticked
    ## name deparses bare and no longer parses, and the text round trip also
    ## put the formula in THIS frame, where the caller's variables are not.
    Z <- stats::model.matrix(ilm_one_sided(b[[2]], fenv), mf)
    nm <- gvar; k <- 1L
    while (nm %in% names(re_list)) { k <- k + 1L; nm <- paste0(gvar, ".", k) }
    re_list[[nm]] <- list(group = g, Z = Z, factor = gvar)
  }

  ## ---- fit ----------------------------------------------------------------
  ## joint defaults to TRUE when the model contains smooths.  Measured: holding
  ## the penalised coefficients fixed breaks their correlation with the smooth's
  ## null-space term and inflates interval width nearly threefold (0.358 vs
  ## 0.123), pushing intervals outside [0, 1).  Joint draws are the correct
  ## default wherever a smooth is present.
  if (is.null(joint)) joint <- length(smsp) > 0L
  w <- if (is.null(wnm)) NULL else as.numeric(mf[[wnm]])
  ## The design's weights enter the likelihood the same way frequency
  ## weights do -- the pseudo-likelihood is the weighted one -- and it is
  ## only the VARIANCE that has to know the difference. The rows the
  ## design describes are the rows the model frame kept, so any dropped
  ## by na.action come out of the weights too.
  if (!is.null(design)) {
    om <- attr(mf, "na.action")
    dw <- design$weights
    if (!is.null(om) && length(om)) {
      if (length(dw) != nrow(mf) + length(om))
        stop("the design describes ", length(dw), " rows and the model ",
             "frame kept ", nrow(mf), " of a different total; build the ",
             "design from the same data.", call. = FALSE)
      design <- ilm_design_subset(design, -as.integer(om))
    } else if (length(dw) != nrow(mf))
      stop("the design describes ", length(dw), " rows and the model ",
           "frame has ", nrow(mf), ".", call. = FALSE)
    w <- design$weights
  }
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
  Zzi <- NULL
  if (!is.null(ziformula)) {
    Zzi <- ilm_zi_design(ziformula, mf, NULL)
    attr(Zzi, "formula") <- ziformula
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
                 Zzi = Zzi, zi_type = zi_type,
                  ylevels = ylevels, weights = w, family = fam, verbose = verbose,
                  restarts = restarts, joint = joint, reml = reml,
                  boundary = boundary)

  ## A covariance that ended at its boundary is said HERE, where the fit was
  ## asked for, and not only in summary(): with verbose = FALSE nothing else
  ## would mention it. Named with its remedy and what the remedy costs, as
  ## every diagnostic in this package is. Refits inside the package go
  ## through ilm_fit() and do not repeat it.
  if (identical(boundary, "hold")) {
    grp <- vapply(fit$re, function(e) !identical(e$kind, "basis"), TRUE)
    small <- names(fit$re)[grp][vapply(which(grp), function(k)
      any(sqrt(pmax(diag(as.matrix(fit$Sigma[[k]])), 0)) < 1e-3), TRUE)]
    at <- union(ilm_boundary_at(fit), small)
    if (length(at))
      message("ilm_model(): the random-effect covariance of ",
              paste(sprintf("`%s`", at), collapse = ", "), " sits at the edge ",
              "of its range -- a variance of zero or a correlation of +/-1 -- ",
              "where the data cannot resolve it. The fixed effects and their ",
              "standard errors are still usable; summary() says what else is. ",
              "If the term belongs in the model, boundary = \"avoid\" keeps it ",
              "inside its range with a small penalty: it is then assumed ",
              "nonzero rather than estimated at zero, so do not test whether ",
              "it is; its variance comes out larger, and for a binary or ",
              "categorical outcome the fixed effects a little further from ",
              "zero -- markedly so when a category is rare.")
  }

  ## ---- everything the ecosystem layer reconstructs a reference grid from ---
  fit$design    <- design
  fit$call      <- cl
  ## why the family was chosen, when it was chosen rather than given
  fit$family_inferred <- fam_why
  fit$formula   <- formula
  fit$fixed_formula <- gp$pf
  fit$terms     <- mt
  fit$xlev      <- xlev
  fit$contrasts <- ctr
  fit$model     <- mf
  fit$data_extra <- data_extra
  fit$smooths   <- sm_store
  fit$bars      <- bars
  fit$na.action <- attr(mf, "na.action")
  fit$n_dropped <- n_drop
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

## The rows a fit was made from, as new rows to predict from: the model frame
## with the variables its transformed terms are built from added back, so
## predict() and model.matrix() can rebuild log(x) or a Fourier basis for
## changed values of x or t. A fit made before those were kept gets its frame.
#' @keywords internal
#' @noRd
ilm_data <- function(object) {
  mf <- object$model
  if (is.null(mf)) return(NULL)
  ex <- object$data_extra
  if (is.null(ex) || !ncol(ex)) return(mf)
  for (v in names(ex)) mf[[v]] <- ex[[v]]
  mf
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

#' Expand a contrast shorthand into what model.matrix() wants
#'
#' `model.matrix()` takes a named list of codings, one entry per factor. Naming
#' every factor to say the same thing about all of them is friction for no
#' gain, so a single string means "this coding, for every unordered factor".
#'
#' Ordered factors are left alone. R codes them with orthogonal polynomials by
#' default, which is a statement about their spacing rather than an arbitrary
#' choice, and silently replacing it because someone asked for sum coding on
#' the unordered ones would change what their coefficients mean.
#'
#' @param contrasts `NULL`, a single string, or a named list.
#' @param mf The model frame.
#' @param mt The terms.
#' @return `NULL` or a named list for `contrasts.arg`.
#' @keywords internal
#' @noRd
ilm_contrasts_arg <- function(contrasts, mf, mt) {
  if (is.null(contrasts) || is.list(contrasts)) return(contrasts)
  if (!is.character(contrasts) || length(contrasts) != 1L)
    stop("`contrasts` must be NULL, a single string such as \"sum\", or a ",
         "named list like lm() takes, not ", class(contrasts)[1],
         " of length ", length(contrasts), ".", call. = FALSE)
  nm <- sub("^contr\\.", "", contrasts)
  ok <- c("treatment", "sum", "helmert", "poly", "SAS")
  if (!nm %in% ok)
    stop("`contrasts` should be one of ", paste(dQuote(ok), collapse = ", "),
         " (with or without the \"contr.\" prefix), or a named list. Got ",
         dQuote(contrasts), ".", call. = FALSE)
  vars <- all.vars(stats::delete.response(mt))
  fac <- vars[vapply(vars, function(v)
    !is.null(mf[[v]]) && is.factor(mf[[v]]) && !is.ordered(mf[[v]]), TRUE)]
  if (!length(fac)) return(NULL)
  stats::setNames(as.list(rep(paste0("contr.", nm), length(fac))), fac)
}
