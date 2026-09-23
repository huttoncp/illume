# Fit a multinomial linear mixed model

The main entry point. Takes an lme4-style formula and a data frame and
fits a mixed-effects model for a nominal categorical outcome with three
or more unordered categories.

Parses the formula, builds the design matrices and calls
[`ilm_fit()`](https://huttoncp.github.io/illume/reference/ilm_fit.md).
Called by `ilm_model()`; documented separately only because
`ilm_model()` dispatches to it.

## Usage

``` r
ilm_model(formula, ...)

ilm_model_formula(
  formula,
  data,
  family = "auto",
  re_struct = NULL,
  ar = NULL,
  weights = NULL,
  contrasts = NULL,
  verbose = TRUE,
  restarts = 3L,
  joint = NULL,
  na.action = stats::na.omit,
  censor = NULL,
  dispformula = NULL,
  rp_df = 3L,
  rp_knots = NULL,
  ziformula = NULL,
  zi_type = c("inflated", "hurdle"),
  design = NULL,
  reml = FALSE
)
```

## Arguments

- formula:

  A formula with random-effect bars and optional smooth terms.

- ...:

  Arguments passed to the formula interface, listed below.

- data:

  A data frame.

- family:

  Response distribution: one of "gaussian", "binomial", "poisson",
  "nbinom", "beta", "multinomial", or one of the ordinal families. See
  [`ilm_family()`](https://huttoncp.github.io/illume/reference/ilm_family.md).
  The default, `"auto"` (or `NULL`), reads the family off the response
  and says which it chose and why: a factor with 2 levels, a logical or
  a 0/1 variable is binomial; an unordered factor with 3 or more levels
  is multinomial and an ordered one ordinal; whole numbers that reach
  down to 0 or 1 are a count, and poisson; whole numbers that never come
  near zero, such as a blood pressure, are read as a measurement, and
  gaussian; values strictly between 0 and 1 are beta; anything else
  numeric is gaussian, as is a response censored at a floor or a
  ceiling. Where the response cannot settle it – a proportion that
  touches 0 or 1 or comes with weights, a numeric variable with two
  values other than 0 and 1, a survival time, a date – the fit stops and
  asks rather than guessing. The choice is a starting point, not a
  verdict: counts are often overdispersed (`"nbinom"`), and ratings on a
  short scale are often better read as ordinal. The chosen family is
  written into `fit$call`, so a refit uses it rather than guessing
  again.

- re_struct:

  Optional named list of category covariance structures, named by
  grouping variable. See
  [`ilm_fit()`](https://huttoncp.github.io/illume/reference/ilm_fit.md).

- ar:

  Optional correlation over time, from
  [`ilm_ar1()`](https://huttoncp.github.io/illume/reference/ilm_ar1.md)
  or
  [`ilm_car1()`](https://huttoncp.github.io/illume/reference/ilm_car1.md).

- weights:

  Optional **frequency** weights: the number of replicate observations
  each row stands for. Evaluated inside `data`. See
  [`ilm_fit()`](https://huttoncp.github.io/illume/reference/ilm_fit.md)
  for when this is valid, and why survey weights are not.

- contrasts:

  How to code factor predictors. `NULL` (the default) uses R's own
  setting, which is treatment coding for unordered factors and
  polynomial for ordered ones – the same as
  [`stats::lm()`](https://rdrr.io/r/stats/lm.html), so a coefficient
  means what it means everywhere else in R. A single string
  (`"treatment"`, `"sum"`, `"helmert"` or `"poly"`, with or without the
  `contr.` prefix) applies that coding to every unordered factor at
  once. A named list sets them one factor at a time, exactly as
  [`stats::lm()`](https://rdrr.io/r/stats/lm.html) takes it.

  Type III tests of a main effect that is also in an interaction need
  orthogonal coding such as
  [`stats::contr.sum()`](https://rdrr.io/r/stats/contrast.html);
  [`ilm_anova()`](https://huttoncp.github.io/illume/reference/ilm_anova.md)
  checks for that and says so rather than reporting a test that is not
  the one it claims. Note that this is a separate matter from the
  sum-to-zero coding a multinomial fit uses across its outcome
  CATEGORIES, which is not a predictor contrast and is not affected by
  this argument.

- verbose:

  Logical. Print checks while fitting.

- restarts:

  Integer. Optimiser restarts.

- joint:

  Logical or `NULL`. Compute the joint precision over fixed and random
  parameters. `NULL` (the default) switches it on when the model
  contains smooths, which is when
  [`predict.ilm_model()`](https://huttoncp.github.io/illume/reference/predict.ilm_model.md)
  needs it.

- na.action:

  How to handle missing values; default
  [`stats::na.omit()`](https://rdrr.io/r/stats/na.fail.html).

- censor:

  Optional censoring specification from
  [`ilm_censor()`](https://huttoncp.github.io/illume/reference/ilm_censor.md),
  for a response with a floor, a ceiling or a detection limit.

- dispformula:

  Optional one-sided formula for the dispersion, modelling its
  logarithm: `~ group` for a separate spread per level, `~ x` for one
  that changes with a covariate, `~ mu` for a power of the fitted mean.
  `mu` is a reserved name. This is the remedy for what
  [`ilm_check_variance()`](https://huttoncp.github.io/illume/reference/ilm_check_variance.md)
  diagnoses.

- rp_df:

  Degrees of freedom for a flexible parametric baseline, used by
  `family = "rp"`, `"rp_odds"` and `"rp_normal"`. `1` is a straight line
  in log time, and so the corresponding parametric model; `3` is the
  usual default and allows two interior knots.

- rp_knots:

  Knot positions on the log-time scale, given directly in place of
  `rp_df`.

- ziformula:

  Optional one-sided formula for the zero part of a count model, on the
  logit scale: `~ 1` for a constant excess-zero probability, `~ x` for
  one that depends on a predictor. This is the remedy for what
  [`ilm_check_zeros()`](https://huttoncp.github.io/illume/reference/ilm_check_zeros.md)
  diagnoses. Fixed effects only – a random effect in the zero part is
  not supported, and a `|` here is an error rather than something
  quietly dropped.

- zi_type:

  `"inflated"` or `"hurdle"`, and they are different models.
  `"inflated"` is a mixture: some rows are structural zeros and the rest
  come from a count that may itself be zero, so a zero in the data could
  have come from either. `"hurdle"` is two processes: whether the
  response clears zero, and how far past it goes, the latter fitted to a
  count that cannot be zero. Choose by what the zeros mean, not by fit –
  a structural zero is a unit that was never at risk. See
  [`ilm_zi_coef()`](https://huttoncp.github.io/illume/reference/ilm_zi_coef.md).

- design:

  An
  [`ilm_design()`](https://huttoncp.github.io/illume/reference/ilm_design.md)
  describing how a complex sample was drawn. Supplying one fits with the
  sampling weights and attaches the design, so
  [`ilm_svy_coef()`](https://huttoncp.github.io/illume/reference/ilm_svy_coef.md)
  can report a variance that reflects the clustering and stratification.
  Do not also pass `weights`: a sampling weight and a replicate count
  are different things and the design already carries one.

- reml:

  Logical. Estimate the variance components by RESTRICTED maximum
  likelihood instead of maximum likelihood. Defaults to `FALSE`, and the
  reason is the order you work in: maximum likelihood is what lets you
  compare fixed-effect structures, because a restricted likelihood
  belongs to contrasts orthogonal to the design matrix and changing that
  matrix changes which data it is the likelihood of. Settle the fixed
  effects under the default, then refit with `reml = TRUE` for the
  estimates you report – maximum likelihood biases the variance
  components downward, and with few clusters that carries through to
  standard errors and to the degrees of freedom from
  [`ilm_denom_df()`](https://huttoncp.github.io/illume/reference/ilm_denom_df.md).
  Once set, any likelihood-ratio test refuses rather than quietly
  comparing things that are not comparable, and so does
  [`ilm_robust()`](https://huttoncp.github.io/illume/reference/ilm_robust.md),
  whose sandwich needs per-observation scores that a restricted
  likelihood does not have. Available for every family, but it delivers
  different amounts depending on the family – see the section below.
  [`ilm_dag_model()`](https://huttoncp.github.io/illume/reference/ilm_dag_model.md)
  defaults to `TRUE`, because there the graph fixed the adjustment set
  before any data were seen.

## Value

An object of class `"ilm_model"`. Beyond the elements listed in
[`ilm_fit()`](https://huttoncp.github.io/illume/reference/ilm_fit.md), a
formula fit also stores `call`, `terms`, `xlev`, `contrasts`, the model
frame and the smooth objects – everything needed to rebuild a reference
grid for
[`predict.ilm_model()`](https://huttoncp.github.io/illume/reference/predict.ilm_model.md)
and for `emmeans` or `marginaleffects`.

An object of class `"ilm_model"`.

## Details

Random-effect terms use lme4 syntax and smooths use mgcv syntax, so a
model can mix them freely:


      y ~ 1                                      intercept only
      y ~ x1 + grp                               ordinary regression
      y ~ x1 + grp + (1 | subj)                  random intercept
      y ~ x1 + time + (1 + time | subj)          random intercept and slope
      y ~ x1 + s(xs, k = 10) + (1 | subj)        penalised smooth of xs
      y ~ x1 + t2(lon, lat) + (1 | site)         2-D surface (spatial)

Two rules for smooth terms. Use `t2()` rather than `te()` for tensor
products, because `te()` cannot be converted to the mixed-model form
this package relies on. And write smooths **unqualified** – `s(x)`, not
`mgcv::s(x)` – because mgcv identifies them by name, so a namespaced
call would be mistaken for an ordinary predictor. This matches
[`mgcv::gam()`](https://rdrr.io/pkg/mgcv/man/gam.html).

The outcome may be a factor or a character vector. Its levels set the
category labels used throughout the output.

Random-effect bars are extracted with `lme4::findbars()` and smooths
with
[`mgcv::interpret.gam()`](https://rdrr.io/pkg/mgcv/man/interpret.gam.html),
then reparameterised by `ilm_smooth()` so that the unpenalised part of
each smooth joins the fixed effects and the penalised part becomes a
random term.

A smooth is centred on the observed data, which changes what the
INTERCEPT means: it becomes the mean response at the sample average of
the smooth, not at a fixed point. That average moves from sample to
sample, and the standard error on the intercept does not carry that
movement, so an interval for it is narrower than its sampling spread.
Measured on 2000 replicates of a 600-row design, the intercept covered a
fixed population value 88% of the time and the sample-specific value
94.3%; the slopes were unaffected, which is what the centring is for.
Read the intercept of a model with a smooth as a property of the sample,
and take conclusions from the slopes and from
[`ilm_plot_model()`](https://huttoncp.github.io/illume/reference/ilm_plot_model.md)'s
effect curves.

## Simple models get exact inference

A gaussian model with no random or smooth terms is an ordinary linear
model. In that case there is nothing to integrate out, so illume reports
**exact** t tests on `n - p` degrees of freedom and F tests in
[`ilm_anova()`](https://huttoncp.github.io/illume/reference/ilm_anova.md),
rather than the large-sample normal and chi-square approximations it
must use when random effects are present. Coefficients, standard errors,
the residual standard deviation and the F tests then agree with
[`stats::lm()`](https://rdrr.io/r/stats/lm.html) and
[`car::Anova()`](https://rdrr.io/pkg/car/man/Anova.html) to numerical
precision.

This applies only where an exact reference genuinely exists. A Poisson
or binomial model without random effects is still a GLM, with no exact
small-sample analogue, so it keeps z and chi-square.

## Reading the coefficients

Categories are coded **sum-to-zero**, so a coefficient is that
category's deviation from the average across categories, *not* a
contrast against a baseline.
[`summary()`](https://rdrr.io/r/base/summary.html) prints a reminder,
because this is easy to misread if you are used to
[`nnet::multinom()`](https://rdrr.io/pkg/nnet/man/multinom.html).

## Always read the checks

This model class fails quietly: a fit can return sensible-looking
coefficients while its covariance matrix is unusable, making the
standard errors meaningless.
[`summary()`](https://rdrr.io/r/base/summary.html) prints the check
verdicts for that reason, and `fit$checks` holds the full table with a
reason and a suggested remedy for anything that is not `"OK"`.

## See also

[`summary.ilm_model()`](https://huttoncp.github.io/illume/reference/summary.ilm_model.md),
[`ilm_anova()`](https://huttoncp.github.io/illume/reference/ilm_anova.md),
[`predict.ilm_model()`](https://huttoncp.github.io/illume/reference/predict.ilm_model.md),
[`ilm_appraise()`](https://huttoncp.github.io/illume/reference/ilm_appraise.md),
[`ilm_pb_lrt()`](https://huttoncp.github.io/illume/reference/ilm_pb_lrt.md).

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(1)
n <- 600
dd <- data.frame(
  subj = factor(sample(40, n, TRUE)),
  x1   = rnorm(n),
  grp  = factor(sample(c("a", "b", "c"), n, TRUE))
)
dd$y <- factor(sample(c("low", "mid", "high"), n, TRUE),
               levels = c("low", "mid", "high"))

fit <- ilm_model(y ~ x1 + grp + (1 | subj), data = dd)
summary(fit)
ilm_anova(fit, type = 3)
head(predict(fit))
} # }
```
