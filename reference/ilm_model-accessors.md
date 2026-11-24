# Standard accessors for a fitted model

Supply the pieces that other packages expect from a model object, so
that [`car::Anova()`](https://rdrr.io/pkg/car/man/Anova.html), `emmeans`
and `marginaleffects` can work with a fit.

Supply the pieces that other packages expect from a model object, so
that [`car::Anova()`](https://rdrr.io/pkg/car/man/Anova.html), `emmeans`
and `marginaleffects` can work with a fit.

Supply the pieces that other packages expect from a model object, so
that [`car::Anova()`](https://rdrr.io/pkg/car/man/Anova.html), `emmeans`
and `marginaleffects` can work with a fit.

Supply the pieces that other packages expect from a model object, so
that [`car::Anova()`](https://rdrr.io/pkg/car/man/Anova.html), `emmeans`
and `marginaleffects` can work with a fit.

Supply the pieces that other packages expect from a model object, so
that [`car::Anova()`](https://rdrr.io/pkg/car/man/Anova.html), `emmeans`
and `marginaleffects` can work with a fit.

## Usage

``` r
# S3 method for class 'ilm_model'
terms(x, ...)

# S3 method for class 'ilm_model'
formula(x, ...)

# S3 method for class 'ilm_model'
model.frame(formula, ...)

# S3 method for class 'ilm_model'
model.matrix(object, ...)

# S3 method for class 'ilm_model'
getCall(x, ...)
```

## Arguments

- x, object, formula:

  A fitted `"ilm_model"` object.

- ...:

  Unused.

## Value

The model terms, formula, model frame, model matrix or original call.

The model terms, formula, model frame, model matrix or original call.

The model terms, formula, model frame, model matrix or original call.

The model terms, formula, model frame, model matrix or original call.

The model terms, formula, model frame, model matrix or original call.

## Details

**The fixed-effects formula.** `formula(fit)` is the formula as it was
written, with its random-effect terms and its smooths. `terms(fit)`
holds the fixed-effects part of it, so `formula(terms(fit))` is the
fixed-effects formula: the random-effect terms and the smooths are left
out, and the ordinary terms, any
[`offset()`](https://rdrr.io/r/stats/offset.html) and a removed
intercept are kept, as in `y ~ g + poly(t, 2) + offset(log(e))` or
`y ~ g + t - 1`. A dispersion formula is a separate formula and is not
part of it. A smooth's unpenalised columns are not terms of it: they are
among the coefficients, as `s(x).f1`.
