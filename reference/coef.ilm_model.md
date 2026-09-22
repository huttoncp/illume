# Coefficients and their covariance matrix

[`coef()`](https://rdrr.io/r/stats/coef.html) returns the fixed effects
as a named vector and [`vcov()`](https://rdrr.io/r/stats/vcov.html) the
matching covariance matrix. Their names line up exactly, which is what
[`car::Anova()`](https://rdrr.io/pkg/car/man/Anova.html), `emmeans` and
`marginaleffects` rely on.

[`coef()`](https://rdrr.io/r/stats/coef.html) returns the fixed effects
as a named vector and [`vcov()`](https://rdrr.io/r/stats/vcov.html) the
matching covariance matrix. Their names line up exactly, which is what
[`car::Anova()`](https://rdrr.io/pkg/car/man/Anova.html), `emmeans` and
`marginaleffects` rely on.

## Usage

``` r
# S3 method for class 'ilm_model'
coef(object, full = FALSE, ...)

# S3 method for class 'ilm_model'
vcov(object, full = FALSE, ...)
```

## Arguments

- object:

  A fitted `"ilm_model"` object.

- full:

  Logical. Include covariance parameters as well as fixed effects.

- ...:

  Unused.

## Value

A named numeric vector ([`coef()`](https://rdrr.io/r/stats/coef.html))
or a named matrix ([`vcov()`](https://rdrr.io/r/stats/vcov.html)).

A named numeric vector ([`coef()`](https://rdrr.io/r/stats/coef.html))
or a named matrix ([`vcov()`](https://rdrr.io/r/stats/vcov.html)).

## Details

Names take the form `category:predictor`, for example `low:x1`. Because
categories are coded sum-to-zero, `low:x1` is the effect of `x1` on the
`low` category *relative to the average across categories*, not relative
to a baseline category.

`full = TRUE` appends the covariance parameters (entries of the Cholesky
factors) to the fixed effects. This matters for population-averaged
predictions, which depend on the random-effect covariance as well as the
fixed effects: a delta-method calculation given only the fixed effects
would silently report no uncertainty at all from the covariance
parameters.

Names take the form `category:predictor`, for example `low:x1`. Because
categories are coded sum-to-zero, `low:x1` is the effect of `x1` on the
`low` category *relative to the average across categories*, not relative
to a baseline category.

`full = TRUE` appends the covariance parameters (entries of the Cholesky
factors) to the fixed effects. This matters for population-averaged
predictions, which depend on the random-effect covariance as well as the
fixed effects: a delta-method calculation given only the fixed effects
would silently report no uncertainty at all from the covariance
parameters.

## See also

[`ilm_coef_table()`](https://craig-hutton.github.io/illume/reference/ilm_coef_table.md)
for a formatted table with tests.

[`ilm_coef_table()`](https://craig-hutton.github.io/illume/reference/ilm_coef_table.md)
for a formatted table with tests.
