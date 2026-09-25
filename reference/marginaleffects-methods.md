# marginaleffects interface

The four methods `marginaleffects` needs. It computes derivatives
**numerically**: it asks for the parameters, nudges them, sets them
back, and re-predicts.

The four methods `marginaleffects` needs. It computes derivatives
**numerically**: it asks for the parameters, nudges them, sets them
back, and re-predicts.

The four methods `marginaleffects` needs. It computes derivatives
**numerically**: it asks for the parameters, nudges them, sets them
back, and re-predicts.

The four methods `marginaleffects` needs. It computes derivatives
**numerically**: it asks for the parameters, nudges them, sets them
back, and re-predicts.

## Usage

``` r
get_coef.ilm_model(model, ...)

set_coef.ilm_model(model, coefs, ...)

get_vcov.ilm_model(model, vcov = NULL, ...)

get_predict.ilm_model(
  model,
  newdata = NULL,
  type = "response",
  groups = NULL,
  ndraw = 100L,
  marginal = NULL,
  ...
)
```

## Arguments

- model:

  A fitted `"ilm_model"` object.

- ...:

  Unused.

- coefs:

  Numeric replacement parameter vector.

- vcov:

  Passed through by `marginaleffects`.

- newdata:

  Optional data frame.

- type:

  Prediction type.

- groups:

  Which groups the predictions are for: `"population"`, the average over
  the groups, or `"typical"`, a group with every random effect at zero.
  Defaults to the `ilm_model.groups` option, which is `"population"`.
  See
  [`predict.ilm_model()`](https://huttoncp.github.io/illume/reference/predict.ilm_model.md).
  Given to a `marginaleffects` function, `groups` is passed on to this
  method with a warning from `marginaleffects` that it does not know the
  argument; setting the option instead avoids the warning.

- ndraw:

  Integer. Random-effect draws when `groups = "population"`, for a
  multinomial outcome.

- marginal:

  Deprecated. `TRUE` is `groups = "population"`, and `FALSE` is
  `groups = "typical"`. The option `ilm_model.marginal` is deprecated in
  the same way.

## Value

Parameters, a covariance matrix, a modified model, or a long-format data
frame of predictions with `rowid`, `group` and `estimate`.

Parameters, a covariance matrix, a modified model, or a long-format data
frame of predictions with `rowid`, `group` and `estimate`.

Parameters, a covariance matrix, a modified model, or a long-format data
frame of predictions with `rowid`, `group` and `estimate`.

Parameters, a covariance matrix, a modified model, or a long-format data
frame of predictions with `rowid`, `group` and `estimate`.

## Details

`get_coef()` returns fixed effects **and** covariance parameters
together, because a population-averaged prediction depends on both.
Returning only the fixed effects would silently report no uncertainty
from the covariance parameters. Variance parameters therefore appear
among the "coefficients" in printed output, which looks odd but is
correct.

`get_coef()` returns fixed effects **and** covariance parameters
together, because a population-averaged prediction depends on both.
Returning only the fixed effects would silently report no uncertainty
from the covariance parameters. Variance parameters therefore appear
among the "coefficients" in printed output, which looks odd but is
correct.

`get_coef()` returns fixed effects **and** covariance parameters
together, because a population-averaged prediction depends on both.
Returning only the fixed effects would silently report no uncertainty
from the covariance parameters. Variance parameters therefore appear
among the "coefficients" in printed output, which looks odd but is
correct.

`get_coef()` returns fixed effects **and** covariance parameters
together, because a population-averaged prediction depends on both.
Returning only the fixed effects would silently report no uncertainty
from the covariance parameters. Variance parameters therefore appear
among the "coefficients" in printed output, which looks odd but is
correct.

## References

Arel-Bundock, V., Greifer, N., & Heiss, A. (2024). How to interpret
statistical models using marginaleffects for R and Python. *Journal of
Statistical Software*, 111(9), 1–32.

Arel-Bundock, V., Greifer, N., & Heiss, A. (2024). How to interpret
statistical models using marginaleffects for R and Python. *Journal of
Statistical Software*, 111(9), 1–32.

Arel-Bundock, V., Greifer, N., & Heiss, A. (2024). How to interpret
statistical models using marginaleffects for R and Python. *Journal of
Statistical Software*, 111(9), 1–32.

Arel-Bundock, V., Greifer, N., & Heiss, A. (2024). How to interpret
statistical models using marginaleffects for R and Python. *Journal of
Statistical Software*, 111(9), 1–32.
