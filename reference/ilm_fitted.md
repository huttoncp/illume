# Fitted values

The fitted mean for each observation: the probability of each category
for a multinomial outcome, the fitted value otherwise. For an ordinal
outcome it is the latent linear predictor, which the thresholds cut into
categories;
[`predict.ilm_model()`](https://huttoncp.github.io/illume/reference/predict.ilm_model.md)
gives each category's probability.

## Usage

``` r
ilm_fitted(object, groups = c("fitted", "typical"), conditional = NULL)
```

## Arguments

- object:

  A fitted `"ilm_model"` object.

- groups:

  `"fitted"` (the default) uses each group's own estimated random
  effects, and each AR, CAR or random-walk cell's, giving the in-sample
  fitted values. `"typical"` sets them all to zero, for a group exactly
  at the average. Smooth terms are included either way, because they are
  mean structure rather than a population of groups. The average over
  the groups is
  [`predict.ilm_model()`](https://huttoncp.github.io/illume/reference/predict.ilm_model.md)'s
  `groups = "population"`.

- conditional:

  Deprecated. `TRUE` is `groups = "fitted"`, and `FALSE` is
  `groups = "typical"`.

## Value

A matrix with one row per observation: one column per category for a
multinomial outcome, one column otherwise.
