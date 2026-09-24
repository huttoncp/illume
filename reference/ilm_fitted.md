# Fitted values

The fitted mean for each observation: the probability of each category
for a multinomial or ordinal outcome, the fitted value otherwise.

## Usage

``` r
ilm_fitted(object, conditional = TRUE)
```

## Arguments

- object:

  A fitted `"ilm_model"` object.

- conditional:

  Logical. `TRUE` evaluates the random effects at their fitted values,
  giving in-sample fitted probabilities. `FALSE` sets grouping and AR
  terms to zero. Smooth terms are included either way, because they are
  mean structure rather than a population to average over.

## Value

A matrix with one row per observation: one column per category for a
multinomial or ordinal outcome, one column otherwise.
