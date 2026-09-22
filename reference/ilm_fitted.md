# Fitted category probabilities

Fitted category probabilities

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

A matrix of probabilities with one column per category.
