# Distances of the fitted random effects from zero

The model assumes random effects are multivariate normal. Mahalanobis
distance reduces each group's vector of effects to one number measuring
how far it lies from the centre, allowing for their covariance, so the
set can be compared against a chi-square reference.

## Usage

``` r
ilm_re_mahalanobis(object, term = NULL)
```

## Arguments

- object:

  A fitted `"ilm_model"` object.

- term:

  Optional name of the grouping term; defaults to the first.

## Value

A list with the distances, degrees of freedom and term name, or `NULL`
if the model has no grouping terms.

## An important caveat

Fitted random effects are **shrunk** toward zero – that is what makes
them useful – so their spread is smaller than the true random effects.
Comparing them to the theoretical chi-square distribution is therefore
anti-conservative: a misspecified covariance can still look acceptable.
Treat this panel as indicative rather than as a test. The plot produced
by
[`ilm_appraise()`](https://huttoncp.github.io/illume/reference/ilm_appraise.md)
is annotated accordingly.
