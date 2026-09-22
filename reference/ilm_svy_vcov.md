# Linearization variance for a fit from a complex sample

The design-based (Taylor linearization) covariance of the fixed effects:
score contributions summed within each primary sampling unit, centred
within each stratum, and scaled by that stratum's number of clusters.

## Usage

``` r
ilm_svy_vcov(object, design = NULL, lonely = c("adjust", "certainty", "fail"))
```

## Arguments

- object:

  A fitted
  [`ilm_model()`](https://craig-hutton.github.io/illume/reference/ilm_model.md)
  with no random effects.

- design:

  An
  [`ilm_design()`](https://craig-hutton.github.io/illume/reference/ilm_design.md),
  or `NULL` to use the one the fit carries.

- lonely:

  What to do with a stratum holding a single PSU, which supplies no
  variance of its own. `"adjust"` centres it at the grand mean of all
  PSU totals instead of its own, which is conservative; `"certainty"`
  treats it as contributing nothing, which assumes it was selected with
  certainty; `"fail"` refuses.

## Value

A covariance matrix, with the design degrees of freedom attached.

## Why the model-based variance is not an option here

A likelihood weighted by sampling weights believes it saw `sum(w)`
observations rather than `n`, so its standard errors are too small by at
least `sqrt(n / sum(w))` – and by more once the sample is clustered,
since the design effect is on top of that. On the stratified two-stage
example in
[`ilm_svy_coef()`](https://craig-hutton.github.io/illume/reference/ilm_svy_coef.md)
the weights alone account for a factor of 5.9 and the observed ratios
are 7.3 and 16.6.

The point estimates are unaffected: a weighted likelihood is
design-consistent for the population parameter. That is what makes the
mistake easy to miss – everything looks right except the uncertainty,
and the uncertainty looks better than right.

## References

Binder, D. A. (1983). On the variances of asymptotically normal
estimators from complex surveys. *International Statistical Review* 51,
279-292.

## See also

[`ilm_svy_coef()`](https://craig-hutton.github.io/illume/reference/ilm_svy_coef.md),
[`ilm_design()`](https://craig-hutton.github.io/illume/reference/ilm_design.md).
