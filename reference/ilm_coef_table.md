# Coefficient table with Wald tests

Estimates, standard errors, z statistics and two-sided p-values, in the
layout [`summary()`](https://rdrr.io/r/base/summary.html) prints.

## Usage

``` r
ilm_coef_table(object)
```

## Arguments

- object:

  A fitted `"ilm_model"` object.

## Value

A data frame with columns `Estimate`, `Std. Error`, `z value` and
`Pr(>|z|)`.

## Details

These are Wald tests, which are fast but have a known weakness: for very
strong effects, or when a category is nearly perfectly predicted, the
Wald statistic can *shrink* rather than grow. This is the Hauck-Donner
effect. If a term matters to your conclusions, confirm it with a
likelihood-ratio test via `ilm_anova(test = "LRT")` or, for small
samples,
[`ilm_pb_lrt()`](https://huttoncp.github.io/illume/reference/ilm_pb_lrt.md).

## References

Hauck, W. W., & Donner, A. (1977). Wald's test as applied to hypotheses
in logit analysis. *Journal of the American Statistical Association*,
72(360), 851–853.
