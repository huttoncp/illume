# Coefficient table with Wald tests

Estimates, standard errors, test statistics and two-sided p-values, in
the layout [`summary()`](https://rdrr.io/r/base/summary.html) prints.

## Usage

``` r
ilm_coef_table(object, df = "auto")
```

## Arguments

- object:

  A fitted `"ilm_model"` object.

- df:

  The degrees of freedom for a gaussian mixed model: `"auto"`
  (Satterthwaite), `"satterthwaite"`, `"kenward-roger"`, `"asymptotic"`
  (the normal), or a number. Ignored for other fits.

## Value

A data frame with columns `Estimate`, `Std. Error`, and either `z value`
and `Pr(>|z|)` or, with finite degrees of freedom, `t value` and
`Pr(>|t|)` – and for a gaussian mixed model a `df` column before them.
The attribute `"df_method"` names the reference.

## Details

**The reference distribution.** With nothing integrated out – a linear
model – the test is exactly t on the residual degrees of freedom. For a
gaussian MIXED model the default is t on Satterthwaite's degrees of
freedom, one per coefficient, as `lmerTest` reports: with few groups the
normal reference is too liberal. `df = "kenward-roger"` gives Kenward
and Roger's df and their adjusted standard errors instead (REML fits); a
number fixes the df; `"asymptotic"` gives the normal, as before. Where a
coefficient's Satterthwaite df cannot be formed it is tested against the
normal and the table says which. Every other family is tested against
the normal: neither approximation is derived for it, and
[`ilm_pb_lrt()`](https://huttoncp.github.io/illume/reference/ilm_pb_lrt.md)
simulates the reference instead.

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

Satterthwaite, F. E. (1946). An approximate distribution of estimates of
variance components. *Biometrics Bulletin*, 2(6), 110–114.
