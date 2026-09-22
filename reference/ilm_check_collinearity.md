# Variance inflation among fixed effects

[`ilm_frame_issues()`](https://craig-hutton.github.io/illume/reference/ilm_frame_issues.md)
catches pairs of columns that are nearly identical, but the case that
actually breaks a model is a predictor collinear with a *combination* of
the others, which no pairwise correlation reveals. That is what variance
inflation measures.

## Usage

``` r
ilm_check_collinearity(object, warn = sqrt(5), fail = sqrt(10))
```

## Arguments

- object:

  A fitted `"ilm_model"` object.

- warn, fail:

  Thresholds on the standard-error inflation factor. The defaults
  correspond to the conventional VIF cut-offs of 5 and 10.

## Value

A data frame with one row per fixed-effect term: `gvif`, `df`,
`se_inflation`, and a `status` of `"OK"`, `"WARN"` or `"FAIL"`.

## Details

For a term with more than one degree of freedom – a factor, a spline
basis – the plain VIF is not comparable across terms, so the generalised
form (GVIF) is reported alongside `gvif^(1/(2*df))`, which is on the
scale of a standard error inflation and can be compared directly.

## References

Fox, J. and Monette, G. (1992). Generalized collinearity diagnostics.
Journal of the American Statistical Association, 87(417), 178-183.

## See also

[`ilm_frame_issues()`](https://craig-hutton.github.io/illume/reference/ilm_frame_issues.md)
for collinearity in the raw data.

## Examples

``` r
set.seed(1)
d <- ilm_sim()
f <- ilm_model(score ~ income + grp + (1 | id), data = d,
               family = "gaussian", verbose = FALSE)
ilm_check_collinearity(f)
#>     term df  gvif se_inflation status
#> 1 income  1 1.004        1.002     OK
#> 2    grp  3 1.004        1.001     OK
```
