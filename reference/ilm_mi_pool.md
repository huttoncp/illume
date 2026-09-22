# Fit a model across imputations and pool the results

Fits the model to each completed data set and combines them by Rubin's
rules. The point estimate is the average; the variance is the average
within-imputation variance plus the variance of the estimates *between*
imputations, inflated by `1 + 1/m`. That second term is the price of not
having observed the values, and leaving it out is exactly what makes
single imputation overconfident.

## Usage

``` r
ilm_mi_pool(object, formula = NULL, ...)
```

## Arguments

- object:

  An
  [`ilm_impute()`](https://craig-hutton.github.io/illume/reference/ilm_impute.md)
  result, or a list of fitted models.

- formula:

  Model formula, when `object` holds data sets.

- ...:

  Passed to
  [`ilm_model()`](https://craig-hutton.github.io/illume/reference/ilm_model.md).

## Value

An object of class `"ilm_pooled"`: a coefficient table with `estimate`,
`se`, `df`, `lower`, `upper`, `p_value`, plus `fmi` (the fraction of
information lost to missingness) per coefficient.

## Details

Degrees of freedom follow Barnard and Rubin (1999), which corrects the
original formula when the complete-data degrees of freedom are small –
the uncorrected version can return more degrees of freedom than the data
could possibly supply.

## References

Barnard, J. and Rubin, D. B. (1999). Small-sample degrees of freedom
with multiple imputation. Biometrika 86(4).

## See also

[`ilm_impute()`](https://craig-hutton.github.io/illume/reference/ilm_impute.md).

## Examples

``` r
set.seed(1); n <- 200
d <- data.frame(x = rnorm(n), z = rnorm(n))
d$y <- 0.5 * d$x + 0.3 * d$z + rnorm(n)
d$x[sample(n, 40)] <- NA
imp <- ilm_impute(d, m = 5, seed = 1, verbose = FALSE)
ilm_mi_pool(imp, y ~ x + z, family = "gaussian")
#> <ilm_pooled> 5 imputations, combined by Rubin's rules
#> 
#>   term                estimate        se       df     lower     upper     fmi
#>   (Intercept)          -0.0349    0.0765    190.8   -0.1858    0.1160   0.024
#>   x                     0.5689    0.0862    107.9    0.3980    0.7397   0.133
#>   z                     0.2787    0.0784    143.2    0.1237    0.4336   0.088
#> 
#>   fmi is the fraction of information lost to missingness.
#>   Above about 0.5, more imputations are worth having.
```
