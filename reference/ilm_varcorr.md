# Variance components of a fitted model

Each random term's covariance, a correlation over time's parameters and
the family's dispersion, all on their natural scale: what
[`lme4::VarCorr()`](https://rdrr.io/pkg/nlme/man/VarCorr.html) reports,
for every structure illume fits.

## Usage

``` r
ilm_varcorr(object)

# S3 method for class 'ilm_model'
VarCorr(x, sigma = 1, ...)

# S3 method for class 'ilm_VarCorr'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)
```

## Arguments

- object, x:

  A fitted `"ilm_model"`.

- sigma:

  Unused; for the generic.

- ...:

  Unused.

- row.names, optional:

  Unused; for the generic.

## Value

An object of class `"ilm_VarCorr"`: a list with `re` (one covariance
matrix per random term, with `"stddev"` and `"correlation"` attributes
as `lme4` gives them), `latent` (the correlation over time, or `NULL`),
`dispersion` (`value`, `meaning`, `modelled`, or `NULL`) and `family`.
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) gives
`lme4`'s columns: `grp`, `var1`, `var2`, `vcov`, `sdcor`.

## Details

**Random terms.** For one linear predictor, a term's covariance is that
of its coefficients: an intercept's variance, or an intercept and
slopes' variances and covariances. For a multinomial outcome it is the
covariance across categories, with the within-group structure beside it
as the `"Sigma_d"` attribute, whose `[1, 1]` the fit holds at one
because the category covariance carries the scale.

**Over time.** An AR(1) term's `Sigma` is its stationary covariance and
`rho` the correlation one grid step apart; a CAR(1) term's `rho` is per
unit of time, and `range` its reciprocal decay; a random walk has
`var_per_time`, the variance of a step one time unit long. The words in
`meaning` are those of
[`ilm_cells()`](https://huttoncp.github.io/illume/reference/ilm_cells.md).

**Dispersion.** On the natural scale and named for what it is. For a
gaussian model it is `sigma`, the residual standard DEVIATION, not the
variance. With a dispersion model it varies by row, and this is its
value at the median row.

These are the variance parameters transformed exactly as a draw of them
would be – through the one transform, `ilm_rebuild()` – so that anything
summing over draws and this function cannot disagree at the estimate.

With `nlme` or `lme4` attached, `VarCorr(fit)` gives the same.

## See also

[`ilm_ranef()`](https://huttoncp.github.io/illume/reference/ilm_ranef.md),
[`ilm_cells()`](https://huttoncp.github.io/illume/reference/ilm_cells.md).

## Examples

``` r
set.seed(1)
d <- data.frame(id = factor(rep(1:12, each = 6)), t = rep(0:5, 12))
d$y <- 1 + rnorm(12)[d$id] + rnorm(12, 0, 0.3)[d$id] * d$t + rnorm(72)
fit <- ilm_model(y ~ t + (1 + t | id), data = d, family = "gaussian",
                 verbose = FALSE)
ilm_varcorr(fit)
#> Variance components (natural scale)
#>  Groups   Name        Variance Std.Dev.
#>  id       (Intercept) 0.3968   0.63    
#>           t           0.09065  0.3011  
#>  Residual             0.7757   0.8807  
#> 
#> Correlations
#>  Groups Between           Corr 
#>  id     (Intercept) and t 0.089
#> 
#> Dispersion: sigma = 0.8807 -- sigma, the residual standard deviation (not the variance)  
as.data.frame(ilm_varcorr(fit))
#>        grp        var1 var2       vcov      sdcor
#> 1       id (Intercept) <NA> 0.39684120 0.62995333
#> 2       id           t <NA> 0.09065218 0.30108500
#> 3       id (Intercept)    t 0.01693502 0.08928698
#> 4 Residual        <NA> <NA> 0.77571995 0.88074965
```
