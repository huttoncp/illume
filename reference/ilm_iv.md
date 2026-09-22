# Instrumental-variables regression by two-stage least squares

For a linear model in which one or more regressors are correlated with
the error – a treatment people chose rather than were assigned, a price
set in response to demand, an exposure measured with error. An
instrument moves the offending regressor without acting on the outcome
any other way, and the part of the regressor it explains is used in its
place.

## Usage

``` r
ilm_iv(
  formula,
  data,
  weights = NULL,
  cluster = NULL,
  robust = FALSE,
  na.action = stats::na.omit
)
```

## Arguments

- formula:

  `y ~ x + w | z + w`; see above.

- data:

  A data frame.

- weights:

  Optional frequency weights.

- cluster:

  Optional one-sided formula or vector for cluster-robust standard
  errors; see
  [`ilm_vcov_cluster()`](https://craig-hutton.github.io/illume/reference/ilm_vcov_cluster.md).

- robust:

  Use heteroskedasticity-robust standard errors. Implied by `cluster`.

- na.action:

  How to handle missing values.

## Value

An object of class `"ilm_iv"`: the coefficients, their standard errors
and tests, the diagnostics above, and the pieces they were built from.

## Writing the formula

`y ~ x + w | z + w`. To the left of the bar is the model you want. To
the right is the **full** list of exogenous variables: the instruments,
plus every covariate that appears on both sides. Leaving `w` out of the
right side does not mean "w is fine", it means "w is endogenous too and
has no instrument", which is a different and usually unintended model.
That case is an error here rather than a silent one.

## What it reports without being asked

- **First-stage F**, one per endogenous regressor, on the excluded
  instruments. Above 10 is the familiar rule of thumb; above about 104
  is what a conventional 5% t-test actually needs, after Lee et al.
  (2022). Both are shown, because the gap between them is the point.

- **Durbin-Wu-Hausman**, testing whether the regressor was endogenous at
  all. If it was not, least squares is unbiased and far more precise,
  and instrumenting has thrown that away for nothing.

- **Sargan's J**, when there are more instruments than endogenous
  regressors, testing whether they agree with each other. It cannot test
  whether they are all invalid together, which is the assumption that
  matters and the one no test reaches.

When the instrument is weak the remedy is
[`ilm_iv_ar()`](https://craig-hutton.github.io/illume/reference/ilm_iv_ar.md),
which builds a confidence set that stays valid however weak it is.

## References

Lee, D. S., McCrary, J., Moreira, M. J. and Porter, J. (2022). Valid
t-ratio inference for IV. *American Economic Review* 112, 3260-3290.

## See also

[`ilm_iv_ar()`](https://craig-hutton.github.io/illume/reference/ilm_iv_ar.md)
for weak-instrument-robust intervals,
[`ilm_dag_model()`](https://craig-hutton.github.io/illume/reference/ilm_dag_model.md)
for choosing an adjustment set when no instrument is available.

## Examples

``` r
set.seed(1); n <- 500
z <- rnorm(n); u <- rnorm(n)
x <- 0.8 * z + u + rnorm(n)          # x is correlated with the error u
y <- 1 + 0.5 * x + u + rnorm(n)
d <- data.frame(y = y, x = x, z = z)
coef(lm(y ~ x, data = d))            # biased upwards by the shared u
#> (Intercept)           x 
#>   0.9378286   0.9317703 
ilm_iv(y ~ x | z, data = d)
#> Instrumental-variables regression (2SLS)
#>   endogenous:  x
#>   instruments: z
#>   500 observations; classical standard errors
#> 
#>             Estimate Std. Error t value  Pr(>|t|)    
#> (Intercept) 0.925233   0.067934 13.6197 < 2.2e-16 ***
#> x           0.522985   0.092160  5.6748  2.36e-08 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Diagnostics
#>   first stage, x          F =   137.56 on 1 and 498 df (p = 3.21e-28)
#>     rule of thumb F > 10: met; F > 104.7 for a valid 5% t-test: met
#>   Durbin-Wu-Hausman    F =    32.67 on 1 and 497 df (p = 1.89e-08)
```
