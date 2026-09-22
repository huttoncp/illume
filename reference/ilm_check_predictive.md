# Compare the observed response with data the model would generate

The frequentist counterpart of a posterior predictive check: simulate
datasets from the fitted model and see whether the real response looks
like one of them. It is a blunt instrument, and that is the point – it
catches gross misspecification (the wrong family, unmodelled zero
inflation, a bounded response fitted as gaussian) at a glance, where a
residual plot makes you work for it.

## Usage

``` r
ilm_check_predictive(object, B = 50L, seed = 1L, plot = TRUE)
```

## Arguments

- object:

  A fitted `"ilm_model"` object.

- B:

  Simulated datasets.

- seed:

  Random seed.

- plot:

  Draw the overlay.

## Value

Invisibly, a list with the observed and simulated summaries and a
`status`.

## See also

[`ilm_check_dispersion()`](https://craig-hutton.github.io/illume/reference/ilm_check_dispersion.md),
[`ilm_check_zeros()`](https://craig-hutton.github.io/illume/reference/ilm_check_zeros.md)
for the specific failures this can only hint at.

## Examples

``` r
set.seed(1)
d <- ilm_sim()
f <- ilm_model(visits ~ income + (1 | id), data = d, family = "poisson",
               verbose = FALSE)
#> Warning: NA/NaN function evaluation
#> Warning: NA/NaN function evaluation
#> Warning: NA/NaN function evaluation
#> Warning: NA/NaN function evaluation
#> Warning: NA/NaN function evaluation
#> Warning: NA/NaN function evaluation
#> Warning: NA/NaN function evaluation
#> Warning: NA/NaN function evaluation
#> Warning: NA/NaN function evaluation
ilm_check_predictive(f, B = 30)
```
