# Is the response more variable than the model allows?

The Pearson statistic divided by its degrees of freedom, compared
against the same quantity computed on data simulated from the fitted
model. The reference distribution is simulated rather than assumed, so
the test does not rely on the chi-square approximation that fails for
sparse counts.

## Usage

``` r
ilm_check_dispersion(object, B = 200L, seed = 1L)
```

## Arguments

- object:

  A fitted `"ilm_model"` object.

- B:

  Simulated datasets. A simulated p-value cannot fall below
  `1 / (B + 1)`, so B must exceed 100 for the strongest verdict to be
  reachable at all.

- seed:

  Random seed.

## Value

Invisibly, a list with the observed ratio, the simulated distribution, a
p-value and a `status`.

## See also

[`ilm_check_zeros()`](https://huttoncp.github.io/illume/reference/ilm_check_zeros.md)
when the excess variance is concentrated at zero.

## Examples

``` r
set.seed(1)
d <- ilm_sim()
f <- ilm_model(claims ~ income + (1 | id), data = d, family = "poisson",
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
#> Warning: NA/NaN function evaluation
#> Warning: NA/NaN function evaluation
#> Warning: NA/NaN function evaluation
ilm_check_dispersion(f, B = 200)
#> ilm_check_dispersion: more variable than poisson allows; try family = "nbinom", or check for an omitted predictor
```
