# A first-order autoregressive structure over evenly spaced time

Observations of the same unit at adjacent time steps share a correlated
latent value. Use this when the time index moves in equal steps: study
visits, months, waves. For irregular times use
[`ilm_car1()`](https://craig-hutton.github.io/illume/reference/ilm_car1.md).

## Usage

``` r
ilm_ar1(time, group, verbose = TRUE)
```

## Arguments

- time:

  Time index, one value per observation. Whole numbers, or values on a
  common step such as `c(0, 5, 10, 15)`.

- group:

  Unit identifier, one value per observation.

- verbose:

  Warn when the latent budget is thin.

## Value

A `"ilm_ar1"` specification, to pass as `ilm_model(ar = )`.

## Details

The latent grid spans every time step between the earliest and the
latest observation, so a unit that misses a wave is handled correctly –
the latent for that step simply has no observation attached, and the
chain continues across it. What is *not* allowed is times that do not
sit on a common grid at all, because then "one step" has no meaning.

## See also

[`ilm_car1()`](https://craig-hutton.github.io/illume/reference/ilm_car1.md)
for arbitrary gaps,
[`ilm_check_ar()`](https://craig-hutton.github.io/illume/reference/ilm_check_ar.md)
to test whether the structure is needed.

## Examples

``` r
d <- ilm_sim(n_id = 20, n_period = 8)
a <- ilm_ar1(as.integer(factor(d$date)), d$id)
#> Warning: AR(1) puts 160 latent values under 160 observations (1 per latent). Below about 1.5 the Laplace approximation frequently fails to converge, and a latent value seen once carries no information the residual does not. Round `time` to a coarser grid so observations share a latent.
a
#> AR(1), evenly spaced time 
#>   20 groups, 160 latent values, 160 observations (1.00 per latent)
#>   step: 1, grid of 8 time slots
#>   the latent budget is thin; see ?ilm_car1
```
