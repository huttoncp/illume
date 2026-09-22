# The sample size a target power implies

Interpolates the simulated curve, and reports the range the Monte Carlo
error allows rather than a single number the simulation cannot support.

## Usage

``` r
ilm_power_n(object, target = 0.8)
```

## Arguments

- object:

  An
  [`ilm_power()`](https://craig-hutton.github.io/illume/reference/ilm_power.md)
  result.

- target:

  Power to reach.

## Value

A data frame with one row per effect size: the interpolated `n`, and the
`n_lower`/`n_upper` implied by the Monte Carlo interval.

## See also

[`ilm_power()`](https://craig-hutton.github.io/illume/reference/ilm_power.md).
