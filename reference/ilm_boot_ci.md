# Bootstrap confidence interval for a statistic

Bootstrap confidence interval for a statistic

## Usage

``` r
ilm_boot_ci(
  data,
  y = NULL,
  by = NULL,
  stat = "mean",
  R = 2000L,
  conf = 0.95,
  ci_type = "percentile",
  seed = NULL,
  progress = NULL
)
```

## Arguments

- data:

  A data frame, or a numeric vector when `y` is `NULL`.

- y:

  Name of the numeric column to summarise.

- by:

  Optional grouping columns.

- stat:

  `"mean"`, `"median"`, `"sd"`, `"var"`, or a function taking a numeric
  vector.

- R:

  Bootstrap replicates.

- conf:

  Confidence level.

- ci_type:

  `"percentile"`, `"bca"`, `"normal"` or `"basic"`. On skewed data
  percentile and BCa hold their nominal coverage better than the other
  two.

- seed:

  Random seed.

- progress:

  Show a progress bar. Defaults to
  [`interactive()`](https://rdrr.io/r/base/interactive.html), so a bar
  appears when someone is watching and nothing is written in a script or
  a knitted document. See
  [ilm_progress_arg](https://huttoncp.github.io/illume/reference/ilm_progress_arg.md).

## Value

A one-row data frame per group, with `observed`, `lower`, `upper` and
the settings used.

## References

Efron, B. and Tibshirani, R. J. (1993). An Introduction to the
Bootstrap. Chapman and Hall.

## See also

[`ilm_boot_diff()`](https://huttoncp.github.io/illume/reference/ilm_boot_diff.md)
for a difference between two groups.

## Examples

``` r
d <- ilm_sim()
ilm_boot_ci(d, "score", R = 200, seed = 1)
#>   stat observed    lower    upper conf   R    ci_type   n
#> 1 mean 49.67677 49.26271 50.05195 0.95 200 percentile 900
ilm_boot_ci(d, "score", by = "grp", R = 200, seed = 1)
#>     grp stat observed    lower    upper conf   R    ci_type   n
#> 1 alpha mean 50.01790 49.46450 50.55394 0.95 200 percentile 404
#> 2  beta mean 49.24146 48.56888 49.91962 0.95 200 percentile 329
#> 3 gamma mean 49.64305 48.79232 50.68667 0.95 200 percentile 164
#> 4 delta mean 53.32000 50.63000 55.67000 0.95 200 percentile   3
```
