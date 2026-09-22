# The gap statistic across every k considered

With the chosen `k` marked. Only available when `k` was searched for
rather than given.

## Usage

``` r
ilm_plot_cluster_gap(x, ...)

ilm_plot_cluster_gap_na(x, ...)
```

## Arguments

- x:

  An
  [`ilm_cluster()`](https://craig-hutton.github.io/illume/reference/ilm_cluster.md)
  or
  [`ilm_profile()`](https://craig-hutton.github.io/illume/reference/ilm_profile.md)
  result.

- ...:

  Passed to
  [`tinyplot::tinyplot()`](https://grantmcdermott.com/tinyplot/man/tinyplot.html).

## Value

`NULL`, invisibly.

## See also

[`ilm_cluster()`](https://craig-hutton.github.io/illume/reference/ilm_cluster.md).

## Examples

``` r
ilm_plot_cluster_gap(ilm_cluster(ilm_reduce(mtcars), k_max = 5, B = 25,
                                 seed = 1))
```
