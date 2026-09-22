# Map the clusters

The observations on two dimensions, coloured by which cluster they
landed in.

## Usage

``` r
ilm_plot_cluster(x, dims = c(1, 2), ...)

ilm_plot_cluster_na(x, dims = c(1, 2), ...)
```

## Arguments

- x:

  An
  [`ilm_cluster()`](https://huttoncp.github.io/illume/reference/ilm_cluster.md)
  or
  [`ilm_profile()`](https://huttoncp.github.io/illume/reference/ilm_profile.md)
  result.

- dims:

  Which two coordinates, by number.

- ...:

  Passed to
  [`tinyplot::tinyplot()`](https://grantmcdermott.com/tinyplot/man/tinyplot.html).

## Value

`NULL`, invisibly.

## See also

[`ilm_cluster()`](https://huttoncp.github.io/illume/reference/ilm_cluster.md),
[`ilm_plot_cluster_gap()`](https://huttoncp.github.io/illume/reference/ilm_plot_cluster_gap.md).

## Examples

``` r
ilm_plot_cluster(ilm_cluster(ilm_reduce(mtcars), k_max = 5, B = 25, seed = 1))
```
