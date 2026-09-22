# Map the clusters from a profile

Convenience wrapper for
[`ilm_plot_cluster()`](https://huttoncp.github.io/illume/reference/ilm_plot_cluster.md)
when what you have is a whole
[`ilm_profile()`](https://huttoncp.github.io/illume/reference/ilm_profile.md)
rather than the clustering on its own.

## Usage

``` r
ilm_plot_profile(x, dims = c(1, 2), ...)

ilm_plot_profile_na(x, dims = c(1, 2), ...)
```

## Arguments

- x:

  An
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

[`ilm_profile()`](https://huttoncp.github.io/illume/reference/ilm_profile.md).

## Examples

``` r
ilm_plot_profile(ilm_profile(mtcars, k_max = 5, B = 25, seed = 1))
```
