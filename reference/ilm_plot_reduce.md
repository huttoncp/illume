# Map the observations from a reduction

Each row placed on two of the dimensions – the usual individuals map of
a PCA, MCA or mixed analysis.

## Usage

``` r
ilm_plot_reduce(x, dims = c(1, 2), by = NULL, ...)

ilm_plot_reduce_na(x, dims = c(1, 2), by = NULL, ...)
```

## Arguments

- x:

  An
  [`ilm_reduce()`](https://craig-hutton.github.io/illume/reference/ilm_reduce.md)
  result.

- dims:

  Which two dimensions, by number.

- by:

  Optional vector to colour by, the same length as the data the
  reduction was built from – a cluster assignment, or a column of the
  original data.

- ...:

  Passed to
  [`tinyplot::tinyplot()`](https://grantmcdermott.com/tinyplot/man/tinyplot.html).

## Value

`NULL`, invisibly.

## See also

[`ilm_plot_cluster()`](https://craig-hutton.github.io/illume/reference/ilm_plot_cluster.md),
[`ilm_plot_reduce_scree()`](https://craig-hutton.github.io/illume/reference/ilm_plot_reduce_scree.md).

## Examples

``` r
r <- ilm_reduce(mtcars)
ilm_plot_reduce(r)

ilm_plot_reduce(r, by = factor(mtcars$cyl))
```
