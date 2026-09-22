# Scree plot for a reduction

How much variance each dimension accounts for, which is what decides how
many are worth keeping.

## Usage

``` r
ilm_plot_reduce_scree(x, ...)

ilm_plot_reduce_scree_na(x, ...)
```

## Arguments

- x:

  An
  [`ilm_reduce()`](https://craig-hutton.github.io/illume/reference/ilm_reduce.md)
  result.

- ...:

  Passed to
  [`tinyplot::tinyplot()`](https://grantmcdermott.com/tinyplot/man/tinyplot.html).

## Value

`NULL`, invisibly.

## See also

[`ilm_reduce()`](https://craig-hutton.github.io/illume/reference/ilm_reduce.md).

## Examples

``` r
ilm_plot_reduce_scree(ilm_reduce(mtcars))
```
