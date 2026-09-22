# Which variables a dimension is made of

A sorted dot plot of squared loadings. Deliberately not the classic
correlation circle: that plot asks more visual literacy of a reader than
a ranked list does, and gives no more here.

## Usage

``` r
ilm_plot_reduce_contrib(x, dim = 1, top_n = 10, ...)

ilm_plot_reduce_contrib_na(x, dim = 1, top_n = 10, ...)
```

## Arguments

- x:

  An
  [`ilm_reduce()`](https://craig-hutton.github.io/illume/reference/ilm_reduce.md)
  result.

- dim:

  Which dimension, by number.

- top_n:

  How many variables to show.

- ...:

  Passed to
  [`tinyplot::tinyplot()`](https://grantmcdermott.com/tinyplot/man/tinyplot.html).

## Value

`NULL`, invisibly.

## See also

[`ilm_reduce()`](https://craig-hutton.github.io/illume/reference/ilm_reduce.md).

## Examples

``` r
ilm_plot_reduce_contrib(ilm_reduce(mtcars), dim = 1)
```
