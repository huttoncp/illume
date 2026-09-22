# Density plot

Density plot

## Usage

``` r
ilm_plot_density(data, x, by = NULL, ...)
```

## Arguments

- data:

  A data frame.

- x:

  Name of the numeric column to plot.

- by:

  Optional grouping column, overlaying one histogram per level.

- ...:

  Passed to
  [`tinyplot::tinyplot()`](https://grantmcdermott.com/tinyplot/man/tinyplot.html).

## Value

`NULL`, invisibly.

## See also

[`ilm_plot_histogram()`](https://craig-hutton.github.io/illume/reference/ilm_plot_histogram.md).

## Examples

``` r
ilm_plot_density(mtcars, "mpg", by = "cyl")
```
