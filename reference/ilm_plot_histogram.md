# Histogram

Histogram

## Usage

``` r
ilm_plot_histogram(data, x, by = NULL, breaks = "Sturges", ...)
```

## Arguments

- data:

  A data frame.

- x:

  Name of the numeric column to plot.

- by:

  Optional grouping column, overlaying one histogram per level.

- breaks:

  Passed to
  [`tinyplot::type_histogram()`](https://grantmcdermott.com/tinyplot/man/type_histogram.html).

- ...:

  Passed to
  [`tinyplot::tinyplot()`](https://grantmcdermott.com/tinyplot/man/tinyplot.html).

## Value

`NULL`, invisibly.

## See also

[`ilm_plot()`](https://craig-hutton.github.io/illume/reference/ilm_plot.md),
which picks a geometry for you.

## Examples

``` r
ilm_plot_histogram(mtcars, "mpg")

ilm_plot_histogram(mtcars, "mpg", by = "cyl")
```
