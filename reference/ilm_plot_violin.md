# Violin plot

Violin plot

## Usage

``` r
ilm_plot_violin(data, y, x = NULL, by = NULL, ...)
```

## Arguments

- data:

  A data frame.

- y:

  Name of the numeric column.

- x:

  Optional categorical column to split by.

- by:

  Optional grouping column.

- ...:

  Passed to
  [`tinyplot::tinyplot()`](https://grantmcdermott.com/tinyplot/man/tinyplot.html).

## Value

`NULL`, invisibly.

## See also

[`ilm_plot_box()`](https://craig-hutton.github.io/illume/reference/ilm_plot_box.md),
which shows the quartiles rather than the shape.

## Examples

``` r
ilm_plot_violin(mtcars, "mpg", x = "cyl")
```
