# Violin plot

Violin plot

## Usage

``` r
ilm_plot_violin(data, y, x = NULL, by = NULL, ..., pch = NULL)
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

- pch:

  Plotting character. Takes a NAME as well as a number:
  `"filled circle"` is 16, and every code from 0 to 25 has one. Case,
  spaces, underscores and hyphens are ignored. A single character is
  drawn literally, so `pch = "x"` is still the letter x.

## Value

`NULL`, invisibly.

## See also

[`ilm_plot_box()`](https://huttoncp.github.io/illume/reference/ilm_plot_box.md),
which shows the quartiles rather than the shape.

## Examples

``` r
ilm_plot_violin(mtcars, "mpg", x = "cyl")
```
