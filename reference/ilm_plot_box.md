# Boxplot

The whiskers use the same 1.5 interquartile-range fence that
[`ilm_outliers()`](https://huttoncp.github.io/illume/reference/ilm_outliers.md)
scores against, so the points beyond them are the values that function
flags under `method = "iqr"`.

## Usage

``` r
ilm_plot_box(data, y, x = NULL, by = NULL, ..., pch = NULL)
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

[`ilm_outliers()`](https://huttoncp.github.io/illume/reference/ilm_outliers.md),
[`ilm_plot_violin()`](https://huttoncp.github.io/illume/reference/ilm_plot_violin.md).

## Examples

``` r
ilm_plot_box(mtcars, "mpg", x = "cyl")
```
