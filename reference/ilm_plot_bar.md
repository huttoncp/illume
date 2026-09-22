# Bar plot

Bar plot

## Usage

``` r
ilm_plot_bar(data, x, by = NULL, ...)
```

## Arguments

- data:

  A data frame.

- x:

  Name of the categorical column.

- by:

  Optional grouping column.

- ...:

  Passed to
  [`tinyplot::tinyplot()`](https://grantmcdermott.com/tinyplot/man/tinyplot.html).

## Value

`NULL`, invisibly.

## See also

[`ilm_counts()`](https://huttoncp.github.io/illume/reference/ilm_counts.md)
for the same information as a table.

## Examples

``` r
ilm_plot_bar(mtcars, "cyl", by = "am")
```
