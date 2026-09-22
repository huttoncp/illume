# Pairwise plots

Every pair of columns at once, handling mixed numeric and categorical
columns rather than only numeric ones as a classic scatterplot matrix
does.

## Usage

``` r
ilm_plot_var_pairs(data, cols = NULL, by = NULL, ...)
```

## Arguments

- data:

  A data frame.

- cols:

  Columns to use. A character vector of names, a regular expression, a
  predicate function such as `is.numeric`, or `NULL` for all of them –
  see
  [ilm_selection](https://huttoncp.github.io/illume/reference/ilm_selection.md).

- by:

  Optional grouping column.

- ...:

  Passed to
  [`tinyplot::tinypairs()`](https://grantmcdermott.com/tinyplot/man/tinyplot.data.frame.html),
  which needs tinyplot 0.7.0 or later. On an earlier tinyplot this is
  the one plot in the package that cannot be drawn, and it says so.

## Value

`NULL`, invisibly.

## See also

[`ilm_plot_scatter()`](https://huttoncp.github.io/illume/reference/ilm_plot_scatter.md)
for one pair,
[`ilm_check_collinearity()`](https://huttoncp.github.io/illume/reference/ilm_check_collinearity.md)
for what a pairs plot cannot show about a fit.

## Examples

``` r
ilm_plot_var_pairs(mtcars, cols = c("mpg", "wt", "hp"), by = "cyl")
```
