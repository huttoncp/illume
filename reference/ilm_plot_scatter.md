# Scatter plot

Scatter plot

## Usage

``` r
ilm_plot_scatter(data, y, x, by = NULL, trend = c("none", "lm", "loess"), ...)
```

## Arguments

- data:

  A data frame.

- y, x:

  Names of the numeric columns.

- by:

  Optional grouping column.

- trend:

  `"none"`, `"lm"` or `"loess"`. A trend line here is a description of
  the two columns shown and nothing more – it holds nothing else fixed,
  so it is not the effect
  [`ilm_model()`](https://craig-hutton.github.io/illume/reference/ilm_model.md)
  would estimate.

- ...:

  Passed to
  [`tinyplot::tinyplot()`](https://grantmcdermott.com/tinyplot/man/tinyplot.html).

## Value

`NULL`, invisibly.

## See also

[`ilm_plot_var_pairs()`](https://craig-hutton.github.io/illume/reference/ilm_plot_var_pairs.md)
for every pair at once.

## Examples

``` r
ilm_plot_scatter(mtcars, "mpg", "wt", by = "cyl", trend = "lm")
#> Warning: 
#> Continuous legends not supported for this plot type. Reverting to discrete legend.
```
