# Line plot

Line plot

## Usage

``` r
ilm_plot_line(data, y, x, by = NULL, ...)
```

## Arguments

- data:

  A data frame.

- y, x:

  Column names; `x` is usually a date or a sequence.

- by:

  Optional grouping column.

- ...:

  Passed to
  [`tinyplot::tinyplot()`](https://grantmcdermott.com/tinyplot/man/tinyplot.html).

## Value

`NULL`, invisibly.

## See also

[`ilm_plot_acf()`](https://craig-hutton.github.io/illume/reference/ilm_plot_acf.md)
for what a line plot of residuals cannot show.

## Examples

``` r
d <- data.frame(t = 1:40, v = cumsum(rnorm(40)))
ilm_plot_line(d, "v", "t")
```
