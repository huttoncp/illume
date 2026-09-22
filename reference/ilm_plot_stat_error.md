# Group means or medians with an error bar

A summary per group with a measure of spread around it.

## Usage

``` r
ilm_plot_stat_error(data, y, x, by = NULL, stat = c("mean", "median"), ...)
```

## Arguments

- data:

  A data frame.

- y:

  Name of the numeric column to summarise.

- x:

  Name of the column defining the groups.

- by:

  Optional secondary grouping column.

- stat:

  `"mean"` with its standard error, or `"median"` with the quartiles.

- ...:

  Passed to
  [`tinyplot::tinyplot()`](https://grantmcdermott.com/tinyplot/man/tinyplot.html).

## Value

`NULL`, invisibly.

## Details

The two `stat` options show different things and are not
interchangeable. `"mean"` draws the standard error, which is about where
the *mean* is; `"median"` draws the quartiles, which is about where the
*data* are. The first shrinks as the sample grows and the second does
not.

Overlapping error bars are a poor test of a difference – they are
conservative and lossy.
[`ilm_boot_diff()`](https://craig-hutton.github.io/illume/reference/ilm_boot_diff.md)
gives the difference itself with its own interval.

## See also

[`ilm_boot_diff()`](https://craig-hutton.github.io/illume/reference/ilm_boot_diff.md)
for the comparison this plot invites.

## Examples

``` r
ilm_plot_stat_error(mtcars, "mpg", "cyl")
```
