# Missing values in one column, across groups

Where
[`ilm_plot_na_all()`](https://huttoncp.github.io/illume/reference/ilm_plot_na_all.md)
compares columns, this compares groups within one column – which is how
a pattern in *who* is missing shows itself.

## Usage

``` r
ilm_plot_na(data, x, by = NULL, stat = c("p_na", "na", "n"), ...)
```

## Arguments

- data:

  A data frame.

- x:

  Name of the column whose missingness to show.

- by:

  Grouping column(s) to split by, as a character vector. Required:
  without it there is only one bar, which is
  [`ilm_plot_na_all()`](https://huttoncp.github.io/illume/reference/ilm_plot_na_all.md)'s
  job.

- stat:

  `"p_na"`, `"na"` or `"n"`.

- ...:

  Passed to
  [`tinyplot::tinyplot()`](https://grantmcdermott.com/tinyplot/man/tinyplot.html).

## Value

`NULL`, invisibly.

## See also

[`ilm_plot_na_all()`](https://huttoncp.github.io/illume/reference/ilm_plot_na_all.md),
[`ilm_check_missing()`](https://huttoncp.github.io/illume/reference/ilm_check_missing.md).

## Examples

``` r
ilm_plot_na(airquality, "Ozone", by = "Month")
```
