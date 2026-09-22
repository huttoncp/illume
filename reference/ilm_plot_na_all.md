# Missing values by column

A bar per column. The first question to ask of an unfamiliar data set,
and the one
[`ilm_check_missing()`](https://huttoncp.github.io/illume/reference/ilm_check_missing.md)
then turns into advice.

## Usage

``` r
ilm_plot_na_all(data, by = NULL, stat = c("p_na", "na", "n"), ...)
```

## Arguments

- data:

  A data frame.

- by:

  Optional grouping column(s), as a character vector – bars are drawn
  per group, which is how you see whether missingness is concentrated
  somewhere.

- stat:

  `"p_na"` (proportion missing), `"na"` (count missing) or `"n"` (count
  present).

- ...:

  Passed to
  [`tinyplot::tinyplot()`](https://grantmcdermott.com/tinyplot/man/tinyplot.html).

## Value

`NULL`, invisibly.

## See also

[`ilm_check_missing()`](https://huttoncp.github.io/illume/reference/ilm_check_missing.md)
for whether it matters,
[`ilm_profile_na()`](https://huttoncp.github.io/illume/reference/ilm_profile_na.md)
for which columns go missing together.

## Examples

``` r
ilm_plot_na_all(airquality)

ilm_plot_na_all(airquality, by = "Month")
```
