# Flag unusual values across a data frame

Runs
[`ilm_outliers()`](https://huttoncp.github.io/illume/reference/ilm_outliers.md)
over every numeric column, optionally within groups, and returns one row
per flagged value with the row it came from.

## Usage

``` r
ilm_outliers_all(
  data,
  by = NULL,
  cols = NULL,
  method = c("iqr", "mad", "zscore"),
  threshold = NULL,
  flagged_only = TRUE,
  na.rm = TRUE
)
```

## Arguments

- data:

  A data frame.

- by:

  Grouping column(s), as a character vector. Reference statistics, and
  so the flags, are computed separately within each group.

- cols:

  Columns to use. A character vector of names, a regular expression, a
  predicate function such as `is.numeric`, or `NULL` for all of them –
  see
  [ilm_selection](https://huttoncp.github.io/illume/reference/ilm_selection.md).

- method:

  `"iqr"`, `"mad"` or `"zscore"`.

- threshold:

  Flagging threshold. Defaults per rule: 1.5 for `"iqr"` (the boxplot
  convention), 3.5 for `"mad"` (Iglewicz and Hoaglin), 3 for `"zscore"`.

- flagged_only:

  Return only the flagged values. `FALSE` returns every value's score.

- na.rm:

  Compute the reference statistics with missing values removed. The
  result always has one row per element of `y`, with `NA` for `score`
  and `is_outlier` wherever `y` is `NA`.

## Value

A data frame with the `by` columns, `row_id` (the row's position in
`data`), `variable`, `value`, `score` and `is_outlier`.

## Details

Grouping matters more than it looks. A value can be perfectly ordinary
for its own group and extreme against the pooled distribution, so
flagging without `by` on data that has groups mostly rediscovers the
groups.

## See also

[`ilm_outliers()`](https://huttoncp.github.io/illume/reference/ilm_outliers.md),
[`ilm_plot_box()`](https://huttoncp.github.io/illume/reference/ilm_plot_box.md)
to see them.

## Examples

``` r
head(ilm_outliers_all(mtcars))
#>   row_id variable   value score is_outlier
#> 1     31       hp 335.000 1.845       TRUE
#> 2     16       wt   5.424 1.602       TRUE
#> 3     17       wt   5.345 1.530       TRUE
#> 4      9     qsec  22.900 1.985       TRUE
#> 5     31     carb   8.000 2.000       TRUE
head(ilm_outliers_all(mtcars, by = "cyl", method = "mad"))
#>   cyl row_id variable  value  score is_outlier
#> 1   6      4     disp 258.00  8.023       TRUE
#> 2   6      6     disp 225.00  5.094       TRUE
#> 3   6     30       hp 175.00  8.768       TRUE
#> 4   6      4     drat   3.08 27.654       TRUE
#> 5   6      6     drat   2.76 38.446       TRUE
#> 6   6     30     drat   3.62  9.443       TRUE
```
