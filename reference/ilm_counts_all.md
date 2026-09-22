# Frequency counts for every column

Values are stacked across columns of different classes, so they are
rendered as character; the counts stay integer.

## Usage

``` r
ilm_counts_all(
  data,
  by = NULL,
  n = "all",
  order = c("d", "a", "i"),
  na.rm = TRUE
)
```

## Arguments

- data:

  A data frame.

- by:

  Optional grouping columns.

- n:

  Number of rows to return, or `"all"`.

- order:

  `"d"` by count descending, `"a"` by count ascending, or `"i"` by the
  value itself (factor level order, not label order).

- na.rm:

  Drop missing values before counting.

## Value

A data frame with `variable`, `value` and `n`.

## Examples

``` r
ilm_counts_all(ilm_sim()[, c("grp", "site")], n = 2)
#>   variable value   n
#> 1      grp alpha 404
#> 2      grp  beta 329
#> 3     site North 277
#> 4     site South 265
```
