# Missingness in one variable

Missingness in one variable

## Usage

``` r
ilm_describe_na(data, y = NULL, by = NULL, digits = 4)
```

## Arguments

- data:

  A data frame, or a vector when `y` is `NULL`.

- y:

  Name of the column.

- by:

  Optional grouping columns.

- digits:

  Rounding for `p_na`.

## Value

A data frame with `obs`, `n`, `na` and `p_na`.

## See also

[`ilm_describe_na_all()`](https://huttoncp.github.io/illume/reference/ilm_describe_na_all.md),
[`ilm_plot_missing()`](https://huttoncp.github.io/illume/reference/ilm_plot_missing.md).

## Examples

``` r
ilm_describe_na(ilm_sim(), "lab_value")
#>   obs   n  na p_na
#> 1 900 792 108 0.12
```
