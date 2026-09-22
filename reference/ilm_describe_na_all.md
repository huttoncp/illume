# Missingness in every variable

Sorted with the most missing first, since those are the variables worth
looking at.

## Usage

``` r
ilm_describe_na_all(data, by = NULL, digits = 4, sort = TRUE)
```

## Arguments

- data:

  A data frame, or a vector when `y` is `NULL`.

- by:

  Optional grouping columns.

- digits:

  Rounding for `p_na`.

- sort:

  Sort by proportion missing, descending.

## Value

A data frame with `variable`, `obs`, `n`, `na` and `p_na`.

## Examples

``` r
ilm_describe_na_all(ilm_sim())
#>     variable obs   n  na p_na
#> 1  lab_value 900 792 108 0.12
#> 2     claims 900 900   0 0.00
#> 3     cohort 900 900   0 0.00
#> 4  consented 900 900   0 0.00
#> 5       date 900 900   0 0.00
#> 6   downtime 900 900   0 0.00
#> 7       flag 900 900   0 0.00
#> 8        grp 900 900   0 0.00
#> 9         id 900 900   0 0.00
#> 10    income 900 900   0 0.00
#> 11     score 900 900   0 0.00
#> 12      site 900 900   0 0.00
#> 13    visits 900 900   0 0.00
```
