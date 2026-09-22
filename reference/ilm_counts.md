# Frequency counts of a vector's unique values

Frequency counts of a vector's unique values

## Usage

``` r
ilm_counts(y, n = "all", order = c("d", "a", "i"), na.rm = TRUE)
```

## Arguments

- y:

  A vector.

- n:

  Number of rows to return, or `"all"`.

- order:

  `"d"` by count descending, `"a"` by count ascending, or `"i"` by the
  value itself (factor level order, not label order).

- na.rm:

  Drop missing values before counting.

## Value

A data frame with `value` and `n`.

## See also

[`ilm_counts_tb()`](https://craig-hutton.github.io/illume/reference/ilm_counts_tb.md)
for both ends at once,
[`ilm_counts_all()`](https://craig-hutton.github.io/illume/reference/ilm_counts_all.md)
for a whole data frame.

## Examples

``` r
d <- ilm_sim()
ilm_counts(d$grp)
#>   value   n
#> 1 alpha 404
#> 2  beta 329
#> 3 gamma 164
#> 4 delta   3
ilm_counts(d$site, n = 3)
#>   value   n
#> 1 North 277
#> 2 South 265
#> 3  East 217
```
