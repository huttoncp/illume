# The most and least frequent values, side by side

The two ends of a frequency distribution are where the problems live: a
level that dominates, and levels too thin to estimate.

## Usage

``` r
ilm_counts_tb(y, n = 10L, na.rm = TRUE)
```

## Arguments

- y:

  A vector.

- n:

  How many values from each end.

- na.rm:

  Drop missing values before counting.

## Value

A data frame with `top_value`, `top_n`, `bot_value` and `bot_n`.

## Examples

``` r
ilm_counts_tb(ilm_sim()$grp, n = 2)
#>   top_value top_n bot_value bot_n
#> 1     alpha   404     delta     3
#> 2      beta   329     gamma   164
```
