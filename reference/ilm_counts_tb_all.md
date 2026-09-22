# Most and least frequent values for every column

Most and least frequent values for every column

## Usage

``` r
ilm_counts_tb_all(data, n = 10L, na.rm = TRUE)
```

## Arguments

- data:

  A data frame.

- n:

  How many values from each end.

- na.rm:

  Drop missing values before counting.

## Value

A data frame with `variable` and the top/bottom columns.

## Examples

``` r
ilm_counts_tb_all(ilm_sim()[, c("grp", "site")], n = 2)
#>   variable top_value top_n bot_value bot_n
#> 1      grp     alpha   404     delta     3
#> 2      grp      beta   329     gamma   164
#> 3     site     North   277    North     46
#> 4     site     South   265              46
```
