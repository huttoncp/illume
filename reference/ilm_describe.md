# Class-aware description of one variable

Summarises a vector, or one column of a data frame, with statistics
chosen for its class. Numeric variables also get a gaussian agreement
index and, where that is low, a plain-language reason.

## Usage

``` r
ilm_describe(
  data,
  y = NULL,
  by = NULL,
  digits = 3,
  gauss = c("index", "ks_d", "both", "none"),
  probs = c(0, 0.5, 1),
  skew = FALSE,
  kurt = FALSE,
  dispersion = TRUE,
  rare_n = 5L,
  cap = 0.12,
  min_n = 20L
)
```

## Arguments

- data:

  A data frame, or a vector when `y` is `NULL`.

- y:

  Name of the column to summarise, or several names. With none, `data`
  is taken as the thing to describe: a bare vector, or a data frame
  whose columns are all of one kind. For a frame of mixed kinds use
  [`ilm_describe_all()`](https://huttoncp.github.io/illume/reference/ilm_describe_all.md),
  which returns one table per kind.

- by:

  Optional character vector of grouping columns.

- digits:

  Rounding for numeric columns, quantiles included.

- gauss:

  One of `"index"` (the 0-1 agreement index), `"ks_d"` (the raw
  Kolmogorov distance behind it), `"both"`, or `"none"`.

- probs:

  Quantiles to report, as proportions. `c(0, 0.5, 1)` gives `p0`, `p50`
  and `p100`; any vector works and the column names follow it.

- skew, kurt:

  Add moment-based skewness and excess kurtosis.

- dispersion:

  Add the variance-to-mean ratio for non-negative integer variables.
  Meaningful only if you intend to model the variable as a count.

- rare_n:

  Levels with fewer observations than this count as rare.

- cap, min_n:

  Control the `gauss` index; see
  [`ilm_gauss_check()`](https://huttoncp.github.io/illume/reference/ilm_gauss_check.md).

## Value

A one-row data frame, or one row per group when `by` is used.

## Details

The default output is deliberately narrow. Quartiles, skewness and
kurtosis are available on request, but `gauss` and `gauss_note` already
say whether they are worth asking for.

## See also

[`ilm_describe_all()`](https://huttoncp.github.io/illume/reference/ilm_describe_all.md)
for every column,
[`ilm_gauss_check()`](https://huttoncp.github.io/illume/reference/ilm_gauss_check.md)
for the index,
[`ilm_describe_na()`](https://huttoncp.github.io/illume/reference/ilm_describe_na.md)
for missingness.

## Examples

``` r
d <- ilm_sim()
ilm_describe(d, "score")
#>   obs   n na      sum   mean   sd    se    p0    p50  p100 p_zero dispersion
#> 1 900 900  0 44709.09 49.677 5.88 0.196 29.99 49.835 65.34      0         NA
#>   gauss gauss_note
#> 1     1           
ilm_describe(d, "income", probs = c(0.1, 0.5, 0.9))
#>   obs   n na      sum     mean       sd      se      p10      p50      p90
#> 1 900 900  0 33361610 37068.46 22902.68 763.423 14765.58 31137.27 66193.31
#>   p_zero dispersion gauss                 gauss_note
#> 1      0         NA 0.305 right-skewed; heavy-tailed
ilm_describe(d, "score", by = "grp")
#>         grp obs   n na      sum   mean    sd    se    p0    p50  p100 p_zero
#> alpha alpha 404 404  0 20207.23 50.018 5.632 0.280 29.99 50.090 65.34      0
#> beta   beta 329 329  0 16200.44 49.241 6.170 0.340 33.83 49.330 64.94      0
#> gamma gamma 164 164  0  8141.46 49.643 5.887 0.460 31.60 49.335 63.91      0
#> delta delta   3   3  0   159.96 53.320 2.537 1.465 50.63 53.660 55.67      0
#>       dispersion gauss            gauss_note
#> alpha         NA     1                      
#> beta          NA     1                      
#> gamma         NA     1                      
#> delta         NA    NA n too small to assess
```
