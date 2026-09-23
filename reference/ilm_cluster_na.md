# Cluster observations by which values they are missing

The missingness counterpart to
[`ilm_cluster()`](https://huttoncp.github.io/illume/reference/ilm_cluster.md):
groups the dimensions that
[`ilm_reduce_na()`](https://huttoncp.github.io/illume/reference/ilm_reduce_na.md)
produced, so that rows missing the same columns end up together. A
type-checked wrapper – the clustering is identical, only the input
differs.

## Usage

``` r
ilm_cluster_na(x, ...)
```

## Arguments

- x:

  An
  [`ilm_reduce_na()`](https://huttoncp.github.io/illume/reference/ilm_reduce_na.md)
  result.

- ...:

  Passed to
  [`ilm_cluster()`](https://huttoncp.github.io/illume/reference/ilm_cluster.md).

## Value

An object of class `"ilm_cluster_na"`, which is also an `"ilm_cluster"`.

## See also

[`ilm_profile_na()`](https://huttoncp.github.io/illume/reference/ilm_profile_na.md).

## Examples

``` r
cl <- ilm_cluster_na(ilm_reduce_na(airquality), k_max = 4, B = 25, seed = 1)
#> ilm_reduce_na(): dropping column(s) whose missingness never varies (always or never missing): Wind, Temp, Month, Day
#> Warning: k was chosen as 4, which is the largest value searched. The curve had not turned, so this is where the search stopped rather than where the evidence pointed. Raise `k_max`, or set `k` from what the design says. On mixed data a selector can also lock onto the number of category combinations rather than the number of clusters; plot(x) shows the gap curve.
cl
#> <ilm_cluster_na> method = kmeans, k = 4 (chosen by gap statistic)
#> 
#>   cluster   size    pct  jaccard  stability     sil
#>         1     35   22.9    1.000     stable   1.000
#>         2      2    1.3    1.000     stable   1.000  small
#>         3    111   72.5    0.973     stable   1.000
#>         4      5    3.3    1.000     stable   1.000  small
#> 
#>   7 observation(s) in a cluster holding less than 5% of the data
```
