# Progress reporting in illume

Long-running functions take a `progress` argument. It defaults to
[`interactive()`](https://rdrr.io/r/base/interactive.html), so a bar
appears when someone is watching and nothing is written in a script, a
test or a knitted document.

## Details

The bar costs nothing worth measuring. On 2000 bootstrap replicates over
20,000 rows, the loop took 2.48s with a bar and 2.53s without.

Where work is spread over several cores, the bar advances as each
**chunk** of the work returns rather than each replicate: the workers
are separate processes and cannot write to the parent's console. It is
coarser, and it still tells you the run is alive and roughly how far
along.

## Examples

``` r
d <- ilm_sim()
ilm_boot_ci(d, "score", R = 200, progress = FALSE)
#>   stat observed    lower    upper conf   R    ci_type   n
#> 1 mean 49.67677 49.32605 50.04253 0.95 200 percentile 900
```
