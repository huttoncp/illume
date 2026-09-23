# The rows an anomaly scan flagged

Pulls the flagged rows out of an
[`ilm_anomaly()`](https://huttoncp.github.io/illume/reference/ilm_anomaly.md)
result, as a data frame, so they can be explored with anything that
takes one.

## Usage

``` r
ilm_anomalous(x, data = NULL, flagged = TRUE, score = TRUE)
```

## Arguments

- x:

  An
  [`ilm_anomaly()`](https://huttoncp.github.io/illume/reference/ilm_anomaly.md)
  result.

- data:

  The data frame that was scanned. Only needed if the result did not
  keep it – see `keep_data` in
  [`ilm_anomaly()`](https://huttoncp.github.io/illume/reference/ilm_anomaly.md).

- flagged:

  Return only the flagged rows (the default), or all of them with the
  score attached.

- score:

  Add the `score`, `p_adj` and `driver` columns.

## Value

A data frame of the flagged rows, in their original order, with `.row`
giving their position in the input.

## Details

The point of a scan is rarely the score. It is whether the odd rows are
alike: one mistake made repeatedly, a subpopulation the model does not
cover, or several unrelated things. That is a clustering or a
description question, and this is what hands the rows over.

## See also

[`ilm_anomaly()`](https://huttoncp.github.io/illume/reference/ilm_anomaly.md),
[`ilm_profile()`](https://huttoncp.github.io/illume/reference/ilm_profile.md),
[`ilm_cluster()`](https://huttoncp.github.io/illume/reference/ilm_cluster.md),
[`ilm_describe_all()`](https://huttoncp.github.io/illume/reference/ilm_describe_all.md).

## Examples

``` r
# \donttest{
a <- ilm_anomaly(mtcars)
odd <- ilm_anomalous(a)
ilm_describe_all(odd)
#>    variable       class obs n na value
#> 1       mpg     numeric   0 0  0  <NA>
#> 2       cyl     numeric   0 0  0  <NA>
#> 3      disp     numeric   0 0  0  <NA>
#> 4        hp     numeric   0 0  0  <NA>
#> 5      drat     numeric   0 0  0  <NA>
#> 6        wt     numeric   0 0  0  <NA>
#> 7      qsec     numeric   0 0  0  <NA>
#> 8        vs     numeric   0 0  0  <NA>
#> 9        am     numeric   0 0  0  <NA>
#> 10     gear     numeric   0 0  0  <NA>
#> 11     carb     numeric   0 0  0  <NA>
#> 12     .row     numeric   0 0  0  <NA>
#> 13   .score     numeric   0 0  0  <NA>
#> 14   .p_adj     numeric   0 0  0  <NA>
#> 15  .driver categorical   0 0  0  <NA>
# }
```
