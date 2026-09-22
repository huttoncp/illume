# Which geom would be drawn, and why

Returns the choice without drawing anything, so the dispatch can be
checked before committing to a plot of many million rows.

## Usage

``` r
ilm_pick_geom(x, y = NULL, n_max = 5000L)
```

## Arguments

- x:

  A vector.

- y:

  An optional second vector.

- n_max:

  Above this many points a scatter becomes a binned density.

## Value

A list with `geom`, `reason` and `n`.

## Examples

``` r
d <- ilm_sim()
ilm_pick_geom(d$score)
#> $geom
#> [1] "histogram"
#> 
#> $reason
#> [1] "one numeric variable"
#> 
#> $n
#> [1] 900
#> 
ilm_pick_geom(d$score, d$income)
#> $geom
#> [1] "point"
#> 
#> $reason
#> [1] "numeric vs numeric"
#> 
#> $n
#> [1] 900
#> 
```
