# What each geom requires

The table the argument checks read, so the documented requirements and
the enforced ones cannot drift apart.

## Usage

``` r
ilm_geom_spec(geom = NULL)
```

## Arguments

- geom:

  Optional geom name or names; all of them if `NULL`.

## Value

A data frame with `geom`, `needs_y`, `x_class`, `y_class` and
`size_controls`.

## Examples

``` r
ilm_geom_spec()
#>        geom needs_y                         x_class
#> 1 histogram   FALSE                         numeric
#> 2   density   FALSE                         numeric
#> 3       bar   FALSE           categorical / logical
#> 4       box    TRUE categorical / logical / numeric
#> 5    violin    TRUE categorical / logical / numeric
#> 6     point    TRUE                         numeric
#> 7     bin2d    TRUE                         numeric
#> 8      line    TRUE                  time / numeric
#> 9     spine    TRUE                     categorical
#>                           y_class     size_controls
#> 1                        not used line/border width
#> 2                        not used line/border width
#> 3                        not used line/border width
#> 4 numeric / categorical / logical line/border width
#> 5 numeric / categorical / logical line/border width
#> 6                         numeric        point size
#> 7                         numeric                 -
#> 8           numeric / categorical line/border width
#> 9                     categorical line/border width
ilm_geom_spec("point")
#>    geom needs_y x_class y_class size_controls
#> 1 point    TRUE numeric numeric    point size
```
