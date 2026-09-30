# Register the marginaleffects interface

Call once per session before using `marginaleffects` with an
`"ilm_model"` fit.

## Usage

``` r
ilm_register_marginaleffects()
```

## Value

`TRUE` if `marginaleffects` is installed, `FALSE` otherwise, invisibly.

## Details

Registering the S3 methods alone is not enough. `marginaleffects` checks
model classes against a fixed list and rejects anything unfamiliar
*before* dispatch happens, so the methods would never be reached. This
also adds `"ilm_model"` to that list through the option the package
provides for the purpose.

**A model with an offset.** `marginaleffects` builds its own grids from
the data, so its predictions and effects are at each row's own exposure
– the counts, not the rates that
[`ilm_emmeans()`](https://huttoncp.github.io/illume/reference/ilm_emmeans.md),
[`ilm_ame()`](https://huttoncp.github.io/illume/reference/ilm_ame.md)
and
[`ilm_scenario()`](https://huttoncp.github.io/illume/reference/ilm_scenario.md)
report per unit of exposure. For rates, give the exposure column the
value 1 in the grid (`newdata = datagrid(e = 1)`, for an offset
`offset(log(e))`), or use illume's own functions.
