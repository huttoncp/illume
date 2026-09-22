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
