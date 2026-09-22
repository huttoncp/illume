# Make illume models readable by the easystats packages

Registers the `insight` methods an `ilm_model` needs, after which
`parameters::model_parameters()`,
[`performance::model_performance()`](https://easystats.github.io/performance/reference/model_performance.html)
and `report::report()` work on one. `insight` stays in Suggests, so this
is called rather than declared.

## Usage

``` r
ilm_register_insight()
```

## Value

`TRUE` invisibly if the methods were registered, `FALSE` if `insight` is
not installed.

## Details

This is not a substitute for
[`ilm_interpret()`](https://craig-hutton.github.io/illume/reference/ilm_interpret.md),
which says things the easystats stack cannot – what illume's own
diagnostics found, and what that implies for the estimates. It is there
because it costs a user nothing to have the option, and because an
independent implementation of the same quantities is worth having to
check against.

## See also

[`ilm_register_marginaleffects()`](https://craig-hutton.github.io/illume/reference/ilm_register_marginaleffects.md),
[`ilm_interpret()`](https://craig-hutton.github.io/illume/reference/ilm_interpret.md).

## Examples

``` r
ilm_register_insight()
```
