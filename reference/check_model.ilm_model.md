# performance::check_model method

`performance::check_model(fit)` draws
[`ilm_appraise()`](https://huttoncp.github.io/illume/reference/ilm_appraise.md)'s
panels. Its own checks cannot be drawn for these fits: they need support
from the insight package that an `"ilm_model"` does not have, and
without this method `check_model()` stops with an error – after
[`ilm_register_insight()`](https://huttoncp.github.io/illume/reference/ilm_register_insight.md)
too. The panels drawn instead are illume's, built for every family it
fits, including a nominal outcome, where the usual residual plots have
no clear meaning.

## Usage

``` r
# S3 method for class 'ilm_model'
check_model(x, ...)
```

## Arguments

- x:

  A fitted `"ilm_model"` object.

- ...:

  Passed to
  [`ilm_appraise()`](https://huttoncp.github.io/illume/reference/ilm_appraise.md).

## Value

Invisibly, the result of
[`ilm_appraise()`](https://huttoncp.github.io/illume/reference/ilm_appraise.md).

## Details

Registered with performance's generic when performance is loaded.
