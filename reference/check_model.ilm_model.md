# performance::check_model method

Routes
[`performance::check_model()`](https://easystats.github.io/performance/reference/check_model.html)
to
[`ilm_appraise()`](https://huttoncp.github.io/illume/reference/ilm_appraise.md),
so the panels are the ones built for this model rather than generic ones
that would not apply.

## Usage

``` r
check_model.ilm_model(x, ...)
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
