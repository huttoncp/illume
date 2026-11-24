# Was a fit made by restricted maximum likelihood?

`TRUE` for a fit made by restricted maximum likelihood (REML), `FALSE`
for one made by maximum likelihood.

## Usage

``` r
ilm_is_reml(object)
```

## Arguments

- object:

  A fitted `"ilm_model"` object.

## Value

`TRUE` or `FALSE`. A fit saved by a version of illume that didn't record
it gives `FALSE`.

## Details

When `reml` is not given to
[`ilm_model()`](https://huttoncp.github.io/illume/reference/ilm_model.md),
it decides: REML for a gaussian model with an identity link and a
dispersion formula (other than `~ mu`), and nothing else: no random term
or smooth, no correlation over time, no zero part, no censoring and no
flexible baseline. Every other model is fitted by maximum likelihood.
This reports what the fit ended up with, whichever way it got there. It
matters because a likelihood, and so an AIC or a likelihood ratio, from
a REML fit compares only with fits that have the same fixed effects and
were also fitted by REML, and because a REML fit's coefficients and
standard errors are not those of maximum likelihood.

## See also

[`ilm_model()`](https://huttoncp.github.io/illume/reference/ilm_model.md),
whose `reml` argument sets it.

## Examples

``` r
set.seed(1)
d <- data.frame(g = factor(rep(1:10, each = 6)), x = stats::rnorm(60))
d$y <- 1 + 0.5 * d$x + stats::rnorm(10)[d$g] + stats::rnorm(60)
f_ml <- ilm_model(y ~ x + (1 | g), data = d, family = "gaussian",
                  verbose = FALSE)
f_reml <- ilm_model(y ~ x + (1 | g), data = d, family = "gaussian",
                    reml = TRUE, verbose = FALSE)
ilm_is_reml(f_ml)
#> [1] FALSE
ilm_is_reml(f_reml)
#> [1] TRUE
```
