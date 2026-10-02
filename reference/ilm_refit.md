# Refit a model to a new response or to other data

The same model – formula, family, random effects, correlation over time,
dispersion and zero parts, fitting options – fitted again, either to a
new response for the rows it was fitted to or to another data set.

## Usage

``` r
ilm_refit(fit, data = NULL, y = NULL)
```

## Arguments

- fit:

  A fitted `"ilm_model"` from the formula interface.

- data:

  A data frame to fit the same model to.

- y:

  A new response for the rows `fit` was fitted to.

## Value

A fitted `"ilm_model"`, with every method a fit from
[`ilm_model()`](https://huttoncp.github.io/illume/reference/ilm_model.md)
has.

## Details

**`y`** is a new response for the rows the model was fitted to, in the
layout
[`ilm_simulate()`](https://huttoncp.github.io/illume/reference/ilm_simulate.md)
returns: a vector with one value per fitted row, or a single column of
its matrix. The design is not rebuilt, so this is the parametric
bootstrap's refit. A censored response is censored as the data were: a
value at or beyond a row's censoring time is censored there, as
[`ilm_simulate()`](https://huttoncp.github.io/illume/reference/ilm_simulate.md)
makes it. A categorical response is given as category numbers, as
[`ilm_simulate()`](https://huttoncp.github.io/illume/reference/ilm_simulate.md)
gives it.

**`data`** is another data set – more rows, fewer, a later window of a
time series – and the model is built from its formula again, so a
factor's levels, a smooth's basis and a correlation over time's cells
come from the new rows. A correlation over time must then have been
given by name, `ilm_ar1(~ time | group)`, so its columns can be found in
the new data.

With neither, the model is refitted to its own data; with both, it
stops.

## See also

[`ilm_simulate()`](https://huttoncp.github.io/illume/reference/ilm_simulate.md),
which draws responses to refit;
[`ilm_apply_remedy()`](https://huttoncp.github.io/illume/reference/ilm_apply_remedy.md),
which refits with one change.

## Examples

``` r
set.seed(1)
d <- data.frame(id = factor(rep(1:12, each = 6)), x = rnorm(72))
d$y <- 1 + 0.5 * d$x + rnorm(12)[d$id] + rnorm(72)
fit <- ilm_model(y ~ x + (1 | id), data = d, family = "gaussian",
                 verbose = FALSE)
## the parametric bootstrap's refit
ystar <- ilm_simulate(fit, nsim = 1, seed = 2)[, 1]
coef(ilm_refit(fit, y = ystar))
#> (Intercept)           x 
#>   0.8747901   0.3474312 
## the same model on the first eight groups
coef(ilm_refit(fit, data = d[d$id %in% 1:8, ]))
#> (Intercept)           x 
#>   0.8355169   0.3989372 
```
