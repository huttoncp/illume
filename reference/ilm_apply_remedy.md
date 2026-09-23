# Refit a model with one of its remedies

Makes one remedy from
[`ilm_remedies()`](https://huttoncp.github.io/illume/reference/ilm_remedies.md):
refits the model through its own call with only that change, so
everything else about it – the family, the weights, a zero part, a
dispersion model, REML – comes along, and the new fit's call is the one
a person would have written, ready for
[`update()`](https://rdrr.io/r/stats/update.html) or a further remedy.
Says what the checks that called for the remedy say now, and which
checks are not OK after the refit.

## Usage

``` r
ilm_apply_remedy(object, remedies, which, data = NULL, verbose = FALSE)
```

## Arguments

- object:

  A fitted `"ilm_model"` object, from the formula interface.

- remedies:

  The list from
  [`ilm_remedies()`](https://huttoncp.github.io/illume/reference/ilm_remedies.md)
  for `object` that the remedy was chosen from. Required rather than
  recomputed, so the remedy made is always the one that was read – a
  list that includes the standalone checks numbers its remedies
  differently from one that does not.

- which:

  The `id` of the remedy to make.

- data:

  The data the model was fitted to, when it cannot be found.

- verbose:

  Logical. Print the new fit's checks as it is fitted.

## Value

The refitted model. Its `remedy_log` holds every remedy applied to reach
it, in order: the check, status, tier, remedy and change.

## Details

The data are found where the model's formula was created, as
[`update()`](https://rdrr.io/r/stats/update.html) would; pass `data`
when they are not there, e.g. when the model was fitted inside a
function from a formula written outside it.

## See also

[`ilm_remedies()`](https://huttoncp.github.io/illume/reference/ilm_remedies.md).

## Examples

``` r
set.seed(1)
d <- data.frame(g = factor(rep(1:20, each = 10)), x = rnorm(200), e = rnorm(200))
d$x <- d$x - ave(d$x, d$g)
d$y <- 1 + 0.5 * d$x + d$e - ave(d$e, d$g)   # no group effect at all
f <- ilm_model(y ~ x + (1 | g), data = d, verbose = FALSE)
#> ilm_model(): family = "gaussian", inferred from `y`: continuous values from -1.71 to 3.31. Pass `family` to choose another.
#> ilm_model(): the random-effect covariance of `g` sits at the edge of its range -- a variance of zero or a correlation of +/-1 -- where the data cannot resolve it. The fixed effects and their standard errors are still usable; summary() says what else is. If the term belongs in the model, boundary = "avoid" keeps it inside its range with a small penalty: it is then assumed nonzero rather than estimated at zero, so do not test whether it is; its variance comes out larger, and for a binary or categorical outcome the fixed effects a little further from zero -- markedly so when a category is rare.
rem <- ilm_remedies(f)
rem
#> 2 remedies
#> 
#> [1] structural -- variance_boundary (WARN)
#>     drop 'g': its variance is estimated at zero, so the fixed effects are
#>     the same without it and their standard errors barely move. Keep it
#>     instead if the design calls for it -- repeated measures, say -- since a
#>     zero estimate is not evidence of no clustering
#>     change: formula = y ~ x
#> 
#> [2] structural -- variance_boundary (WARN)
#>     refit with boundary = "avoid", a small penalty that keeps every
#>     random-effect covariance inside its range. Each variance is then
#>     assumed nonzero rather than estimated at zero, so do not test whether
#>     it is; it comes out larger, and for a binary or categorical outcome the
#>     fixed effects come out a little further from zero -- markedly so when a
#>     category is rare
#>     change: boundary = "avoid"
#> 
#> Refit with one by ilm_apply_remedy(fit, <this list>, id). Numerical is the
#> same model fitted harder; structural changes the random-effect or variance
#> structure and not what the fixed effects mean; estimand changes what they
#> estimate or what their standard errors account for, so apply one of those
#> only by choice.
f2 <- ilm_apply_remedy(f, rem, 1)
#> ilm_apply_remedy(): refitted with formula = y ~ x.
#>   variance_boundary: WARN -> no longer checked, as the term is gone
#>   every check is OK after the refit
f2$remedy_log
#>               check status       tier
#> 1 variance_boundary   WARN structural
#>                                                                                                                                                                                                                                                            remedy
#> 1 drop 'g': its variance is estimated at zero, so the fixed effects are the same without it and their standard errors barely move. Keep it instead if the design calls for it -- repeated measures, say -- since a zero estimate is not evidence of no clustering
#>            change
#> 1 formula = y ~ x
```
