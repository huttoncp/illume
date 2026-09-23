# Remedies for what a model's checks found

Lists a remedy for every check of a fitted model that is not `"OK"`,
each one written out as the change to
[`ilm_model()`](https://huttoncp.github.io/illume/reference/ilm_model.md)
that makes it, so that it can be made with
[`ilm_apply_remedy()`](https://huttoncp.github.io/illume/reference/ilm_apply_remedy.md)
rather than retyped. The checks made at every fit are read from
`object$checks`. The ones that simulate –
[`ilm_check_dispersion()`](https://huttoncp.github.io/illume/reference/ilm_check_dispersion.md),
[`ilm_check_zeros()`](https://huttoncp.github.io/illume/reference/ilm_check_zeros.md)
and
[`ilm_check_variance()`](https://huttoncp.github.io/illume/reference/ilm_check_variance.md)
– are passed in, since running them is a choice.

## Usage

``` r
ilm_remedies(object, dispersion = NULL, zeros = NULL, variance = NULL)

# S3 method for class 'ilm_remedies'
print(x, ...)
```

## Arguments

- object:

  A fitted `"ilm_model"` object.

- dispersion, zeros, variance:

  Optional results of
  [`ilm_check_dispersion()`](https://huttoncp.github.io/illume/reference/ilm_check_dispersion.md),
  [`ilm_check_zeros()`](https://huttoncp.github.io/illume/reference/ilm_check_zeros.md)
  and
  [`ilm_check_variance()`](https://huttoncp.github.io/illume/reference/ilm_check_variance.md)
  run on `object`, whose remedies are then listed too.

- x:

  An `"ilm_remedies"` object.

- ...:

  Unused.

## Value

A data frame of class `"ilm_remedies"`, one row per remedy, with `id`,
the `check` (or checks) it answers and its `status`, the `tier`, the
`remedy` in words, and the `change`: the
[`ilm_model()`](https://huttoncp.github.io/illume/reference/ilm_model.md)
arguments that make it, as code, or `""` when it is made by hand. An
empty data frame when every check is OK.

## Details

Each remedy has a tier, which says what applying it would change:

- `numerical`:

  The same model, fitted harder: more optimiser restarts.

- `structural`:

  A different random-effect or variance structure, whose fixed effects
  mean what they meant before: a term whose variance is estimated at
  zero removed (the fit is the same without it), a covariance of lower
  rank, a dispersion model, a zero part, a negative binomial in place of
  a poisson, or the boundary-avoiding penalty with its measured costs.

- `estimand`:

  A change to what the fixed effects estimate, or to what their standard
  errors account for: a random effect whose variance is not zero
  removed, categories merged, levels pooled. Apply one only because the
  question calls for it, not because a check did.

Remedies are listed numerical first, then structural, then estimand.
Some can only be made by hand – which categories to merge is a question
about what they mean – and those have no `change`. A remedy is a
candidate, not a cure: refit, and read the checks of the new fit.

## See also

[`ilm_apply_remedy()`](https://huttoncp.github.io/illume/reference/ilm_apply_remedy.md)
to refit with one,
[`summary.ilm_model()`](https://huttoncp.github.io/illume/reference/summary.ilm_model.md)
for the checks themselves.

## Examples

``` r
## x and y vary within groups only, so the groups differ by nothing and
## the random intercept's variance is estimated at zero
set.seed(1)
d <- data.frame(g = factor(rep(1:20, each = 10)), x = rnorm(200), e = rnorm(200))
d$x <- d$x - ave(d$x, d$g)
d$y <- 1 + 0.5 * d$x + d$e - ave(d$e, d$g)
f <- ilm_model(y ~ x + (1 | g), data = d, verbose = FALSE)
#> ilm_model(): family = "gaussian", inferred from `y`: continuous values from -1.71 to 3.31. Pass `family` to choose another.
#> ilm_model(): the random-effect covariance of `g` sits at the edge of its range -- a variance of zero or a correlation of +/-1 -- where the data cannot resolve it. The fixed effects and their standard errors are still usable; summary() says what else is. If the term belongs in the model, boundary = "avoid" keeps it inside its range with a small penalty: it is then assumed nonzero rather than estimated at zero, so do not test whether it is; its variance comes out larger, and for a binary or categorical outcome the fixed effects a little further from zero -- markedly so when a category is rare.
ilm_remedies(f)
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
```
