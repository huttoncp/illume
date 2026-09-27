# A table of remedies from another package's diagnostics

For a package whose diagnostics find something wrong with an illume
model and name a remedy that is a refit of it – a negative binomial for
counts more variable than a poisson allows, a dispersion or zero part –
so that its users apply those remedies as they apply illume's own, with
[`ilm_apply_remedy()`](https://huttoncp.github.io/illume/reference/ilm_apply_remedy.md).
Its
[`ilm_remedies()`](https://huttoncp.github.io/illume/reference/ilm_remedies.md)
method builds the table here:

## Usage

``` r
ilm_remedy_table(object, check, status, tier, remedy, args = NULL)

# S3 method for class 'ilm_remedies'
c(...)
```

## Arguments

- object:

  The fitted `"ilm_model"` the remedies would refit.

- check:

  Character: the name of the check each remedy answers. Not one of the
  fit's own checks, nor an `ilm_check_` name – those are
  [`ilm_remedies()`](https://huttoncp.github.io/illume/reference/ilm_remedies.md)'s,
  and the report after a refit reads them from the fit.

- status:

  Character: what the check found, `"WARN"`, `"FAIL"`, `"BOUNDARY"` or
  `"INCONCLUSIVE"`.

- tier:

  Character: `"numerical"`, `"structural"` or `"estimand"`, as in
  [`ilm_remedies()`](https://huttoncp.github.io/illume/reference/ilm_remedies.md).

- remedy:

  Character: the remedy in words, a sentence a person can act on.

- args:

  A list with one element per row: a named list of
  [`ilm_model()`](https://huttoncp.github.io/illume/reference/ilm_model.md)
  arguments that make the remedy, or `NULL` for one made by hand. `NULL`
  makes every row by hand. The data and `verbose` are not arguments a
  remedy sets. A formula is evaluated where the model's own formula was
  written, so names in it mean what they mean there.

- ...:

  For [`c()`](https://rdrr.io/r/base/c.html): remedy tables for the same
  fit, from
  [`ilm_remedies()`](https://huttoncp.github.io/illume/reference/ilm_remedies.md)
  or `ilm_remedy_table()`.

## Value

A data frame of class `"ilm_remedies"`, as from
[`ilm_remedies()`](https://huttoncp.github.io/illume/reference/ilm_remedies.md).

## Details

    ilm_remedies.my_calibration <- function(object, ...)
      ilm_remedy_table(object$fit, check = "pit_shape", status = "FAIL",
                       tier = "structural", remedy = "...",
                       args = list(list(family = "nbinom")))

The table is tied to `object` as illume's own are: listed, merged when
several checks name the same change, and ordered by tier the same way,
and refused by
[`ilm_apply_remedy()`](https://huttoncp.github.io/illume/reference/ilm_apply_remedy.md)
for any other fit. The `change` column is written from `args`, never
passed in, so the change a person reads is the refit that is made.

A remedy that is not a refit of the model – more simulations, another
kind of fold, a recalibration – is not something
[`ilm_apply_remedy()`](https://huttoncp.github.io/illume/reference/ilm_apply_remedy.md)
can make. List it by hand, with `args = NULL` for its row, or leave it
to the package's own output.

[`c()`](https://rdrr.io/r/base/c.html) combines remedy tables for the
same fit – illume's own and other packages' – into one list, numbered
afresh: a change several lists name is listed once, with every check
that named it, and the remedies are ordered by tier as in each list
alone. Where two lists put one change in different tiers it takes the
more cautious. Tables for different fits are not combined.

## See also

[`ilm_remedies()`](https://huttoncp.github.io/illume/reference/ilm_remedies.md),
[`ilm_apply_remedy()`](https://huttoncp.github.io/illume/reference/ilm_apply_remedy.md).

## Examples

``` r
set.seed(2)
d <- data.frame(x = rnorm(300))
d$y <- rnbinom(300, mu = exp(1 + 0.5 * d$x), size = 1.5)
f <- ilm_model(y ~ x, data = d, family = "poisson", verbose = FALSE)
## what another package's diagnostic of this fit would hand back
rem <- ilm_remedy_table(f, check = "pit_shape", status = "FAIL",
  tier = "structural",
  remedy = "the intervals are too narrow: refit as a negative binomial",
  args = list(list(family = "nbinom")))
rem
#> 1 remedy
#> 
#> [1] structural -- pit_shape (FAIL)
#>     the intervals are too narrow: refit as a negative binomial
#>     change: family = "nbinom"
#> 
#> Refit with one by ilm_apply_remedy(fit, <this list>, id). Numerical is the
#> same model fitted harder; structural changes the random-effect or variance
#> structure, or how the variances are estimated, and not what the fixed
#> effects mean; estimand changes what they estimate or what their standard
#> errors account for, so apply one of those only by choice.
f2 <- ilm_apply_remedy(f, rem, 1)
#> ilm_apply_remedy(): refitted with family = "nbinom".
#>   pit_shape: FAIL before; run it again on the new fit
#>   every check is OK after the refit
## one list, illume's remedies and the other package's together
c(ilm_remedies(f), rem)
#> 1 remedy
#> 
#> [1] structural -- pit_shape (FAIL)
#>     the intervals are too narrow: refit as a negative binomial
#>     change: family = "nbinom"
#> 
#> Refit with one by ilm_apply_remedy(fit, <this list>, id). Numerical is the
#> same model fitted harder; structural changes the random-effect or variance
#> structure, or how the variances are estimated, and not what the fixed
#> effects mean; estimand changes what they estimate or what their standard
#> errors account for, so apply one of those only by choice.
```
