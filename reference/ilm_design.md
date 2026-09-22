# Describe how a sample was drawn

Bundles the weights, clustering and stratification of a complex sample
so that
[`ilm_model()`](https://craig-hutton.github.io/illume/reference/ilm_model.md)
can fit with them and report a variance that reflects them. The
arguments mirror
[`survey::svydesign()`](https://rdrr.io/pkg/survey/man/svydesign.html).

## Usage

``` r
ilm_design(data, weights, ids = NULL, strata = NULL, fpc = NULL)
```

## Arguments

- data:

  The data frame the sample lives in.

- weights:

  A one-sided formula naming the sampling weight, or a numeric vector.
  These are inverse probabilities of selection, not replicate counts.

- ids:

  A one-sided formula naming the primary sampling unit, or `NULL` for an
  unclustered sample, in which case every row is its own PSU.

- strata:

  A one-sided formula naming the stratum, or `NULL`.

- fpc:

  A one-sided formula naming the finite population correction – either
  the population size in each stratum or the sampling fraction – or
  `NULL` to ignore it, which is conservative.

## Value

An object of class `"ilm_design"`.

## References

Lumley, T. (2010). *Complex Surveys: A Guide to Analysis Using R*.
Wiley.

## See also

[`ilm_model()`](https://craig-hutton.github.io/illume/reference/ilm_model.md)
with `design =`, and
[`ilm_svy_coef()`](https://craig-hutton.github.io/illume/reference/ilm_svy_coef.md).

## Examples

``` r
set.seed(1); n <- 500
d <- data.frame(x = rnorm(n), psu = rep(1:50, each = 10),
                st = rep(1:2, each = 250))
d$w <- ifelse(d$st == 1, 40, 10)
d$y <- 1 + 0.5 * d$x + rnorm(n)
des <- ilm_design(d, weights = ~ w, ids = ~ psu, strata = ~ st)
des
#> <ilm_design>
#>   500 rows in 50 primary sampling unit(s), 2 stratum(a)
#>   weights 10 to 40, summing to 12500
#>   design degrees of freedom: 48
```
