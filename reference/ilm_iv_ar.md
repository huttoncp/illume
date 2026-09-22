# Anderson-Rubin confidence set for an instrumented coefficient

A confidence set that is valid however weak the instrument is, which the
usual estimate-plus-or-minus-two-standard-errors is not. It works by
testing each candidate value directly: if the coefficient really is `b`,
then `y - x b` has the instrument's influence removed from it, so the
instruments should predict nothing of what is left. The set is every `b`
that survives.

## Usage

``` r
ilm_iv_ar(object, level = 0.95, range = NULL, n_grid = 2001L)
```

## Arguments

- object:

  An [`ilm_iv()`](https://huttoncp.github.io/illume/reference/ilm_iv.md)
  fit with exactly one endogenous regressor.

- level:

  Confidence level.

- range:

  Optional `c(lo, hi)` to search over; by default a wide band around the
  2SLS estimate.

- n_grid:

  Number of grid points.

## Value

A list with `lower`, `upper`, whether the set is `bounded`, `empty`, and
the grid it was built from.

## It can be unbounded, and that is information

When the instrument is weak the set can run to infinity in one or both
directions, or be empty. An unbounded set is the honest statement that
the data cannot rule out arbitrarily large effects; a Wald interval in
the same situation reports a tidy finite range and is simply wrong. An
empty set means no value of the coefficient reconciles the instruments
with each other, which is evidence against the model rather than against
any particular value.

## References

Anderson, T. W. and Rubin, H. (1949). Estimation of the parameters of a
single equation in a complete system of stochastic equations. *Annals of
Mathematical Statistics* 20, 46-63.

## See also

[`ilm_iv()`](https://huttoncp.github.io/illume/reference/ilm_iv.md).

## Examples

``` r
set.seed(2); n <- 400
z <- rnorm(n); u <- rnorm(n)
x <- 0.15 * z + u + rnorm(n)          # a weak instrument
y <- 1 + 0.5 * x + u + rnorm(n)
fit <- ilm_iv(y ~ x | z, data = data.frame(y = y, x = x, z = z))
ilm_iv_ar(fit)
#> Anderson-Rubin 95% confidence set for x
#>   [-Inf, Inf]
#>   UNBOUNDED. The instrument is too weak to rule out arbitrarily
#>   large effects. A Wald interval would report a tidy finite range
#>   here and would be wrong.
```
