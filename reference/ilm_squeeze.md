# Move proportions off 0 and 1

The beta family is defined on the open interval, so a proportion
recorded as exactly 0 or 1 has no likelihood. When those values are a
rounding of interior ones rather than a separate process, the usual
remedy is the transformation of Smithson and Verkuilen (2006):

## Usage

``` r
ilm_squeeze(y, n = length(y), quiet = FALSE)
```

## Arguments

- y:

  A numeric vector of proportions in `[0, 1]`.

- n:

  Sample size to use in the transformation; the number of observations
  by default.

- quiet:

  Suppress the message reporting the shift.

## Value

A numeric vector strictly inside `(0, 1)`.

## Details

\$\$y' = \frac{y (n - 1) + 1/2}{n}\$\$

which compresses the whole scale towards one half by a factor of
`(n - 1) / n`. Every value moves, not only the offending ones: nudging
the endpoints alone would leave a gap just inside them and change the
shape of the distribution more than this does.

It is a choice, not a repair. The function reports how far it moved
things so the choice is visible, and with `n` small the shift is not
negligible – at `n = 50` an observed 0 becomes 0.01.

## References

Smithson, M. and Verkuilen, J. (2006). A better lemon squeezer?
Maximum-likelihood regression with beta-distributed dependent variables.
*Psychological Methods* 11, 54-71.

## See also

[`ilm_model()`](https://craig-hutton.github.io/illume/reference/ilm_model.md)
with `family = "beta"`, and `ziformula` for the case where the boundary
values are a separate process rather than rounding.

## Examples

``` r
y <- c(0, 0.2, 0.5, 0.9, 1)
ilm_squeeze(y)
#> ilm_squeeze: 2 values at a boundary; every value moved towards 1/2 by up to 0.1. An observed 0 is now 0.1 and an observed 1 is 0.9.
#> [1] 0.10 0.26 0.50 0.82 0.90
```
