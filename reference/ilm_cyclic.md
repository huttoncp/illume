# Cyclic cubic spline basis for a cycle of known period

A smooth periodic basis that joins up at the wrap-around: the value and
the first two derivatives match at the end of one cycle and the start of
the next, so December runs into January without a step.

## Usage

``` r
ilm_cyclic(x, period, df = 5L)
```

## Arguments

- x:

  Time, as a number, integer or `Date`. Values are reduced modulo
  `period`, so the raw time works as well as the phase.

- period:

  Length of one full cycle, in the units of `x`.

- df:

  Number of columns, and so the degrees of freedom spent. Must be at
  least 3.

## Value

A numeric matrix with `df` columns.

## Details

The alternative to
[`ilm_fourier()`](https://craig-hutton.github.io/illume/reference/ilm_fourier.md),
and the one whose flexibility is *local*: changing the curve near one
knot leaves the rest of the cycle alone, where adding a Fourier harmonic
changes the shape everywhere.

Do not expect a large gain from that. Compared against
[`ilm_fourier()`](https://craig-hutton.github.io/illume/reference/ilm_fourier.md)
at matched degrees of freedom on a twelve-phase cycle, the two were
within a few AIC of each other for every smooth shape tried, and where
the shape really jumped a factor beat both by a hundred or more. On a
long cycle – day of year with one narrow summer feature – Fourier led
below about six degrees of freedom and this led above it, by at most 23
AIC.

Reach for it when the cycle is long, when you want more resolution in
one part of it without disturbing the rest, or when the alternative is a
factor with more levels than the data can pay for. Reach for
[`ilm_fourier()`](https://craig-hutton.github.io/illume/reference/ilm_fourier.md)
on a short cycle or a plain rise and fall, and for a factor when the
shape has a genuine step in it.

`df` is the number of columns returned, and so the degrees of freedom
the term spends – the same convention as
[`splines::ns()`](https://rdrr.io/r/splines/ns.html). One basis function
is dropped, because a cyclic basis sums to one at every point and would
otherwise be collinear with the intercept; dropping it leaves the fitted
curve and the predictions unchanged, since the remaining columns and the
intercept span exactly the same space.

Do not use an ordinary
[`splines::bs()`](https://rdrr.io/r/splines/bs.html) or
[`splines::ns()`](https://rdrr.io/r/splines/ns.html) on a phase variable
for this. Neither knows the two ends of the cycle are the same point, so
both leave a discontinuity at the wrap-around and waste degrees of
freedom estimating a jump that is not there.

## See also

[`ilm_fourier()`](https://craig-hutton.github.io/illume/reference/ilm_fourier.md)
for the cheaper smooth alternative,
[`ilm_check_ar()`](https://craig-hutton.github.io/illume/reference/ilm_check_ar.md),
which names the period when the residuals contain a cycle.

## Examples

``` r
b <- ilm_cyclic(1:12, period = 12, df = 4)
dim(b)
#> [1] 12  4
# the basis is periodic: one cycle on is the same point
all.equal(ilm_cyclic(1:6, 12, 4), ilm_cyclic(13:18, 12, 4),
          check.attributes = FALSE)
#> [1] TRUE
```
