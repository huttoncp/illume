# Fourier terms for a cycle of known period

Sine and cosine columns for use directly in a model formula, when the
response has a repeating cycle – a month-of-year effect in monthly data,
a day-of-week effect in daily data, a tide, a shift pattern.

## Usage

``` r
ilm_fourier(x, period, K = 1L)
```

## Arguments

- x:

  Time, as a number, integer or `Date`. One unit is one step, so the
  period is expressed in the same units.

- period:

  Length of one full cycle: 12 for months of a year, 7 for days of a
  week, 365.25 for days of a year.

- K:

  Number of harmonics.

## Value

A numeric matrix with one row per element of `x` and `2 * K` columns
(one fewer when `2 * K` equals `period`, where the last sine column is
identically zero).

## Details

`K` sets how much shape the cycle is allowed to have. `K = 1` is a
single smooth peak and trough per cycle; higher `K` adds harmonics and
can bend into sharper or double-peaked shapes, at two more columns each.
A cycle of length `period` supports at most `floor(period / 2)`
harmonics, beyond which the columns are not distinguishable at the
observed times.

Unlike [`stats::poly()`](https://rdrr.io/r/stats/poly.html) this depends
on nothing but its arguments, so a model containing it predicts
correctly on new data, including times outside the range that was
fitted.

## Choosing a basis

A factor for the phase – `factor(month)` – fits any shape at all, but
spends `period - 1` degrees of freedom doing it, and gives no reason for
neighbouring phases to resemble each other. Fourier terms spend `2K` and
are smooth by construction. For a 12-month cycle, `K = 2` (four columns)
covers most real seasonal shapes.

Go to a factor when the cycle has a genuine step in it, such as a policy
that starts in April: measured on a narrow two-month spike, a factor
beat both smooth bases by more than 140 AIC at matched sample size,
because no smooth basis can make a step cheaply. See
[`ilm_cyclic()`](https://craig-hutton.github.io/illume/reference/ilm_cyclic.md)
for the local alternative.

## See also

[`ilm_check_ar()`](https://craig-hutton.github.io/illume/reference/ilm_check_ar.md),
which names the period when the residuals contain a cycle.

## Examples

``` r
head(ilm_fourier(1:24, period = 12))
#>              sin1          cos1
#> [1,] 5.000000e-01  8.660254e-01
#> [2,] 8.660254e-01  5.000000e-01
#> [3,] 1.000000e+00  6.123234e-17
#> [4,] 8.660254e-01 -5.000000e-01
#> [5,] 5.000000e-01 -8.660254e-01
#> [6,] 1.224647e-16 -1.000000e+00

set.seed(1)
d <- data.frame(id = factor(rep(1:20, each = 24)), t = rep(1:24, 20))
d$y <- 2 * sin(2 * pi * d$t / 12) + rnorm(480)
f <- ilm_model(y ~ ilm_fourier(t, 12) + (1 | id), data = d,
               family = "gaussian", verbose = FALSE)
coef(f)
#>            (Intercept) ilm_fourier(t, 12)sin1 ilm_fourier(t, 12)cos1 
#>             0.02037751             1.92025693             0.01996844 
```
