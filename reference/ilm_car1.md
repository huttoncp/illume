# A continuous-time autoregressive structure over irregular time

The correlation between two observations of the same unit is `rho`
raised to the power of the gap between them, so observations close in
time are similar and the spacing need not be regular. This is the
structure `nlme`'s `corCAR1()` fits, and the AR(1) special case when
every gap is one.

## Usage

``` r
ilm_car1(time, group, verbose = TRUE)
```

## Arguments

- time:

  Time, one value per observation. Any numeric scale or `Date`; `rho` is
  the correlation one unit apart, so the units matter. Or a one-sided
  formula `~ time | group` naming two columns of the model's data.

- group:

  Unit identifier, one value per observation. Omitted when `time` is a
  formula.

- verbose:

  Warn when the latent budget is thin.

## Value

A `"ilm_car1"` specification, to pass as `ilm_model(ar = )`.

## Details

Written as a Markov chain the transition from one observation to the
next is `phi_k = rho ^ d_k` with innovation variance `1 - phi_k^2`,
which is exact for the Ornstein-Uhlenbeck process rather than an
approximation to it. Irregular spacing therefore costs nothing in
accuracy.

What it does cost is latent values. One exists per distinct (group,
time) pair, so genuinely continuous times with no repeats give one
latent per observation – the regime where the Laplace approximation
stops converging. The constructor reports the budget and warns when it
is thin; rounding `time` to a coarser grid, so that observations share a
latent, is the usual fix and is what makes the difference between a
model that fits and one that does not.

`ilm_car1(~ time | group)` names two columns of the model's data instead
of passing them; see
[`ilm_ar1()`](https://huttoncp.github.io/illume/reference/ilm_ar1.md)
for why that is the safer form.

## See also

[`ilm_ar1()`](https://huttoncp.github.io/illume/reference/ilm_ar1.md)
for evenly spaced time,
[`ilm_rw1()`](https://huttoncp.github.io/illume/reference/ilm_rw1.md)
for a level that drifts rather than reverting,
[`ilm_check_ar()`](https://huttoncp.github.io/illume/reference/ilm_check_ar.md)
and
[`ilm_plot_acf()`](https://huttoncp.github.io/illume/reference/ilm_plot_acf.md)
to test whether the structure is needed,
[`ilm_cells()`](https://huttoncp.github.io/illume/reference/ilm_cells.md)
for the fitted cells.

## Examples

``` r
set.seed(1)
d <- data.frame(id = rep(1:20, each = 6),
                day = as.vector(replicate(20, sort(sample(1:60, 6)))))
ilm_car1(d$day, d$id)
#> Warning: CAR(1) puts 120 latent values under 120 observations (1 per latent). Below about 1.5 the Laplace approximation frequently fails to converge, and a latent value seen once carries no information the residual does not. Round `time` to a coarser grid so observations share a latent.
#> CAR(1), continuous time 
#>   20 groups, 120 latent values, 120 observations (1.00 per latent)
#>   gaps: min 1, median 7, max 24
#>   the latent budget is thin; see ?ilm_car1
```
