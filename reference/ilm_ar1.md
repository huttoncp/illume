# A first-order autoregressive structure over evenly spaced time

Observations of the same unit at adjacent time steps share a correlated
latent value. Use this when the time index moves in equal steps: study
visits, months, waves. For irregular times use
[`ilm_car1()`](https://huttoncp.github.io/illume/reference/ilm_car1.md).

## Usage

``` r
ilm_ar1(time, group, verbose = TRUE)
```

## Arguments

- time:

  Time index, one value per observation. Whole numbers, or values on a
  common step such as `c(0, 5, 10, 15)`. Or a one-sided formula
  `~ time | group` naming two columns of the model's data; see Details.

- group:

  Unit identifier, one value per observation. Omitted when `time` is a
  formula.

- verbose:

  Warn when the latent budget is thin.

## Value

A `"ilm_ar1"` specification, to pass as `ilm_model(ar = )`.

## Details

The latent grid spans every time step between the earliest and the
latest observation, so a unit that misses a wave is handled correctly –
the latent for that step simply has no observation attached, and the
chain continues across it. What is *not* allowed is times that do not
sit on a common grid at all, because then "one step" has no meaning.

**By name.** `ilm_ar1(~ time | group)` names two columns of the model's
data instead of passing them.
[`ilm_model()`](https://huttoncp.github.io/illume/reference/ilm_model.md)
reads them from its model frame, after rows with missing values are
dropped, so the structure cannot come out of step with the response;
vectors taken from the whole data frame beforehand can. The fitted model
remembers the names, so
[`ilm_cells()`](https://huttoncp.github.io/illume/reference/ilm_cells.md)
and predictions for new rows know where to find the time and the group.

## See also

[`ilm_car1()`](https://huttoncp.github.io/illume/reference/ilm_car1.md)
for arbitrary gaps,
[`ilm_rw1()`](https://huttoncp.github.io/illume/reference/ilm_rw1.md)
for a level that drifts,
[`ilm_check_ar()`](https://huttoncp.github.io/illume/reference/ilm_check_ar.md)
to test whether the structure is needed,
[`ilm_cells()`](https://huttoncp.github.io/illume/reference/ilm_cells.md)
for the fitted cells.

## Examples

``` r
d <- ilm_sim(n_id = 20, n_period = 8)
a <- ilm_ar1(as.integer(factor(d$date)), d$id)
#> Warning: AR(1) puts 160 latent values under 160 observations (1 per latent). Below about 1.5 the Laplace approximation frequently fails to converge, and a latent value seen once carries no information the residual does not. Round `time` to a coarser grid so observations share a latent.
a
#> AR(1), evenly spaced time 
#>   20 groups, 160 latent values, 160 observations (1.00 per latent)
#>   step: 1, grid of 8 time slots
#>   the latent budget is thin; see ?ilm_car1
```
