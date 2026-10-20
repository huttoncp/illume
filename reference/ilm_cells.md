# The cells of a correlation over time

One row per latent value of a fitted model's
[`ilm_ar1()`](https://huttoncp.github.io/illume/reference/ilm_ar1.md),
[`ilm_car1()`](https://huttoncp.github.io/illume/reference/ilm_car1.md)
or [`ilm_rw1()`](https://huttoncp.github.io/illume/reference/ilm_rw1.md)
term: which group and time it belongs to, how many observations sit on
it, and where the fit holds it. Code that works with the latent values –
forecasts, draws, checks – reads the layout here rather than rebuilding
the index arithmetic, which is how the two come apart.

## Usage

``` r
ilm_cells(object)
```

## Arguments

- object:

  A fitted `"ilm_model"`.

## Value

`NULL` when the model has no correlation over time. Otherwise a data
frame with one row per cell, ordered by group and then time, and columns

- `term`:

  `"ar"`, the name the term has in `Sigma` and in the parameter names.

- `type`:

  `"ar1"`, `"car1"` or `"rw1"`.

- `group`:

  the unit, a factor with the fitted levels.

- `time`:

  the cell's time, on the scale it was given: numeric, `Date` or
  date-time.

- `index`:

  where the fit holds the cell's value: the `index`-th of the entries
  named `"B_ar"` in `object$sdr$par.random` and in the joint precision.
  `NA` at a random walk's anchors, which are held at zero rather than
  estimated. With a multinomial outcome a cell has one value per linear
  predictor, stored predictor after predictor, so predictor `c`'s value
  is at `index + (c - 1) * attr(, "n_latent")`.

- `n_obs`:

  how many fitted rows sit on the cell; 0 at an AR(1) grid step no one
  was observed at.

- `anchor`:

  `TRUE` at each group's first cell of a random walk.

- `last`:

  `TRUE` at each group's final cell, where a forecast starts. An AR(1)
  grid runs to the latest time in the data for every group, so a group
  seen only early ends in cells without observations, whose values are
  the chain's own forecast.

Attributes: `obs_cell`, each fitted row's cell (a row number of this
table) in the order of the model frame; `vars`, the time and group
columns when the term was given by name, `~ time | group`, and `NULL`
when it was given vectors; `n_latent`, the number of estimated cells;
`parameterisation`, what the fitted parameters mean, in words; and for
AR(1), `step` and `origin`, the grid's spacing and first time.

## See also

[`ilm_ar1()`](https://huttoncp.github.io/illume/reference/ilm_ar1.md),
[`ilm_car1()`](https://huttoncp.github.io/illume/reference/ilm_car1.md),
[`ilm_rw1()`](https://huttoncp.github.io/illume/reference/ilm_rw1.md).

## Examples

``` r
set.seed(2)
d <- data.frame(id = rep(c("a", "b"), each = 12), day = rep(1:12, 2))
d$y <- rnorm(24) + rep(cumsum(rnorm(12, 0, 0.5)), 2)
fit <- ilm_model(y ~ 1, data = d, ar = ilm_rw1(~ day | id), verbose = FALSE)
#> ilm_model(): family = "gaussian", inferred from `y`: continuous values from -3.59 to 2.42. Pass `family` to choose another.
cl <- ilm_cells(fit)
head(cl)
#>   term type group time index n_obs anchor  last
#> 1   ar  rw1     a    1    NA     1   TRUE FALSE
#> 2   ar  rw1     a    2     1     1  FALSE FALSE
#> 3   ar  rw1     a    3     2     1  FALSE FALSE
#> 4   ar  rw1     a    4     3     1  FALSE FALSE
#> 5   ar  rw1     a    5     4     1  FALSE FALSE
#> 6   ar  rw1     a    6     5     1  FALSE FALSE
attr(cl, "parameterisation")
#> [1] "Random walk: each group's walk is zero at its first cell; a step across d time units adds an independent change with covariance d * Sigma$ar, so Sigma$ar is the variance per unit of time and the variance grows with the time since the group's first cell."
```
