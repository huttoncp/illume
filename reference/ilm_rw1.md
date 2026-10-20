# A random walk over time

Observations of the same unit share a latent level that wanders: each
step adds an independent change whose variance is proportional to the
time that passed, so the level has no mean to return to. This is the
local-level model of structural time series. Use it when a unit's level
drifts rather than fluctuating about a fixed value – where
[`ilm_car1()`](https://huttoncp.github.io/illume/reference/ilm_car1.md)
would put its correlation at the edge of its range.

## Usage

``` r
ilm_rw1(time, group, verbose = TRUE)
```

## Arguments

- time:

  Time, one value per observation: any numeric scale, a `Date` or a
  date-time. The variance is per unit of it. Or a one-sided formula
  `~ time | group` naming two columns of the model's data.

- group:

  Unit identifier, one value per observation. Omitted when `time` is a
  formula.

- verbose:

  Warn about nearly coincident times.

## Value

A `"ilm_rw1"` specification, to pass as `ilm_model(ar = )`.

## Details

**The walk starts at zero.** Each group's walk is held at zero at its
first time, so its level there is the one the fixed effects give, and
every later value is measured from it. Without that anchor the walk and
the intercept would describe the same thing and neither could be
estimated. It also means the groups are taken to start from a common
level, give or take the fixed effects; where they do not, add a random
intercept, `(1 | group)`, and each starts where it does.

**The variance is per unit of time.** Between two times `d` units apart
the walk moves by a normal change with variance `d` times the fitted
variance, which the fitted model holds as `Sigma$ar`. So the units of
`time` matter: days and weeks give variances seven times apart for the
same walk. Irregular gaps enter as they are, and cost nothing in
accuracy.

**What is estimated.** One latent value per distinct (group, time) pair
after each group's first. With a gaussian response the Laplace
approximation is exact, so the fit is the exact likelihood of the
local-level model, by maximum likelihood or by REML (`reml = TRUE`).
With other families it is the Laplace approximation, and the fit-time
checks judge the number of observations per latent value as they do for
[`ilm_car1()`](https://huttoncp.github.io/illume/reference/ilm_car1.md).

`ilm_rw1(~ time | group)` names two columns of the model's data instead
of passing them; see
[`ilm_ar1()`](https://huttoncp.github.io/illume/reference/ilm_ar1.md)
for why that is the safer form. It is also what lets
`predict(groups = "population")` place new rows on the walk, whose
variance depends on how long after its group's first time each one
falls.

## References

Durbin, J., & Koopman, S. J. (2012). *Time Series Analysis by State
Space Methods* (2nd ed.). Oxford University Press. (Chapter 2, the local
level model.)

## See also

[`ilm_car1()`](https://huttoncp.github.io/illume/reference/ilm_car1.md)
for a process that returns to its mean,
[`ilm_cells()`](https://huttoncp.github.io/illume/reference/ilm_cells.md)
for the fitted cells.

## Examples

``` r
set.seed(1)
d <- data.frame(t = 1:60, unit = "a")
d$y <- 2 + cumsum(c(0, rnorm(59, 0, 0.7))) + rnorm(60)
fit <- ilm_model(y ~ 1, data = d, ar = ilm_rw1(~ t | unit),
                 verbose = FALSE)
#> ilm_model(): family = "gaussian", inferred from `y`: continuous values from 0.414 to 7.2. Pass `family` to choose another.
fit$Sigma$ar                 # variance of the walk per unit of time
#>           [,1]
#> [1,] 0.2615382
head(ilm_cells(fit))
#>   term type group time index n_obs anchor  last
#> 1   ar  rw1     a    1    NA     1   TRUE FALSE
#> 2   ar  rw1     a    2     1     1  FALSE FALSE
#> 3   ar  rw1     a    3     2     1  FALSE FALSE
#> 4   ar  rw1     a    4     3     1  FALSE FALSE
#> 5   ar  rw1     a    5     4     1  FALSE FALSE
#> 6   ar  rw1     a    6     5     1  FALSE FALSE
```
