# Power for a term, by simulating from a fitted model

Treats a fitted model as the truth, generates datasets of the requested
sizes, refits each, and counts how often the named term is detected. It
works for any family
[`ilm_model()`](https://craig-hutton.github.io/illume/reference/ilm_model.md)
fits, because it never needs a closed-form variance.

## Usage

``` r
ilm_power(
  object,
  n = NULL,
  term = NULL,
  effect = NULL,
  sims = 200L,
  alpha = 0.05,
  seed = 1L,
  progress = NULL
)
```

## Arguments

- object:

  A fitted
  [`ilm_model()`](https://craig-hutton.github.io/illume/reference/ilm_model.md)
  to treat as the truth.

- n:

  Sample sizes to try. Defaults to a spread around the fitted size. For
  a model with a grouping factor this is the number of ROWS, and the
  number of clusters moves with it.

- term:

  The fixed-effect coefficient to test, named either as the coefficient
  (`"armtreatment"`) or as the variable that produced it (`"arm"`). The
  variable form is resolved when it maps to a single coefficient, which
  covers continuous predictors and two-level factors; a term with more
  columns is refused, with its coefficients listed, because a power
  curve follows one effect size at a time. Defaults to the first
  coefficient that is not an intercept.

- effect:

  Values for that coefficient on the LINK scale. Defaults to the fitted
  value. Supplying several traces power across effect sizes.

- sims:

  Replicates per cell.

- alpha:

  Two-sided level.

- seed:

  Random seed.

- progress:

  Show a progress bar; see
  [ilm_progress_arg](https://craig-hutton.github.io/illume/reference/ilm_progress_arg.md).

## Value

An object of class `"ilm_power"`: one row per `n` by `effect` cell with
`power`, its Monte Carlo interval, `power_converged` and `converged`.

## The estimate has a standard error

Power from `sims` replicates is a proportion, so it carries a standard
error of `sqrt(p (1 - p) / sims)` – 2.8% at 0.80 from 200 replicates.
The returned interval is a Wilson interval on that, and the sample size
implied by a target power is a range rather than a number. Raising
`sims` narrows it, and nothing else does.

## Fits that do not converge

A replicate that fails to fit has not detected anything, so `power`
counts it as a non-detection. `power_converged` divides by the
replicates that worked, which is the number most software reports and is
**power conditional on convergence** – a different and more flattering
quantity. When `converged` is below one the two differ, and the gap is
the size of the problem rather than something to smooth over.

## References

Arnold, B. F., Hogan, D. R., Colford, J. M. and Hubbard, A. E. (2011).
Simulation methods to estimate design power. *BMC Medical Research
Methodology* 11, 94.

## See also

[`ilm_power_n()`](https://craig-hutton.github.io/illume/reference/ilm_power_n.md)
to read off the size for a target power,
[`plot.ilm_power()`](https://craig-hutton.github.io/illume/reference/plot.ilm_power.md)
for the curve,
[`ilm_simulate()`](https://craig-hutton.github.io/illume/reference/ilm_simulate.md)
for the generator.

## Examples

``` r
set.seed(1); n <- 200
d <- data.frame(x = rnorm(n))
d$y <- rbinom(n, 1, plogis(-0.5 + 0.5 * d$x))
f <- ilm_model(y ~ x, data = d, family = "binomial", verbose = FALSE)
ilm_power(f, n = c(200, 400), sims = 50)
#> Simulated power for x (binomial family, alpha = 0.05)
#>   50 replicates per cell; the fitted study had 200 rows
#>    n effect power mc_lower mc_upper power_converged converged
#>  200 0.3856  0.68    0.542    0.792            0.68         1
#>  400 0.3856  0.96    0.865    0.989            0.96         1
#> 
#>   mc_lower/mc_upper is a Wilson interval on the power ESTIMATE: at
#>   0.80 from 50 replicates the standard error is 0.057, so a
#>   sample size read off this curve is a range. More replicates narrow it;
#>   nothing else does.
```
