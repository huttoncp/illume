# Power for a term, by simulating from a fitted model

Treats a fitted model as the truth, generates datasets of the requested
sizes, refits each, and counts how often the named term is detected. It
works for any family
[`ilm_model()`](https://huttoncp.github.io/illume/reference/ilm_model.md)
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
  [`ilm_model()`](https://huttoncp.github.io/illume/reference/ilm_model.md),
  or an
  [`ilm_scaffold()`](https://huttoncp.github.io/illume/reference/ilm_scaffold.md),
  to treat as the truth.

- n:

  Sample sizes to try. Defaults to a spread around the fitted size. For
  a model with a grouping factor this is the number of ROWS, and the
  number of clusters moves with it.

- term:

  What to test: one coefficient, named as the coefficient
  (`"armtreatment"`, or `"b:x"` in a multinomial model), or a whole
  term, named as the variable or term that produced it (`"arm"`, `"x"`).
  A term that produced a single coefficient is that coefficient; one
  that produced several is tested jointly. Defaults to the first
  coefficient that is not an intercept, and in a multinomial model to
  the first term.

- effect:

  For a single coefficient, values for it on the LINK scale; defaults to
  the fitted value. For a term tested jointly, MULTIPLES of its fitted
  coefficients – `c(0.5, 1)` asks what happens if the effect is half
  what was assumed; defaults to 1. Several values trace power across
  effect sizes.

- sims:

  Replicates per cell.

- alpha:

  Two-sided level.

- seed:

  Random seed.

- progress:

  Show a progress bar; see
  [illumex::ilm_progress_arg](https://huttoncp.github.io/illumex/reference/ilm_progress_arg.html).

## Value

An object of class `"ilm_power"`: one row per `n` by `effect` cell with
`power`, its Monte Carlo interval, `power_converged` and `converged`.

## The studies that are simulated

From a model fitted to data, each replicate resamples the fitted rows –
whole clusters at a time when the model has a grouping factor, so the
number of groups moves with `n` – and draws a new response from the
model: fresh random effects for the new groups, fresh serial
correlation, the fitted smooths, the zero part, the dispersion model and
the censoring. It is then refitted as the same model, contrasts and all.

From an
[`ilm_scaffold()`](https://huttoncp.github.io/illume/reference/ilm_scaffold.md),
which has no data behind it, each replicate is a fresh draw of the
planned design instead: the allocation balanced as the protocol would
balance it, and any covariate given as a function drawn again, because a
new study recruits new people.

## The test that is counted

The one the analysis will report. A single coefficient gets
[`summary()`](https://rdrr.io/r/base/summary.html)'s test: t on the
residual degrees of freedom when nothing is integrated out, the Wald z
otherwise. A term with several coefficients – a factor with three
levels, or any term of a multinomial model, which has one coefficient
per category – is tested jointly, F or chi-square, as
[`ilm_anova()`](https://huttoncp.github.io/illume/reference/ilm_anova.md)
tests it. Name one coefficient (`"b:x"`) to follow a single category
instead.

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

A check the analysis makes BEFORE fitting – too few observations per
random effect, say – is a verdict on the design, the same for every
study drawn from it, and is not counted against each replicate. The
result reports it instead, with its remedy, because every analysis of
such a study will print it as a failure and the plan should hear it
first.

## References

Arnold, B. F., Hogan, D. R., Colford, J. M. and Hubbard, A. E. (2011).
Simulation methods to estimate design power. *BMC Medical Research
Methodology* 11, 94.

## See also

[`ilm_power_n()`](https://huttoncp.github.io/illume/reference/ilm_power_n.md)
to read off the size for a target power,
[`plot.ilm_power()`](https://huttoncp.github.io/illume/reference/plot.ilm_power.md)
for the curve,
[`ilm_power_design()`](https://huttoncp.github.io/illume/reference/ilm_power_design.md)
for a study with no data yet,
[`ilm_simulate()`](https://huttoncp.github.io/illume/reference/ilm_simulate.md)
for the generator.

## Examples

``` r
set.seed(1); n <- 200
d <- data.frame(x = rnorm(n))
d$y <- rbinom(n, 1, plogis(-0.5 + 0.5 * d$x))
f <- ilm_model(y ~ x, data = d, family = "binomial", verbose = FALSE)
ilm_power(f, n = c(200, 400), sims = 50)
#> Simulated power for x (binomial family, alpha = 0.05)
#>   counting the Wald z test, as the analysis reports it
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
