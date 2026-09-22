# Difference in differences

Estimates the effect of a treatment that switches on for some units at
some time, by comparing the change in the treated group with the change
in the control group. The estimate is the interaction in a mixed model,
so the fit is an ordinary
[`ilm_model()`](https://craig-hutton.github.io/illume/reference/ilm_model.md)
and every diagnostic and method applies to it.

## Usage

``` r
ilm_did(
  data,
  y,
  unit,
  time,
  treated = NULL,
  post = NULL,
  treat_time = NULL,
  treatment = NULL,
  covariates = NULL,
  family = NULL,
  ar = FALSE,
  allow_staggered = FALSE,
  verbose = TRUE,
  ...
)
```

## Arguments

- data:

  A data frame, one row per unit-period.

- y:

  Outcome column.

- unit:

  Unit identifier.

- time:

  Period. Numeric or anything that orders.

- treated:

  Column marking units ever treated. Omit when `treatment` is given.

- post:

  Column marking periods after treatment starts. Omit when `treat_time`
  or `treatment` is given.

- treat_time:

  First treated period, when `post` is not supplied.

- treatment:

  Per-row treatment indicator, as an alternative to `treated` and
  `post`. Staggered timing is detected from it.

- covariates:

  Further columns for the mean structure. These must be things the
  treatment cannot have affected; a covariate the treatment changes is a
  mediator and adjusting for it removes part of the effect.

- family:

  Response distribution; inferred from `y` when `NULL`.

- ar:

  Add AR(1) within unit on top of the random intercept.

- allow_staggered:

  Report a pooled estimate even when adoption is staggered. Off by
  default, and see the section above.

- verbose:

  Narrate each step.

- ...:

  Passed to
  [`ilm_model()`](https://craig-hutton.github.io/illume/reference/ilm_model.md).

## Value

An object of class `"ilm_did"`: `att`, `fit`, `parallel`, `event` and
the settings used.

## What identifies the estimate

Difference in differences rests on **parallel trends**: without the
treatment, the two groups' outcomes would have moved together. That
cannot be checked where it matters, because it concerns a world in which
the treatment did not happen, but it can be checked before the
treatment, and that is worth doing. Two things here do it:

- a formal test of whether the groups' pre-treatment slopes differ,
  reported as `parallel`;

- the event study, `$event`, one coefficient per period relative to
  treatment. Pre-treatment coefficients should sit near zero. This is
  the more informative of the two, because it shows *how* a violation
  looks, and
  [`ilm_plot_did()`](https://craig-hutton.github.io/illume/reference/ilm_plot_did.md)
  draws it.

## Serial correlation

Outcomes from the same unit over time are correlated, and difference in
differences with naive standard errors over-rejects badly as a result
(Bertrand, Duflo and Mullainathan 2004). A unit random intercept is
fitted by default. With enough periods, `ar = TRUE` adds an AR(1)
process on top, which is the part that matters when the within-unit
correlation decays rather than being constant.

## Staggered adoption

When units are treated at different times, the two-way fixed effects
regression is **not** an average treatment effect if the effect varies
across units or over time: already-treated units end up serving as
controls for later-treated ones, and those comparisons can enter with
negative weight. This function detects staggered timing and refuses to
report a single pooled number for it unless `allow_staggered = TRUE`, in
which case the result carries the caveat. The honest remedy is a
cohort-based estimator – the `did` and `fixest` packages implement them
– and illume does not yet have one.

## What this is calibrated for

On 800 replicates of a 40-unit, 8-period panel with unit random effects
and a true effect of 0.800: the estimate was unbiased (+0.0029), its
interval covered 0.946 of the time, and the parallel-trends check raised
a false alarm 0.054 of the time under trends that really were parallel,
with null p-values uniform by Kolmogorov-Smirnov (p = 0.111) rather than
merely correct at the 5% threshold.

Getting that check calibrated took choosing the right reference. The
comparison is between two groups of UNITS' pre-treatment slopes, so
units are the independent replicates and the large-sample normal is too
generous when there are few: it rejected 0.060, 0.068 and 0.050 of the
time at 40, 20 and 80 units, where t on `units - 2` gave 0.048, 0.055
and 0.045. A check that cries wolf is worse than no check here, because
the remedy for a failed parallel-trends test is to abandon the design.

## References

Bertrand, M., Duflo, E. and Mullainathan, S. (2004). How much should we
trust differences-in-differences estimates? Quarterly Journal of
Economics 119(1).

Goodman-Bacon, A. (2021). Difference-in-differences with variation in
treatment timing. Journal of Econometrics 225(2).

## See also

[`ilm_plot_did()`](https://craig-hutton.github.io/illume/reference/ilm_plot_did.md),
[`ilm_rdd()`](https://craig-hutton.github.io/illume/reference/ilm_rdd.md),
[`ilm_model()`](https://craig-hutton.github.io/illume/reference/ilm_model.md).

## Examples

``` r
set.seed(1)
d <- expand.grid(unit = 1:40, time = 1:8)
d$treated <- as.integer(d$unit <= 20)
d$post <- as.integer(d$time >= 5)
d$y <- 1 + 0.3 * d$time + rnorm(40)[d$unit] +
       0.8 * d$treated * d$post + rnorm(nrow(d))
ilm_did(d, "y", "unit", "time", treated = "treated", post = "post",
        verbose = FALSE)
#> <ilm_did> y ~ treatment  | 20 treated, 20 control units
#>   family: gaussian   random intercept: unit  
#> 
#>   ATT    0.7254  ( 0.3000,  1.1507)  p = 0.000831
#> 
#>   parallel trends before treatment: OK  (slope difference 0.2003, p = 0.15)
#>   event study: 8 periods, 0 pre-treatment coefficient(s) excluding zero  (ilm_plot_did)
```
