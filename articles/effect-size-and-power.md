# Effect size and power

Three questions that bracket an analysis. Before it: how large does the
study need to be. After it: how large is the effect, on a scale someone
can act on. And then: what would happen if we changed something.

``` r

library(illume)
```

## Before: how large a study

[`ilm_power()`](https://huttoncp.github.io/illume/reference/ilm_power.md)
treats a fitted model as the truth, simulates studies at the sizes you
are considering, refits each one, and counts how often the term is
detected. Nothing in it needs a closed-form variance, so it covers mixed
models, zero-inflated counts and ordinal outcomes as readily as a linear
regression.

``` r

pilot <- ilm_model(y ~ dose, data = pilot_data, family = "binomial")
pw <- ilm_power(pilot, n = c(200, 400, 800), effect = c(0.3, 0.5), sims = 500)
pw
plot(pw)
```

Two things it reports that most software does not, and both change what
the number means.

### The power estimate is itself an estimate

A power of 0.80 from 200 replicates is a proportion, so it carries a
standard error of `sqrt(0.8 * 0.2 / 200)`, about 2.8%. The power is 0.74
to 0.86, and a sample size read off that curve is a range rather than a
number.

``` r

ilm_power_n(pw, target = 0.8)
#>   effect     n  n_lower  n_upper
#>      0.3   240      207      272
```

Quoting “you need 240 participants” from that is false precision.
Quoting “somewhere between 210 and 270, and more replicates will narrow
it” is the answer the simulation supports. Raising `sims` narrows the
interval and nothing else does – collecting more pilot data changes the
estimate, not its Monte Carlo error.

The interval is Wilson rather than Wald, because power estimates live
near 1, where a Wald interval reaches above it.

### Some fits will not converge

A replicate that fails to fit has not detected anything. `power` counts
it as a non-detection; `power_converged` divides only by the replicates
that worked, which is **power conditional on convergence** – a different
and more flattering quantity, and the one most software reports.

``` r

#>   n effect power mc_lower mc_upper power_converged converged
#>  40  2.072 0.850    0.756    0.912           0.872     0.975
```

When `converged` is below one the two differ, and the gap is the size of
the problem rather than something to smooth over. A design where a fifth
of studies fail to fit is a design worth reconsidering, and averaging
over only the successes hides exactly that.

### Mixed designs

For a model with a grouping factor, `n` is the number of rows but the
number of **clusters** moves with it, because that is what a mixed
design’s power actually depends on. Resampling rows would hold the
number of groups fixed while n grew, which is the wrong thing to hold
fixed. A cluster drawn twice becomes two clusters, or the design has
fewer independent groups than it appears to.

A random **slope** is drawn as a random slope, not as an intercept
shift. A `(1 + time | id)` design whose participants were all given the
same slope would have them moving in near-parallel, which understates
how far they differ and so overstates power for anything interacting
with `time`.

### When there is no pilot study

[`ilm_power()`](https://huttoncp.github.io/illume/reference/ilm_power.md)
treats a fitted model as the truth, and a study being planned has
nothing to fit.
[`ilm_scaffold()`](https://huttoncp.github.io/illume/reference/ilm_scaffold.md)
builds the model from assumptions instead: a design grid, the parameters
you are prepared to defend, and an ordinary `ilm_model` object that
everything else in the package will accept.

Assumptions go in as cell means on the response scale, which is the form
most people actually hold one in.

``` r

d <- list(arm = c("control", "treatment"), time = c("pre", "post"))
cells <- c(control.pre = 20, control.post = 21,
           treatment.pre = 20, treatment.post = 25)

s <- ilm_scaffold(change ~ arm * time + (1 | id), design = d,
                  within = "time", n_unit = 60,
                  cells = cells, sd = 6, icc = 0.5)
#> <ilm_scaffold>  parameters ASSUMED, not estimated
#>   formula    : change ~ arm * time + (1 | id)
#>   family     : gaussian
#>   size       : 60 ids x 2 = 120 rows
#>   between    : arm
#>   within     : time
#>
#>   assumed coefficients
#>           (Intercept)          armtreatment              timepost
#>                    20                     0                     1
#> armtreatment:timepost
#>                     4
#>
#>   residual sd : 6
#>   id sd     : 6
```

`n_unit` counts participants, not rows, so the number you give is the
number that goes in a protocol. `arm` varies between participants and
`time` within, and the grid is built accordingly – allocating a
between-participants variable to rows would put the same person in both
arms.

The point of returning a model rather than a number is that the
assumptions can be interrogated before anything is built on them. The
cell means went in; here they are coming back out through an entirely
different route.

``` r

ilm_emmeans(s, c("arm", "time"))
#>        arm time estimate    se lower upper
#>    control  pre       20 1.595 16.87 23.13
#>  treatment  pre       20 1.595 16.87 23.13
#>    control post       21 1.595 17.87 24.13
#>  treatment post       25 1.595 21.87 28.13
```

Cell means the formula cannot produce are refused rather than
approximated. A least-squares solve always returns *something*, and the
something it returns for crossed means under an additive formula is a
scaffold for a different study:

``` r

ilm_scaffold(change ~ arm + time, design = d, within = "time",
             n_unit = 60, cells = cells, sd = 6, icc = 0.5)
#> Error: this formula cannot produce the cell means given: cell
#> 'control.pre' is off by 1 on the link scale. Cell means that differ by
#> more than the terms in the formula allow need the interaction between
#> them -- arm * time rather than arm + time.
```

[`ilm_power_design()`](https://huttoncp.github.io/illume/reference/ilm_power_design.md)
is the scaffold plus the power curve in one call.

``` r

pw <- ilm_power_design(change ~ arm * time + (1 | id), design = d,
                       within = "time", n_unit = c(60, 120, 240, 400),
                       cells = cells, sd = 6, icc = 0.5,
                       term = "arm:time", sims = 300)
#>  n_unit   n effect power mc_lower mc_upper power_converged converged
#>      60 120      4 0.470    0.414    0.526           0.470         1
#>     120 240      4 0.700    0.646    0.749           0.700         1
#>     240 480      4 0.967    0.940    0.982           0.967         1
#>     400 800      4 1.000    0.987    1.000           1.000         1

ilm_power_n(pw, target = 0.8)
#>  effect   n  n_lower n_upper n_unit n_unit_lower n_unit_upper
#>       4 330 292.5298 365.872    165     146.2649      182.936
```

Read `n_unit`, not `n`: 165 participants measured twice, not 330 people.

The standard errors printed on a scaffold come from **one realisation**
of the design at its own size. They are not a property of the
assumptions, and power read off them would be one coin flip. That is
what
[`ilm_power_design()`](https://huttoncp.github.io/illume/reference/ilm_power_design.md)
is for.

None of this makes an assumption true. It makes the assumption explicit,
and separates “we need 165 participants” from the four numbers that
claim implies.

## After: how large the effect

[`ilm_effects()`](https://huttoncp.github.io/illume/reference/ilm_effects.md)
reports each coefficient on the scale its family puts it on.

``` r

fit <- ilm_model(relapse ~ dose + age, data = d, family = "binomial")
ilm_effects(fit)
```

| family | scale |
|----|----|
| binomial, beta | odds ratio |
| poisson, negative binomial | incidence rate ratio |
| ordinal | proportional odds ratio |
| weibull, lognormal, loglogistic | time ratio (an acceleration factor, **not** a hazard ratio) |
| Royston-Parmar | hazard ratio |
| gaussian | the coefficient, plus partial variance explained |

Three details that are easy to get wrong.

**The interval is built on the link scale and transformed.**
[`exp()`](https://rdrr.io/r/base/Log.html) of the endpoints is exact and
asymmetric, which is the shape a ratio’s uncertainty really has. A
standard error computed on the ratio scale and doubled either way is
symmetric and can reach below zero.

**A ratio from a mixed model is conditional.** `exp(beta)` there is the
odds ratio for a *given cluster* – two patients in the same hospital,
one exposed – not for the population, and the gap widens as the random
effects grow. The print says so and names
[`ilm_ame()`](https://huttoncp.github.io/illume/reference/ilm_ame.md)
for the population-averaged quantity:

``` r

ilm_effects(fit)      # conditional odds ratio, e.g. 2.58
ilm_ame(fit, "dose")  # marginal effect on the probability scale, e.g. 0.193
```

**Standardising a binary predictor is meaningless.** Dividing by the
standard deviation of a 0/1 variable answers what happens per standard
deviation of treatment, which is not a question anyone asks. Those are
left alone, and the reason is printed rather than the number quietly
omitted. `standardise = "gelman"` divides numeric predictors by *two*
standard deviations, which is what makes them comparable with a binary
one.

``` r

ilm_effects(fit, standardise = "gelman")
```

A robust covariance can be substituted, so the same table can be built
on clustered or design-based standard errors:

``` r

ilm_effects(fit, vcov = ilm_vcov_cluster(fit, ~ clinic))
```

## Then: what would happen if

[`ilm_scenario()`](https://huttoncp.github.io/illume/reference/ilm_scenario.md)
sets predictors to values someone is considering and reports what the
model implies.

``` r

ilm_scenario(fit, dose = c(0, 10, 20), contrast = "first")
```

The default averages over **the observed units**: set the intervention
for everyone, keep each unit’s own covariates, predict, average. That is
standardisation, and when the adjustment set is valid it is the causal
estimand.

This is not the same as building one row with the predictors set and
predicting from it. That row is a *representative unit* – what happens
to someone who looks like this – and under any non-linear link it is a
different number, because the average of a prediction is not the
prediction at the average. On a logistic fit at the same dose: 0.657
averaged over the units, 0.683 at the average covariate.
`over = "reference"` gives the second, and the print says which was
used.

``` r

ilm_scenario(fit, dose = 10)                        # the population
ilm_scenario(fit, dose = 10, over = "reference")    # a unit like this
```

### It checks whether the scenario is inside the data

A model returns a number for a dose nobody received, and the interval
around it looks like any other, because the interval carries uncertainty
about the coefficients rather than about whether the functional form
survives out there. Each scenario is checked twice.

**By value**, against each predictor’s observed range – with a 2%
tolerance, so asking for `dose = 0` when the smallest observed is 0.013
is not flagged. Flagging that every time is how a warning gets ignored.

**By combination**, which no range check catches. Every value can sit
inside its own range while the point sits in a corner nobody occupies.
On personnel data where years of service is roughly age minus 22, a
25-year-old with 35 years of service has both values well inside their
ranges, and sits 2.39 standardised units from the nearest real person
against 0.25 for a typical one.

``` r

ilm_scenario(fit, age = c(45, 25), service = c(23, 35))
#>  age service estimate lower upper extrapolating
#>   45      23    3.095 3.003 3.180         FALSE
#>   25      35    2.820 1.292 4.192          TRUE
```

A flagged scenario is reported rather than refused. Sometimes
extrapolating is the point – but the reader should know they are doing
it, and nothing in the number tells them.

## See also

[`vignette("workflow")`](https://huttoncp.github.io/illume/articles/workflow.md)
for where these sit in an analysis,
[`?ilm_power`](https://huttoncp.github.io/illume/reference/ilm_power.md),
[`?ilm_scaffold`](https://huttoncp.github.io/illume/reference/ilm_scaffold.md)
and
[`?ilm_power_design`](https://huttoncp.github.io/illume/reference/ilm_power_design.md)
for planning without a pilot,
[`?ilm_effects`](https://huttoncp.github.io/illume/reference/ilm_effects.md),
[`?ilm_scenario`](https://huttoncp.github.io/illume/reference/ilm_scenario.md),
and [`?ilm_ame`](https://huttoncp.github.io/illume/reference/ilm_ame.md)
for the marginal effect of a one-unit change rather than a named
scenario.
