# Predicted survival curve from a fitted model

The probability of surviving past each time, for one or more covariate
patterns, with confidence intervals.

## Usage

``` r
ilm_survival(
  object,
  newdata = NULL,
  times = NULL,
  conf = 0.95,
  groups = c("fitted", "typical", "population"),
  nsim = 200L,
  seed = 1L
)
```

## Arguments

- object:

  A fitted `"ilm_model"` with a time-to-event family.

- newdata:

  Covariate patterns, one row each. With none, the curve is drawn at the
  median of each numeric predictor and the commonest level of each
  factor, for a typical group.

- times:

  Times at which to evaluate. With none, a grid spanning the observed
  follow-up.

- conf:

  Confidence level.

- groups:

  `"fitted"` (the default), `"typical"` or `"population"`; see "Which
  groups". Only a model with random effects or a correlation over time
  has a choice to make.

- nsim:

  Draws behind the intervals of `"fitted"` and `"population"`. A typical
  group's interval uses none.

- seed:

  Random seed for those draws.

## Value

A data frame with `row`, `time`, `surv`, `lower` and `upper`.

## Details

All three accelerated failure time families have a closed-form survivor
function in `z = (log t - eta) / scale`, and a flexible (Royston-Parmar)
model one in its baseline spline plus `eta`, so the curve comes straight
from the linear predictor.

## Which groups

In a model with random effects a curve is for some group, and `groups`
says which, as it does for
[`predict.ilm_model()`](https://huttoncp.github.io/illume/reference/predict.ilm_model.md):

- `"fitted"` (the default) draws each row at **its own group's**
  estimated effects. A row whose group the fit has not seen, or
  `newdata` without the grouping columns, is an error naming the two
  choices below.

- `"typical"` draws each row with every random effect at zero: the curve
  of a **typical group**. With no `newdata` the curve is at a typical
  covariate pattern, which belongs to no group, and this is what it
  gives.

- `"population"` averages the **survival curve itself** over the
  distribution of random effects, by quadrature at each time: the curve
  for the population of groups as a whole. It is flatter than the
  typical group's curve, not the curve at zero, since survival is not
  linear in the linear predictor.

Intervals for a typical group use no draws. The curve is worked out at
`eta - c * se` and `eta + c * se`, where `se` is the standard error of
the linear predictor `eta` and `c` the normal quantile for `conf`; this
is the delta method on `eta`, and it gives the same interval on every
run. The standard error comes from the covariance of the fixed effects
or, when the model has a smooth other than a factor-smooth interaction
(`bs = "fs"`), from the joint covariance of the fixed effects and those
smooths' penalised coefficients. Every random effect is at zero, and so
is a factor-smooth interaction. For `"fitted"` and `"population"` the
intervals are percentile intervals from `nsim` joint draws of the
coefficients – and, for `"fitted"`, of the groups' own effects, as
[`predict.ilm_model()`](https://huttoncp.github.io/illume/reference/predict.ilm_model.md)
draws them. Uncertainty in an accelerated failure time family's scale
parameter is not included in any of these intervals, so they are
slightly narrow, most in the tail beyond the last observed event.

## See also

[`ilm_plot_survival()`](https://huttoncp.github.io/illume/reference/ilm_plot_survival.md),
which draws it against the Kaplan-Meier estimate,
[`ilm_surv()`](https://huttoncp.github.io/illume/reference/ilm_surv.md).

## Examples

``` r
set.seed(1)
d <- data.frame(x = rnorm(200))
tt <- exp(1.5 + 0.8 * d$x + 0.7 * log(rexp(200)))
ct <- rexp(200, rate = 1 / (2 * median(tt)))
d$time <- pmin(tt, ct); d$event <- as.integer(tt <= ct)
f <- ilm_model(time ~ x, data = d, family = "weibull",
               censor = ilm_surv(d$time, d$event), verbose = FALSE)
head(ilm_survival(f, newdata = data.frame(x = c(-1, 1))))
#>   row       time      surv     lower     upper
#> 1   1 0.02367112 0.9969355 0.9960635 0.9976145
#> 2   1 0.23229399 0.9385701 0.9217585 0.9518640
#> 3   1 0.44091686 0.8621862 0.8264971 0.8910207
#> 4   1 0.64953972 0.7804836 0.7272392 0.8245985
#> 5   1 0.85816259 0.6986832 0.6307910 0.7565302
#> 6   1 1.06678545 0.6197232 0.5406977 0.6891251
```
