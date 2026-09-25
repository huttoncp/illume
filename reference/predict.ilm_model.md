# Predictions from a fitted model

Returns the fitted mean on the response scale – one probability per
category for a multinomial or ordinal outcome, as
`nnet::multinom(type = "probs")` does – or the linear predictor,
optionally with standard errors and intervals.

## Usage

``` r
# S3 method for class 'ilm_model'
predict(
  object,
  newdata = NULL,
  type = c("response", "link", "class"),
  groups = c("typical", "population"),
  se.fit = FALSE,
  interval = c("none", "confidence"),
  level = 0.95,
  nsim = 200L,
  ndraw = 200L,
  seed = 1L,
  marginal = NULL,
  ...
)
```

## Arguments

- object:

  A fitted `"ilm_model"` object.

- newdata:

  Optional data frame. If omitted, predictions are for the data used to
  fit the model.

- type:

  `"response"` for probabilities (the default), `"link"` for linear
  predictors, or `"class"` for the most likely category.

- groups:

  `"typical"` (the default) for a group with every random effect at
  zero, or `"population"` for the average over the groups. See "Which
  groups".

- se.fit:

  Logical. Return standard errors.

- interval:

  `"none"` or `"confidence"`.

- level:

  Numeric. Interval coverage, default 0.95.

- nsim:

  Integer. Parameter draws used for uncertainty.

- ndraw:

  Integer. Random-effect draws used when `groups = "population"` for a
  multinomial outcome; a single linear predictor is averaged by
  quadrature and does not use it.

- seed:

  Integer. Random seed, so results are reproducible.

- marginal:

  Deprecated. `TRUE` is `groups = "population"`, and `FALSE` is
  `groups = "typical"`.

- ...:

  Unused.

## Value

A matrix of probabilities (or linear predictors), or a factor for
`type = "class"`. When standard errors or intervals are requested, a
list with `fit`, `se.fit`, `lower`, `upper`, `level` and `joint`.

## Which groups

In a model with random effects every prediction is for some group, and
this is the choice that matters most. There is no safe default that
suits everyone. `groups` says which:

- `"typical"` (the default) sets every random effect to zero, giving the
  prediction for a **typical group**: one exactly at the average.

- `"population"` averages the prediction over the distribution of random
  effects, giving it for the **population of groups as a whole**.

These differ, sometimes substantially, because averaging and a nonlinear
inverse link do not commute: the average of the transformed values is
not the transform of the average. Through a softmax, the population's
probabilities are pulled toward being more even across categories. Which
you want depends on the question – "what do I expect for an average
subject?" or "what proportion of the population falls in each category?"

Each fitted group's own effects are what
[`ilm_fitted()`](https://huttoncp.github.io/illume/reference/ilm_fitted.md)
uses, for the rows the model was fitted to, and
[`ilm_ranef()`](https://huttoncp.github.io/illume/reference/ilm_ranef.md)
returns them.

`marginal` is the old name for this choice: `marginal = FALSE` is
`groups = "typical"`, and `marginal = TRUE` is `groups = "population"`.
It still works, with a warning, and will be removed after the next
release.

Under an **identity link** the two coincide exactly, because the random
effects have mean zero and nothing nonlinear stands between.
`"population"` is then answered in closed form rather than by
simulation, so the result does not depend on `ndraw` and carries no
Monte Carlo noise.

A **random slope** is averaged over as a slope. The amount being
integrated over then depends on the row – it grows with distance from
wherever the slope is centred – so the population's curve and the
typical group's separate by more at the ends of the range than in the
middle. Averaging such a term as if it were an intercept understates
that, and the error grows with the slope variance and with distance from
centre.

An **AR(1) or CAR(1) term** is averaged over too: at any one row its
latent value has the stationary distribution, whatever the time. A
**random walk**
([`ilm_rw1()`](https://huttoncp.github.io/illume/reference/ilm_rw1.md))
is not stationary: at a row its variance is the variance per unit of
time multiplied by the time since the row's group started, so each row
is averaged over its own spread, and the population mean moves away from
the typical group's as time passes. New rows need their time and group
for that, which the fit can find only when the walk was given by name,
`ilm_rw1(~ time | group)`; a group the fit has not seen is taken to
start at its earliest time among the new rows, as each fitted group
started at its first.

With **one linear predictor** – every family but the multinomial – a
row's whole latent contribution is a single normal, and the average is
taken by Gauss-Hermite quadrature, with more nodes as the latent SD
grows: better than 1e-12 through a logit up to an SD of 6 on the link
scale, and better than 1e-9 through a complementary log-log up to 5. It
is the same on every call, and free of `ndraw`. A multinomial outcome
has one dimension per category, and is averaged over `ndraw` draws with
common random numbers.

## Uncertainty

Standard errors and intervals come from simulation rather than a
formula, because the softmax makes the quantity nonlinear in the
parameters. Intervals are **percentile** intervals from the simulated
draws, so they always lie within 0 and 1; a symmetric interval on the
probability scale would not.

If the fit was made with `joint = TRUE` the draws include the penalised
smooth coefficients. Without it only the fixed effects vary, which
breaks the correlation described in `ilm_joint_draws()` and distorts
intervals around smooths; a warning says so.
[`ilm_model()`](https://huttoncp.github.io/illume/reference/ilm_model.md)
enables it automatically when the model contains smooths.

Random-effect draws are held fixed across rows and across parameter
draws ("common random numbers"). Without that, Monte Carlo noise would
swamp comparisons between grid points, and `marginaleffects` would be
unable to compute stable numerical derivatives.

## References

Skrondal, A., & Rabe-Hesketh, S. (2009). Prediction in multilevel
generalized linear models. *Journal of the Royal Statistical Society,
Series A*, 172(3), 659–687. (On the distinction between conditional and
marginal prediction.)

## See also

[`ilm_model()`](https://huttoncp.github.io/illume/reference/ilm_model.md),
[`ilm_fitted()`](https://huttoncp.github.io/illume/reference/ilm_fitted.md).
