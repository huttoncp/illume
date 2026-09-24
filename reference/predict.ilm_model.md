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
  marginal = FALSE,
  se.fit = FALSE,
  interval = c("none", "confidence"),
  level = 0.95,
  nsim = 200L,
  ndraw = 200L,
  seed = 1L,
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

- marginal:

  Logical. Average over the random-effect distribution
  (population-averaged) rather than setting it to zero (conditional).

- se.fit:

  Logical. Return standard errors.

- interval:

  `"none"` or `"confidence"`.

- level:

  Numeric. Interval coverage, default 0.95.

- nsim:

  Integer. Parameter draws used for uncertainty.

- ndraw:

  Integer. Random-effect draws used when `marginal = TRUE` for a
  multinomial outcome; a single linear predictor is averaged by
  quadrature and does not use it.

- seed:

  Integer. Random seed, so results are reproducible.

- ...:

  Unused.

## Value

A matrix of probabilities (or linear predictors), or a factor for
`type = "class"`. When standard errors or intervals are requested, a
list with `fit`, `se.fit`, `lower`, `upper`, `level` and `joint`.

## Conditional versus population-averaged

This is the choice that matters most, and there is no safe default that
suits everyone.

With `marginal = FALSE` the random effects are set to zero, giving the
probabilities for a **typical** group – one exactly at the population
average. With `marginal = TRUE` the prediction is averaged over the
distribution of random effects, giving the probabilities for the
**population as a whole**.

These differ, sometimes substantially, because averaging and the softmax
transform do not commute: the average of the transformed values is not
the transform of the average. The population-averaged probabilities are
pulled toward being more even across categories. Which you want depends
on the question – "what do I expect for an average subject?" or "what
proportion of the population falls in each category?"

Under an **identity link** the two coincide exactly, because the random
effects have mean zero and nothing nonlinear stands between. `marginal`
is then answered in closed form rather than by simulation, so the result
does not depend on `ndraw` and carries no Monte Carlo noise.

A **random slope** is averaged over as a slope. The amount being
integrated over then depends on the row – it grows with distance from
wherever the slope is centred – so the marginal and conditional curves
separate by more at the ends of the range than in the middle. Averaging
such a term as if it were an intercept understates that, and the error
grows with the slope variance and with distance from centre.

An **AR(1) or CAR(1) term** is averaged over too: at any one row its
latent value has the stationary distribution, whatever the time.

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
