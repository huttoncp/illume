# How much of an effect runs through a mediator

Splits the effect of a treatment on an outcome into the part that
operates through a mediator and the part that does not, using the
counterfactual definitions rather than a product of coefficients.

## Usage

``` r
ilm_mediate(
  model_m,
  model_y,
  treat,
  mediator,
  control_value = NULL,
  treat_value = NULL,
  sims = 1000L,
  level = 0.95,
  seed = 1L,
  progress = NULL
)
```

## Arguments

- model_m:

  A fitted
  [`ilm_model()`](https://huttoncp.github.io/illume/reference/ilm_model.md)
  for the mediator.

- model_y:

  A fitted
  [`ilm_model()`](https://huttoncp.github.io/illume/reference/ilm_model.md)
  for the outcome, including both the treatment and the mediator among
  its predictors.

- treat:

  Name of the treatment variable, in both models.

- mediator:

  Name of the mediator: the response of `model_m` and a predictor in
  `model_y`.

- control_value, treat_value:

  The two treatment values to contrast. Defaults are 0 and 1, or the
  first two levels of a factor.

- sims:

  Simulation draws.

- level:

  Confidence level for the percentile intervals.

- seed:

  Random seed.

- progress:

  Show a progress bar; see
  [illumex::ilm_progress_arg](https://rdrr.io/pkg/illumex/man/ilm_progress_arg.html).

## Value

An object of class `"ilm_mediate"`: a data frame of estimates and
intervals, with the draws attached.

## What it estimates

For each simulation draw, parameters are drawn from both models'
sampling distributions, the mediator is predicted under each treatment
value, and the outcome is predicted under each of the four combinations
of treatment and mediator. Averaging over the observed units gives

- `ACME` – the average causal mediation (indirect) effect, reported at
  each treatment value because they differ whenever the outcome model
  has a treatment-by-mediator interaction;

- `ADE` – the average direct effect, likewise;

- `total` – their sum, which is the same whichever pairing is used;

- `prop_mediated` – ACME over total, which is a ratio of estimates and
  behaves badly when the total is near zero. It is reported and should
  be read with that in mind.

With a linear outcome model and no interaction this reduces exactly to
the product of the treatment-to-mediator and mediator-to-outcome
coefficients. The tests check that it does, to 1e-10.

## The assumption that cannot be checked

A decomposition needs no unmeasured confounding of the treatment-outcome
relationship, of the treatment-mediator relationship, and of the
**mediator-outcome** relationship. The first two can be addressed by
design; the third rarely can, because the mediator was not randomised –
whatever makes someone's mediator high may also make their outcome high
for reasons that have nothing to do with the treatment.

Nothing in the data tests this.
[`ilm_mediate_sens()`](https://huttoncp.github.io/illume/reference/ilm_mediate_sens.md)
asks the answerable question instead: how strong would such confounding
have to be before the indirect effect went away.

## References

Imai, K., Keele, L. and Tingley, D. (2010). A general approach to causal
mediation analysis. *Psychological Methods* 15, 309-334.

## See also

[`ilm_mediate_sens()`](https://huttoncp.github.io/illume/reference/ilm_mediate_sens.md)
for the untestable assumption,
[`ilm_dag_model()`](https://huttoncp.github.io/illume/reference/ilm_dag_model.md)
for whether the adjustment sets are available at all.

## Examples

``` r
set.seed(1); n <- 400
d <- data.frame(x = rbinom(n, 1, 0.5), c = rnorm(n))
d$m <- 0.3 + 0.7 * d$x + 0.2 * d$c + rnorm(n)
d$y <- 1 + 0.4 * d$x + 0.6 * d$m + 0.1 * d$c + rnorm(n)
fm <- ilm_model(m ~ x + c, data = d, family = "gaussian", verbose = FALSE)
fy <- ilm_model(y ~ x + m + c, data = d, family = "gaussian",
                verbose = FALSE)
ilm_mediate(fm, fy, treat = "x", mediator = "m", sims = 200)
#> Causal mediation: x -> m -> outcome
#>   x = 1 against 0, 400 units, 200 simulation draws
#>               effect estimate   lower  upper    p
#>       ACME (control)   0.3101 0.16360 0.4927 0.00
#>       ACME (treated)   0.3101 0.16360 0.4927 0.00
#>        ADE (control)   0.2616 0.07091 0.4922 0.01
#>        ADE (treated)   0.2616 0.07091 0.4922 0.01
#>         Total effect   0.5717 0.33440 0.7945 0.00
#>  Proportion mediated   0.5513 0.33540 0.8466   NA
#> 
#>   The outcome model has no x:m interaction, so the two ACMEs are equal by
#>   construction rather than by evidence. Adding one lets the indirect
#>   effect differ by arm, and is worth doing before concluding it does not.
#> 
#>   All of this rests on there being no unmeasured confounding of the
#>   MEDIATOR and the outcome, which the data cannot check because the
#>   mediator was not randomised. ilm_mediate_sens() asks how strong such
#>   confounding would have to be to remove the indirect effect.
```
