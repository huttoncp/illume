# What the model says would happen under specified scenarios

Sets named predictors to chosen values and reports the outcome the model
implies, with an interval, for each combination. Built for the question
a client asks before scaling something: what would this do to our
population.

## Usage

``` r
ilm_scenario(
  object,
  ...,
  over = c("sample", "reference"),
  sims = 1000L,
  level = 0.95,
  contrast = FALSE,
  seed = 1L,
  progress = NULL
)
```

## Arguments

- object:

  A fitted
  [`ilm_model()`](https://craig-hutton.github.io/illume/reference/ilm_model.md).

- ...:

  Named predictors and the values to set them to, for instance
  `dose = c(0, 10, 20), age = 40`. Every combination is used.

- over:

  `"sample"` or `"reference"`; see above.

- sims:

  Parameter draws for the interval.

- level:

  Confidence level.

- contrast:

  Compare scenarios against one another: `FALSE`, `"first"` (each
  against the first) or `"pairwise"`.

- seed:

  Random seed.

- progress:

  Show a progress bar; see
  [ilm_progress_arg](https://craig-hutton.github.io/illume/reference/ilm_progress_arg.md).

## Value

An object of class `"ilm_scenario"`: a data frame of the grid with
`estimate`, `lower` and `upper`, plus any contrasts.

## Which population

`over = "sample"` (the default) sets the named predictors for **every
observed unit**, leaves each unit's other covariates as they are,
predicts, and averages. That is standardisation, and it answers "what
would happen to this population". When the adjustment set is valid it is
also the causal estimand.

`over = "reference"` builds a single row with the named predictors at
their scenario values and everything else at its mean or modal value,
and predicts once. That answers "what would happen to a unit that looks
like this".

They are different questions and they give different answers under any
non-linear link, because the average of a prediction is not the
prediction at the average. The printed output says which was used.

## Extrapolation

A model returns a number for a dose nobody received, and the interval
around it looks like any other, because the interval carries uncertainty
about the coefficients rather than about whether the functional form
holds out there. Each scenario is therefore checked twice: each value
against that predictor's observed range, and the whole combination
against how far it sits from the nearest real unit, in standardised
units, compared with how far real units sit from each other. A flagged
scenario is not refused – it is reported, because sometimes
extrapolating is the point and the reader should know they are.

## References

Hernán, M. A. and Robins, J. M. (2020). *Causal Inference: What If*.
Chapman & Hall/CRC.

## See also

[`ilm_ame()`](https://craig-hutton.github.io/illume/reference/ilm_ame.md)
for the effect of a one-unit change rather than a named scenario,
[`ilm_emmeans()`](https://craig-hutton.github.io/illume/reference/ilm_emmeans.md)
for group means,
[`ilm_dag_model()`](https://craig-hutton.github.io/illume/reference/ilm_dag_model.md)
for whether the adjustment set licenses a causal reading.

## Examples

``` r
set.seed(1); n <- 400
d <- data.frame(dose = runif(n, 0, 20), age = rnorm(n, 50, 10))
d$y <- rbinom(n, 1, plogis(-3 + 0.12 * d$dose + 0.03 * d$age))
f <- ilm_model(y ~ dose + age, data = d, family = "binomial",
               verbose = FALSE)
ilm_scenario(f, dose = c(0, 10, 20), sims = 200, contrast = "first")
#> Scenario projection (averaged over the observed units)
#>   every unit set to the scenario's values, its own covariates kept
#>   held as observed: age
#>  dose estimate  lower  upper extrapolating
#>     0   0.2067 0.1417 0.2813         FALSE
#>    10   0.4066 0.3610 0.4602         FALSE
#>    20   0.6434 0.5528 0.7217         FALSE
#> 
#>   Differences between scenarios
#>          contrast estimate  lower  upper p
#>  dose=10 - dose=0   0.1999 0.1333 0.2548 0
#>  dose=20 - dose=0   0.4367 0.2839 0.5587 0
```
