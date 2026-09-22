# Effect sizes for the fixed effects

Reports each coefficient on the scale its family puts it on – an odds
ratio, an incidence rate ratio, a time ratio, a hazard ratio – with an
interval built on the link scale and transformed, which is exact and
asymmetric rather than symmetric and possibly negative.

## Usage

``` r
ilm_effects(object, level = 0.95, standardise = FALSE, vcov = NULL)
```

## Arguments

- object:

  A fitted
  [`ilm_model()`](https://craig-hutton.github.io/illume/reference/ilm_model.md).

- level:

  Confidence level.

- standardise:

  `FALSE`, `TRUE` (one standard deviation) or `"gelman"` (two), for
  numeric predictors only.

- vcov:

  Optional covariance matrix to use instead of the model's, for instance
  from
  [`ilm_vcov_cluster()`](https://craig-hutton.github.io/illume/reference/ilm_vcov_cluster.md)
  or
  [`ilm_svy_vcov()`](https://craig-hutton.github.io/illume/reference/ilm_svy_vcov.md).

## Value

A data frame of `term`, `estimate` (on the link scale), `se`, `effect`
(on the family's scale), `lower`, `upper`, `z` and `p`.

## Details

For a gaussian response there is no ratio to take, so it reports the
coefficient, the standardised coefficient, and the partial variance each
term explains.

## Conditional against marginal

In a model with random effects, `exp(beta)` is the ratio for a **given
cluster**: two patients in the same hospital, one exposed. It is not the
ratio for the population, and the gap widens as the random effects grow.
Which one a reader wants is almost always the population one, and which
one they are usually given is this one. It is labelled, and
[`ilm_ame()`](https://craig-hutton.github.io/illume/reference/ilm_ame.md)
computes the population-averaged quantity.

## Standardised coefficients

`standardise = TRUE` divides each numeric predictor's coefficient by
that predictor's standard deviation, so it reads per standard deviation.
Binary and categorical predictors are left alone: dividing by the
standard deviation of a 0/1 variable answers what happens per standard
deviation of treatment, which is not a question anyone asks. Gelman's
suggestion of dividing by **two** standard deviations makes a numeric
predictor comparable with a binary one; `standardise = "gelman"` does
that.

## See also

[`ilm_ame()`](https://craig-hutton.github.io/illume/reference/ilm_ame.md)
for the population-averaged effect,
[`ilm_emmeans()`](https://craig-hutton.github.io/illume/reference/ilm_emmeans.md)
for group means,
[`ilm_vcov_cluster()`](https://craig-hutton.github.io/illume/reference/ilm_vcov_cluster.md)
for a robust covariance to pass in.

## Examples

``` r
set.seed(1); n <- 400
d <- data.frame(x = rnorm(n), g = factor(sample(c("a", "b"), n, TRUE)))
d$y <- rbinom(n, 1, plogis(-0.4 + 0.8 * d$x + 0.5 * (d$g == "b")))
f <- ilm_model(y ~ x + g, data = d, family = "binomial", verbose = FALSE)
ilm_effects(f)
#> Fixed effects as odds ratios (95% intervals)
#>   exp(coefficient): a change of one unit MULTIPLIES the odds of the modelled outcome
#>         term estimate     se odds_ratio lower  upper      z         p
#>  (Intercept)  -0.3221 0.1567     0.7246 0.533 0.9852 -2.055 3.984e-02
#>            x   0.7808 0.1235     2.1830 1.714 2.7810  6.324 2.554e-10
#>           gb   0.4409 0.2143     1.5540 1.021 2.3650  2.058 3.961e-02
```
