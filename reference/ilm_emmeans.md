# Estimated marginal means

The model's predicted mean for each level of the variables asked for,
with the other predictors averaged over rather than left at whatever
value happens to be first.

## Usage

``` r
ilm_emmeans(
  object,
  specs,
  at = NULL,
  weights = c("equal", "proportional", "cells"),
  type = c("link", "response"),
  level = 0.95
)
```

## Arguments

- object:

  A fitted
  [`ilm_model()`](https://huttoncp.github.io/illume/reference/ilm_model.md).

- specs:

  Variables to keep, as a character vector. Everything else is averaged
  over or held at its mean.

- at:

  Named list fixing particular values, for instance
  `list(age = c(40, 60))`.

- weights:

  `"equal"`, `"proportional"` or `"cells"`; see above.

- type:

  `"link"` or `"response"`.

- level:

  Confidence level.

## Value

An object of class `"ilm_emm"`: a data frame of the grid with
`estimate`, `se`, `lower`, `upper`, plus the contrast machinery it
carries.

## What is averaged, and how

A marginal mean is a prediction on a grid of every factor combination,
with numeric predictors held at their means. The factors *not* named in
`specs` have to be averaged over, and how they are averaged changes the
answer whenever the design is unbalanced:

- `weights = "equal"` treats every cell alike, answering "what would
  these groups look like in a balanced design". This is the usual
  default and matches what most software calls an estimated marginal
  mean.

- `weights = "proportional"` weights each cell by the product of the
  *marginal* frequencies of the variables being averaged over, answering
  "what do these groups look like in a population with this sample's
  composition".

- `weights = "cells"` weights by the *joint* frequency actually
  observed, so each group is averaged over its own covariate
  distribution.

The first two keep a comparison clean: the averaged-over part is
identical for every level of `specs`, so a difference between two
marginal means is the model's effect and nothing else. `"cells"` does
not, because each group gets its own mix. When the averaged-over
variables are associated with `specs`, a difference of cell-weighted
means carries that difference in composition as well as the effect, and
[`ilm_contrast()`](https://huttoncp.github.io/illume/reference/ilm_contrast.md)
says so. That is sometimes the quantity you want – it describes the
groups as they are – but it is not an adjusted comparison.

## An ordered response

For an ordinal fit the link-scale values are marginal means of the
**latent scale** – the linear predictor the thresholds cut up – and not
of the categories, which have no mean to take. A contrast between two of
them is a difference in log odds of being in a higher category, constant
across cuts by the same assumption
[`ilm_check_proportional()`](https://huttoncp.github.io/illume/reference/ilm_check_proportional.md)
tests.

`type = "response"` gives each category's PROBABILITY instead, computed
in every cell of the grid and averaged over the cells, with a
delta-method standard error that carries the thresholds' uncertainty as
well as the slopes'.
[`ilm_contrast()`](https://huttoncp.github.io/illume/reference/ilm_contrast.md)
then compares groups within each category as differences in probability.
This agrees with `emmeans` on a
[`MASS::polr()`](https://rdrr.io/pkg/MASS/man/polr.html) fit with
`mode = "prob"`.

## A multinomial response

Every category gets its own row for every level of `specs`.

On the **link** scale the value is the category's centred log-odds: its
log-probability less the average log-probability over all the
categories, which is what the sum-to-zero coefficients describe. These
are exact linear combinations, as for any other family, and a contrast
between two groups within a category is a difference of log-odds against
the same average.

On the **response** scale the value is the category's PROBABILITY,
computed in every cell of the grid and then averaged over the cells with
the chosen weights – the probabilities of each group sum to 1 – with a
delta-method standard error and an interval formed on the logit scale,
so it stays inside 0 and 1.
[`ilm_contrast()`](https://huttoncp.github.io/illume/reference/ilm_contrast.md)
then compares groups within each category as differences in probability.
This is what `emmeans` computes for an
[`nnet::multinom()`](https://rdrr.io/pkg/nnet/man/multinom.html) fit
with `mode = "prob"`, and the two agree; the coefficients differ between
the packages, because `nnet` codes against a baseline category, but the
probabilities do not.

## Which scale

Marginal means are computed on the **link** scale, where they are exact
linear combinations of the coefficients and their variance is `L V L'`
with no approximation. `type = "response"` back-transforms the
endpoints, which keeps the interval's coverage but means the reported
centre is a median rather than a mean on that scale. For a model with a
random effect that centre is also conditional on the group rather than
population-averaged –
[`ilm_ame()`](https://huttoncp.github.io/illume/reference/ilm_ame.md) is
the population-averaged quantity.

## See also

[`ilm_contrast()`](https://huttoncp.github.io/illume/reference/ilm_contrast.md)
to compare them,
[`ilm_ame()`](https://huttoncp.github.io/illume/reference/ilm_ame.md)
for the average marginal effect of a predictor.

## Examples

``` r
set.seed(1); n <- 200
d <- data.frame(g = factor(sample(c("a", "b", "c"), n, TRUE)), x = rnorm(n))
d$y <- 1 + 0.5 * (d$g == "b") + 0.2 * d$x + rnorm(n)
fit <- ilm_model(y ~ g + x, data = d, family = "gaussian", verbose = FALSE)
ilm_emmeans(fit, "g")
#> <ilm_emm> marginal means of g (link scale, equal weights)
#> 
#>  g estimate     se  lower upper
#>  a    1.253 0.1257 1.0052 1.501
#>  b    1.263 0.1239 1.0186 1.507
#>  c    1.027 0.1259 0.7787 1.275
#>   Compare them with ilm_contrast().
```
