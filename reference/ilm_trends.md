# Estimated marginal slopes

The slope of a continuous predictor within each level of a factor, with
the other predictors averaged over. This is what takes a significant
interaction between a factor and a covariate apart.

## Usage

``` r
ilm_trends(
  object,
  specs,
  var,
  at = NULL,
  weights = c("equal", "proportional", "cells"),
  delta = NULL,
  level = 0.95,
  df = "auto"
)
```

## Arguments

- object:

  A fitted `"ilm_model"` object.

- specs:

  Factor(s) to compute the slope within, as a character vector.

- var:

  The continuous predictor whose slope is wanted.

- at:

  Named list fixing other predictors, as in
  [`ilm_emmeans()`](https://huttoncp.github.io/illume/reference/ilm_emmeans.md).
  The slope of `var` is constant across these unless `var` enters
  non-linearly or interacts with them, in which case where you evaluate
  it matters and the result says so.

- weights:

  How to average over factors not in `specs`; see
  [`ilm_emmeans()`](https://huttoncp.github.io/illume/reference/ilm_emmeans.md).

- delta:

  Step for the numerical derivative. The default is 1/1000 of the
  observed range of `var`, which is exact for any linear term and small
  enough elsewhere.

- level:

  Confidence level.

- df:

  Denominator degrees of freedom for the within-level tests: one of
  `"auto"`, `"satterthwaite"`, `"kenward-roger"`, `"residual"`,
  `"asymptotic"`, or a single number. See
  [`ilm_denom_df()`](https://huttoncp.github.io/illume/reference/ilm_denom_df.md).
  `"auto"` gives an exact t where nothing was integrated out,
  Satterthwaite for a gaussian mixed model, and a z test otherwise.

## Value

A data frame of class `"ilm_emm"`, one row per level of `specs`, with
`estimate`, `se`, `df`, `statistic`, `p.value`, `lower` and `upper`.
Being an `"ilm_emm"` it can be passed straight to
[`ilm_contrast()`](https://huttoncp.github.io/illume/reference/ilm_contrast.md).

## The two questions, which are not the same question

A treatment-by-time interaction says the slope on time depends on the
arm. Two different follow-ups are then available, and they answer
different things:

- **Do the slopes differ between arms?** Pass the result to
  [`ilm_contrast()`](https://huttoncp.github.io/illume/reference/ilm_contrast.md),
  which differences the rows. This is the interaction restated one pair
  at a time.

- **Is the slope different from zero within an arm?** Read the
  `statistic` and `p.value` columns of this table, which test each row
  against zero.

They can disagree in both directions: two arms can have slopes that
differ significantly while neither is distinguishable from zero, and
both can be strongly non-zero while not differing from one another.
Reporting one and describing it as the other is the common way a
within-arm result gets written up as a between-arm claim.

## Why a marginal mean cannot do this

[`ilm_emmeans()`](https://huttoncp.github.io/illume/reference/ilm_emmeans.md)
averages predictions with the covariate held at its mean, which
collapses the very thing the interaction is about. The slope has to be
estimated as a slope. The two share all their machinery – the same
reference grid, the same weighting, the same covariance – and differ
only in what fills the contrast matrix.

## See also

[`ilm_emmeans()`](https://huttoncp.github.io/illume/reference/ilm_emmeans.md)
for means,
[`ilm_contrast()`](https://huttoncp.github.io/illume/reference/ilm_contrast.md)
for differences between the slopes,
[`ilm_anova()`](https://huttoncp.github.io/illume/reference/ilm_anova.md)
for the omnibus interaction this takes apart.

## Examples

``` r
set.seed(1)
d <- data.frame(arm = factor(rep(c("ctl", "trt"), each = 60)),
                week = rep(0:5, 20))
d$y <- 2 + 0.1 * d$week + 0.4 * d$week * (d$arm == "trt") + rnorm(120)
f <- ilm_model(y ~ arm * week, data = d, family = "gaussian",
               verbose = FALSE)
tr <- ilm_trends(f, "arm", var = "week")
tr                      # is each arm's slope different from zero?
#> Estimated marginal slopes of `week` within arm
#>   averaged with equal weights; t tests on residual df
#>  arm estimate     se  df statistic p.value   lower  upper
#>  ctl  -0.0058 0.0667 116   -0.0865  0.9312 -0.1379 0.1264
#>  trt   0.5263 0.0667 116    7.8871  0.0000  0.3941 0.6585
#> 
#>   p-values above test each slope against ZERO. To test whether the
#>   slopes DIFFER from one another, pass this to ilm_contrast().
ilm_contrast(tr)        # do the two arms' slopes differ from each other?
#> <ilm_contrast> 1 comparison(s), adjust = none 
#> 
#>   contrast estimate      se  lower upper     p_adj
#>  trt - ctl   0.5321 0.09437 0.3452 0.719 1.228e-07
```
