# Fit indices for a multinomial mixed model

A
[`performance::model_performance()`](https://easystats.github.io/performance/reference/model_performance.html)
method. Returns AIC, AICc, BIC (both sample-size conventions), the
log-likelihood, degrees of freedom, McFadden's R-squared, and the
scoring rules from
[`ilm_scores()`](https://huttoncp.github.io/illume/reference/ilm_scores.md).

## Usage

``` r
model_performance.ilm_model(model, metrics = "all", ..., verbose = TRUE)
```

## Arguments

- model:

  A fitted `"ilm_model"` object.

- metrics:

  Unused; present for compatibility with the generic.

- ...:

  Unused.

- verbose:

  Logical. Message when the fit failed its checks.

## Value

A one-row data frame of class `"performance_model"`.

## What is deliberately absent

Nakagawa's marginal and conditional R-squared, and the ICC, are **not**
reported. Both require a distribution-specific residual variance, and a
nominal multinomial outcome does not have one. Reporting them would mean
inventing a quantity rather than estimating it.

McFadden's R-squared is relative to an intercept-only model with the
same grouping random effects, so it answers "how much do the predictors
add?" rather than "how much does the whole model explain?". Unlike
R-squared in linear regression, values around 0.2 to 0.4 already
indicate a good fit.

The scoring rules are computed **in sample** and conditional on the
fitted random effects, so both are optimistic. Use cross-validation or a
held-out set if you need an honest predictive comparison.

## References

McFadden, D. (1974). Conditional logit analysis of qualitative choice
behavior. In P. Zarembka (Ed.), *Frontiers in Econometrics*. Academic
Press.

Nakagawa, S., & Schielzeth, H. (2013). A general and simple method for
obtaining R-squared from generalized linear mixed-effects models.
*Methods in Ecology and Evolution*, 4(2), 133–142. (The approach not
used here, and why: it needs a distribution-specific variance.)
