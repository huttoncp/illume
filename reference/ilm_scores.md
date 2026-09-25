# Proper scoring rules for a categorical outcome

R-squared has no agreed definition for a nominal outcome, but **proper
scoring rules** do. A scoring rule is "proper" if it is optimised by
reporting your true beliefs, which makes it a fair way to score
probabilistic predictions.

## Usage

``` r
ilm_scores(object, groups = c("fitted", "typical"), conditional = NULL)
```

## Arguments

- object:

  A fitted `"ilm_model"` object.

- groups:

  `"fitted"` (the default) or `"typical"`: which random effects the
  probabilities being scored use, as in
  [`ilm_fitted()`](https://huttoncp.github.io/illume/reference/ilm_fitted.md).

- conditional:

  Deprecated. `TRUE` is `groups = "fitted"`, and `FALSE` is
  `groups = "typical"`.

## Value

A named numeric vector: `log_score`, `brier`, `accuracy`, and for an
ordinal fit `rps`, on the 0 to 1 scale.

## Details

The log score is the average of `-log(p)` at the observed category:
heavily penalises confident mistakes. The Brier score is the average
squared distance between the predicted probability vector and the
observed one-hot outcome, ranging from 0 to 2. Lower is better for both.

Accuracy is also returned, though it ignores how confident the
predictions were and is the least informative of the three.

Defined for a multinomial fit, an ordinal one, and a binomial one with a
0/1 response, whose single fitted probability is scored as the two
categories it implies. Other families have no category probabilities to
score; compare those fits with
[`AIC()`](https://rdrr.io/r/stats/AIC.html) or
[`ilm_anova()`](https://huttoncp.github.io/illume/reference/ilm_anova.md),
and check their predictions with
[`ilm_check_predictive()`](https://huttoncp.github.io/illume/reference/ilm_check_predictive.md).

An ordinal fit also gets the **ranked probability score**, which the
other two do not provide: it compares the cumulative distributions, so a
prediction that puts its weight one category away from the truth scores
better than one that puts it three away. The log and Brier scores treat
every miss alike, as they should for categories with no order.

## References

Gneiting, T., & Raftery, A. E. (2007). Strictly proper scoring rules,
prediction, and estimation. *Journal of the American Statistical
Association*, 102(477), 359–378.

Brier, G. W. (1950). Verification of forecasts expressed in terms of
probability. *Monthly Weather Review*, 78(1), 1–3.
