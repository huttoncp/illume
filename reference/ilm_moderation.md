# Does the effect of X depend on something else?

Searches a set of candidate moderators for evidence that the effect of a
treatment or exposure varies, testing each by the joint test of its
interaction block and adjusting for having looked at several.

## Usage

``` r
ilm_moderation(
  object,
  x = NULL,
  moderators = NULL,
  adjust = "holm",
  split = FALSE,
  data = NULL,
  n_keep = 3L,
  seed = 1L,
  progress = NULL
)
```

## Arguments

- object:

  A fitted
  [`ilm_model()`](https://huttoncp.github.io/illume/reference/ilm_model.md),
  or an
  [`ilm_dag_model()`](https://huttoncp.github.io/illume/reference/ilm_dag_model.md).

- x:

  The treatment, exposure or intervention whose effect might be
  moderated. Required, **except** for an
  [`ilm_dag_model()`](https://huttoncp.github.io/illume/reference/ilm_dag_model.md),
  where the graph has already named the exposure. There is no sensible
  default otherwise: `x:m` is the same term whichever of the two you
  call the moderator, and only the interpretation distinguishes them.

- moderators:

  Candidate moderators. Defaults to the model's other fixed-effect terms
  – variables you already judged worth adjusting for – which keeps the
  multiplicity burden honest. Naming them explicitly can reach anything
  in the model frame, or in `data`.

- adjust:

  Passed to
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html); any
  method it supports, including `"holm"` (the default), `"BH"`,
  `"bonferroni"`, `"BY"`, `"hochberg"`, `"hommel"` and `"none"`.

- split:

  Use an honest sample split instead of an adjustment. See above: better
  calibrated, roughly half the power.

- data:

  Needed only for a moderator that is not in the model frame.

- n_keep:

  Refits to retain on the result, strongest first, so
  [`ilm_plot_moderation()`](https://huttoncp.github.io/illume/reference/ilm_plot_moderation.md)
  can draw them without refitting. The rest are refitted on demand.

- seed:

  Random seed, used only when `split = TRUE`.

- progress:

  Show a progress bar. Defaults to
  [`interactive()`](https://rdrr.io/r/base/interactive.html).

## Value

An object of class `"ilm_moderation"`: a data frame with one row per
candidate (`moderator`, `df`, `statistic`, `p`, `p_adj`), carrying the
retained refits and the exposure.

## What this is for, and what it is not

This is the **exploratory** tool. A moderator you hypothesised a priori
belongs in the model, where
[`ilm_anova()`](https://huttoncp.github.io/illume/reference/ilm_anova.md)
tests it and
[`ilm_emmeans()`](https://huttoncp.github.io/illume/reference/ilm_emmeans.md)
and
[`ilm_contrast()`](https://huttoncp.github.io/illume/reference/ilm_contrast.md)
describe it – and it owes no multiplicity penalty, because you did not
go looking. Any such interaction already present in `object` is
therefore excluded from the search and said so, rather than being
charged for a search it was not part of.

## Why the adjustment is not optional

Trying several moderators and reporting the strongest is a search. With
six candidates and no moderation present at all, that procedure rejected
29.3% of the time at a nominal 5%. Every adjustment method brings it
back: Holm, BH and Bonferroni all landed at 0.027 in the same
simulation.

## Holm, or an honest split

`split = TRUE` picks the strongest moderator on half the data and tests
it on the other half, which needs no adjustment because the test half
was never used to choose. It is better calibrated – 0.053 against Holm's
conservative 0.027 – and costs about half the power (0.320 against 0.547
for a numeric moderator, 0.200 against 0.380 for a three-level factor).
The gap widens with degrees of freedom, so splitting is weakest exactly
where experimental designs live. Holm is the default for that reason.

Splitting is genuinely *required* only when the hypotheses cannot be
counted, as in an open-ended tree search over covariates and split
points. This function searches a named set, so they can be.

## See also

[`ilm_plot_moderation()`](https://huttoncp.github.io/illume/reference/ilm_plot_moderation.md)
to see one,
[`ilm_anova()`](https://huttoncp.github.io/illume/reference/ilm_anova.md)
for a moderator you specified in advance,
[`ilm_trends()`](https://huttoncp.github.io/illume/reference/ilm_trends.md)
for slopes within levels.

## Examples

``` r
# \donttest{
set.seed(1)
d <- data.frame(tx = rbinom(300, 1, 0.5), age = rnorm(300),
                site = factor(sample(c("a", "b", "c"), 300, TRUE)))
d$y <- 0.4 * d$tx + 0.3 * d$age + 0.8 * d$tx * (d$site == "c") + rnorm(300)
fit <- ilm_model(y ~ tx + age + site, data = d, verbose = FALSE)
#> ilm_model(): family = "gaussian", inferred from `y`: continuous values from -2.6 to 3.97. Pass `family` to choose another.
ilm_moderation(fit, x = "tx")
#> Moderation of `tx`: 2 candidates, holm-adjusted
#>  moderator df statistic       p   p_adj
#>       site  2     6.440 0.00183 0.00366
#>        age  1     0.279 0.59700 0.59700
#> 
#>   `site`: the effect of `tx` is not the same for everyone.
#>   ilm_plot_moderation() to see one. Each test is the joint test of the
#>   interaction block, so `df` is what it cost.
# }
```
