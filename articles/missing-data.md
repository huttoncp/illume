# Missing data

Most models drop incomplete rows without telling you. That is sometimes
exactly right and sometimes ruinous, and which it is depends on a
question that has a testable answer.

``` r

library(illume)
```

## The question is not how much is missing

It is whether the missingness is related to what you are modelling. The
usual framing – missing completely at random, at random, not at random –
is about untestable conditions, and it leaves people either ignoring the
problem or imputing reflexively.

There is a narrower question that the data *can* answer, and it is the
one that decides whether complete-case analysis is safe:

> given the predictors, does the probability of being missing depend on
> the **outcome**?

If it does not, complete cases stay unbiased – even if missingness
depends heavily on the predictors, which is weaker than
missing-completely-at-random and much more often true. If it does,
complete cases are biased and imputation is the remedy.

``` r

ilm_check_missing(d, y ~ x + z)
```

This fits a model for being-missing on the outcome **conditional on the
predictors**, rather than asking whether missingness and the outcome are
marginally associated. The distinction is not academic: on data where
missingness follows a covariate that the outcome also depends on, the
marginal test returns p = 0.273 – a real association, carrying no
information once that covariate is held fixed – while the conditional
test returns 0.454 and 0.039 in the two cases that matter, correctly
separating them.

It also reports how much is missing, per column and per row, and whether
any column is missing in a pattern that makes it unusable.

## Imputing

``` r

im <- ilm_impute(d, m = 20)
```

`m = 20` datasets, not one. Filling in a value and then treating it as
data understates the uncertainty, and that understatement is the whole
reason for multiple imputation. The fits are then pooled by Rubin’s
rules:

``` r

fits <- lapply(im$imputations, function(z)
  ilm_model(y ~ x + z, data = z, family = "gaussian"))
ilm_mi_pool(fits)
```

The pooled table reports `fmi`, the fraction of information lost to
missingness. Above about 0.5, more imputations are worth having.

### Which method

`method = "auto"` picks. The three routes answer different situations.

**`"fcs"`, chained equations.** Each incomplete variable is regressed on
the others, cycling until the imputations settle. It is the default
wherever it can be fitted, and it wins wherever it can be fitted – there
is no design in the ablation study where another method beats it on data
it can handle.

**`"lowrank"`, regularised iterative PCA.** Chained equations needs more
rows than predictors. Past that point there is no regression to fit, and
a low-rank reconstruction takes over: `n = 60` with `p = 80` is the case
it exists for. Bootstrapped so the imputations differ, pooled the same
way, after Josse and Husson.

**`"glrm"`, a generalized low rank model.** A loss chosen per column
type rather than squared error on everything, so it imputes
**categorical** columns as categories – which `"lowrank"` says plainly
it cannot, and leaves as `NA`. A category is sampled from its fitted
probabilities rather than set to the most likely level; taking the
argmax would make the `m` datasets agree far more than the data support,
and the pooled variance would be a single fit’s.

``` r

ilm_impute(d, m = 20, method = "glrm")
```

### When imposing structure is worse than nothing

A low-rank method assumes there is a low-rank structure to find. When
there is not, imposing one is worse than filling in column means. The
ablation study measures it, and
[`ilm_impute()`](https://huttoncp.github.io/illume/reference/ilm_impute.md)
detects it by asking whether a rank-`k` fit predicts held-out cells
better than the column means do – warning when it does not, and naming
`method = "fcs"` as the alternative.

That check replaced an earlier one that warned when cross-validation hit
its rank ceiling, and which missed the case entirely: on pure noise the
search picked rank 5 of an allowed 7 and reported nothing.

## Reconstruction accuracy is not the measure

How close the filled values are to the truth is what everyone reports,
and it is not what an analysis needs. What matters is whether the
interval around a coefficient estimated *after* imputation still
contains the truth as often as it claims.

The two can disagree. Single imputation scores respectably on
reconstruction and covers at 0.79 to 0.89 against a nominal 0.95,
because it fills in values and then treats them as though they had been
observed. The studies directory reports both, and where they disagree
the coverage decides.

## Missingness as data

Sometimes the pattern of missingness is itself the thing worth looking
at – which variables go missing together, and whether those groups of
respondents differ.

``` r

ilm_describe_na_all(d)          # how much, and where
ilm_plot_missing(d)             # the pattern
ilm_reduce_na(d)                # structure IN the missingness
ilm_profile_na(d)               # groups of respondents by what they skipped
```

[`ilm_reduce_na()`](https://rdrr.io/pkg/illumex/man/ilm_reduce_na.html)
runs the same dimension reduction on the missingness indicators rather
than the values. A questionnaire where income and savings go missing
together, and separately from the health block, has structure in its
non-response that is worth knowing before deciding what to do about it.

## See also

[`vignette("workflow")`](https://huttoncp.github.io/illume/articles/workflow.md)
for where this sits,
[`?ilm_check_missing`](https://rdrr.io/pkg/illumex/man/ilm_check_missing.html),
[`?ilm_impute`](https://huttoncp.github.io/illume/reference/ilm_impute.md),
[`?ilm_mi_pool`](https://huttoncp.github.io/illume/reference/ilm_mi_pool.md),
and [`?ilm_glrm`](https://rdrr.io/pkg/illumex/man/ilm_glrm.html) for the
low-rank machinery underneath two of the three routes.
