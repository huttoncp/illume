# The illume workflow

This is a map rather than a tutorial. It walks one analysis from a
sample-size calculation through to reporting, eleven stages, with a
short paragraph and a minimal block for each, and a pointer to the
vignette or help page that covers it properly.

Nothing here is compulsory. Most analyses skip several stages, and a few
use only one. The point of the map is that when you do need a stage, the
function for it is in the same package, speaks the same object, and
names its own remedy when it finds a problem.

``` r

library(illume)
```

------------------------------------------------------------------------

## 1. Power

Before the data exist.
[`ilm_power()`](https://huttoncp.github.io/illume/reference/ilm_power.md)
treats a model as the truth, simulates studies at the sizes you are
considering, refits each, and counts how often the effect is found. It
works for any family the package fits, because it never needs a
closed-form variance.

``` r

pilot <- ilm_model(y ~ dose, data = pilot_data, family = "binomial")
pw <- ilm_power(pilot, n = c(200, 400, 800), effect = 0.4, sims = 500)
ilm_power_n(pw, target = 0.8)
plot(pw)
```

Two numbers it reports that most software does not: the **Monte Carlo
interval** around the power estimate, because 0.80 from 200 replicates
is really 0.74 to 0.86 and a required sample size read off that curve is
a range; and the **convergence rate**, because a replicate that failed
to fit has not detected anything, and dividing only by the ones that
worked gives power conditional on convergence, which is a different and
more flattering quantity.

With no pilot to fit,
[`ilm_power_design()`](https://huttoncp.github.io/illume/reference/ilm_power_design.md)
takes the assumptions instead – the design, the expected cell means, the
residual SD and the ICC – and
[`ilm_scaffold()`](https://huttoncp.github.io/illume/reference/ilm_scaffold.md)
returns the assumed model on its own, so what it implies can be checked
before a sample size is quoted from it.

``` r

d <- list(arm = c("control", "treatment"), time = c("pre", "post"))
pw <- ilm_power_design(change ~ arm * time + (1 | id), design = d,
                       within = "time", n_unit = c(60, 120, 240, 400),
                       cells = c(control.pre = 20, control.post = 21,
                                 treatment.pre = 20, treatment.post = 25),
                       sd = 6, icc = 0.5, term = "arm:time", sims = 300)
ilm_power_n(pw, target = 0.8)
#>  effect   n  n_lower n_upper n_unit n_unit_lower n_unit_upper
#>       4 330 292.5298 365.872    165     146.2649      182.936
```

`n_unit` is participants and `n` is rows: 165 people measured twice.

See
[`?ilm_power`](https://huttoncp.github.io/illume/reference/ilm_power.md),
[`?ilm_scaffold`](https://huttoncp.github.io/illume/reference/ilm_scaffold.md)
and the *Effect size and power* vignette.

------------------------------------------------------------------------

## 2. Exploration

Stages 2 to 4 use `illumex`, the exploration package `illume` attaches,
so nothing extra needs loading. Descriptive statistics, plots, and
bootstrap intervals for the quantities you are about to model.
[`ilm_describe_all()`](https://huttoncp.github.io/illumex/reference/ilm_describe_all.html)
summarises every column at once;
[`ilm_plot_all()`](https://huttoncp.github.io/illumex/reference/ilm_plot_all.html)
picks a sensible plot per variable;
[`ilm_boot_ci()`](https://huttoncp.github.io/illumex/reference/ilm_boot_ci.html)
and
[`ilm_boot_diff()`](https://huttoncp.github.io/illumex/reference/ilm_boot_diff.html)
put intervals on means, medians and differences between groups without
assuming normality.

``` r

ilm_describe_all(d)
ilm_plot_all(d)
ilm_boot_diff(score ~ arm, data = d)        # every pair, simultaneous by default
```

See the *Exploring data* vignette.

------------------------------------------------------------------------

## 3. Structure and profiling

When the columns are many and correlated,
[`ilm_reduce()`](https://huttoncp.github.io/illumex/reference/ilm_reduce.html)
finds the few directions they share – PCA, MCA or FAMD depending on the
types, or a generalized low rank model with a loss per column type –
[`ilm_cluster()`](https://huttoncp.github.io/illumex/reference/ilm_cluster.html)
groups the rows in that space, and
[`ilm_profile()`](https://huttoncp.github.io/illumex/reference/ilm_profile.html)
does both and describes what distinguishes each group.

``` r

pr <- ilm_profile(d, ndim = 5)               # closed-form FAMD
pr <- ilm_profile(d, ndim = 5, method = "glrm")   # a loss per column type
ilm_plot_profile(pr)
```

See the *Profiling* vignette and
[`?ilm_glrm`](https://huttoncp.github.io/illumex/reference/ilm_glrm.html).

------------------------------------------------------------------------

## 4. Outliers and anomalies

Two different questions.
[`ilm_outliers()`](https://huttoncp.github.io/illumex/reference/ilm_outliers.html)
asks whether a value is extreme in its own column.
[`ilm_anomaly()`](https://huttoncp.github.io/illumex/reference/ilm_anomaly.html)
asks whether a **row** is a plausible combination: someone 150 cm tall
is unremarkable, someone weighing 110 kg is unremarkable, and someone
who is both is not, and nothing in either column’s distribution says so.

``` r

ilm_outliers_all(d)                          # one column at a time
ilm_anomaly(d)                               # the combination
```

[`ilm_anomaly()`](https://huttoncp.github.io/illumex/reference/ilm_anomaly.html)
also names the column driving each flagged row, because an anomaly
nobody can explain is not actionable.

See the *Anomaly detection* vignette.

------------------------------------------------------------------------

## 5. Missing data

[`ilm_check_missing()`](https://huttoncp.github.io/illumex/reference/ilm_check_missing.html)
asks the question that matters, which is not how much is missing but
whether it is missing in a way that biases the answer: complete cases
stay unbiased when missingness is unrelated to the **outcome given the
predictors**, which is weaker than missing-completely-at-random and is
testable. When it is not,
[`ilm_impute()`](https://huttoncp.github.io/illume/reference/ilm_impute.md)
is the remedy.

``` r

ilm_check_missing(d, y ~ x + z)
im <- ilm_impute(d, m = 20)                  # chained equations by default
ilm_mi_pool(im, y ~ x + z)                   # fit each; Rubin's rules
```

`method = "glrm"` imputes categorical columns as categories;
`method = "lowrank"` handles data wider than it is long.

See the *Missing data* vignette.

------------------------------------------------------------------------

## 6. A causal model, if you have one

Optional, and the stage that decides whether the next one can be read
causally.
[`ilm_dag()`](https://huttoncp.github.io/illume/reference/ilm_dag.md)
records what you believe causes what;
[`ilm_adjust_sets()`](https://huttoncp.github.io/illume/reference/ilm_adjust_sets.md)
finds the minimal sets of variables that close the back-door paths;
[`ilm_dag_test()`](https://huttoncp.github.io/illume/reference/ilm_dag_test.md)
checks the conditional independencies the graph implies against the
data.

``` r

g <- ilm_dag("treat -> outcome; conf -> treat; conf -> outcome")
ilm_adjust_sets(g, exposure = "treat", outcome = "outcome")
ilm_dag_test(g, data = d)                    # does the data agree with it?
fit <- ilm_dag_model(g, exposure = "treat", outcome = "outcome", data = d)
```

[`ilm_dag_model()`](https://huttoncp.github.io/illume/reference/ilm_dag_model.md)
fits the model the adjustment set implies, so the mean structure is not
chosen by looking at p-values.

See the *Causal models* vignette.

------------------------------------------------------------------------

## 7. The model

One function. The formula is lme4’s, the families run from gaussian
through counts, proportions, ordered and unordered categories, and
survival times, and the same object comes back whatever you fit.

``` r

ilm_model(y ~ x + (1 | subj), data = d, family = "gaussian")
ilm_model(y ~ x + (1 | subj), data = d, family = "multinomial")
ilm_model(y ~ x, data = d, family = "ordinal")
ilm_model(y ~ x, data = d, family = "poisson", ziformula = ~ z)
ilm_model(y ~ x, data = d, family = "beta")
ilm_model(y ~ x + t2(time), data = d, family = "gaussian", ar = ilm_ar1(subj, time))
ilm_iv(y ~ x + w | z + w, data = d)          # an instrument instead
```

Designs the model cannot absorb get their own function:
[`ilm_did()`](https://huttoncp.github.io/illume/reference/ilm_did.md)
for difference-in-differences,
[`ilm_rdd()`](https://huttoncp.github.io/illume/reference/ilm_rdd.md)
for a regression discontinuity,
[`ilm_design()`](https://huttoncp.github.io/illume/reference/ilm_design.md)
for a complex sample.

See the *Regression models* vignette.

------------------------------------------------------------------------

## 8. Diagnostics

Every check says whether an assumption holds **and names a remedy that
exists in this package**. That is the rule the diagnostics layer is
built around; a check with no remedy is a complaint.

``` r

ilm_appraise(fit)                            # the whole battery
ilm_check_variance(fit)                      # -> dispformula =
ilm_check_zeros(fit)                         # -> ziformula =
ilm_check_proportional(fit)                  # -> family = "multinomial"
ilm_check_ar(fit)                            # -> ar = ilm_ar1(...)
ilm_variogram(fit)                           # spatial or temporal correlation
ilm_rqr_test(fit)                            # residuals, calibrated by simulation
```

Most remedies are a change to one argument of
[`ilm_model()`](https://huttoncp.github.io/illume/reference/ilm_model.md),
and
[`ilm_remedies()`](https://huttoncp.github.io/illume/reference/ilm_remedies.md)
writes each one out as that change. Each carries a tier:

- `numerical`: the same model fitted harder;
- `structural`: a different random-effect or variance structure, with
  fixed effects that mean what they did;
- `estimand`: a change to what the model estimates, to be made by
  choice, not because a check asked.

[`ilm_apply_remedy()`](https://huttoncp.github.io/illume/reference/ilm_apply_remedy.md)
refits with the one you pick.

``` r

rem <- ilm_remedies(fit, zeros = ilm_check_zeros(fit))
rem                                          # tier, remedy, and the change
fit2 <- ilm_apply_remedy(fit, rem, 1)        # refit with remedy 1
fit2$remedy_log                              # how fit2 was reached
```

If the variance structure is the worry rather than the mean structure,
[`ilm_robust()`](https://huttoncp.github.io/illume/reference/ilm_robust.md)
gives cluster-robust standard errors without needing to say what the
right structure would be.

------------------------------------------------------------------------

## 9. What the model says

[`ilm_anova()`](https://huttoncp.github.io/illume/reference/ilm_anova.md)
for which terms matter,
[`ilm_effects()`](https://huttoncp.github.io/illume/reference/ilm_effects.md)
for each on its family’s own scale,
[`ilm_emmeans()`](https://huttoncp.github.io/illume/reference/ilm_emmeans.md)
and
[`ilm_contrast()`](https://huttoncp.github.io/illume/reference/ilm_contrast.md)
for group means and the comparisons between them,
[`ilm_ame()`](https://huttoncp.github.io/illume/reference/ilm_ame.md)
for the average effect of a one-unit change, and
[`ilm_plot_model()`](https://huttoncp.github.io/illume/reference/ilm_plot_model.md)
for partial dependence.

``` r

ilm_anova(fit)                               # Type II by default
ilm_effects(fit)                             # odds ratios, rate ratios, ...
em <- ilm_emmeans(fit, "group")
ilm_contrast(em)                             # simultaneous by default
ilm_ame(fit, "x")                            # population-averaged
ilm_plot_model(fit)
```

A ratio from a mixed model is **conditional** on the cluster, and
[`ilm_effects()`](https://huttoncp.github.io/illume/reference/ilm_effects.md)
says so rather than letting it be read as a population quantity.

When the significant term is an interaction between a factor and a
*continuous* predictor – treatment by time, most often – a marginal mean
cannot take it apart, because averaging at the mean of the covariate
collapses the thing being asked about.
[`ilm_trends()`](https://huttoncp.github.io/illume/reference/ilm_trends.md)
estimates the slope within each level instead, and answers the two
questions separately: whether each slope differs from zero, and whether
the slopes differ from each other.

``` r

tr <- ilm_trends(fit, "group", var = "time")
tr                                           # is each group's slope non-zero?
ilm_contrast(tr)                             # do the groups' slopes differ?
```

For a mixed model the reference is a large-sample one unless you ask
otherwise.
[`ilm_denom_df()`](https://huttoncp.github.io/illume/reference/ilm_denom_df.md)
supplies Satterthwaite or Kenward-Roger denominator degrees of freedom,
which matter when the clusters are few.

If the design is a factorial or repeated-measures one and you want it
reported the way those designs are usually reported – F, mean squared
error, generalized eta squared, a sphericity correction –
[`ilm_aov_ez()`](https://huttoncp.github.io/illume/reference/ilm_aov_ez.md)
does the whole of stages 7 to 9 from a description of the columns. See
[`vignette("anova")`](https://huttoncp.github.io/illume/articles/anova.md).

------------------------------------------------------------------------

## 10. Scenario projection

Optional, and the stage a decision usually needs.
[`ilm_scenario()`](https://huttoncp.github.io/illume/reference/ilm_scenario.md)
sets predictors to values a stakeholder is considering and reports what
the model says would happen, with intervals and with differences between
scenarios.

``` r

ilm_scenario(fit, dose = c(0, 10, 20), contrast = "first")
```

The default averages over **the observed units**: set the intervention
for everyone, keep each unit’s own covariates, predict, average. That is
a different question from predicting for a unit at the average
covariate, and under any non-linear link it is a different answer. It
also checks whether each scenario sits inside the data – by value, and
by combination, since every value can be ordinary while the combination
is one nobody occupies.

------------------------------------------------------------------------

## 11. Reporting

[`ilm_interpret()`](https://huttoncp.github.io/illume/reference/ilm_interpret.md)
writes the findings as prose from deterministic templates, so the same
fit always produces the same words. Causal language appears only when a
DAG has licensed it.

``` r

ilm_interpret(fit)
```

------------------------------------------------------------------------

## What to read next

| Vignette | Covers |
|----|----|
| *Introduction to illume* | orientation and one worked example |
| *Regression models* | [`ilm_model()`](https://huttoncp.github.io/illume/reference/ilm_model.md) in depth, every family, mixed multinomial |
| *Exploring data* (illumex) | descriptives, plots, bootstrap intervals |
| *Profiling* (illumex) | [`ilm_reduce()`](https://huttoncp.github.io/illumex/reference/ilm_reduce.html), [`ilm_cluster()`](https://huttoncp.github.io/illumex/reference/ilm_cluster.html), [`ilm_profile()`](https://huttoncp.github.io/illumex/reference/ilm_profile.html) |
| *Anomaly detection* (illumex) | [`ilm_anomaly()`](https://huttoncp.github.io/illumex/reference/ilm_anomaly.html) against [`ilm_outliers()`](https://huttoncp.github.io/illumex/reference/ilm_outliers.html); [`ilm_plot_anomaly()`](https://huttoncp.github.io/illumex/reference/ilm_plot_anomaly.html) to see whether the flagged rows stand apart |
| *Missing data* | diagnosing it, imputing it, pooling |
| *Causal models* | DAGs, difference-in-differences, discontinuities |
| *Effect size and power* | [`ilm_effects()`](https://huttoncp.github.io/illume/reference/ilm_effects.md), [`ilm_power()`](https://huttoncp.github.io/illume/reference/ilm_power.md), [`ilm_power_design()`](https://huttoncp.github.io/illume/reference/ilm_power_design.md), [`ilm_scaffold()`](https://huttoncp.github.io/illume/reference/ilm_scaffold.md), [`ilm_scenario()`](https://huttoncp.github.io/illume/reference/ilm_scenario.md) |
| *Moderation* | [`ilm_moderation()`](https://huttoncp.github.io/illume/reference/ilm_moderation.md), [`ilm_plot_moderation()`](https://huttoncp.github.io/illume/reference/ilm_plot_moderation.md) |
