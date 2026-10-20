# Introduction to illume

## What this package is for

`illume` is a frequentist analysis workflow aimed at **inference**
rather than prediction. It covers the path from a sample-size
calculation, through exploration and missing data, to a model, its
diagnostics, and a scenario a stakeholder can act on – with one object
and one vocabulary the whole way.

Two commitments shape it.

**One engine fits everything.**
[`ilm_model()`](https://huttoncp.github.io/illume/reference/ilm_model.md)
goes from a simple regression to a mixed model with smooths and
correlated errors, across families from gaussian through counts,
proportions, ordered and unordered categories and survival times.
Because a single engine fits them all,
[`summary()`](https://rdrr.io/r/base/summary.html),
[`ilm_anova()`](https://huttoncp.github.io/illume/reference/ilm_anova.md),
[`predict()`](https://rdrr.io/r/stats/predict.html) and every diagnostic
stay the same as the model gets harder.

**Every check names a remedy that exists in this package.** A diagnostic
that tells you an assumption fails and leaves you to find the fix
elsewhere is a complaint.
[`ilm_check_variance()`](https://huttoncp.github.io/illume/reference/ilm_check_variance.md)
names `dispformula`;
[`ilm_check_zeros()`](https://huttoncp.github.io/illume/reference/ilm_check_zeros.md)
names `ziformula`;
[`ilm_check_proportional()`](https://huttoncp.github.io/illume/reference/ilm_check_proportional.md)
names `family = "multinomial"`. Where a remedy did not exist, it was
built – which is most of what this package is. The remedies are code as
well as words:
[`ilm_remedies()`](https://huttoncp.github.io/illume/reference/ilm_remedies.md)
writes each one out as the change to
[`ilm_model()`](https://huttoncp.github.io/illume/reference/ilm_model.md)
that makes it, with a tier saying what applying it would change, and
[`ilm_apply_remedy()`](https://huttoncp.github.io/illume/reference/ilm_apply_remedy.md)
refits with the one you choose.

``` r

library(illume)
```

## A worked example, end to end

``` r

d <- ilm_sim()                       # a grouped dataset with a known structure

## what is in it
ilm_describe_all(d)
ilm_check_missing(d, downtime ~ income + region)

## a model: a count outcome, repeated within worker
fit <- ilm_model(downtime ~ income + region + (1 | id), data = d,
                 family = "poisson")

## does it hold up, and if not, what is the fix?
ilm_appraise(fit)                    # the whole battery at once
zz  <- ilm_check_zeros(fit)          # more zeros than poisson allows?
rem <- ilm_remedies(fit, zeros = zz) # each remedy, as the change to make
fit <- ilm_apply_remedy(fit, rem, 1) # refit with remedy 1, or your pick

## what does it say?
ilm_anova(fit)
ilm_effects(fit)                     # incidence rate ratios
em <- ilm_emmeans(fit, "region")
ilm_contrast(em)                     # every pair, simultaneously

## what would happen if we changed something?
ilm_scenario(fit, income = c(20, 40, 60), contrast = "first")

## say it in words
ilm_interpret(fit)
```

If
[`ilm_check_zeros()`](https://huttoncp.github.io/illume/reference/ilm_check_zeros.md)
flags excess zeros,
[`ilm_remedies()`](https://huttoncp.github.io/illume/reference/ilm_remedies.md)
names the fix – `ziformula = ~ 1`, as a mixture or as a hurdle –
[`ilm_apply_remedy()`](https://huttoncp.github.io/illume/reference/ilm_apply_remedy.md)
makes it, and everything downstream is unchanged. That is the shape of
the whole package.

## What is in it

|  |  |
|----|----|
| **Models** | [`ilm_model()`](https://huttoncp.github.io/illume/reference/ilm_model.md) – gaussian, binomial, poisson, negative binomial, beta, multinomial, three ordinal links, three accelerated failure time families, Royston-Parmar survival; random effects, smooths, AR(1)/CAR(1) and random walks over time, dispersion models, zero-inflation and hurdles |
| **Other designs** | [`ilm_iv()`](https://huttoncp.github.io/illume/reference/ilm_iv.md) instrumental variables, [`ilm_did()`](https://huttoncp.github.io/illume/reference/ilm_did.md) difference in differences, [`ilm_rdd()`](https://huttoncp.github.io/illume/reference/ilm_rdd.md) regression discontinuity, [`ilm_design()`](https://huttoncp.github.io/illume/reference/ilm_design.md) complex samples |
| **Diagnostics** | [`ilm_appraise()`](https://huttoncp.github.io/illume/reference/ilm_appraise.md) and around twenty individual checks, each naming its remedy; [`ilm_remedies()`](https://huttoncp.github.io/illume/reference/ilm_remedies.md) writes the remedies out as code and [`ilm_apply_remedy()`](https://huttoncp.github.io/illume/reference/ilm_apply_remedy.md) makes one |
| **Inference** | [`ilm_anova()`](https://huttoncp.github.io/illume/reference/ilm_anova.md), [`ilm_effects()`](https://huttoncp.github.io/illume/reference/ilm_effects.md), [`ilm_emmeans()`](https://huttoncp.github.io/illume/reference/ilm_emmeans.md)/[`ilm_contrast()`](https://huttoncp.github.io/illume/reference/ilm_contrast.md), [`ilm_trends()`](https://huttoncp.github.io/illume/reference/ilm_trends.md), [`ilm_moderation()`](https://huttoncp.github.io/illume/reference/ilm_moderation.md), [`ilm_ame()`](https://huttoncp.github.io/illume/reference/ilm_ame.md), [`ilm_robust()`](https://huttoncp.github.io/illume/reference/ilm_robust.md), [`ilm_denom_df()`](https://huttoncp.github.io/illume/reference/ilm_denom_df.md), [`ilm_pb_lrt()`](https://huttoncp.github.io/illume/reference/ilm_pb_lrt.md) |
| **ANOVA** | [`ilm_aov_ez()`](https://huttoncp.github.io/illume/reference/ilm_aov_ez.md) – factorial and repeated measures by naming columns, with generalized eta squared and sphericity corrections |
| **Exploration** (illumex) | [`ilm_describe_all()`](https://huttoncp.github.io/illumex/reference/ilm_describe_all.html), the data plots, [`ilm_boot_ci()`](https://huttoncp.github.io/illumex/reference/ilm_boot_ci.html), [`ilm_outliers()`](https://huttoncp.github.io/illumex/reference/ilm_outliers.html), [`ilm_anomaly()`](https://huttoncp.github.io/illumex/reference/ilm_anomaly.html) and [`ilm_plot_anomaly()`](https://huttoncp.github.io/illumex/reference/ilm_plot_anomaly.html) |
| **Structure** (illumex) | [`ilm_reduce()`](https://huttoncp.github.io/illumex/reference/ilm_reduce.html), [`ilm_cluster()`](https://huttoncp.github.io/illumex/reference/ilm_cluster.html), [`ilm_profile()`](https://huttoncp.github.io/illumex/reference/ilm_profile.html), [`ilm_glrm()`](https://huttoncp.github.io/illumex/reference/ilm_glrm.html) |
| **Missing data** | [`ilm_check_missing()`](https://huttoncp.github.io/illumex/reference/ilm_check_missing.html) (illumex) to describe it; [`ilm_impute()`](https://huttoncp.github.io/illume/reference/ilm_impute.md) and [`ilm_mi_pool()`](https://huttoncp.github.io/illume/reference/ilm_mi_pool.md) to fill it in and pool |
| **Causal** | [`ilm_dag()`](https://huttoncp.github.io/illume/reference/ilm_dag.md), [`ilm_adjust_sets()`](https://huttoncp.github.io/illume/reference/ilm_adjust_sets.md), [`ilm_dag_test()`](https://huttoncp.github.io/illume/reference/ilm_dag_test.md), [`ilm_dag_model()`](https://huttoncp.github.io/illume/reference/ilm_dag_model.md), [`ilm_mediate()`](https://huttoncp.github.io/illume/reference/ilm_mediate.md) |
| **Design and decision** | [`ilm_power()`](https://huttoncp.github.io/illume/reference/ilm_power.md), [`ilm_scaffold()`](https://huttoncp.github.io/illume/reference/ilm_scaffold.md), [`ilm_power_design()`](https://huttoncp.github.io/illume/reference/ilm_power_design.md), [`ilm_scenario()`](https://huttoncp.github.io/illume/reference/ilm_scenario.md), [`ilm_interpret()`](https://huttoncp.github.io/illume/reference/ilm_interpret.md) |

## The multinomial claim, stated carefully

If an outcome has three or more categories with no natural order – a
diagnosis, a habitat type, a choice among options – and the data are
grouped, the frequentist options are thin.
[`nnet::multinom()`](https://rdrr.io/pkg/nnet/man/multinom.html) has no
random effects;
[`lme4::glmer()`](https://rdrr.io/pkg/lme4/man/glmer.html) and `glmmTMB`
have no multinomial family; `brms` is Bayesian and much slower.

``` r

ilm_model(choice ~ price + (1 | household), data = dd, family = "multinomial")
```

**They are not absent, though.** `mclogit::mblogit()` fits this model,
by penalised quasi-likelihood rather than by maximising a
Laplace-approximated marginal likelihood. The difference is measurable,
and it is not uniform:

                           illume    mclogit
      40 clusters x 25
        coverage              0.953     0.949      <- no meaningful difference
        attenuation           1.005     0.960
      100 clusters x 4
        coverage              0.949     0.891
        attenuation           1.004     0.731      <- 27% shrunk toward zero
        bias                  0.035     0.188
        RMSE                  0.409     0.336      <- mclogit WINS here

With well-populated clusters the two agree and there is nothing to
choose between them. PQL degrades where the clusters are small, because
it shrinks the fixed effects toward zero – and its standard errors do
not know it, so the interval is centred in the wrong place at close to
the right width.

The last row is the honest caveat: mclogit’s **RMSE is lower**, because
shrinkage buys a variance reduction that more than pays for the bias. If
the goal is prediction that is a defensible trade. If the goal is a
coefficient with an interval you intend to interpret, it is not, and
coverage is the measure that says so. The full table is in
`studies/findings/mclogit.md` in the source repository, which is not
shipped with the installed package:
<https://github.com/huttoncp/illume/blob/main/studies/findings/mclogit.md>.

## Working with other packages

`illume` fits into the usual ecosystem.
[`car::Anova()`](https://rdrr.io/pkg/car/man/Anova.html) dispatches to
the correct joint test,
[`performance::model_performance()`](https://easystats.github.io/performance/reference/model_performance.html)
works,
[`performance::check_model()`](https://easystats.github.io/performance/reference/check_model.html)
draws the panels of
[`ilm_appraise()`](https://huttoncp.github.io/illume/reference/ilm_appraise.md),
and `marginaleffects` and the easystats stack can each be enabled with
one call:

``` r

ilm_register_marginaleffects()
marginaleffects::avg_slopes(fit, variables = "x1")

ilm_register_insight()          # parameters, performance, report
parameters::model_parameters(fit)
```

[`ilm_ame()`](https://huttoncp.github.io/illume/reference/ilm_ame.md)
agrees with
[`marginaleffects::avg_slopes()`](https://rdrr.io/pkg/marginaleffects/man/slopes.html)
to 1e-6 on the estimates and 1e-7 on the standard errors, which is the
point of having both: an independent implementation is the only thing
that has ever found a real bug here.

## Where to read more

[`vignette("workflow")`](https://huttoncp.github.io/illume/articles/workflow.md)
is the map: eleven stages from power to reporting, with a pointer for
each.

| Vignette | Covers |
|----|----|
| `workflow` | the whole path, stage by stage |
| `benchmarking` | what was measured, against what, and where illume comes off worse |
| `regression-models` | [`ilm_model()`](https://huttoncp.github.io/illume/reference/ilm_model.md) in depth, every family, mixed multinomial |
| `exploring-data` (illumex) | descriptives, plots, bootstrap intervals |
| `profiling` (illumex) | [`ilm_reduce()`](https://huttoncp.github.io/illumex/reference/ilm_reduce.html), [`ilm_cluster()`](https://huttoncp.github.io/illumex/reference/ilm_cluster.html), [`ilm_profile()`](https://huttoncp.github.io/illumex/reference/ilm_profile.html) |
| `anomaly-detection` (illumex) | rows that are implausible as combinations |
| `missing-data` | diagnosing it, imputing it, pooling |
| `causal-models` | DAGs, difference in differences, discontinuities |
| `effect-size-and-power` | effect sizes, power with or without a pilot, scenario projection |
| `moderation` | whether an effect holds for everyone, without p-hacking the subgroups |
| `anova` | factorial and repeated-measures ANOVA, in the language those designs are taught in |

The three marked (illumex) belong to `illumex`, the exploration package
`illume` attaches:
[`vignette("exploring-data", package = "illumex")`](https://huttoncp.github.io/illumex/articles/exploring-data.html).

## How this package was built

`illume` was developed in collaboration with [Claude Opus
5.5](https://www.anthropic.com) under a human-in-the-loop model: the
model prototyped and drafted the implementation, the documentation and
the simulation studies; the author specified what the package should be,
set the standing constraints, made the decisions the design turned on,
and rejected proposals. The author’s independent review and testing is
*in progress at the time of writing* and will be complete before any
release.

It is disclosed because the provenance of a statistical tool is a
reasonable thing for its users to know, and because venues increasingly
require it – the Journal of Open Source Software has mandated an AI
usage disclosure since January 2026.

But the disclosure is less load-bearing than it looks. The question that
matters about a statistical package is not who typed it; it is whether
the claims are checkable and whether anyone checked them. Every validity
claim here is a comparison against an independent implementation
somebody else wrote, and the scripts and retained runs are in
`studies/`. Every real defect this package has had was found that way –
none by reading the code, none by the test suite. That is the standard
it asks to be judged by, and it is the same standard whoever wrote it.
The README has the longer version.

## Standing on other people’s shoulders

Almost nothing here is new mathematics. What this package mostly does is
put existing methods behind one object with one vocabulary, and it can
only do that because other people wrote the hard parts first and wrote
them in the open.

`RTMB` and `TMB` are the engine – every model here is an
automatic-differentiation tape and a Laplace approximation, and the
multinomial mixed model is feasible at all only because Kristensen and
colleagues made writing one a matter of declaring what to integrate out.
`Matrix` carries the sparse algebra, `collapse` makes the exploratory
half fast enough to use, `tinyplot` draws everything, and `mgcv`
supplies the smooth bases.

Every validity claim in this package is a comparison against somebody
else’s independent implementation – `lme4`, `glmmTMB`, `emmeans`,
`afex`, `car`, `lmerTest`, `pbkrtest`, `survey`, `sandwich`, `ordinal`,
`pscl`, `mice`, `nlme`, `survival`, `mclogit`, `brms`. Every real defect
this package has had was found that way and by nothing else.

The exploratory half is a direct descendant of
[`elucidate`](https://github.com/bcgov/elucidate), by the same author:
`describe()`, `counts()`, `dupes()`, `wash_df()` and the `plot_*()`
family are where
[`ilm_describe()`](https://huttoncp.github.io/illumex/reference/ilm_describe.html),
[`ilm_counts()`](https://huttoncp.github.io/illumex/reference/ilm_counts.html),
[`ilm_dupes()`](https://huttoncp.github.io/illumex/reference/ilm_dupes.html)
and the rest came from, along with the idea that the routine work before
a model should be one call with a consistent interface.

And some of it was read straight out of their source: the sphericity
arithmetic from `car`, the generalized eta squared formula from `afex`,
the way REML is obtained from `glmmTMB`. Where a method has a name
attached to it in these pages, that is deliberate – it should be
possible to find out whose idea any part of this was. The README has the
longer version.

## References

The multinomial model is a baseline-category logit with random effects;
Agresti (2013), *Categorical Data Analysis*, chapter 8, is the standard
reference. The fitting method is described in Kristensen et al. (2016),
*Journal of Statistical Software* 70(5). The mixed-model representation
of smooths follows Wood (2017), *Generalized Additive Models*, section
5.4. The parametric bootstrap follows Halekoh and Højsgaard (2014),
*Journal of Statistical Software* 59(9).
