# illume

**A Unified Engine for Exploration and Frequentist Inference.**

`illume` is an analysis workflow for R aimed at **inference** rather
than prediction. It covers the path from a sample-size calculation,
through exploration and missing data, to a model, its diagnostics, and a
scenario a stakeholder can act on – with one object and one vocabulary
the whole way.

The exploratory half – describing, cleaning, plotting and profiling data
before a model – is the companion package
[`illumex`](https://github.com/huttoncp/illumex). `illume` attaches it,
so [`library(illume)`](https://github.com/huttoncp/illume) gives you
both.

> Not on CRAN yet. This is pre-release software under active
> development.

## Two commitments

**One engine fits everything.**
[`ilm_model()`](https://huttoncp.github.io/illume/reference/ilm_model.md)
goes from a simple regression to a mixed model with smooths and
correlated errors, across families from gaussian through counts,
proportions, ordered and unordered categories, and survival times. All
of them are fitted by approximate maximum likelihood using the Laplace
approximation via [RTMB](https://github.com/kaskr/RTMB). Because a
single engine fits them all,
[`summary()`](https://rdrr.io/r/base/summary.html),
[`ilm_anova()`](https://huttoncp.github.io/illume/reference/ilm_anova.md),
[`predict()`](https://rdrr.io/r/stats/predict.html) and every diagnostic
stay the same as the model gets harder. Leave `family` off and it is
read off the response, with the choice and the evidence for it stated;
what the response cannot settle is asked about rather than guessed.

**Every check names a remedy that exists in this package.** A diagnostic
that tells you an assumption fails and leaves you to find the fix
elsewhere is a complaint, not a tool.
[`ilm_check_variance()`](https://huttoncp.github.io/illume/reference/ilm_check_variance.md)
names `dispformula`;
[`ilm_check_zeros()`](https://huttoncp.github.io/illume/reference/ilm_check_zeros.md)
names `ziformula`;
[`ilm_check_proportional()`](https://huttoncp.github.io/illume/reference/ilm_check_proportional.md)
names `family = "multinomial"`. Where a remedy did not exist, it was
built – which is most of what this package is.

## Installation

``` r

# install.packages("remotes")
remotes::install_github("huttoncp/illume")
```

This installs `illumex` from GitHub along with it. R \>= 4.1. The only
compiled dependency is RTMB/TMB, which is on CRAN.

## A worked example, end to end

``` r

library(illume)

d <- ilm_sim()                       # a grouped dataset with a known structure

## what is in it
ilm_describe_all(d)
ilm_check_missing(d, downtime ~ income + region)

## a model: a count outcome, repeated within worker
fit <- ilm_model(downtime ~ income + region + (1 | id), data = d,
                 family = "poisson")

## does it hold up?
ilm_appraise(fit)                    # the whole battery at once
ilm_check_zeros(fit)                 # more zeros than poisson allows?

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
flags excess zeros, the next line is `ilm_model(..., ziformula = ~ 1)`
and everything downstream is unchanged. That is the shape of the whole
package.

## What is in it

|  |  |
|----|----|
| **Models** | [`ilm_model()`](https://huttoncp.github.io/illume/reference/ilm_model.md) – gaussian, binomial, poisson, negative binomial, beta, multinomial, three ordinal links, three accelerated failure time families, Royston-Parmar survival; random intercepts and slopes, penalised smooths, AR(1)/CAR(1), dispersion models, zero-inflation and hurdles |
| **Other designs** | [`ilm_iv()`](https://huttoncp.github.io/illume/reference/ilm_iv.md) instrumental variables, [`ilm_did()`](https://huttoncp.github.io/illume/reference/ilm_did.md) difference in differences, [`ilm_rdd()`](https://huttoncp.github.io/illume/reference/ilm_rdd.md) regression discontinuity, [`ilm_design()`](https://huttoncp.github.io/illume/reference/ilm_design.md) complex samples |
| **Diagnostics** | [`ilm_appraise()`](https://huttoncp.github.io/illume/reference/ilm_appraise.md) and around twenty individual checks, each naming its remedy |
| **Inference** | [`ilm_anova()`](https://huttoncp.github.io/illume/reference/ilm_anova.md), [`ilm_effects()`](https://huttoncp.github.io/illume/reference/ilm_effects.md), [`ilm_emmeans()`](https://huttoncp.github.io/illume/reference/ilm_emmeans.md)/[`ilm_contrast()`](https://huttoncp.github.io/illume/reference/ilm_contrast.md), [`ilm_trends()`](https://huttoncp.github.io/illume/reference/ilm_trends.md), [`ilm_ame()`](https://huttoncp.github.io/illume/reference/ilm_ame.md), [`ilm_robust()`](https://huttoncp.github.io/illume/reference/ilm_robust.md), [`ilm_denom_df()`](https://huttoncp.github.io/illume/reference/ilm_denom_df.md), [`ilm_pb_lrt()`](https://huttoncp.github.io/illume/reference/ilm_pb_lrt.md) |
| **ANOVA** | [`ilm_aov_ez()`](https://huttoncp.github.io/illume/reference/ilm_aov_ez.md) – factorial and repeated measures by naming columns; F, generalized eta squared, sphericity corrections |
| **Exploration** (illumex) | [`ilm_describe_all()`](https://rdrr.io/pkg/illumex/man/ilm_describe_all.html), the data plots, [`ilm_boot_ci()`](https://rdrr.io/pkg/illumex/man/ilm_boot_ci.html), [`ilm_outliers()`](https://rdrr.io/pkg/illumex/man/ilm_outliers.html), [`ilm_anomaly()`](https://rdrr.io/pkg/illumex/man/ilm_anomaly.html) and [`ilm_plot_anomaly()`](https://rdrr.io/pkg/illumex/man/ilm_plot_anomaly.html) |
| **Structure** (illumex) | [`ilm_reduce()`](https://rdrr.io/pkg/illumex/man/ilm_reduce.html), [`ilm_cluster()`](https://rdrr.io/pkg/illumex/man/ilm_cluster.html), [`ilm_profile()`](https://rdrr.io/pkg/illumex/man/ilm_profile.html), [`ilm_glrm()`](https://rdrr.io/pkg/illumex/man/ilm_glrm.html) |
| **Missing data** | [`ilm_check_missing()`](https://rdrr.io/pkg/illumex/man/ilm_check_missing.html) (illumex) to describe it; [`ilm_impute()`](https://huttoncp.github.io/illume/reference/ilm_impute.md) and [`ilm_mi_pool()`](https://huttoncp.github.io/illume/reference/ilm_mi_pool.md) to fill it in and pool |
| **Causal** | [`ilm_dag()`](https://huttoncp.github.io/illume/reference/ilm_dag.md), [`ilm_adjust_sets()`](https://huttoncp.github.io/illume/reference/ilm_adjust_sets.md), [`ilm_dag_test()`](https://huttoncp.github.io/illume/reference/ilm_dag_test.md), [`ilm_dag_model()`](https://huttoncp.github.io/illume/reference/ilm_dag_model.md), [`ilm_mediate()`](https://huttoncp.github.io/illume/reference/ilm_mediate.md) |
| **Design and decision** | [`ilm_power()`](https://huttoncp.github.io/illume/reference/ilm_power.md) from a fit, [`ilm_power_design()`](https://huttoncp.github.io/illume/reference/ilm_power_design.md) and [`ilm_scaffold()`](https://huttoncp.github.io/illume/reference/ilm_scaffold.md) from assumptions alone – for every family, multinomial included, each simulated study analysed with the test the analysis will report; [`ilm_scenario()`](https://huttoncp.github.io/illume/reference/ilm_scenario.md), [`ilm_interpret()`](https://huttoncp.github.io/illume/reference/ilm_interpret.md) |

## The multinomial claim, stated carefully

If an outcome has three or more categories with no natural order and the
data are grouped, the frequentist options are thin.
[`nnet::multinom()`](https://rdrr.io/pkg/nnet/man/multinom.html) has no
random effects;
[`lme4::glmer()`](https://rdrr.io/pkg/lme4/man/glmer.html) and `glmmTMB`
have no multinomial family; `brms` is Bayesian and much slower.

``` r

ilm_model(choice ~ price + (1 | household), data = dd, family = "multinomial")
```

They are not absent, though. `mclogit::mblogit()` fits this model by
penalised quasi-likelihood. At nominal 0.95:

                           illume    mclogit
      40 clusters x 25
        coverage              0.953     0.949      <- no meaningful difference
      100 clusters x 4
        coverage              0.949     0.891
        attenuation           1.004     0.731      <- 27% shrunk toward zero
        RMSE                  0.409     0.336      <- mclogit wins here

With well-populated clusters the two agree. PQL degrades where the
clusters are small, and its standard errors do not know it – the
interval is centred in the wrong place at close to the right width. The
last row is the caveat that travels with the result: mclogit’s RMSE is
*lower*, because shrinkage buys a variance reduction that more than pays
for the bias. For prediction that is a defensible trade; for a
coefficient you intend to interpret it is not.

On data that misbehave – a sparse category, a variance at its boundary,
lopsided clusters, random effects that are not normal – the comparison
is less tidy, and not all of it favours illume.
[`vignette("benchmarking")`](https://huttoncp.github.io/illume/articles/benchmarking.md)
has the whole account, the results against illume included.

## How it is validated

Every claim in the documentation is a measurement. The
[`studies/`](https://github.com/huttoncp/illume/tree/main/studies)
directory of this repository (not shipped with the installed package)
holds the scripts, the retained runs, and a generated findings log per
study:

| Study | Against | Reports |
|----|----|----|
| `coverage` | 32 model configurations | interval coverage per family, structure and design |
| *(validated against)* | `afex`, `lmerTest`, `pbkrtest`, `emmeans` | ANOVA tables, degrees of freedom, marginal slopes |
| `power` | closed form where one exists | 13 designs |
| `bench` | `lme4`, `glmmTMB`, `nlme`, `survreg` | agreement and timing |
| `mclogit` | `mclogit::mblogit` | the table above |
| `messy` | `mclogit::mblogit`, on data that misbehave | coverage, bias, RMSE and convergence in six regimes |
| `brms` | `brms` | the Bayesian comparison |
| `imputation` | 5 methods over 6 designs | reconstruction *and* downstream coverage |

[`vignette("benchmarking")`](https://huttoncp.github.io/illume/articles/benchmarking.md)
collects them, with the results that go against illume stated as plainly
as the ones that do not.

Individual features are checked against outside implementations rather
than against expectations: `emmeans`, `pscl`, `ordinal`, `sandwich`,
`survey`, `glmmTMB`, closed-form 2SLS, and the SVD. **Every real defect
found in this package came from one of those comparisons, or from
widening a simulation. None came from the test suite** – which is worth
knowing about test suites.

## Working with other packages

[`car::Anova()`](https://rdrr.io/pkg/car/man/Anova.html) dispatches to
the correct joint test, and
[`performance::model_performance()`](https://easystats.github.io/performance/reference/model_performance.html)
and `check_model()` work directly. `marginaleffects` and the easystats
stack each need one call:

``` r

ilm_register_marginaleffects()
marginaleffects::avg_slopes(fit, variables = "x1")

ilm_register_insight()
parameters::model_parameters(fit)
```

## Documentation

| Vignette | Covers |
|----|----|
| `workflow` | the whole path, eleven stages |
| `benchmarking` | what was measured, against what, and where illume comes off worse |
| `regression-models` | [`ilm_model()`](https://huttoncp.github.io/illume/reference/ilm_model.md) in depth, every family, mixed multinomial |
| `missing-data` | diagnosing it, imputing it, pooling |
| `causal-models` | DAGs, difference in differences, discontinuities |
| `effect-size-and-power` | effect sizes, power, scenario projection |
| `moderation` | whether an effect holds for everyone, without p-hacking the subgroups |
| `anova` | factorial and repeated-measures ANOVA, for a psychology audience |

``` r

vignette("workflow", package = "illume")
```

Exploring data, profiling and anomaly detection are `illumex`’s
vignettes: `vignette("exploring-data", package = "illumex")`.

## How this package was built

`illume` was developed in collaboration with [Claude Opus
5.5](https://www.anthropic.com), run at maximum reasoning effort, under
a human-in-the-loop model. The division of labour was consistent
throughout and is worth stating plainly rather than leaving to be
inferred.

**What the model did.** Prototyped and drafted the implementation, the
documentation and the simulation studies. Proposed designs, and argued
for them when it thought a choice was wrong.

**What the author did.** Specified what the package should be and what
it should refuse to do. Set the standing constraints – no new hard
dependencies, every diagnostic names a remedy that exists in the
package, no claim without a measurement behind it. Made the decisions
the design turned on: which methods to include, which defaults to set,
what to do when two defensible options existed. Directed the ordering of
the work, and rejected proposals.

**Status of independent review.** The author is reviewing the package
and testing it independently of the development process, and that work
is *in progress at the time of writing*. It will be complete before any
release, and this paragraph will say so when it is. Until then the
package is pre-release and should be treated as such.

### Why this is disclosed rather than mentioned quietly

Partly because it is true and the provenance of a statistical tool is a
reasonable thing for its users to know. Partly because it is
increasingly required: the Journal of Open Source Software has mandated
an AI usage disclosure since January 2026, and other venues are moving
the same way.

But mostly because the disclosure is less load-bearing than it looks,
and saying so is the honest position. The question that matters about a
statistical package is not who typed it. It is whether its claims are
checkable and whether anyone checked them. Every validity claim here is
a comparison against an independent implementation that somebody else
wrote – `lme4`, `glmmTMB`, `emmeans`, `afex`, `car`, `lmerTest`,
`survey`, `sandwich`, `mice` and the rest – and the scripts, the
retained runs and the generated findings are in
[`studies/`](https://github.com/huttoncp/illume/tree/main/studies) for
anyone who wants to re-run them. Every real defect this package has had
was found that way. None was found by reading the code, and none by the
test suite.

That is the standard the package asks to be judged by, and it is the
same standard whoever wrote it. \## Standing on other people’s shoulders

Almost nothing here is new mathematics. What `illume` mostly does is put
existing methods behind one object with one vocabulary, and it can only
do that because other people wrote the hard parts first and wrote them
in the open. Three different debts are worth separating, because the
third is the one that usually goes unsaid.

**Without these there would be no package.** `RTMB` and `TMB` are the
engine: every model here is an automatic-differentiation tape and a
Laplace approximation, and the reason the multinomial mixed model is
feasible at all is that Kristensen, Nielsen, Berg, Skaug and Bell made
writing one a matter of declaring which parameters to integrate out.
`Matrix` carries the sparse algebra underneath it. `collapse` does the
grouped operations fast enough that the exploratory half is usable on
real data. `tinyplot` draws every plot. `mgcv` supplies the smooth basis
construction, and Simon Wood’s mixed-model representation of a penalised
smooth is what lets one engine fit both.

**Without these it could not be trusted.** Every validity claim in this
package is a comparison against an independent implementation, and each
of those is someone else’s careful work being used as a measuring stick:
`lme4` and `glmmTMB` for mixed models, `emmeans` for marginal means and
slopes, `afex` and `car` for ANOVA tables and sphericity, `lmerTest` and
`pbkrtest` for degrees of freedom, `survey` for design-based inference,
`sandwich` for robust covariance, `ordinal`, `pscl`, `mice`, `nlme`,
`survival`, `nnet`, `mclogit` and `brms` for their respective families
and methods. Every real defect this package has ever had was found by
one of those comparisons and by nothing else, which is a debt of a
fairly literal kind.

**And some of it was read straight out of their source.** The
Greenhouse-Geisser, Huynh-Feldt and Mauchly arithmetic in
[`ilm_aov_ez()`](https://huttoncp.github.io/illume/reference/ilm_aov_ez.md)
was taken from `car`, because there are several non-equivalent forms in
print and the one that matters is the one people compare against. The
generalized eta squared formula came from `afex` for the same reason.
The way REML is obtained – by adding the fixed effects to the block that
gets integrated out – was taken from `glmmTMB`, and so was the guard
that a sandwich cannot be computed from a REML fit. The parametric
bootstrap follows Halekoh and Højsgaard. The multiple-imputation-by-PCA
route follows Josse and Husson; the chained-equations default is the
approach `mice` established; the profiling workflow is a lighter-weight
rebuild of `FactoMineR`’s `HCPC()` and `catdes()`.

**And one predecessor.** The exploratory half, now the companion package
`illumex`, is a direct descendant of
[`elucidate`](https://github.com/bcgov/elucidate), written by the same
author for the BC Public Service. The lineage is visible in the function
names: `describe()`, `describe_all()`, `counts()`, `counts_tb()`,
`dupes()`, `copies()`, `wash_df()`, `recode_errors()`, `translate()` and
the whole `plot_*()` family became
[`ilm_describe()`](https://rdrr.io/pkg/illumex/man/ilm_describe.html),
[`ilm_counts()`](https://rdrr.io/pkg/illumex/man/ilm_counts.html),
[`ilm_dupes()`](https://rdrr.io/pkg/illumex/man/ilm_dupes.html) and the
rest, with the same idea behind them – that the routine work before a
model should be one call with a consistent interface, rather than six
lines of [`sapply()`](https://rdrr.io/r/base/lapply.html) reassembled
from memory every time. `illumex` reimplements that on a different
backend and extends it – elucidate is built on `data.table`, `dplyr` and
`ggplot2`, this on `collapse` and `tinyplot`, and the two share an
interface rather than an implementation – but the design is elucidate’s
and it would be poor form to pretend otherwise.

Where a method has a name attached to it in the documentation, that is
deliberate. It should be possible to find out whose idea any part of
this was. \## References

Kristensen, Nielsen, Berg, Skaug and Bell (2016). TMB: Automatic
Differentiation and Laplace Approximation. *Journal of Statistical
Software* 70(5).

Halekoh and Hojsgaard (2014). A Kenward-Roger Approximation and
Parametric Bootstrap Methods for Tests in Linear Mixed Models. *Journal
of Statistical Software* 59(9).

Agresti (2013). *Categorical Data Analysis*, 3rd ed., chapter 8.

Wood (2017). *Generalized Additive Models: An Introduction with R*, 2nd
ed., section 5.4.

## License

MIT. See `LICENSE`.

## Contributing

See
[CONTRIBUTING.md](https://huttoncp.github.io/illume/CONTRIBUTING.md).
Two things get asked of anything new: that it adds no hard dependency,
and that any diagnostic names a remedy which exists in this package.

Please note that this project is released with a [Contributor Code of
Conduct](https://huttoncp.github.io/illume/CODE_OF_CONDUCT.md). By
contributing, you agree to abide by its terms.
