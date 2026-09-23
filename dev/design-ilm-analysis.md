# `ilm_analysis()`: a conductor for the analysis workflow

Design notes, not a plan of work. Written on 2026-09-23 at Craig's request, to
return to once the package is more mature. Nothing here is scheduled.

The prerequisite has been built: remedies the code can act on,
`ilm_remedies()` and `ilm_apply_remedy()`. The design below uses them.

---

## What was asked

A single call that takes a formula (or a DAG) and a data frame, and then works
through the analysis with sensible defaults and little further input:

1. Run the exploratory functions to find problems in the data that should shape
   the fitting strategy.
2. Fit the model by maximum likelihood.
3. Run the diagnostics that suit that model.
4. Apply the remedies the diagnostics call for, rather than only suggesting
   them, and repeat until the fit is sound or the remedies run out.
5. Refit the final model by REML.
6. Report the results with `summary()` and `anova()`.
7. Break down the significant fixed effects and interactions with
   partial-dependence plots, marginal means and contrasts.
8. Write a `.qmd` or `.Rmd` report from the prose templates: the steps, notes,
   methods, plots and findings.
9. Render it to a self-contained HTML file or a PDF.

## What already exists for each step

| Step | Existing functions | What is missing |
|---|---|---|
| 1 | illumex: `ilm_describe()`, `ilm_frame_issues()`, `ilm_check_missing()`, `ilm_gauss_check()` | A rule for how each finding changes the fit (see the plan below) |
| 2 | `ilm_model()`, and `ilm_dag_model()` for the DAG route | -- |
| 3 | The 17 checks in `fit$checks`, made at every fit, plus `ilm_check_dispersion()`, `ilm_check_zeros()`, `ilm_check_variance()`, `ilm_check_ar()`, `ilm_check_covariate()`, `ilm_check_omitted()`, `ilm_check_proportional()` and `ilm_appraise()` | Which standalone checks run for which family: a table, fixed in advance |
| 4 | `ilm_remedies()` and `ilm_apply_remedy()` (built), with tiers and `fit$remedy_log` | The loop, its stopping rule and a revert rule (below) |
| 5 | `ilm_model(reml = TRUE)` | The rule for when it applies (below) |
| 6 | `summary()`, `ilm_anova()`, `ilm_pb_lrt()` | -- |
| 7 | `ilm_plot_model()`, `ilm_emmeans()`, `ilm_contrast(adjust = "max_t")`, `ilm_trends()`, `ilm_ame()` | Follow-ups declared in advance (below) |
| 8 | `ilm_interpret()` writes the prose | `ilm_report()`, which writes the document |
| 9 | -- | Rendering, with rmarkdown or quarto in Suggests |

Most of the pieces exist. The conductor is a thin layer over them. What makes it
hard is deciding what it may do without asking.

## Three steps that cannot be automated as asked

### Step 4: "try remedies until the diagnostics pass"

Remedies are not all alike. `ilm_remedies()` sorts them into three tiers:

- **numerical**: the same model, fitted harder (more restarts). Safe to
  automate.
- **structural**: a different random-effect or variance structure whose fixed
  effects mean what they meant before. Examples are a negative binomial in
  place of a poisson, a zero part, a dispersion model, a covariance of lower
  rank, a term at a variance of zero removed, or the boundary-avoiding penalty.
  These are standard, but they are chosen by looking at the data, and the
  final p-values do not account for that choice.
- **estimand**: changes what the fixed effects estimate or what their standard
  errors account for. Examples are removing a random effect whose variance is
  not zero, merging categories, pooling levels, transforming the response,
  adding a spline, or dropping a predictor or outliers. These should never run
  automatically.

"Iterate until the diagnostics pass" also treats the diagnostics as stopping
rules, and they are not: with a large n every check fails, and with a small n
none can. A search that runs until the checks pass and then writes a polished
report is exactly what naming the remedies, rather than applying them, was
meant to prevent.

### Step 5: "refit the final model by REML"

- For a gaussian response, REML is exact and is what the variance components
  should be reported from.
- For any other family, `reml = TRUE` gives an approximately-restricted
  likelihood (`fit$reml_exact` is `FALSE`). It reduces the downward bias of
  the variance components but has none of REML's exact properties.
- Likelihood-ratio tests of fixed effects need ML fits under either family.

So the rule is:

1. Select by ML.
2. Refit by REML for reporting when the family is gaussian.
3. For any other family, report the ML fit unless the plan says otherwise.

It is a rule, not a blanket step.

### Step 7: "break down the significant effects"

Choosing follow-ups because they are significant biases them: the effects that
cleared the bar are, on average, the ones overestimated. So:

- Follow-ups are declared in the plan and adjusted for multiplicity
  (`ilm_contrast()` already offers `adjust = "max_t"`).
- On the DAG route only the exposure's effect is causal. Breaking down the
  adjustment covariates as if they were effects is the Table 2 fallacy
  (Westreich & Greenland 2013). The conductor reports the exposure and says why
  it reports nothing else.

## The design

### 1. A plan fixed before fitting: `ilm_plan()`

It holds everything the conductor may decide, decided before it sees an
outcome:

- the formula or DAG, and the focal terms;
- which remedy tiers may run automatically, **numerical only by default**;
- which standalone checks run, by family, from a table the package owns;
- the follow-ups: marginal means, contrasts, the adjustment;
- the estimation rule above (ML to select, REML to report where it is exact);
- the report format.

A plan can come from `ilm_scaffold()` or `ilm_power()`, so the path is plan,
then power, then analyse, then report. A pre-specified plan with decision rules
is defensible in a way that automated tinkering is not.

Step 1's exploratory findings go into the plan as warnings a person reads. They
do not change the model silently. For example, a category with 12 rows becomes
"expect the category_counts check to fail; its remedy is estimand-tier".

### 2. The remedy loop

This uses what is already built. Each round:

1. List the remedies with
   `rem <- ilm_remedies(fit, dispersion = , zeros = , variance = )`, passing
   the standalone results the plan calls for.
2. Keep the remedies whose tier the plan allows and that a refit can make
   (`nzchar(rem$change)`).
3. Take the first. `ilm_remedies()` already orders them by tier, and within a
   tier the no-change simplification comes before the penalty.
4. Refit with `ilm_apply_remedy(fit, rem, id)`.
5. **Revert rule.** If the new fit fails more checks than the old one, discard
   it, log the attempt, and do not try that remedy again.
6. **Stopping rule.** Stop when no allowed remedy remains or after a fixed
   budget of refits (5, say), whichever comes first. The loop never stops
   because every check passed. A check still failing when it stops is reported
   as failing.

The decision log starts as `fit$remedy_log`, which `ilm_apply_remedy()`
already keeps: the check, status, tier, remedy and change for every step.

### 3. A decision log and a sensitivity table

The log records every check, every remedy tried, why, and the model it led to,
including reverted attempts. Next to it goes a table of the focal estimates
from the pre-specified model and from the final one. If a remedy moved a focal
estimate, the reader sees by how much.

### 4. `ilm_report()`

`ilm_report()` is separate from the conductor, and useful without it. It writes
a `.qmd` or `.Rmd` holding the actual code for every step, so the report is a
script someone can rerun and edit, not a black box.

- The prose comes from `ilm_interpret()`, the figures from `ilm_plot_model()`,
  and the tables from `summary()` and `ilm_anova()`.
- It renders to HTML by default, and to PDF only when LaTeX is installed.
- rmarkdown and quarto go in Suggests, behind `requireNamespace()`.

### 5. Validation before release

This is the package's usual standard. Run the conductor on the messy-data
scenarios, and measure the error rate and interval coverage of its final
inference against those of the pre-specified model. If structural remedies
inflate the error rate, the default stays numerical-only and the vignette says
so, with the numbers.

## What `ilm_remedies()` covers now, and what the conductor would still need

**Covered:**

- every check made at every fit, 17 of them across 15 names, each with a rule;
  `test-remedies.R` fails if a check has none;
- the three simulation checks whose remedies are single `ilm_model()`
  arguments: dispersion (`family = "nbinom"`), zeros (`ziformula`, both zero
  types) and variance (`dispformula`);
- refits through the fit's own call, keeping the family, weights, zero part,
  dispersion model and REML, with a call a person could have written;
- `fit$remedy_log`.

**Still needed:**

- Remedies for the other standalone checks:
  - `ilm_check_ar()`: adding `ar = ilm_ar1(time, group)` needs the time and
    group the check was given, which its result does not yet carry.
  - `ilm_check_covariate()`: a smooth of the covariate, which is estimand-tier.
  - `ilm_check_omitted()`: adding the variable, also estimand-tier.
  - `ilm_check_proportional()` and `ilm_check_collinearity()`.
- An order within the structural tier that the validation study supports. The
  current order is the rules' own and has not been measured.
- The revert rule and the budget, both above.
- Fixed-seed snapshot tests of the decision log.

## Size

A first version with numerical remedies only is a moderate project, about the
size of the power work. The expensive parts are the standalone-check table and
the validation study. The conductor touches everything, so it should stay a
thin layer over functions that are each useful alone. Whatever it drives
defines the package's core.

## Open questions for Craig

- Should structural remedies ever run by default, or only when the plan opts
  in?
- The refit budget, and whether a reverted remedy may be retried after a later
  one.
- The name: `ilm_analysis()`, or something that says it follows a plan.
- Whether `ilm_report()` comes first, since it is useful on its own and the
  conductor needs it.

## Reference

Westreich, D., & Greenland, S. (2013). The Table 2 fallacy: presenting and
interpreting confounder and modifier coefficients. *American Journal of
Epidemiology*, 177(4), 292--298.
