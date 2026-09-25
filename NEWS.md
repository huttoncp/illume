# illume 0.0.8.9000

* illume now requires illumex 0.0.8.9000.
* **Which groups** a prediction, a fitted value or a residual is for is now
  one argument, `groups`, with the same words everywhere:
  - `"fitted"`: each group's own estimated effects;
  - `"typical"`: every random effect at zero, a group exactly at the average;
  - `"population"`: averaged over the groups;
  - `"new"`: a group the fit has not seen. No function takes it yet: such a
    group's prediction is a spread, not one value.

  Each function takes the words that apply to it, with its default unchanged.
  `predict()` takes `"typical"` (the default), `"population"` or `"fitted"`,
  and `ilm_ame()` `"typical"` (the default) or `"population"`.
  `ilm_fitted()`, `ilm_scores()` and `ilm_rqr()` take `"fitted"` (the
  default) or `"typical"`. The `marginaleffects` bridge takes `"population"`
  (the default, from the option `ilm_model.groups`) or `"typical"`. Any other
  word is an error that says where it is answered.
  - `predict(groups = "fitted")` is new. It gives each row its own group's
    estimated effects, as lme4's `predict()` does by default, and with a
    correlation over time the fitted value of the row's cell. New rows are
    matched to the fitted groups by label and to the fitted cells by time.
    A group the fit has not seen, or a time off its group's fitted cells,
    is an error that names `"typical"` and `"population"`, not a silent zero.
    Its intervals draw the random effects with everything else, through
    `ilm_draws()`.
  - The words replace `marginal` and `conditional`, which still work for one
    release, with a warning once per session: `marginal = TRUE` is
    `groups = "population"`, and `conditional = TRUE` is `groups = "fitted"`.
    The same holds for the option `ilm_model.marginal`.
  - "Conditional" had meant two things: a fitted group's effects in
    `ilm_fitted()`, and effects at zero in `predict()`'s help.
* Fixed: the prints of `ilm_effects()` and of `ilm_emmeans()` for categories,
  and the "Effect size and power" and "Workflow" vignettes, said `ilm_ame()`
  gives the population-averaged effect. By default it gives a typical
  group's. They now name `ilm_ame(groups = "population")`.
* Fixed: `ilm_fitted()`'s help said it gives each category's probability for
  an ordinal outcome. It gives the latent linear predictor, which the
  thresholds cut into categories, as the checks and scores built on it
  expect. `predict()` gives the probabilities.
* `ilm_ranef()` lists a fitted model's random effects: the conditional modes,
  labelled by grouping variable, level and coefficient, with their conditional
  SDs. That is lme4's `condVar`, taken from the Laplace approximation's inner
  Hessian. There is a row per basis function of a smooth and per cell of a
  correlation over time. `ilm_varcorr()` gives the variance components on
  their natural scale: each term's covariance with lme4's `stddev` and
  `correlation` attributes, the correlation over time's parameters, and the
  dispersion, named for what it is. For a gaussian model that is `sigma`, the
  residual SD, not the variance.
  - Both are registered as `ranef()` and `VarCorr()` methods for nlme's
    generics, which lme4 re-exports.
  - On a model both packages fit, the modes, SDs and covariances agree with
    lme4's to 1e-3.
  - Each random effect's `row` is its position in the fit's full parameter
    vector, where it matches the fit's own value exactly, by ML and by REML.
    That is the check that would have caught the REML defect below.
* `ilm_draws()` draws every parameter of a fit jointly from its joint
  precision: the fixed effects, the variance and dispersion parameters, each
  group's random effects and each cell of a correlation over time.
  - **Layout.** The rows are named by block, with a map that labels each
    one: coefficient, group and dimension, or cell.
  - **Natural scale.** Each draw's variance components are given on the
    natural scale, through the same transform as `ilm_varcorr()`, so at the
    estimate the two agree exactly.
  - **Boundaries.** Where a fit holds a boundary direction, the draws hold
    it exactly, by conditioning on it. The fit now keeps the directions it
    held and the Hessian its covariance came from.
  - **Options.** `given = "theta"` holds the variance parameters at their
    estimates. `blocks` returns part of the vector. A fit made without
    `joint = TRUE` has its precision formed on demand.
* `ilm_matrices()` gives the designs for new rows: the fixed design, each
  random term's design with each row's group matched by label, each
  smooth's penalised basis, and the zero part's and dispersion model's
  designs. It also places each new row among a correlation over time's
  cells: the cell it falls on, the cells either side and the time to each.
  A prediction assembled from these matrices is `predict()`'s. The smooth
  basis comes from the same code, now shared.
* `ilm_dist()` gives a fit's response distribution as functions: the
  density, distribution and quantile functions, a random generator and the
  mean, in illume's own parameterisation.
  - **Inputs.** They take linear predictors and internal parameters, a draw
    from `ilm_draws()` say, rather than means. The inverse link, the
    dispersion's transform, the zero part's mixing or truncation, a censored
    row's interval and an ordered response's thresholds are all handled
    inside.
  - **Checked against the fit.** On a fit with nothing integrated out, the
    density summed over the rows is the objective the optimiser minimised,
    to 1e-9. That holds for every family, for censored gaussian and survival
    responses, for frequency weights, for zero-inflated and hurdle models
    and for dispersion models.
  - **Not covered.** The flexible parametric survival families, whose linear
    predictor depends on time itself.
* `ilm_normal_expect()` exports the Gauss-Hermite quadrature behind
  `predict(groups = "population")`, so code built on a fit averages over a
  latent spread as its predictions do.
* Fixed: `fixef()` did not reach illume's method. It was registered on a
  generic of illume's own, so with nlme or lme4 attached, `fixef(fit)`
  stopped with "no applicable method". It is now registered on nlme's
  generic.
* A fit keeps each grouping term's level labels, and a refit keeps them too.
* Fixed: `ilm_denom_df(method = "kenward-roger")` was not Kenward-Roger.
  - **The covariance shrank.** Its "inflated" covariance was built from
    differences of the fixed-effect covariance with the wrong algebra: the Q
    term had the wrong sign and the P-Phi-P term was missing. So the
    covariance came out smaller than the unadjusted one, by a factor of
    1 - 4/n for a single variance component.
  - **The df were Satterthwaite's.** Its df were Satterthwaite's formula
    applied to that matrix, not Kenward and Roger's. Its one test asserted
    only that the df were below Satterthwaite's, which the defect
    guaranteed.
  - **The gate did not refuse.** It claimed to refuse a correlated random
    slope but tested the category width, which is always 1.

  It is now computed as `pbkrtest` computes it: from the covariance of the
  observations, with Kenward and Roger's own df and F scaling. It agrees
  with `pbkrtest` to 2e-6 in df and 3e-8 in the F scaling, and in the
  adjusted covariance to within the two fits' REML difference. That holds
  for a random intercept, correlated and uncorrelated random slopes, and
  crossed factors. In a balanced design it reproduces the exact test: 6 df,
  to 1e-8.

  A correlated slope is linear in its covariance's elements, so it is now
  covered. Kenward-Roger now needs a REML fit (`reml = TRUE`), because it
  is derived for REML estimates. It stops for smooths, weights and more
  than 4000 rows. `ilm_trends(df = "kenward-roger")` now uses the adjusted
  standard errors as well as the df. Found in a review.
* `ilm_rw1()` fits a random walk over time, the local-level model of
  time-series analysis. Each group's walk is held at zero at its first time,
  so the fixed effects give its level there, and each step after that adds an
  independent change with variance proportional to the time elapsed, so
  irregular gaps enter as they are. With a gaussian response the fit is the
  exact local-level likelihood. On five reference series (one to three
  groups, unit steps and irregular gaps), by maximum likelihood and by REML,
  the estimates agree with an independent maximisation of the exact
  likelihood to within 5e-6 and the log-likelihood to 1e-10. A forecast of
  a group's last cell from the fit (its mode, and its variance from the
  joint precision) agrees with the closed form to 1e-8. For other families it is
  the Laplace approximation, and the fit-time check on observations per
  latent value applies to it as it does to CAR(1).
  `predict(groups = "population")` averages each row over its own spread,
  which grows with the time since its group started.
* `ilm_ar1()`, `ilm_car1()` and `ilm_rw1()` can name their two columns:
  `ilm_car1(~ day | id)`. `ilm_model()` reads them from its model frame, after
  rows with missing values are dropped, so the structure cannot fall out of
  step with the response. Vectors taken from the whole data frame could, and
  stopped the fit with a length mismatch. The fit keeps the names, which is
  what lets `predict()` place new rows on a walk.
* `ilm_cells()` lists the latent cells of a fitted AR(1), CAR(1) or
  random-walk term, one row per cell. Each row gives the group and time (as
  numbers, dates or date-times, as given), how many observations sit on the
  cell, where the fit holds its value, and whether it is its group's first or
  last cell. Code that works with the latent values, such as a forecast, reads
  the layout from it instead of rebuilding it.
* `summary()` names the correlation over time it fitted and gives its standard
  deviation. It used to call a CAR(1) term "ar1".
* `ilm_interpret()` says each effect in the response's own units, over a change
  a reader can picture: across the middle half of a numeric predictor, level by
  level for a factor, with the model's predictions at both ends, averaged over
  the rows the model was fitted to as `ilm_ame()`'s effects are. Where it said
  "a higher age is associated with a higher value of income_k (estimate 0.610,
  95% interval 0.493 to 0.726, p = <1e-04). In the units of the response: an
  increase of 0.610 on average", it now says "Across the middle half of age,
  predicted income_k is 49.4 at 29 against 64.7 at 54: 15.2 higher (95%
  interval 12.3 to 18.1). That is 0.61 per unit of age." A probability is given
  as one ("57% at 29 against 36% at 54: 20 percentage points lower"), with the
  odds ratio after it. Causal language, where a design licenses it, reads as a
  change: "moving x across its middle half ... raises predicted y". In a mixed
  model with a nonlinear link, it says that the predictions hold the random
  effects at zero -- a typical group -- rather than averaging over groups,
  which on the logit or log scale is not the same.
* The evidence for a term is its joint test, so a factor with several levels
  gets one verdict, and so does a predictor in a multinomial model, which used
  to get one per category, each against the average of the categories -- and
  could then report lower odds of a category beside a higher probability of it.
  Its effect is now the predicted share of every category at both ends.
* P-values in the prose are given to three significant figures, or as
  "p < 0.001"; a factor level is named as a value ("tenure 'rent' rather than
  'own'", not "tenurerent: being rent rather than own").
* Fixed: a number shorter than four characters came into `ilm_interpret()`'s
  sentences with spaces in front of it: "(95% interval  13 to  32)".
  Printing re-wraps the text, which hid it. The strings themselves, which a
  report or a paste takes, carried it.
* `ilm_interpret()` writes up a cluster profile from `illumex::ilm_profile()`
  or `ilm_profile_na()`: the clustering, a paragraph per cluster, the rows
  between clusters and the variables that only add distance, and why none of
  it is a test.
* Fixed: every REML fit predicted `NA`. Under REML the fixed effects sit in the
  random block and the stored `beta` was built from a parameter vector that did
  not hold them, so `predict()` -- and with it `ilm_ame()`, `ilm_scenario()`
  and the effects in `ilm_interpret()` -- returned `NA` for every
  `ilm_model(reml = TRUE)` and every `ilm_dag_model()`, which fits by REML.
  `coef()` and the standard errors were never affected.
* Fixed: a REML fit's random effects were read from the wrong entries.
  - **Cause.** TMB orders its random block by the parameter list, and under
    REML the fixed effects are in that block and come first. Read by
    position, every random effect was the entry p * C places before its
    own, and the first few were the fixed effects themselves.
  - **What it broke.** Everything conditional on the random effects, on
    every REML mixed model: fitted values with each group's own effects, the
    quantile residuals, `ilm_check_predictive()`, the dispersion and zeros
    checks, and `ilm_plot_model(what = "random")`. On a gaussian random
    intercept, the conditional predictor averaged 0.208 against a response
    mean of 0.163.
  - **Smooths.** A penalised smooth's coefficients are random effects too,
    so `predict()` on a REML fit with a smooth was wrong outright: a sin()
    curve came out as -2.08 where it was 0.48. `ilm_ame()` and
    `ilm_scenario()` on such a fit were wrong with it.
  - **What was not affected.** The fit itself, `coef()`, `vcov()` and the
    intervals from joint draws were right all along.
  - **Now.** The random effects are read by name, as the AR block already
    was. They agree with lme4's REML modes, and a REML smooth with mgcv's.
    Found independently by two other agents.
* Fixed: `predict()`'s intervals for a model with a smooth were centred on the
  wrong parameters whenever the family has one the model declares after its
  random effects -- a dispersion (gaussian, negative binomial, beta, the
  survival families), cut points (ordinal), a zero part, a dispersion
  formula, an AR term. The draws behind them had the right spread around the
  wrong centre: for a gaussian `y ~ s(x)`, a fit of 0.245 had a "95% interval"
  of -3.30 to -2.98, and every interval missed its own fit. Through a
  nonlinear link the standard errors were wrong too -- a negative binomial
  smooth's ran from 0.02 to 16 times mgcv's -- and so were the bands of
  `ilm_plot_model(what = "effect")`, which are built from them. Binomial and
  Poisson smooths, which have no such parameter, were unaffected. The draws
  are now centred parameter by parameter, and `test-joint-draws.R` holds the
  intervals to mgcv's.
* Fixed: `predict(groups = "population")` left an AR(1) or CAR(1) term out of
  the average over the random effects. At any one row its latent value has the
  stationary distribution, and it is now averaged over with the grouping
  terms: a Poisson model with a stationary AR variance of 0.59 averaged 1.40
  where its own simulations averaged 1.86, and now agrees with them and with
  the closed form.
* With one linear predictor -- every family but the multinomial -- a row's
  whole latent contribution is a single normal, so
  `predict(groups = "population")` now averages over it by Gauss-Hermite
  quadrature rather than by draws, with more nodes as the latent SD grows. Against numerical integration it is
  better than 1e-12 through a logit up to a latent SD of 6, and better than
  1e-9 through a complementary log-log up to 5; a fixed 40 nodes, as first
  written, drifted to 2e-5 at a logit SD of 3.5. It gives the same answer on
  every call and does not use `ndraw`; the effects and scenarios built on it
  no longer move with the seed. A multinomial outcome is still averaged by draws, now
  including an AR term.
* Fixed: `ilm_emmeans(type = "response")` on a beta model returned the
  link-scale means: 0.37 where the proportion was 0.59. The inverse link was
  chosen by switching on the family's name, and beta fell through to the
  identity. It now comes from the family object itself, which closes that for
  any family added later. The link name illume gives the easystats packages
  had the same fault for beta and the ordinal links, which it reported as
  "log", and now reads the family's own link. `logLik()`'s help now says that
  for a binomial with trials it leaves out the binomial coefficients, which
  `glm()`, `lme4` and `glmmTMB` include: a constant of the data, so no
  comparison of models moves, but the values differ by exactly
  `sum(lchoose(trials, successes))`. Found by another agent.
* Fixed: `ilm_scenario()` computed its means from the fixed-effect design and
  each random term's intercept variance alone. A random slope was averaged as
  if it were an intercept, an AR term was left out, and a smooth lost its
  penalised part entirely: a gaussian `sin()` curve came out as -1.31 at its
  peak of +1. Each mean now comes from `predict(groups = "population")`, and the
  interval draws the whole parameter vector, variance components included,
  where it drew the fixed effects alone. The estimate is the mean at the
  fitted parameters, which matches `predict()` exactly and does not move with
  `seed` or `sims`; it was the mean of the draws, which through a logit was
  pulled towards the middle by up to half a percentage point. An ordinal fit is refused, as a
  multinomial one was, since it too has no single number to report.
* `ilm_ame()` says which effect it reports: by default the one for a typical
  group, with the random effects at zero -- the effect `ilm_interpret()`
  describes -- and with `marginal = TRUE` the one averaged over them, which
  through a logit is the flatter. Its help used to justify its standard errors
  by a population average it did not take. The `marginaleffects` bridge
  averages over the random effects by default, as `marginal = TRUE` does.
* `car::Anova()`, `performance::model_performance()` and
  `performance::check_model()` now reach illume's methods, as the README and
  the introduction always said they did. The methods were exported as
  ordinary functions rather than registered with the generics, and R's
  method lookup does not find those: `car::Anova()` quietly fell back to its
  own tests -- on a multinomial fit, 1 and 2 degrees of freedom where the
  joint tests have 2 and 4 -- `model_performance()` returned `NULL` with a
  warning, and `check_model()` stopped with an error, because performance
  cannot read these fits (not even after `ilm_register_insight()`). They are
  now registered for whenever car and performance are loaded.
  `check_model()` draws `ilm_appraise()`'s panels, and its help page says
  whose they are.
* The documentation says what the package fits. `ilm_model()` is "Fit
  generalized linear and additive mixed models" rather than "Fit a
  multinomial linear mixed model", its `family` argument lists all fifteen
  families, and `summary()`, `predict()`, the fitted values, the quantile
  residuals, `ilm_appraise()`, `ilm_simulate()`, `ilm_anova()`, the fit
  indices and `ilm_fit()` no longer describe themselves as multinomial-only.
  Where something IS about categories, it now says so.
* `re_struct` is documented with examples in `ilm_model()`, `ilm_fit()` and
  the regression-models vignette: one element per random term, named by its
  grouping variable, holding `type` (and `rank`) for the covariance across a
  multinomial outcome's categories, and `d_cor = FALSE` for an uncorrelated
  random slope in any family. An element may leave `type` to its default
  now: `list(subj = list(d_cor = FALSE))`, which is all an uncorrelated slope
  needs outside a multinomial model, used to stop and ask for one.
  `ilm_model()`'s example, which never ran, is replaced by ones that do.
* `ilm_fit()` defaults to `family = "gaussian"`, as `glm.fit()` does, and to
  no random terms (`re_list = list()`), so `ilm_fit(X, y)` is a linear model.
  Its default was `"multinomial"`, where the package began. A call that gives
  `J` three or more categories and no family is one written for that default:
  it stops and asks for `family = "multinomial"` rather than fitting a
  gaussian model to the category codes. `ilm_model()` is unchanged: it reads
  the family off the response.
* The parameter-aliasing check no longer starts a sentence with a capitalised
  parameter name ("Retired:(Intercept) <-> retired:age cannot be separated");
  it reads "These data cannot separate retired:(Intercept) from retired:age".
* A covariance term flagged at its boundary but curved in every direction
  always takes the recomputed Hessian now, as a fit with nothing at a
  boundary does. It still went to TMB's own Hessian whenever TMB called that
  positive definite, so its standard errors depended on the verdict the
  change below was meant to take out of the decision: 2 to 10% apart in one
  fit of the boundary study's re-run.
* `ilm_consistency()` counts a refit that lands on a covariance boundary.
  Such a refit has that direction held, so `pdHess` is `FALSE` by design and
  its fixed effects are usable, but the check asked for `pdHess` and dropped
  it. Since every boundary is held, that meant every refit at one: on a
  multinomial model whose subjects vary along one direction, 11 of 30 refits
  were dropped and the check reported the model "UNSTABLE", unable to recover
  itself. It refits all 30.
* `ilm_consistency()` calls a variance fitted at zero a boundary whenever the
  fit does -- below 1e-3 -- as well as when it is small beside the model's
  largest. With one random term, that term was compared with itself, so a
  random intercept fitted at 0.0001 was reported as biased, and the Laplace
  approximation blamed, in a gaussian model, which has none.
* The study set was re-run for this version and reproduces 0.0.7.9000
  wherever nothing changed (`studies/findings/`). Two of its scripts had
  fallen behind the package, and both are fixed. Four counted a fit as usable
  only when its Hessian was positive definite, which a fit held at a
  covariance boundary is not, by design: the coverage study's "convergence
  rate" for a thin five-category model read 0.088 where 95% of its fits were
  usable. And the imputation study still looked for a function that moved to
  illumex. Counted as the package counts them, the multinomial cells cover at
  0.946 to 0.954, and the fits held at a boundary, on their own, at 0.942 to
  0.955.

# illume 0.0.7.9000

Ten features, each validated against an outside implementation where one
exists, and then a second round that arrived as defect reports from a separate
project built against the package. It now covers the path from a sample-size
calculation -- with or without a pilot to base one on -- to a scenario a
stakeholder can act on.

Every one of them produced at least one real bug, and none of those came from
the test suite. They came from comparing against something outside, from
widening a simulation, or from someone calling illume from outside it, which
remain the only things that have ever found a defect here.

## The exploratory half is its own package: illumex

* Describing, cleaning, plotting and profiling data before a model -- 65 of
  the 147 exported functions -- moved to a new package, `illumex`, which
  `illume` attaches (`Depends`). Every function keeps its name, arguments and
  behaviour, and `library(illume)` makes them available as before, so no code
  that uses them changes. What moved never fitted or read a model: the call
  graph was checked before the line was drawn, and the only call across it is
  `ilm_impute()`'s use of `ilm_glrm()`.
* What moved: `ilm_describe*()`, `ilm_counts*()`, `ilm_dupes()`, `ilm_copies()`,
  `ilm_wash_df()`, `ilm_recode_errors*()`, `ilm_translate()`,
  `ilm_frame_issues()`, `ilm_gauss_check()`, `ilm_boot_ci()`, `ilm_boot_diff()`,
  outliers and anomalies, `ilm_reduce()`, `ilm_cluster()`, `ilm_profile()`,
  `ilm_glrm()`, `ilm_check_missing()` and the `*_na()` functions, the plots of
  data (as against those of a model), `ilm_sim()`, and the `progress` argument's
  help page. So did the vignettes on exploring data, profiling and anomaly
  detection. Imputing and pooling stay here, because pooling needs a model.
* Why: the package had grown to span five fields, each a package elsewhere.
  The split moves 44% of the functions but only about 20% of the code -- the
  engine, where the hard maintenance is, stays -- and it gives the exploratory
  half an install with no compiled dependency, a faster check, and a release
  cycle of its own. `cluster`, `isotree`, `PCAmixdata`, `boot` and `tibble`
  are no longer suggested here.
* A handful of small helpers exist in both packages rather than one reaching
  into the other's internals; `test-shared-helpers.R` fails if the copies
  differ. Until `illumex` is on CRAN, `illume` finds it on GitHub through
  `Remotes:`.

## Remedies the code can make

* `ilm_remedies(fit)` lists a remedy for every check that is not OK. Each one
  is written as the change to `ilm_model()` that makes it --
  `re_struct = list(site = list(type = "rr", rank = 1L))`,
  `formula = y ~ x + (1 | g)`, `family = "nbinom"` -- and
  `ilm_apply_remedy(fit, remedies, id)` makes it. The refit goes through the
  model's own call, so the family, weights, zero part, dispersion model and
  REML all come along. It then says what the checks that asked for the remedy
  say now. The new fit's call is the one a person would have written, ready
  for `update()` or a further remedy, and `fit$remedy_log` records every
  remedy that led to it.
* The checks that simulate -- `ilm_check_dispersion()`, `ilm_check_zeros()`
  and `ilm_check_variance()` -- are passed in rather than run, since running
  them is a choice:
  `ilm_remedies(fit, dispersion = ilm_check_dispersion(fit))`.
  `ilm_check_variance()` now returns the column its groups came from, so the
  dispersion formula can be written in it.
* Every remedy has a tier.
  * `numerical`: the same model fitted harder.
  * `structural`: a different random-effect or variance structure, whose fixed
    effects mean what they meant before. Examples are a term whose variance is
    estimated at zero removed, a covariance of lower rank, a dispersion model,
    a zero part, or the boundary-avoiding penalty with its measured costs.
  * `estimand`: changes what the fixed effects estimate or what their standard
    errors account for. Apply one only because the question calls for it.

  Whether a variance is at zero is read off the fit, so "drop the term" is
  structural where the fit is the same without it, and estimand everywhere
  else. Some remedies can only be made by hand, such as which categories to
  merge or which levels to pool; they are listed with no change, and applying
  one is refused.
* A remedy is a candidate, not a cure. The tests hold it to what it claims:
  * dropping a random intercept at zero leaves the fixed effects at lm's to
    1e-5;
  * a smooth shrunk to a line, replaced by the line, reaches the same
    log-likelihood with one parameter fewer;
  * a category covariance at its edge given rank 1 passes the check it failed.

  `test-remedies.R` also fails if any check the package can report has no rule.
* The README, the *Introduction* and the *Regression models* vignette show
  them. The vignette works through a five-category covariance too rich for
  twenty subjects: the checks name a rank, `ilm_remedies()` writes it out, and
  `ilm_apply_remedy()` takes the level check from FAIL to WARN.
* Writing the rules out found five remedies that were wrong:
  * `re_struct` naming only some of a model's terms stopped with an error about
    `re_struct$NA`. The argument reads as if a term left out keeps the default,
    and now it does. A name that is not a term of the model now stops, listing
    the terms, instead of being ignored.
  * For any family but the multinomial, a random-effect term with too few
    levels was told to use `rr(1)` -- "1 parameters instead of 1", which is
    the structure already fitted. It now names what can come down: the
    intercept-slope correlation, the random slope, the levels, the term. The
    latent-budget check had the same fault and the same fix.
  * A smooth penalised to nothing was told "drop term 's(x)'", but its
    unpenalised part -- for `s(x)`, a straight line in `x` -- is still fitted,
    and dropping the term takes the line out too. The remedy is now to put `x`
    in its place.
  * Two remedies named things the package does not have: "fix one of the two
    parameters" for aliased parameters, and "a reduced-rank AR" for a thin AR
    grid. Both are gone.

## The family is read off the response

* `ilm_model()` defaulted to `family = "gaussian"`, so leaving the family off
  a model of a three-category outcome stopped with "gaussian family needs a
  numeric response", and leaving it off a count fitted a straight line to it
  without comment. The default is now `"auto"`: the family is read off the
  response, and the fit says which it chose and on what evidence --
  ``family = "poisson", inferred from `visits`: whole numbers from 0 to 14, a
  count`` -- and names the usual competitor, `"nbinom"` for counts and
  `"ordinal"` for categories with an order. The choice is written into
  `fit$call`, so anything that refits through the call keeps it instead of
  guessing again on simulated data, and `print()` and `summary()` mark the
  family as inferred.
* Whole numbers are the hard case, because a count and a measurement recorded
  to the nearest unit look alike. Whole numbers that reach down to 0 or 1 are
  read as a count; ones that never come near zero, such as a blood pressure,
  as a measurement. The range is in the message either way, so a misreading
  shows at a glance.
* What the response cannot settle is asked about, not guessed: a proportion
  that touches 0 or 1 or comes with weights, a two-valued variable coded other
  than 0/1, a survival time -- which survival family is a modelling decision --
  and a date.
* `ilm_dag_model()`, `ilm_did()` and `ilm_rdd()` infer through the same rules.
  Their own guesser sent an ordered factor to the multinomial, proportions to
  the gaussian, and a variable coded 1/2 to the binomial, which then refused
  it.
* `ilm_mi_pool()` infers the family once, on the first imputation, and holds
  it, so every imputation is fitted with the same likelihood. An imputation
  that fails to fit is now reported; it used to be dropped without a word, and
  the draw that made the model hardest to fit is not a random one to lose.

## Power that is the power of the analysis you will run

* `ilm_power()` refitted each simulated study with `ilm_model(formula, data,
  family)` rebuilt by hand. That dropped the contrasts -- a model fitted with
  sum-to-zero coding came back with differently named coefficients and a
  power of exactly **zero** -- along with the zero part, the dispersion model
  and the censoring, and it failed outright on a transformed response or
  predictor, which the model frame holds only transformed. Each replicate is
  now the fitted model's own design restricted to the rows the study drew,
  refitted through the machinery every other test in the package uses.
* It counts the test the analysis will report: t on the residual degrees of
  freedom where nothing is integrated out, the Wald z otherwise, and for a
  term with several coefficients the joint F or chi-square `ilm_anova()`
  reports. A z test in place of the t put the power of a 20-person two-arm
  trial at 0.43 where the exact value is 0.395; it is now 0.387 (0.374 to
  0.399 over 6,000 replicates).
* A multinomial model used to fail with "non-conformable arguments". A term
  is now tested across all its categories -- or one category at a time, by
  naming the coefficient -- and `effect` is then a multiple of the assumed
  coefficients. Against the asymptotic noncentral chi-square it gives 0.818
  where theory says 0.813 at 400 rows. At 100 and 200 rows it runs 3 to 4
  points lower, as the Wald test's conservative size there (4.2% and 4.5% at
  a nominal 5%) predicts: that is the test the analysis reports.
* Random effects are drawn as the matrix normal the model fits, with a
  multinomial term's covariance across categories and a random slope's own
  variance, and serial correlation, dispersion models, zero parts, censoring
  and flexible survival baselines are simulated and refitted as themselves.
  Frequency weights are drawn as the observations they stand for.
* A check made before fitting is a verdict on the design, the same for every
  study drawn from it. It is now reported with its remedy rather than counted
  against each replicate: a three-category multinomial with four visits per
  participant fails the latent-budget check whatever the data, and counting
  it put the power at exactly zero.
* Recorded rather than smoothed: with 60 clusters, the Wald test of a
  between-cluster effect in a multinomial mixed model rejected a true null
  5.8% of the time with four visits per cluster and 7.8% with eight. The
  power reported is the power of that test, and `ilm_interpret()` says so.
  `ilm_pb_lrt()` is the calibrated alternative for the analysis itself.

## Planning a study with no data behind it

* `ilm_scaffold()` and `ilm_power_design()` take multinomial and ordinal
  outcomes, stated as the probability of each category in each cell or as
  coefficients (with thresholds, for an ordinal one). Probabilities no
  proportional-odds model can produce are refused, and the model that can
  produce them is named. Both used to fail, one of them with a message
  calling it a bug in `ilm_scaffold()`.
* An ICC for a binary or ordinal outcome is taken on the latent scale, the
  usual convention (`pi^2 / 3` for a logit link). It used to demand an `sd`
  such a model does not have.
* `reml = TRUE` plans for an analysis fitted by REML. On a 2 x 2 design with
  an ICC of 0.5 the maximum-likelihood analysis's power ran about 1.5 points
  above the REML one (0.358 against 0.342 at 20 per arm), because its
  variance components run small; the REML figure is within Monte Carlo error
  of `simr`, which also uses REML.
* Each simulated study draws the planned design afresh, instead of
  resampling the scaffold's own grid: the allocation is balanced the way a
  protocol randomises a factorial, and a covariate given as a function is
  drawn again. Resampling left cells of a small factorial empty, so 18
  participants in a 3 x 3 design converged 28.5% of the time. They now
  converge every time.
* Allocation itself was balanced factor by factor, which left the CELLS to
  chance, even in the scaffold's own grid. It is now balanced by cell, with
  any remainder spread so that each factor's own levels stay within one of
  each other.

## Marginal means for categorical outcomes

* `ilm_emmeans()` gives a multinomial fit a row for every category. On the
  link scale that row is the category's centred log-odds. On the response
  scale it is the probability, averaged over the grid, with a delta-method
  standard error. `ilm_contrast()` compares groups within each category, as
  differences in probability. Estimates and standard errors agree with
  `emmeans` on `nnet::multinom()` to 1e-6 under all three weightings, and so
  do the contrasts.
* An ordinal fit's `type = "response"` gives each category's probability,
  with standard errors that carry the thresholds' uncertainty as well as the
  slopes'. It agrees with `emmeans` on `MASS::polr()` to 1e-6, logit and
  probit.
* `ilm_trends()` gives a multinomial fit a slope for every category, and
  agrees with `emtrends`.
* Naming the response among the variables to average over now says so. It
  used to fail on "undefined columns selected".

## A survival model's refits had the wrong baseline

* A flexible parametric (`rp`) model's baseline is a spline in log time, so
  its columns are built from the response. Every refit to a simulated
  response kept the columns built from the observed times: the parametric
  bootstrap, the simulation envelopes and the consistency check all refitted
  a model whose baseline did not belong to its data. The columns are now
  rebuilt from the new times at the same knots, and such a refit matches a
  fresh fit to the same times exactly.

## Writing it up

* `ilm_interpret()` says when the family was read off the response, and
  describes an ordinal fit in categories and cumulative odds ratios rather
  than as an identity-link model. A covariance at its boundary gets the same
  verdict `summary()` gives; the interpretation used to say that all fitting
  checks had passed. With fewer than 100 groups it notes that the tests are
  large-sample ones and names `ilm_pb_lrt()`.
* It also writes up an `ilm_power()` result: the power at each size with its
  Monte Carlo interval, the size that reaches 80% with the range the
  simulation supports, and any verdict on the design. An `ilm_contrast()`
  result is written up as the comparisons it makes.

## A variance at its boundary no longer takes the fixed effects with it

* When a random-effect variance sits at zero, or a correlation at +/-1, the
  likelihood is flat in that direction and the Hessian over every parameter is
  singular however well the fixed effects are determined. The whole fit used to
  be graded FAIL, standard errors and all. Now the covariance of that term is
  held at its estimate in the direction the data cannot resolve and the rest
  of the Hessian inverted (the next section says why only that direction),
  and the fit says, in `summary()` and `print()`, what can be trusted: the
  fixed effects, their standard errors and tests, yes; that covariance, no.
  The check reads BOUNDARY, not FAIL.
* Measured on the messy-data regime with a true random-effect SD of 0.05, on
  the study's own 400 seeds: fits with usable fixed effects rose from **63.3%
  to 97.5%**. The 128 newly usable ones cover at **0.944**, and their standard
  errors are 1.011 times those of the same model without the random term,
  which is what a variance of zero says they should be. In the combined regime
  the share rose from 91.0% to 99.0%, and the 30 fits recovered there cover at
  **0.908** -- a little under nominal, and recorded rather than smoothed. The
  clean, unbalanced and heavy-tailed regimes are unchanged.
* The study's own tables still report what was measured at the time; the next
  run of `messy_compare.R` will count fits the new way.
* glmmTMB recovers some of these by recomputing the Hessian more accurately
  before giving up, and illume now does the same without its dependency. On
  these regimes that alone recovered **no** fit: every failure was a genuine
  boundary, which glmmTMB lets through only because its threshold is machine
  epsilon, reporting a standard error in the billions for the variance.
* At a boundary nlminb stops with "false" or "singular convergence" as a log
  standard deviation drifts towards minus infinity. The boundary term is now
  held where it got to and the rest allowed to finish, as a variance
  constrained at zero would in lme4.
* A covariance fitted at a correlation of -1, with a positive definite Hessian
  and standard errors within 4% of the model without the term, used to be
  graded FAIL with "the standard errors above are not usable" beneath them.
  It is a BOUNDARY now too.
* Fixed effects that are not identified are not rescued. Aliased columns passed
  a `chol()` test on a finite-difference Hessian, whose noise makes a singular
  matrix slightly positive; the test is now scale-free and the fit stays FAIL.
* On R release, with RTMB 2.0, a multinomial covariance heading for a
  correlation of -1 kept going past the point where the Laplace
  approximation's arithmetic holds. A correlation of +/-1 is a log-Cholesky
  diagonal at minus infinity, and the random effects' density is evaluated
  through the inverse of that factor, so near -15 the inner Hessian carries a
  term of size e^30 that its factorisation has to cancel back down to the
  data's curvature. The objective there is noise, and the optimiser found a
  "maximum" 0.3 better than the true one, stopped with false convergence at a
  gradient of 1.08, and printed standard errors of **0**. It was not only
  RTMB 2.0: with RTMB 1.9, 2 fits in 400 on the flat regime above did the
  same and had been counted unusable. Those diagonals now have a floor of -10, an SD of
  4.5e-5 on the logit scale -- zero for every purpose, and still read as a
  boundary. It applies only where the arithmetic needs it (an unstructured
  category covariance, and the AR innovation covariance, which is built the
  same way), and only to a fit that went past it. On the messy-data study's
  2400 fits it engaged in 144 and moved no coefficient of a usable fit by
  more than 2e-5. The two zero-SE fits, whose coefficients had come from the
  noise and moved by 0.012 and 0.015, now have standard errors within 1% and
  8% of the model without the term.
* A Hessian is no longer recomputed, nor a term held, at a point that is not
  stationary: the inverse of a Hessian differenced where the gradient is 1.08
  is not a covariance, and it was what gave those standard errors of 0. The
  line is the gradient check's own FAIL line, so no fit it calls usable is
  refused.

## At a boundary, only the flat direction is held

* A covariance at its boundary is flat only in the direction that reaches
  it. At a correlation of 1, both standard deviations are still estimated,
  with real uncertainty. The fixed effects' standard errors now hold only
  that direction and estimate everything else. The direction is the
  eigenvectors of the term's block of the Hessian with curvature below 1e-3
  of its largest. The standard errors are then those of the reduced model
  the boundary implies: a covariance of lower rank, or the term dropped.
* Before, TMB's verdict on its own Hessian decided. If it called the Hessian
  positive definite, everything was inverted; if not, the whole term was
  held, as lme4 does. At a boundary that verdict is noise, so the same data
  took different routes on different platforms: R 4.4 held a fit in
  `test-boundary.R` that R 4.6.1 did not, and standard errors could differ by
  up to 10%. The old rule also passed four fits with a quasi-separated
  category as usable, with standard errors of up to 2.9e6. Both problems are
  gone: a boundary is held the same way whatever TMB said.
* Measured on 4,000 simulated fits across ten regimes:
  * the messy-data study's six regimes;
  * a rank-one multinomial covariance;
  * binomial random slopes with a slope SD of 0, a slope SD of 0.25, and an
    intercept-slope correlation of 1.

  Each dataset was fitted once, and the standard errors were computed all
  three ways at the same optimum (`studies/scripts/boundary_se.R`). Among the
  1,433 fits at a boundary:
  * Holding only the flat direction reproduced the reduced model, refitted,
    to a median ratio of **1.000**.
  * Holding the whole term understated the standard errors, by 0.3-2%
    typically and by up to 45% in the worst fit. It covered least in every
    regime. With a rank-one truth, coverage was 0.938 for the whole-term
    hold, 0.943 for the old rule and **0.947** for the flat direction. With an
    intercept-slope correlation of 1 it was 0.943, 0.945 and **0.951**.
  * The flat-direction hold covered 33 intervals that the whole-term hold
    missed, and never the reverse.

  The differences are about a point of coverage at most. The regime with
  everything wrong at once under-covers under all three rules (0.904 to
  0.910). That comes from the rare category's bias, which no choice of
  standard error touches.
* The threshold is not delicate. Every fit at a boundary had a gap of at
  least 2.3 orders of magnitude between its flat and its curved directions,
  and thresholds of 1e-3 and 1e-4 gave the same coverage.
* For a random slope at a variance of zero, the standard errors now come out
  about 1% above those of the random-intercept model, rather than equal to
  them. The intercept-slope covariance sits at zero but is curved, so its
  uncertainty is kept.
* AR(1) and CAR(1) terms at a correlation of +/-1 are held the same way, but
  they were not part of the measurement.

## Or keep it off the boundary

* `ilm_model(boundary = "avoid")` adds the boundary-avoiding penalty of Chung
  et al. (2013, 2015): half the log-determinant of each grouping term's
  covariance, which in one dimension is a gamma(2) prior on the standard
  deviation. No variance can reach zero, and no correlation +/-1. Checked
  against blme, to 3e-7 for a gaussian random intercept, 5e-6 for a binomial
  one and 5e-4 for a binomial random slope; for a gaussian random slope,
  where blme's route failed inside blme, against the exact penalised
  likelihood maximised by brute force, which it matches to every printed
  digit.
* Measured against the default on the messy-data regimes, 400 datasets each.
  With a true SD of 0.05 every penalised fit was usable, against 395 of 400,
  and coverage was no worse (0.949 against 0.946). The cost is in the
  estimates: that SD came out at a median of 0.21 rather than 0.10, and --
  because a larger between-group variance means larger within-group effects
  on a logit scale -- the fixed effects moved further from zero, by about 1%
  in four regimes of six, 11% where an outcome category was rare (coverage
  0.930 to 0.916) and 20% in the combined regime. So `"hold"`, maximum
  likelihood, stays the default.
* A fit that ends at a boundary under `"hold"` now says so as it is fitted --
  with `verbose = FALSE` nothing else would -- that its fixed effects are
  still usable, and names `"avoid"` with what it costs. `summary()`,
  `ilm_interpret()` and the check table name it the same way.
* `logLik()` of a penalised fit is the likelihood of the data at its
  estimate, not the penalised objective, and every refit -- a reduced model,
  a bootstrap replicate, a simulation envelope -- keeps the penalty.

## Refits are the same model

* Eight places refit a model -- a reduced model for a test, a bootstrap
  replicate, a simulation envelope, a consistency check, the null model behind
  an R-squared -- and seven of them wrote out their own list of what "the same
  model" means. Every one of those lists had gone out of date, and none carried
  the **zero part** of a zero-inflated or hurdle model, so
  `ilm_anova(test = "LRT")` compared it
  against a model without one: a predictor with no effect came out at
  chi-square **68.8, p < 2e-16**, where the right answer was 0.02. Type II
  tests of a term inside an interaction, `ilm_consistency()` and the
  simulation-based checks were affected the same way, and `ilm_pb_lrt()`'s
  bootstrap refits also lacked censoring, a dispersion model and a flexible
  survival baseline.
* All of them are now built from one list, next to the one function that reads
  it, and a test refits every kind of model on its own response and requires
  its own likelihood back.
* Refits of a REML fit are REML too. They were maximum likelihood, so a Type
  III recode or a consistency check quietly used a different estimator.
* McFadden's R-squared is NA for a REML fit, with a note: restricted
  likelihoods of models with different fixed effects are not comparable.
  `model_performance()` also failed outright on every fit whose outcome is not
  a set of categories, on "subscript out of bounds" inside the scoring rules;
  those scores are now NA there, and a 0/1 binomial fit is scored as the two
  categories it is.

## The flagship model, category by category

* `ilm_ame()` on a multinomial or ordinal fit took the LAST column of the
  predicted probabilities and reported it as the effect, unlabelled. It now
  gives one row per category, in a `category` column. They sum to zero, and
  the slopes agree with `marginaleffects::avg_slopes()` on
  `nnet::multinom()` to 1e-3, the level at which the two fits agree.
* `ilm_interpret()` described only the first category's coefficients of a
  multinomial fit, and quoted that last-category effect beside each. It now
  describes every category, in the language of a categorical outcome, each
  with its own effect in percentage points.
* `ilm_emmeans()` failed on every multinomial fit with an error blaming "a
  smooth or a matrix column". It now says it does not yet average a
  multinomial fit, and names `predict()` and `ilm_ame()` for what it would
  have given.

## Column names that need backticks

* ``ilm_model(`my y` ~ `x 1` + (1 | `site id`))`` failed with "unexpected
  symbol". Names turned into text and parsed back lose their backticks --
  in the formula front end, in mgcv, in the prediction code and in half a
  dozen functions that take column names as strings. Each place now quotes a
  name when it becomes code and not when it is looked up. mgcv cannot take
  such a name at all, so a smooth of one is built on a stand-in and evaluated
  through the same map.
* Checked against the only oracle that matters here, the same data under
  ordinary names: identical likelihoods, coefficients and standard errors for
  random intercepts and slopes, nested terms, smooths, a dispersion model, a
  zero part and weights, and identical results from `ilm_iv()`, `ilm_did()`,
  `ilm_rdd()`, `ilm_aov_ez()`, `ilm_impute()`, `ilm_check_missing()`,
  `ilm_moderation()`, `ilm_pb_lrt()` and `ilm_power()`.
* Functions that take a term by name take it the way a person writes it,
  `"x 1"`, as well as the way R does.
* Two defects found on the way, both independent of names. A smooth's `by`
  variable never reached the model frame unless it appeared elsewhere in the
  formula, and mgcv stopped with "Can't find by variable"; `s(x, by = z)` now
  matches `mgcv::gam()`. And a random slope's design on new data was built
  from text in a `tryCatch()`, so a failure dropped the slope from the draws
  without a word.

## Seeing an anomaly scan

* `ilm_plot_anomaly()`: `"scores"`, the default, plots every row's score
  against its rank beside the band the scan simulated. It is the one view that
  tells five genuine outliers from the top 5% of a smooth continuum.
  `"drivers"` counts which column drives the flags, `"map"` places the rows on
  the first two dimensions of `ilm_reduce()`, and `"row"` shows one row's
  z-scores beside its residuals.
* The band **restarts at the line**. Setting the flagged rows aside moves every
  other row up that many ranks, so against the band as simulated a clean
  remainder sits above it for a long stretch: over 60 scans, 0.69 of the ranks
  after the line with 2% planted anomalies, against 0.82 to 1.00 for noise that
  really has heavy tails.
* The verdict is a reading, not a test, and its rates are in the help page:
  planted anomalies were said to stand clear in 89 and 97 scans of 100 and
  never to run on; t-tailed noise was said to run on in half to two thirds of
  the scans that flagged anything. Nothing is read off the band when nothing is
  flagged -- in data with no anomalies the top twenty scores sat half above it
  in 21 scans of 100.
* An isolation forest gets no band and no row view, and says why.

## Smaller fixes

* `summary()` called a REML fit "maximum likelihood" and a Poisson or
  multinomial mixed model a "linear mixed model". `print()` called a gaussian
  fit "2 categories", and never showed `[CHECKS FAILED]` for an ordinal one.
* `ilm_moderation()` refitted with `all.vars(formula)[1]` as the response, so
  `log(y) ~ ...` was refitted on `y`.
* `ilm_aov_ez()` matched an `observed` factor as a regular expression.
* Asking for censoring on a family without a censored form said only the
  gaussian had one; the survival families have one too, and it says so.
* `ilm_scores()` and `ilm_calibration()` take an ordinal fit, which predicts
  category probabilities as a multinomial one does, and `ilm_scores()` adds
  the ranked probability score for it, which charges a near miss less than a
  far one.
* `ilm_pb_lrt()` called the chi-square reference "MISCALIBRATED" on however
  few replicates were run -- five were enough. The verdict now needs the
  Monte Carlo interval on the implied size to exclude 0.05, and says when
  there are too few replicates to judge.
* The dispersion and zero-inflation checks gave a reason for declining a
  categorical family that was wrong for it ("the residual variance is
  estimated"); they now point to the check that does apply.
* `set_coef()` takes the fixed effects alone, as `coef()` returns them. It
  used to stop on an unexplained "'names' attribute" error.
* Refits inside a simulation no longer pass each replicate's optimiser
  warnings through; the replicates that fail are counted, as before.
* The pkgdown reference index is generated, and hand edits to `_pkgdown.yml`
  had drifted from the generator; they are in the generator now.
  `ilm_anomalous()` is indexed, and the `benchmarking` and `moderation`
  articles are listed, without which a pkgdown build stops.
* The benchmarking vignette said 1000 replications where the agreement study
  ran 160 (1000 was each dataset's size), gave mclogit's RMSE advantage as 8 to
  13% where it is 3 to 14%, and its reproduction commands passed the wrong
  arguments to two scripts. The messy-data findings are now generated by
  `summarise_run.R` like every other study's, and its run lives in
  `studies/runs/0.0.7.9000/` with the rest.
* Random-effect bars are parsed with reformulas, where lme4 2.0 moved its
  parser; reached through lme4 2.0 it warned on every fit. reformulas is
  suggested, and an older lme4 still serves without it.
* The GitHub Actions workflows use the Node 24 releases of `actions/checkout`
  and the Pages deploy action, which the runners had started warning about.
* illume has a hex sticker, shown in the README and on the pkgdown site, which
  also has favicons made from it.

## iml_*(), for the transposition

* Every exported `ilm_*()` also answers to `iml_*()`, an easy transposition to
  type. Each alias is the function itself rather than a wrapper, so arguments,
  defaults, autocompletion and help pages are the same, and a test fails when
  a new export has none.

## Estimated marginal means and contrasts

* `ilm_emmeans()` and `ilm_contrast()`. Marginal means are exact linear
  combinations `L %*% beta` with variance `L V L'`, computed on the link scale
  where that is exact. Pairwise, against-a-control and polynomial contrasts,
  with a simultaneous adjustment.
* Against **emmeans** on the same fit: means to 2.0e-06 and standard errors to
  3.5e-08 across all three weighting schemes; pairwise contrasts to 3.6e-06;
  unadjusted and Bonferroni p-values to 4.3e-06 and 2.6e-05.
* The joint adjustment simulates the maximum of a multivariate **t**, not a
  normal. Each contrast is divided by an estimated standard error, so the
  reference carries that estimate's uncertainty too. At n = 30 a normal
  reference gives intervals that are too narrow.
* What was called `"proportional"` weighting was the joint cell frequency,
  which is emmeans' `"cells"` and a different estimand: it averages each group
  over ITS OWN mix of the other variables, so a difference carries composition
  as well as effect. On a design where the two factors are associated that
  moved a contrast by more than 0.5 against a coefficient of 0.4. Proportional
  is now the product of the one-way margins; `"cells"` remains available, and
  `ilm_contrast()` warns when asked to difference them.

## Zero-inflated and hurdle counts

* `ilm_model(ziformula = , zi_type = )`. A mixture and a hurdle are different
  models: under a mixture a zero has two possible origins and `p` is the
  structural share, under a hurdle `p` is every zero there is. Which applies is
  a question about the subject rather than about fit, though AIC separates them.
* Against **pscl**: count coefficients to 4.4e-07, zero coefficients to
  1.5e-06, zero standard errors to 1.3e-07, and equal log-likelihoods for
  zero-inflated Poisson, hurdle Poisson and zero-inflated negative binomial.
* **The residuals had to learn about it, and skipping that is silent.** With the
  zero part in the quantile residuals, 30 fits gave a median KS p of 0.80 and
  nothing below 0.05; scoring the same fits against the count part alone gave a
  median of 1.4e-38 and every one below 0.05.
* `ilm_sim_cond()`, which five diagnostics build their reference from, had to
  learn it too -- otherwise `ilm_check_zeros()` compares a zero-inflated fit
  against draws from the count part and reports the inflation it was told about
  as a failure. It did: FAIL at 1482 observed zeros against 484 expected. Now
  1042 against 1041, and 0 of 60 false alarms.

## Ordered outcomes

* `family = "ordinal"`, `"ordinal_probit"` and `"ordinal_cloglog"`: cumulative
  link models, with the thresholds fitted as a first value and log increments
  so the ordering holds by construction.
* Against **ordinal::clm** on all three links: coefficients, thresholds and both
  sets of standard errors to 1e-7, log-likelihood to 1e-10. The mixed model
  matches **ordinal::clmm** exactly.
* The intercept is dropped as a COLUMN, after `model.matrix()` has coded the
  factors. A `- 1` in the formula instead makes `model.matrix()` expand the
  first factor to all its levels: the same likelihood, a singular Hessian, and
  every standard error `NaN`.
* `ilm_check_proportional()` tests the assumption that makes one coefficient
  per predictor enough. Its reference is **simulated, not chi-squared** -- the
  per-cut fits share their data, and a chi-squared reference flagged a
  perfectly proportional binary predictor 31.5% of the time. Simulating gives
  0.032 to 0.072 against a nominal 0.05, p-values essentially uniform, and 400
  of 400 detections of a genuine violation.

## Type II by default

* `ilm_anova()` now defaults to `type = 2`, which does not depend on how the
  factors are coded and coincides with Type III whenever there is no
  interaction.
* `type = 3` refits with `contr.sum` for the offending factors -- and **centres
  uncentred numerics in interactions**, which the warning had always claimed
  mattered and had never looked for -- then says so, as afex does. The fit
  passed in is untouched. `recode = FALSE` keeps the old behaviour.
* On `y ~ g * h + x`, `g` reads F = 4.55 under treatment coding and F = 7.66
  under sum coding: the first tests it at `h = "p"`, only the second averages
  over `h`. On `y ~ x * z` with both uncentred, `x` reads F = 187 against 642
  centred.
* `contrasts` also takes a single string now (`"sum"`, `"treatment"`,
  `"helmert"`, `"poly"`) applying to every unordered factor. Ordered factors
  keep their polynomial coding, which is a statement about spacing rather than
  a default.

## Cluster-robust standard errors

* `ilm_vcov_cluster()` and `ilm_robust()`, with CR0, CR1 and CR2 and a t
  reference on Bell-McCaffrey degrees of freedom.
* Against **sandwich**: HC0/HC1/HC2 unclustered to 1.3e-09, CR0/CR1 clustered to
  1.0e-08 for gaussian, 3.0e-08 for poisson and 1.1e-06 for binomial.
* The defaults are what cover. For a cluster-level predictor, over 600
  replicates:

  ```
                         G=10    G=20    G=40   G=20 unbalanced
    model-based         0.402   0.427   0.427   0.440
    CR0 + normal        0.798   0.877   0.932   0.897
    CR1 + t(G-1)        0.860   0.907   0.940   0.923
    CR2 + t(G-1)        0.912   0.922   0.940   0.937
    CR2 + Bell-McCaffrey 0.950  0.948   0.952   0.965
  ```

  Bell-McCaffrey degrees of freedom come out far below `G - 1`: median 4.0
  where `G - 1` is 9, and 15.2 where it is 39.

## Beta regression, and a zero part for it

* `family = "beta"`, parameterised by a mean and a PRECISION -- larger means
  less spread, the opposite of every other dispersion parameter here.
  `dispformula` models the precision.
* Against **glmmTMB**: coefficients to 2.4e-06, standard errors to 1.1e-08,
  log-likelihood to 1.5e-09, and the same with a precision model.
* Boundary values have no likelihood, and the error names both honest
  responses: a separate process, which `ziformula` now really does model for a
  beta response, or rounding, which `ilm_squeeze()` handles by the
  Smithson-Verkuilen shift while reporting how far it moved things.
* Zero-inflated beta against glmmTMB: mean coefficients to 9.9e-07, zero
  coefficients to 5.2e-06, log-likelihood to 3.1e-09. It is a hurdle and only a
  hurdle -- a continuous density has no probability of producing a zero, so a
  mixture has nothing to add, and asking for one is an error rather than a
  worse approximation.

## Instrumental variables

* `ilm_iv(y ~ x + w | z + w)`, matching a closed-form 2SLS to 1.1e-14 on
  coefficients and 9.1e-16 on standard errors, with bias falling from +0.0034
  at n = 800 to +0.0003 at n = 50,000 and 0.955 coverage.
* Running `lm()` twice gives the same point estimate and the wrong standard
  errors, in a direction that depends on the sign of the coefficient and of the
  endogeneity -- 24% too LARGE in the documented example.
* The first-stage F is reported against both thresholds: the familiar 10, and
  the 104.7 a conventional 5% t-test actually needs (Lee et al. 2022).
  Durbin-Wu-Hausman and Sargan's J print unasked.
* `ilm_iv_ar()` is the remedy the F diagnoses: an Anderson-Rubin set that stays
  valid however weak the instrument is. With F = 0.235 the Wald interval reports
  a tidy [-2.62, 5.82] and the AR set is unbounded.

## Multivariate anomaly detection

* `ilm_anomaly()` finds rows implausible as a COMBINATION. On test data with ten
  rows pushed off the correlation structure, the column-at-a-time scan caught
  none of them.
* The rank comes from parallel analysis, **not** the cross-validation
  `ilm_impute()` uses: that chose 6 or 7 on a rank-2 structure, on clean data as
  well as contaminated, and detection fell from 0.975 to 0.560.
* The fit is trimmed rather than held out in folds. Folds did nothing --
  matching in-sample scoring to three decimals at four contamination levels --
  because most anomalies remain in every training fold.
* The reference is simulated with PER-COLUMN noise, calibrated to the observed
  median. Getting either wrong flagged **38.9%** of the rows of clean data. As
  it stands, 11 rows in 40,000 across 100 clean datasets, and 98% detection at
  2% contamination.

## Generalized low rank models

* `ilm_glrm()`, and `method = "glrm"` on `ilm_reduce()`, `ilm_profile()` and
  `ilm_impute()`. A loss per column type rather than squared error on one-hot
  indicators, so a category is reconstructed AS a category.
* With quadratic loss and no penalty it IS principal components, verified
  against the SVD to 1e-16 with principal angles of zero.
* `ilm_impute(method = "glrm")` imputes CATEGORICAL columns, which the low-rank
  route says plainly it cannot. A category is sampled from its fitted
  probabilities rather than set to the most likely level.
* The ridge penalty is chosen by cross-validation, because a fixed one is wrong
  across sizes: held-out error ran 0.82 at 0.1, 0.68 at 2 and 1.16 at 25,
  against 0.78 for an iterative SVD and 1.18 for column means.

## Survey weights

* `ilm_design()` and `ilm_svy_coef()`: Taylor linearization with clustering,
  stratification and finite population correction.
* Against **survey::svyglm**: coefficients to 8.1e-09, standard errors to
  1.6e-10, and the same design degrees of freedom -- PSUs less strata, 56 rather
  than 720 rows.
* Sampling weights handed to a model that treats weights as frequencies give
  standard errors too small by at least `sqrt(n / sum(w))`, and by more once
  clustered: 7.3 and 16.6 times in the documented example.
* The weights check gained two things it was missing. It only flagged
  NON-INTEGER weights, and sampling weights are routinely rounded; weights
  totalling more than ten times the number of rows are now flagged too. And it
  flagged constant weights above 1, which for a binomial response are trial
  counts and entirely ordinary; neither rule applies to binomial now.

## Causal mediation

* `ilm_mediate()` estimates the counterfactual ACME and ADE, reproducing the
  Baron-Kenny product exactly where the product is right (to 0.0004) and
  differing where it is not: with an interaction the two ACMEs are 0.416 and
  0.782, and a single product cannot be both.
* `ilm_mediate_sens()` addresses the assumption nothing can test. Over 30
  replicates the predicted bias matched the realised one to three decimals at
  every confounder strength, recovering a true ACME of 0.42 as 0.4205, 0.4188
  and 0.4167 from observed values of 0.52 to 0.82.

## Effect sizes, scenarios and power

* `ilm_effects()` reports each coefficient on its family's own scale -- odds
  ratio, incidence rate ratio, time ratio, hazard ratio, proportional odds
  ratio -- with intervals built on the link scale and transformed. Odds ratios
  match `glm()` to 1e-6 and their intervals `confint.default()` to 1e-4.
  A ratio from a MIXED model is labelled conditional, with `ilm_ame()` named for
  the marginal one.
* `ilm_scenario()` projects named settings. The default standardises over the
  observed units, which is not the same as predicting for a unit at the average
  covariate -- 0.657 against 0.683 on a logistic fit -- and says which it used.
  It flags extrapolation by value AND by combination: a 25-year-old with 35
  years of service has both values in range and sits 2.39 standardised units
  from the nearest real person against 0.25 for a typical one.
* `ilm_power()` simulates studies at each size and effect. Against the closed
  form for a linear model, 0.450/0.737/0.956/0.999 predicted against
  0.450/0.762/0.952/0.997 simulated, and a required n of 240 (207 to 272)
  against 234. It reports the **Monte Carlo interval**, because 0.80 from 200
  replicates is 0.74 to 0.86, and the **convergence rate**, because power
  conditional on convergence is not power.

## Two defects, found from outside

* **`ilm_check_covariate()` and `ilm_check_omitted()` failed on any missing
  data.** The fit drops rows through `na.omit`, so its residuals are shorter
  than the data frame a covariate is passed from, and both died with
  `arguments must have same length` -- an error that names no remedy and does
  not say which two lengths disagreed. They now align the variable to the rows
  the fit kept, which the fit already records. This mattered more than a
  crash: the package's own documentation argues that targeted covariate checks
  are the only residual checks with real power here (z = 4.04 against z = 0.24
  for a broad one), and they were unreachable for anyone whose data had a gap
  in it.
* **Three vignettes called `ilm_check_missing(d, y ~ x + z)` while the
  signature took a column name.** Every chunk was `eval = FALSE`, so nothing
  caught it. Rather than correct three call sites, the function now takes a
  formula as well: `outcome ~ predictors` says exactly what it asks, and three
  vignettes reaching for it independently is reasonable evidence about what
  the call looks like to someone who has not read the signature.

Both were found by calling illume from outside it, while building a separate
package against it. Neither was found by this package's own test suite, which
is the same story every defect here has had.

The regression test for the first one checks that automatic alignment gives
the **same answer** as subsetting by hand, not merely that it stops erroring:
a wrong alignment also runs cleanly, and dropping the last rows rather than
the right ones turns z = 23.5 into z = 0.41.

* **`ilm_check_ar()` and `ilm_variogram()` had the same defect and did not get
  the same fix.** They take `time` and `group` straight from the data frame
  while the fit kept only the complete rows, so every call on data with a gap
  in it stopped. On a 320-row example with 21 missing outcomes both were
  unusable. They now align the same way, and `ilm_plot_acf()` picks it up too
  by going through `ilm_ar_envelope()`.
* The alignment is attempted and falls back to the column as given, rather
  than letting `ilm_align_rows()` raise. A length that matches neither the
  data nor the fit therefore still produces the original message, which two
  existing tests pin deliberately -- fixing one error by replacing another
  one's wording is not a fix.

## A tibble is a data frame and has to behave like one

* **`ilm_reduce()` and `ilm_profile()` failed on any tibble with a numeric
  column**, which is to say on most real data: on `gapminder::gapminder` they
  stopped with `All variables in X.quanti must be numeric`. `ilm_reduce()`
  handed its numeric half straight to `PCAmixdata::PCAmix()`, which checks
  columns with `is.numeric(X.quanti[, j])` -- and `[` on a tibble does not drop
  to a vector, so every numeric column looked non-numeric. The categorical
  half was already coerced, which is why only the numeric one broke.
* Nothing in that message mentions tibbles, and a tibble is what anyone gets
  from readr, dplyr or gapminder, so this closed the most ordinary route into
  the function while looking like a complaint about the data.
* The mixed-data handling itself was never the problem. `ilm_reduce()` picks
  PCA, MCA or FAMD by column type and `gapminder` takes the FAMD branch, as
  intended, once the frame reaches PCAmix in a shape it accepts.
* The rest of the package was swept for the same fault and is clean:
  `ilm_model()`, `ilm_aov_ez()`, `ilm_describe_all()`, `ilm_anomaly()`,
  `ilm_impute()`, `ilm_glrm()`, `ilm_check_missing()` and the plotting
  functions all give a tibble and a data frame the same answer.
* **`ilm_cluster()` refused mixed columns without naming the remedy.** It
  clusters coordinates -- a k-means centroid is not defined on a factor -- so
  refusing is right, but `the coordinates to cluster must be numeric` left the
  caller to discover `ilm_reduce()` for themselves, which is not how anything
  else in the package behaves. It now names the offending columns and the
  call that fixes it, as `ilm_anomaly()` already did for the same situation.

The regression test builds its tibble with `tibble::as_tibble()`. Setting the
class by hand does **not** reproduce this: `[.tbl_df` is only dispatched to
when tibble's namespace is loaded, so without it `[` falls through to
`[.data.frame`, which drops correctly and the bug vanishes. The first version
of the test did it by hand and passed against the broken code.

## Does the effect hold for everyone?

`ilm_moderation()` searches a set of candidate moderators for evidence that the
effect of a treatment varies, and charges for the search.

* **Trying several moderators and reporting the strongest is a search, and the
  p-value from a search is not the p-value from a test.** With six candidates
  and no moderation present at all, that procedure rejected **29.3%** of the
  time at a nominal 5%. Every adjustment brings it back: Holm, BH and
  Bonferroni all landed at 0.027.
* **The test is the JOINT test of the interaction block**, from `ilm_anova()`,
  with its degrees of freedom reported. A single interaction coefficient tests
  moderation only when the moderator has 1 df -- on a four-level factor the
  three coefficients gave p = 0.62, 0.0018 and 0.00002 for the same variable.
* `x` must be named and has no default. `x:m` is the same term whichever of
  the two is called the moderator; only the interpretation distinguishes them.
  An `ilm_dag_model()` is the exception, since the graph has named the
  exposure already.
* **An interaction already in the model is excluded from the search**, because
  it was specified a priori and owes no multiplicity penalty. Charging it one
  would make a pre-registered hypothesis weaker for having been tested beside
  exploratory ones. `ilm_anova()` is named as where it belongs.
* `adjust` passes straight to `stats::p.adjust()`, so every method it supports
  is available; Holm is the default because it controls the family-wise rate
  and dominates Bonferroni.
* Each candidate is refit **through the fit's own call**, so family,
  `ziformula`, `dispformula`, `ar`, weights, contrasts and `reml` all survive.
  Rebuilding the formula by hand drops every one of them silently -- a
  zero-inflated model would be tested without its zero part and nothing would
  say so. Candidates are tested on equal footing: each refit starts from the
  user's model minus any *other* candidate's interaction.
* `split = TRUE` picks the moderator on half the data and tests it on the
  other half, needing no adjustment. It is better calibrated -- 0.053 against
  Holm's conservative 0.027 -- and costs about half the power (0.320 against
  0.547 for a numeric moderator, 0.200 against 0.380 for a three-level
  factor). **The gap widens with degrees of freedom**, 1.71x to 1.90x, so
  splitting is weakest precisely in the three-level-factor case experimental
  work is full of. It is genuinely required only when the hypotheses cannot be
  counted, as in an open-ended tree search; this searches a named set.

`ilm_plot_moderation()` draws one: cell means for a categorical exposure, and
the **slope** of a continuous one within each level of the moderator, from
`ilm_trends()` -- so the picture and the p-value are the same estimator on the
same fit. It is deliberately not `ilm_plot_model()`'s effect plot, which holds
the other predictors at typical values; pinning the moderator is exactly what
would hide the moderation.

A 0/1 treatment stored as a **number** is categorical in meaning and numeric in
type, so `ilm_emmeans()` would hold it at its mean and draw one curve where two
are wanted. Its values are named explicitly.

One bug worth recording because the tests missed it: the plot builds its call
with `do.call()`, and tinyplot deparses its arguments to title a legend, so a
`by` factor arrives as a deparsed vector and the width computation throws
`invalid graphics state`. Whether it fires depends on how long the level NAMES
are -- a test using `"a"`/`"b"`/`"c"` passed while `"north"`/`"central"`/
`"south"` did not. Fixed by titling the legend explicitly, and the test now
uses realistic labels. This is the second time this exact mechanism has bitten.

## Mixed data, measured rather than assumed

`ilm_cluster()` and `ilm_anomaly()` both refused to use categorical columns.
Deciding what to do about that meant comparing seven approaches against data
with known structure, and three of the results overturned the plan they were
meant to confirm.

* **`ilm_cluster(data)` now takes raw mixed data**, reducing it through
  `ilm_reduce()` (FAMD) and saying so. `ilm_profile()` remains the front door:
  it is the same pipeline plus a description of each cluster, in one call.
* **FAMD then k-means is the default because it is the only one that survives
  correlated variables.** Adjusted Rand index against three known clusters:

  ```
                              FAMD+kmeans   Gower+PAM   RF proximity   VarSelLCM
    balanced                     0.433        0.447        0.257         0.440
    within-cluster r = 0.40      0.395          --           --          0.286
    within-cluster r = 0.75      0.265          --           --          0.047
    two irrelevant factors       0.006          --         0.093         0.397
  ```

  VarSelLCM is remarkable in its own regime -- it selects variables, and on the
  irrelevant-factor case it scored 0.397 where everything else collapsed -- and
  it assumes the variables are independent within a cluster. When that fails it
  fails the other way, to 0.047. Real variables are correlated, so it is not
  the default. Random-forest proximity was measured and is not offered: 0.257
  against 0.433 on well-behaved data.
* **Gaussian mixtures are not suitable here, for a structural reason.** On
  elliptical clusters a GMM beats k-means on the numeric columns alone (0.547
  against 0.452) and collapses to 0.253 on the FAMD coordinates, where k-means
  reaches 0.620. The categorical part of an embedding takes only as many
  distinct values as there are level combinations, so a Gaussian mixture spends
  components modelling the lattice. Its strength is what breaks it.
* **Every `k`-selector over-selects on mixed data**, and `ilm_cluster()` now
  warns when the chosen `k` is the largest searched -- the curve had not
  turned, so that is where the search stopped rather than where the evidence
  pointed. With truth `k = 3`, the gap statistic recovered it in 25% of runs at
  best, BIC in 0%, and average silhouette in 0% -- silhouette on a Gower
  dissimilarity peaked at exactly the number of level combinations, finding the
  lattice rather than the clusters.
* **Both functions now report what the missing values will cost before they
  cost it**, through `ilm_describe_na_all()`. Twelve per cent missing across
  three columns leaves 68% of rows complete, and the other 32% are dropped
  silently by any complete-case method. The note names `ilm_impute()` and says
  which single column, if dropped, would recover the most rows.

## Anomalies that live in the categories

* **`ilm_anomaly(method = "iforest")`**, an isolation forest over every column.
  AUC against planted anomalies:

  ```
                            reconstruction   iForest   MCD     LOF+Gower
    extreme value                1.000        0.975    1.000     0.998
    implausible combination      0.999        0.950    1.000     0.997
    category contradicts numbers 0.771        0.961    0.937     0.990
    category pairing never seen  0.505        0.886    0.487     0.363
  ```

  The reconstruction stays the default: unbeaten at what it was built for, and
  the only method with a calibrated null and FDR-adjusted p-values. It is a
  coin toss on the bottom row, which is what the forest is for.
* `ndim = 1` deliberately. The extended isolation forest is better on numeric
  combinations (0.965 against 0.950) and worse on rare category pairings (0.796
  against 0.886) -- and the categories are the reason it is here.
* A forest returns a score with no null behind it, so `p` and `p_adj` are `NA`
  and `alpha` is the share of rows being called anomalous rather than an error
  rate being controlled. The print method says so rather than letting two
  columns of `NA` imply an oversight.
* It keeps the `driver`. Replacing one column at a time with its median or
  commonest level and taking the largest score drop recovered the planted
  driver 100% of the time for an extreme value and 94% for a category
  contradicting the numbers, at one extra prediction per column.
* Not added: MCD (ties the reconstruction on numeric anomalies, blind to
  categorical ones, loses the calibration) and LOF, which is **worse than
  random** on rare category pairings at 0.363 -- a rare combination forms its
  own small tight group, which looks locally dense.

## Which variables are carrying a clustering

`ilm_var_contrib()`. Nothing in the package selects variables, and an
irrelevant one is not neutral: on three known clusters with two informative
columns, adding two pure-noise factors took recovery from 0.301 to **0.006**.
That is the largest single effect measured anywhere in this comparison --
larger than the choice of method, of distance, or of `k`.

It reports the between-cluster share of variance for each variable against a
permutation reference, so a variable the clustering ignored can be seen and
dropped.

The interesting case is the other one. When a variable separates the clusters
*almost perfectly* while nothing else does, the partition **is** its levels
under another name -- and every informative variable then scores low against
it, so the obvious reading of the table ("drop the ones at the bottom") is
exactly backwards. That case is detected and reported separately, because from
inside a clustering a variable that defines it looks like the best variable.
The first draft of this function gave the backwards advice.

What it cannot do is settle relevance: the clustering was fitted to these
variables, so one it used separates the clusters it helped make whether or not
it means anything. A pure-noise column scored 0.441 against 0.576 for a real
one. The documentation says this rather than implying otherwise.

**`ilm_profile()` runs it, rather than leaving it for the user to find.** A
diagnostic in a function nobody calls produces exactly the analyses this
package exists to prevent, so the check is in the path: the result is stored
on the object, the verdict is printed under the cluster descriptions, and the
dominance case additionally **warns**, because otherwise the cluster
descriptions get read at face value and they are all true and all about one
variable. `var_contrib = FALSE` turns it off and `var_contrib_B` tunes it; it
costs about 70% on top of the clustering at 500 rows and five columns, which
is the price of the check being run rather than merely available.

On thirty columns with two informative ones, the two rank first at 0.45 and
0.40 with everything else at 0.13 or below. The verdicts are less useful there
than the ranking -- a noise column the clustering happened to split on still
beats a shuffled label -- which is why the output says to read the ranking.

## Parameters that were set and then ignored

Same source, same method, one round later.

* **`set_coef()` rebuilt the mean and the covariance structure and nothing
  else.** Setting a dispersion changed the stored parameter vector and left
  `$dispersion` at its fitted value, so `ilm_simulate()` went on drawing from
  the fitted residual SD: `log_sigma` set to `log(4)` produced draws with SD
  1.00. Nothing errored, and the number that came back was a perfectly
  ordinary standard deviation. The same gap left ordinal thresholds,
  zero-inflation parameters and a dispersion model stale, and it reached
  further than `ilm_simulate()` -- `ilm_power()` draws its responses through
  the same field, so a power curve run on an object with assumed parameters
  used the fitted dispersion instead of the assumed one.
* This is the failure the function's own documentation warns about, in the
  places it was not actually guarded against: `marginaleffects` works by
  nudging a parameter and re-predicting, so a parameter the rebuild ignores
  contributes exactly zero uncertainty rather than an error.
* **A rebuilt CAR(1) correlation used AR(1)'s transform.** AR(1) stores
  `atanh(rho)`; CAR(1) stores `log(range)` and `rho` is `exp(-1/range)`. Both
  return something inside (-1, 1), so the wrong one is invisible in the value:
  on the test fit it gave 0.703 where the fit reported 0.659, and `$ar_range`
  was not updated at all.
* `get_coef()` deliberately returns the variance parameters beside the fixed
  effects, so setting one and having it silently ignored was the worst of both
  conventions.

Every test in `test-set-coef.R` sets a parameter to a value whose correct
consequence is known in closed form, because "it changed" is not "it changed
to the right thing". Nine of them fail against the previous code.

* **`ilm_power(term = )` now takes a variable name as well as a coefficient
  name.** A researcher planning a trial thinks `"arm"`, not `"armtreatment"`.
  The variable is resolved through the model matrix's `assign` attribute, not
  by matching a name prefix, so an exposure `x` does not collect a covariate
  called `xray`. A term with more than one column -- a three-level factor, a
  spline -- is refused with its coefficients listed and `ilm_anova()` named,
  because a power curve follows one effect size at a time and choosing a
  column silently would report power for a comparison nobody asked for.

## Planning a study that has not been run

`ilm_scaffold()` builds a design grid from a study specification, attaches the
parameters you assume, and returns an ordinary `"ilm_model"`. Everything
downstream then works unchanged, which is the point: the assumptions can be
interrogated before anything is built on them. `ilm_power_design()` is the
same thing plus a power curve, in one call.

* Assumptions go in **either** as coefficients on the link scale, or as
  expected **cell means on the response scale** -- "controls average 12, the
  treated group averages 14.5" -- which is the form most people actually hold
  one in. illume solves back to coefficients.
* A least-squares solve always returns something, so cell means that the
  formula cannot produce would otherwise give a scaffold for a different study
  than the one being planned. Crossed means under an additive formula stop,
  name the cell that is off and by how much, and name `arm * time` as the
  remedy.
* `n_unit` counts participants when the formula has a grouping bar and rows
  when it does not, so the number given is the number that goes in a protocol.
  A between-unit variable is allocated to UNITS and repeated down their rows;
  allocating it to rows would put a participant in both arms.
* `icc` is accepted in place of `re_sd`, since an intraclass correlation is
  what a repeated-measures literature reports.
* Against the closed form for a two-sample **t** test: 0.858/0.995/1.000
  simulated against 0.872/0.993/1.000 predicted at 100/200/300 participants.
  Against a treatment-by-time mixed model hand-rolled in plain R and fitted
  with **lme4/lmerTest**, at 40/80/140 participants: 0.284/0.452/0.654 against
  0.242/0.440/0.632. All within Monte Carlo error, but all three lean the same
  way, which is what a **z** reference against a Satterthwaite **t** does.
* `ilm_emmeans()` on a scaffold returns the cell means that were assumed, to
  1e-6 -- the solve goes one way and that comes back the other, so the check
  does not share an error with the thing it checks.
* The print method says the parameters were assumed, not estimated, and that
  the standard errors are one realisation of the design rather than a property
  of the assumptions. An `"ilm_model"` that was never fitted to anything is
  exactly the sort of object that gets mistaken for a result.

## A random slope was being averaged as a random intercept

In both places it could be. Found while building the above, since a
treatment-by-time design is what that is mostly for.

* **`ilm_power()` drew only the intercept variance.** The draw took
  `Sigma[[k]][1, 1]` and added it as a per-group shift, so a fit with
  `(1 + time | id)` produced simulated participants whose individual
  trajectories varied no more than residual noise allows: a per-subject slope
  SD of 0.50 where the fit said 1.14 and the generating value was 1.20. Every
  simulated participant moved in near-parallel, which understates how much
  they differ and so **overstates** power for anything interacting with the
  within-subject variable -- 0.450 against 0.333 on the test case, a third too
  high.
* The scale of a term lives in `Sigma` and its shape in `Sigma_d`, so the
  covariance of a group's whole random-effect vector is their product. The
  draw now uses it, with the bar re-evaluated on the simulated study's own
  rows. Drawn data now has a slope SD of 1.248 against the real data's 1.250.
* **`predict(marginal = TRUE)` made the same simplification, and is fixed too.**
  It drew `Sigma` alone and added it to every row, which drops the slope
  variance and also holds the shift constant across rows when the spread being
  integrated over grows with distance from where the slope is centred. On a
  logistic model with a fitted slope SD of 0.81 the population-averaged
  probability at the far end of the range came back as 0.703 where the correct
  average is 0.618.
* The draw is now the matrix normal the objective actually specifies: one
  group's effect is a `dk x C` matrix with row covariance `Sigma_d` and column
  covariance `Sigma`, so `U = A Z B` with `AA' = Sigma_d` and `B'B = Sigma`,
  and the contribution at a row is `z_row %*% U`. `ilm_power()` and
  `predict()` now read that factorisation from one place, having already gone
  out of step once.
* When `dk == 1` this is the previous draw exactly, down to the random number
  stream, so nothing that was already right moved.
* Checked against a brute-force average over 400,000 draws computed
  independently, at every point of the prediction grid, for a random intercept
  and for a random slope.
* **Under an identity link the average is now exact rather than simulated.**
  `E[eta + z'u] = eta`, so there was a closed form all along and simulating it
  returned a noisy estimate of a number already known -- 0.14 on the response
  scale with 200 draws and a random slope. The result no longer depends on
  `ndraw`.
* A bar that varies over a column the prediction data does not carry now
  warns and says it integrated the intercept only. Averaging over part of a
  term is a different quantity from averaging over the term.

## Plotting characters by name

* `pch = "filled circle"` works wherever `pch = 16` did. Base R's `pch` is 26
  integers nobody remembers, nothing in the argument says which is which, and
  the difference between 16, 19, 20 and 21 is not guessable -- so a plot gets
  whichever code the author recalled and a reader comparing two plots cannot
  tell whether a difference in the markers was meant.
* Names are also checkable in a way numbers are not: `pch = 26` is accepted by
  `graphics` and quietly draws nothing, while a name that is not in the table
  stops and lists the ones that are. All 26 codes have a name, most have
  aliases, and case, spaces, underscores, hyphens and dots are all ignored.
* A single character is left alone, because base R draws it literally --
  translating `pch = "x"` would silently turn the plot into crosses.
* `pch` is now a documented argument of `ilm_plot()`, `ilm_plot_scatter()`,
  `ilm_plot_line()`, `ilm_plot_stat_error()`, `ilm_plot_box()` and
  `ilm_plot_violin()` rather than something to be discovered inside `...`. It
  sits after `...`, so no existing positional call can be matched to it.

The obvious implementation -- route every call through one `do.call()` funnel
so the translation lives in a single place -- **is wrong here, and passed the
whole test suite before `R CMD check` caught it.** `do.call()` puts the
evaluated arguments into the call, and tinyplot deparses its arguments to
title a legend, so a numeric `by` column arrived as a deparsed 32-element
vector and the legend width computation threw `invalid graphics state`. The
failing call was `ilm_plot_scatter(mtcars, "mpg", "wt", by = "cyl", trend =
"lm")` -- an example that had been in the package, working, for months.

## Finite degrees of freedom for a mixed model

* `ilm_denom_df()` adds Satterthwaite and Kenward-Roger. A Wald statistic
  treats the variance components as known; they are not, and with few clusters
  the chi-square reference is anti-conservative -- which is the regime a
  repeated-measures design lives in.
* Satterthwaite is the default and works for **every structure this package
  fits**, because it asks the objective rather than reimplementing `V(theta)`:
  the fixed-effect covariance with `theta` held fixed is the inverse of the
  beta block of the Hessian of the Laplace objective, recovered by differencing
  `obj$gr()` since `obj$he()` is unavailable once anything is integrated out.
  Checked against `vcov()`, that route reproduces it to **8e-11**.
* Against `lmerTest::lmer(REML = FALSE)`, the same estimator illume uses,
  per-coefficient df agree to between **1e-06 and 5e-06** across four designs,
  and multi-row F denominators match exactly.
* **Kenward-Roger does not come free from the same route and says so.** It
  needs `Q_ij` and `R_ij` separately, and differentiating the beta block twice
  yields only `Q_ij + Q_ji - R_ij`. The pieces separate when the marginal
  covariance is LINEAR in the variance parameters -- random intercepts and
  slopes with one residual variance, where Kenward and Roger derived it and
  where repeated-measures designs live. Once a correlation parameter enters it
  does not, so `ilm_denom_df()` checks first and refuses with the reason,
  naming Satterthwaite instead.

## Restricted maximum likelihood

* `ilm_model(reml = TRUE)`. Adding the fixed effects to the block TMB
  integrates out is not an approximation to REML, it **is** REML: under a flat
  prior that integral is the restricted likelihood, and for a linear-gaussian
  model the Laplace approximation to it is exact. Against `lme4`: variance
  components to 1e-05, coefficients to 1.9e-14, covariance to 1.6e-07.
* On a fixed-effects gaussian model REML reproduces `lm()`'s standard errors to
  **3.2e-11**, against maximum likelihood's 3.7e-08 -- as it should, since
  restricting the likelihood to contrasts orthogonal to `X` is where the
  `n - p` divisor comes from. `vcov()` no longer applies its own correction on
  top, which would have inflated every standard error twice.
* **The default stays maximum likelihood, and the reason is the order you work
  in.** A restricted likelihood belongs to contrasts orthogonal to the design
  matrix, so changing the fixed effects changes which data it is the likelihood
  of. Settle the mean structure under ML, then refit with `reml = TRUE` for what
  you report. `ilm_dag_model()` defaults to REML instead, because the graph
  fixed the adjustment set before any data were seen.
* Available for **every family**, not only gaussian. `glmmTMB` does exactly
  this -- `if (REML) randomArg <- c(randomArg, "beta")` -- so refusing was
  stricter than the reference implementation without being more correct.
  Against glmmTMB with `REML = TRUE`: poisson coefficients to 3.1e-09, binomial
  to 2.8e-16, variance components to 1.3e-08, with REML lifting them 5.1% and
  8.4% off ML's downward bias. What differs by family is how much it delivers,
  which `reml_exact` records: exact for a linear model, approximately
  restricted elsewhere.
* Worth recording because it looked like a defect: illume's REML standard
  errors differ from glmmTMB's by about 0.3%. On a gaussian fit, where `lmer`
  is the canonical reference, illume matches `lmer` to **7.2e-09** and glmmTMB
  differs from it by 9.1e-04. The gap is glmmTMB's convention.
* Guards, because every one of these fails silently otherwise. A
  likelihood-ratio `ilm_anova()` and `ilm_pb_lrt()` refuse on a REML fit and
  say what to do. `ilm_robust()` refuses too -- integrating the coefficients
  out leaves no per-observation score for a sandwich to sum, which is the same
  combination `glmmTMB`'s `estfun` rejects.
* The help pages lagged behind the change above. `ilm_fit()`'s still said
  non-gaussian REML was refused, `ilm_model()`'s pointed to a section that did
  not exist, and `ilm_dag_model()`'s said a restricted likelihood has no
  meaning outside a gaussian model. All three now say what the code does:
  exact REML for a gaussian response, and an approximately restricted
  likelihood for every other family, recorded in `fit$reml_exact`.
  `ilm_dag_model()` still uses REML for a gaussian response only, now for the
  stated reason that it is exact there. The *Regression models* vignette says
  the same.

## Marginal slopes

* `ilm_trends()` takes an interaction apart when one side of it is continuous
  -- the treatment-by-time case, where the omnibus test says an interaction
  exists and says nothing about what is driving it.
* Two questions follow and they are **not the same question**: whether the
  slopes DIFFER between arms, and whether each slope differs from ZERO within
  an arm. They disagree in both directions -- two arms can have slopes that
  differ significantly while neither is distinguishable from zero, and both can
  be strongly non-zero while not differing from one another.
* A marginal mean cannot answer either, because averaging at the mean of the
  covariate collapses the thing being asked about. But a slope is the same
  `L %*% beta` with a different `L`, so `ilm_trends()` returns the class
  `ilm_emmeans()` does and `ilm_contrast()` differences the rows without
  knowing the difference.
* Against `emmeans::emtrends`: slopes to 3e-07 and standard errors to 7e-09 on
  a fixed-effects fit; on a REML mixed fit slopes to 5.8e-10 with Satterthwaite
  df matching to three decimals. Contrasts between slopes match
  `emmeans::contrast` exactly.

## ANOVA, specified by naming columns

* `ilm_aov_ez()`. Name the participant, the outcome, and which factors vary
  between or within participants. No formula with an error term, no reshaping.
* The omnibus table is computed **classically** rather than read off the mixed
  model, and that is the whole reason it exists: `ilm_model()` reports a Wald
  chi-square, because once random effects are integrated out there is no exact
  residual degrees of freedom to divide by. An F, a mean squared error, a
  generalized eta squared and a Greenhouse-Geisser correction do not fall out
  of that fit.
* Against **afex** across five designs -- between only, within only, between by
  within, with a GG correction, and with a covariate -- F agrees to 1.4e-12,
  generalized eta squared to 1.7e-14, p to 9.5e-14, and the fractional
  corrected degrees of freedom to 4.4e-16. Mauchly's W, its p-value and the GG
  epsilon match `car`'s arithmetic.
* Two things had to be right and neither announced itself. The univariate F is
  a ratio of TRACES, and a trace is basis-free only in an orthonormal basis;
  `contr.sum` columns are neither unit-length nor mutually orthogonal, so both
  traces are taken in the `(P'P)^-1` metric. Without it every
  within-participant F is wrong while every between-participant one stays
  right -- the pattern that looks like a modelling disagreement rather than a
  bug. And generalized eta squared needs every stratum's error in its
  denominator including the participant stratum, which in a purely
  within-participants design never reaches the table; without putting it back,
  every `ges` came out more than twice too large and each value looked
  plausible.
* Beyond afex, each check names a remedy that exists here. Sphericity fails ->
  an unstructured mixed model does not assume it. Participants have missing
  cells -> they are dropped, counted, and the mixed model uses them. A
  covariate varies within participant -> refused, because one coefficient for a
  time-varying covariate blends how participants who score higher on average
  differ with what happens when a participant scores higher than usual, and
  those can have opposite signs.
* Interactions are followed up with **simple effects in both directions**
  rather than every cell pair: all nine cells of a 3-by-3 give thirty-six
  differences, most of which move both factors at once. Where the interaction
  is with a covariate the follow-up is `ilm_trends()`.

## The studies, re-run in full

* All five standing studies re-run at this version alongside the new imputation
  ablation, at the parameters used at 0.0.2.9000 so the two are comparable.
  Every cell that existed then reproduces it.
* Coverage grew from 23 cells to **32**. Everything added since 0.0.2.9000 had
  agreement against an outside implementation and no coverage of its own, and
  those are different claims: agreement says the point estimate is right,
  coverage says the interval is. The new cells hold -- beta 0.949,
  zero-inflated beta 0.951, hurdle Poisson 0.949, zero-inflated Poisson mixed
  0.950, ordered logit 0.946 fixed and 0.953 mixed, 2SLS 0.946.
* Two do not sit at nominal and both are recorded. `mediate` **over**-covers at
  0.967 with a standard-error ratio of 1.107: the ACME interval is a percentile
  interval from simulation draws and runs about 11% wider than the sampling
  variability of the estimate. `svy_strat` covers at 0.939 with a ratio of
  0.956, about two Monte Carlo standard errors below nominal -- small, but in
  the anti-conservative direction.

## Where illume loses to mclogit

* A second comparison against `mclogit::mblogit()`, this time on data that
  misbehaves: unbalanced clusters with many singletons, a 3.6% outcome
  category, a variance component at the boundary, non-Gaussian random
  effects, and all four together. 400 replications per regime.
  `studies/findings/messy.md`, and the new **Benchmarking and validation**
  vignette.
* illume's intervals cover better in all six regimes, and the gap widens
  where theory says PQL should struggle -- 0.931 against 0.853 with a sparse
  category, 0.945 against 0.894 with everything combined. illume is 6 to 13
  times faster throughout.
* **Three results go the other way and are recorded rather than smoothed.**
  illume *inflates* coefficients when a category is sparse (attenuation 1.46),
  and its mean absolute bias there is worse than mclogit's, 0.394 against
  0.271; the intervals are wide enough to cover anyway, but the point estimate
  should not be read at face value. illume converged on **63.3%** of
  replications when the true random-effect sd was 0.05, against mclogit's
  100%, so its coverage in that cell describes only those 63.3%. And
  mclogit's RMSE is lower in four of the six regimes.
* The two methods fail differently rather than one dominating. PQL always
  converges, shrinks, and does not say so. The Laplace approximation means
  what its intervals claim and declines more often. Declining loudly is the
  intent, but it is still a cost.

## What the imputation ablation says

* Five methods across six designs, scored on reconstruction AND on coverage of
  a downstream coefficient after Rubin pooling -- the two disagree, and only
  the second is what an analysis needs. `studies/findings/imputation.md`.
* **Chained equations wins wherever it can be fitted**: best or joint-best
  coverage in all five such cells (0.900, 0.900, 0.975, 0.900, 0.950), the
  smallest bias, the best reconstruction. It stays the default.
* It also recovers CATEGORIES better than the generalized low rank model does
  -- 0.634 against 0.548, 0.615 against 0.539 -- which was not the expected
  result. GLRM's claim on imputation is that it can do categorical columns at
  all where the low-rank route returns `NaN`, not that it does them better.
* **The low-rank route earns its place in exactly one case and holds it
  there**: at n = 60, p = 80 chained equations cannot be fitted, and the
  low-rank reconstruction halves the error against column means (0.559 against
  1.007) at 0.975 coverage. Elsewhere the low-rank methods under-cover badly,
  0.50 to 0.65 in four of six cells.
* **The rank selector is NOT changed, and the earlier finding did not
  generalise.** One design had suggested `ilm_lowrank_ncp()` simply picks
  badly. Across six, cross-validation is worse than parallel analysis in four
  cells and materially better in one -- 0.950 against 0.475 on the design with
  the most data and the clearest structure, where there are enough held-out
  cells to choose well. Neither criterion dominates, so neither is imposed.
  `ilm_anomaly()` continues to use parallel analysis for its own reason: it
  needs the directions that are real shared structure, not the rank that best
  predicts a cell.
* Everything under-covers somewhat, mean-fill included (0.775 to 0.925 against
  a nominal 0.95). At 40 replicates the Monte Carlo error is 0.034, so the
  table separates acceptable from broken rather than ranking 0.90 against 0.95.

## Dependency minimums, measured rather than assumed

* Every `Imports:` entry that needed a floor now has one, and each was
  established by test rather than by reading a changelog: `RTMB (>= 1.7)`,
  `TMB (>= 1.9.7)`, `collapse (>= 2.0.2)`, `tinyplot (>= 0.6.1)`. Previously
  there were none at all, so a resolver was free to pick a version that would
  fail somewhere unhelpful.
* **RTMB 1.7 was built and every family fitted against it**: gaussian, gaussian
  mixed, binomial, Poisson mixed, negative binomial, beta, ordinal,
  zero-inflated and multinomial mixed, plus `summary()`, `ilm_anova()`,
  `predict()`, `ilm_rqr()` and `ilm_emmeans()`. Fourteen of fourteen.
* The AD objective touches eight RTMB symbols -- `AD`, `ADREPORT`, `diag`,
  `getAll`, `logspace_add`, `logspace_sub`, `matrix`, `solve` -- and all eight
  are in 1.7. Six other RTMB exports illume uses are newer than 1.7
  (`findInterval`, `order`, `pchisq`, `pnbinom`, `qchisq`, `sort`), but every
  one of those is called on ordinary numerics, where base and stats answer.
* Two RTMB changelog entries looked like they might bite and do not. The 1.8
  `dgamma()` rate fix does not reach the only `dgamma()` here, which is
  `stats::dgamma()` on a plain vector. The 1.9 addition of `log.p` and
  `lower.tail` to `pnorm()` is unused, because the censored and lognormal
  likelihoods write `pnorm(mu - y)` rather than `1 - pnorm(y - mu)` -- a choice
  made for tail precision that happens to keep the API two versions older.
* `collapse`'s floor is `fmatch`, absent at 1.9.6 and present from 2.0.2.
  `Matrix` needs no floor: only `solve()` and `chol()` are used.

## Plotting no longer requires the newest tinyplot

* `tinytheme_list()` arrived in tinyplot 0.7.0 and was being used to enumerate
  valid theme names for an error message -- which meant `theme =` was refused
  outright on older tinyplot, although `tinytheme()` and the `theme` argument
  have both been there since 0.6.x. The shape of `theme` is now always checked,
  the names are enumerated when the helper exists, and otherwise tinyplot is
  left to object. Themes work again on 0.6.1.
* `ilm_plot_var_pairs()` genuinely needs `tinypairs()`, which is 0.7.0-only. It
  now says so and names `ilm_plot_scatter()` and `ilm_plot_var_all()` as the
  alternatives, rather than failing with "not an exported object".
* Measured on tinyplot 0.6.1: twelve of thirteen plotting paths work, the
  thirteenth being the pairs plot reporting its own requirement. On 0.7.0 the
  suite is unchanged.

## Documentation

* Twelve vignettes. `workflow` is new and is the map: eleven stages from a
  power calculation through to reporting, with scenario projection at stage 10.
  The introduction is now an orientation rather than a tutorial, and the
  modelling material it used to carry has become `regression-models`. Also
  new: `profiling`, `anomaly-detection`, `missing-data`,
  `effect-size-and-power`, `anova`, `moderation`, and `benchmarking`, which is
  the evidence -- what was measured, against what, and where illume comes off
  worse.
* `regression-models` now opens with why the multinomial mixed model is fitted
  the way it is, against `mclogit` and `brms`, with the results that go against
  illume beside the ones that do not.

# illume 0.0.6.9000

Faster where it was slow, a way to choose columns, progress where it is worth
showing, and an imputation route for data too wide to regress.

## Choosing columns

* `cols` and `by` accept a character vector, a **regular expression**, a
  **predicate function** such as `is.numeric`, or `NULL` -- see [ilm_selection].
  All four are ordinary values, so a selection can be held in a variable and
  passed on, which a bare-name interface gives up. No non-standard evaluation
  and no new dependency.
* A lone string is ambiguous between a name and a pattern. A name wins, and
  when neither works the error says both readings were tried.

## Progress

* `progress` on the bootstrap functions, `ilm_impute()`, `ilm_cluster()` and
  the refit-based diagnostics. It defaults to [interactive()], so a bar appears
  when a person is watching and nothing is written into a script, a test or a
  knitted document.
* Where work is spread over cores the bar advances per chunk, since workers are
  separate processes and cannot write to the parent's console.

## Imputation for data too wide to regress

* `ilm_impute(method = )` takes `"auto"`, `"fcs"` or `"lowrank"`. Chained
  equations needs more rows than predictors; past that point there is no
  regression to fit, and `"auto"` falls back to a regularised low-rank
  reconstruction, bootstrapped so the imputations differ and pooled by
  `ilm_mi_pool()` -- multiple imputation PCA, after Josse and Husson (2016).
* The rank is chosen by cross-validation over held-out observed cells.

## Findings behind those changes

* **The BCa jackknife was quadratic in n.** It needs every leave-one-out value
  of the statistic, which is n evaluations on vectors of n - 1: 4.31s of a
  6.07s call at n = 20,000, against 2.55s for the resampling itself. The mean,
  variance and standard deviation have exact leave-one-out forms in the running
  sums, which is linear and measured at 0.00s, matching the loop to 1e-16. A
  BCa interval for the default statistic now costs the same as a percentile
  one, 6.07s to 2.52s.
* **The resampler could not scale.** Drawing every index at once into an n by R
  matrix is a 160 MB allocation at n = 20,000 and 8 GB at a million rows, so
  the function stopped working rather than slowing down. Per replicate it is
  also 12% faster.
* **Binding a frame per group was most of the cost of a grouped outlier scan.**
  28 columns by 8 groups is 224 rbind calls, each copying everything
  accumulated: 1.58s to 0.98s by accumulating into vectors and binding once.
* **A low-rank fit is better than chained equations only where chained
  equations cannot be fitted.** Hiding known cells and scoring against them,
  noise floor 0.500:

  ```
    design                          mean-fill    fcs   lowrank
    n=200 p=8  rank 3                   2.289  1.177     1.963
    n=200 p=8  full rank                2.739  2.622     3.581
    n=400 p=12 rank 4, 30% missing      2.010  1.054     1.683
    n=60  p=80 rank 3  (p > n)          1.920      -     0.887
  ```

  The middle row is the warning: with no low-rank structure to find, imposing
  one is worse than filling in column means. `ilm_impute()` detects that by
  asking whether a rank-k fit predicts held-out cells better than the column
  means do, and says so when it does not. A first attempt warned when
  cross-validation hit its rank ceiling and MISSED the case entirely -- on pure
  noise it picked rank 5 of an allowed 7.
* Reconstruction accuracy is not the measure that governs here. It says how
  close the filled values are, not whether inference afterwards is calibrated;
  single imputation scores respectably on the first and covers at 0.79 to 0.89.

## Fixes

* `ilm_reduce()` failed outright on any two-column selection. The number of
  dimensions available from p columns is min(n - 1, p), not p - 1 -- two
  columns have two components, and asking PCAmix for one is an error rather
  than a smaller answer. It also under-counted everywhere else: three columns
  returned two.
* A bootstrap draw leaves some rows with weight zero, and the low-rank fit
  rescaled its reconstruction by the inverse of those weights, producing `Inf`
  for exactly those rows and failing the next decomposition. The shrinkage is
  now a ratio applied to an unweighted projection.
* The test for "too wide for chained equations" counted rows complete across
  every column, where what matters is rows on which each variable was observed:
  8 columns at 20% missing leaves 17% of rows complete and every one of them
  usable.
* The fitted noise scale was indexed by column name while carrying none.

# illume 0.0.5.9000

The rest of the exploration layer: unusual values, and named plots.

## Flagging unusual values

* `ilm_outliers()` scores every value by how far it sits from the centre and
  flags those past a threshold, by Tukey fence (`"iqr"`), modified z-score
  (`"mad"`) or ordinary z-score. `ilm_outliers_all()` runs it over every numeric
  column, optionally within groups, and returns the row each flagged value came
  from.
* Grouping is not a detail. A value can be ordinary for its own group and
  extreme against the pooled distribution, so flagging without `by` on grouped
  data largely rediscovers the groups.
* The documentation says what a flag is **not**. It is not a verdict that a
  value is wrong: a flagged value may be a recording error, a member of another
  population, or an ordinary draw from a heavy tail. Dropping flagged rows
  because they are flagged changes the estimand; if the tail is real the remedy
  is a model that expects it.

## Named plots

* `ilm_plot_histogram()`, `ilm_plot_density()`, `ilm_plot_box()`,
  `ilm_plot_violin()`, `ilm_plot_scatter()`, `ilm_plot_bar()`,
  `ilm_plot_line()` and `ilm_plot_stat_error()`. `ilm_plot()` picks a geometry
  for you; these are the same drawing done by naming the plot you want.
* `ilm_plot_var()` and `ilm_plot_var_all()` choose from the column types;
  `ilm_plot_var_pairs()` plots every pair, handling mixed numeric and
  categorical columns rather than only numeric ones; `ilm_plot_c()` composes
  several plots into one figure.
* `ilm_plot_na_all()` and `ilm_plot_na()` show missingness across columns and
  across groups. `ilm_plot_missing()` now draws through the first of these, so
  there is one implementation.
* Everything is on tinyplot. That was the last base-graphics plot in the
  package, and columns are named as strings throughout, the way the rest of
  illume takes them.

## Bootstrap differences can be seen, not just summarised

* `ilm_boot_diff()` keeps the replicate differences, labelled by comparison,
  and reports `p_superiority` -- the share of replicates above zero.
* `ilm_plot_boot_diff()` draws them for any comparison in the result, with zero
  and the interval marked. The interval says where the difference is; this says
  what the resampling produced, which is whether it is symmetric, skewed, or
  piled against a boundary.

## Findings behind those changes

* **The Tukey fence is drawn from HINGES, not from quantiles, and the
  difference is visible.** `ilm_outliers(method = "iqr")` claimed to reproduce
  the rule a boxplot's whiskers use, and did not: built on `quantile()`'s
  type-7 default it put the upper fence for `mtcars$wt` at 5.153 where
  `fivenum()`'s hinges put it at 5.311, so a value of 5.25 was flagged while
  sitting inside the whisker it was supposed to match. It now uses the hinges
  and agrees with `grDevices::boxplot.stats()` exactly, which a test pins
  across four data sets.
* **A numeric grouping column made tinyplot warn on every call.** Handing it
  one for a discrete geometry makes it attempt a continuous legend, fail, and
  revert with a warning. The grouped geometries now take `by` as a factor,
  which is what a grouping variable is; scatter and line still accept a
  continuous one, where it means something.

# illume 0.0.4.9000

Missing data, and the profiling set that reads its patterns.

## Missing values are no longer dropped in silence

* `ilm_model()` reports how many rows went and which columns took them. Past a
  tenth of the data it also says what that does and does not imply. Dropping
  incomplete rows is usually right; doing it silently is not, because a model
  fitted to 61% of the data with no note of it invites conclusions the data
  cannot carry. The count is kept on the fit as `$n_dropped`.

## Deciding whether it matters

* `ilm_check_missing()` reports how much is missing, which columns go missing
  together, whether the pattern is monotone, and what missingness is related to
  -- then says what that implies.
* **It distinguishes the two cases that need different answers.** For a
  regression, dropping incomplete rows is unbiased whenever missingness is
  independent of the OUTCOME *given* the covariates -- which is far weaker than
  MCAR, and often true. So missingness tied to a covariate gets "complete cases
  stay unbiased", and only missingness tied to the outcome gets sent to
  `ilm_impute()`.
* **MAR versus MNAR is not tested, because it cannot be.** The data that would
  separate them are the data that are missing. What is testable is MCAR, and
  that is what is tested; the result says plainly that a dependence on the
  unseen values themselves cannot be ruled out, and that the remedy for that is
  a sensitivity analysis rather than a test.

## Filling them in, honestly

* `ilm_impute()` does multiple imputation by chained equations, one
  [ilm_model()] per incomplete variable with the family its own type calls for.
  Each value is **drawn** from the predictive distribution -- coefficients from
  their sampling distribution, the residual scale from its own posterior, then
  the response's randomness on top -- rather than set to a fitted mean.
* `ilm_mi_pool()` fits across the imputations and combines by Rubin's rules,
  with Barnard-Rubin degrees of freedom, and reports the fraction of
  information lost to missingness per coefficient.
* `single = TRUE` gives one completed data set and warns, because anything
  computed from it will be overconfident.

## Profiling: dimension reduction, clustering, and what the clusters are

* `ilm_reduce()` reduces a frame's columns to a few dimensions, choosing PCA,
  MCA or a mixed method from the column types rather than making the user name
  it.
* `ilm_cluster()` groups the rows, choosing `k` by the gap statistic, and
  reports two different things a bare assignment does not: per-cluster
  **stability** by bootstrap Jaccard, and per-observation **ambiguity** from
  the silhouette width. A point can sit on a boundary inside a large, stable
  cluster, and a size-based flag alone would never show it.
* `ilm_profile()` runs both and says what each cluster *is* -- "cluster 4 is
  characterised by dim 1 (high disp, cyl)" -- by v-test, the device FactoMineR's
  `catdes()` uses. It is documented as a threshold rather than a test, because
  the clusters were found from the coordinates being tested.
* `ilm_reduce_na()`, `ilm_cluster_na()` and `ilm_profile_na()` do the same to
  the *pattern of missingness*: which columns go missing together, and for
  whom. A block of variables lost as one points at a shared cause, which is a
  different problem from values going one at a time.
* `ilm_plot_reduce()`, `ilm_plot_reduce_scree()`, `ilm_plot_reduce_contrib()`,
  `ilm_plot_cluster()`, `ilm_plot_cluster_gap()` and `ilm_plot_profile()`, each
  with an `_na` counterpart.
* `PCAmixdata` and `cluster` join Suggests behind require-guards. Neither is a
  hard dependency and no new Imports were added.

## Findings behind those changes

* **Complete cases are fine more often than they are given credit for, and the
  simulation says by how much.** 200 replicates, 400 rows, `y = 0.5x + 0.3z`,
  coverage of a nominal 95% interval for the coefficient on `x`:

  ```
                         full   complete   multiple    single
                         data      cases  imputation  imputation
    MCAR, 30% missing
      bias            -0.0028    -0.0012     -0.0035    -0.0069
      coverage          0.955      0.960       0.970      0.890
    MAR on a covariate, 40% missing
      bias            -0.0028     0.0001     -0.0017     0.0002
      coverage          0.955      0.980       0.970      0.835
    MAR on the OUTCOME, 41% missing
      bias            -0.0028    -0.0989     -0.0090    -0.0058
      coverage          0.955      0.615       0.955      0.790
  ```

  Complete cases are unbiased at 40% missing when missingness follows a
  covariate, and fail badly when it follows the outcome: a bias of -0.099 is a
  fifth of the effect and coverage collapses to 0.615. Multiple imputation
  repairs exactly that case.
* **Single imputation is the cautionary column.** Its point estimates are no
  worse than multiple imputation's anywhere in that table, and its coverage runs
  0.790 to 0.890, because nothing in its standard errors knows part of the data
  was invented.
* **The imputation agrees with `mice`.** On 41% MAR-on-outcome missingness at
  n = 500 with m = 20, the pooled estimates differ by 0.34 of a pooled standard
  error and the standard errors by a factor of 1.016, with comparable degrees of
  freedom -- while complete cases sat at 0.397 against a truth of 0.500.
* **A marginal test answers the wrong question, and looks right doing it.** The
  first version of `ilm_check_missing()` asked whether missingness was
  associated with the outcome, full stop. But when missingness follows a
  covariate the outcome also depends on, the two are marginally associated
  while carrying no information about each other once that covariate is held
  fixed -- measured at 0.273 on data where complete cases covered 0.980. The
  verdict now comes from a conditional test, which reads 0.039 there and 0.454
  where the outcome really does drive it. The marginal table is still reported,
  as description.
* **Missingness profiling recovers the structure it is given.** On 300 rows
  where two columns were made to go missing as a block and a third
  independently, the first dimension loaded the block at 0.9997 each and the
  independent column at 0.0013, which landed on the second dimension at 0.9987
  instead.

## Fixes

* `ilm_pool()` was already taken. `R/ilm_parallel.R` defines it to build a
  cluster of worker processes and six diagnostics call it, so the Rubin's-rules
  pooler is `ilm_mi_pool()`. Load order happened to favour the existing
  function, so the new one was unreachable rather than breaking `ilm_anova()`,
  `ilm_check_ar()` and the rest.
* A gaussian imputation drew its noise from `sd(y)` where `sigma()` and
  `residuals()` are unavailable for an `ilm_model`. That is the marginal
  spread, before the predictors explain any of it, so every imputation carried
  more noise than the model said was there. `$dispersion` is the residual scale
  and agrees with `lm()`'s sigma to the printed digits.
* A cluster summary read "1 of its member sits"; "members" is plural whatever
  the count, and only the verb agrees.
* A cluster characterised by two dimensions with the same top-loading variables
  named them twice in one sentence, which happens whenever there are few
  variables to go round -- as in a missingness profile of a frame with two
  incomplete columns.

# illume 0.0.3.9000

Designs that identify an effect, and a way to read a fit back in words.

The inference engine is unchanged from 0.0.2.9000 -- `ilm_model()` returns the
same estimates -- so the coverage and Type I error evidence recorded under that
version still applies to everything here. What is new sits on top of it.

## Causal graphs

* `ilm_dag()` reads a causal graph from a dagitty-style string, an edge frame or
  a `dagitty` object. A bidirected edge is stored as an unobserved common cause,
  so there is one representation rather than two.
* `ilm_adjust_sets()` gives the minimal sets of **measured** variables that
  identify an exposure effect. An empty result is a finding: it says the data
  cannot answer the question, which is the most useful thing a graph can say and
  can only be said before the modelling.
* `ilm_dag_implied()` lists what the graph claims about the data -- one testable
  independence per missing edge -- and `ilm_dag_test()` checks those claims,
  with a verdict per claim and one overall.
* `ilm_dsep()` exposes the d-separation the rest is built on.
* The graph algorithms are written in-package. `dagitty` imports V8, and a
  JavaScript engine is a heavy thing to require of someone who wants to fit a
  regression; it stays in `Suggests` and the tests pin these functions to it.

## The DAG-guided workflow

* `ilm_dag_model()` takes a graph and a data frame and runs the analysis:
  checks the graph against the data, finds what must be adjusted for, picks a
  response distribution, looks for grouping structure, fits, diagnoses, and
  repairs the error structure when a check asks. Verbose by default.
* **The graph fixes the mean structure and nothing searches over it.** Only the
  error structure is adjusted -- a dispersion model, a random effect for
  grouping the graph does not mention. Searching over covariates and reporting
  the winner's p-values is post-selection inference and would undo the nominal
  coverage the studies establish. A DAG is a pre-registration device and this
  treats it as one. Every step is recorded in `$steps`.
* Where several minimal sets are admissible, all are fitted and reported side by
  side. They target the same quantity, so the spread across them is a
  sensitivity analysis that costs only compute.
* Grouping variables are taken only from **outside** the graph. A variable the
  graph mentions has a causal role and belongs in the mean structure; anything
  else can only be structure in how the data were collected.

## Difference in differences

* `ilm_did()` estimates the ATT, tests parallel trends in the pre-period, and
  returns an event study that `ilm_plot_did()` draws.
* Staggered adoption is **refused** rather than quietly estimated: under
  heterogeneous effects a pooled two-way fixed effects estimate is not an
  average treatment effect, because already-treated units serve as controls for
  later-treated ones. `allow_staggered = TRUE` returns it with the caveat
  attached.
* A unit random intercept is fitted by default and `ar = TRUE` adds AR(1), which
  is what the serial correlation in a panel calls for.

## Regression discontinuity

* `ilm_rdd()` gives the jump at the cutoff from a local linear fit with a
  triangular kernel, plus four design checks: density at the cutoff, covariate
  balance, placebo cutoffs, and how far the estimate moves with the bandwidth.
  `ilm_plot_rdd()` draws binned means with the fitted lines.
* Fuzzy assignment is detected and not reported as sharp.
* The bias-corrected robust intervals of Calonico, Cattaneo and Titiunik are
  `rdrobust`'s and are not reimplemented. What is reported is an honest local
  fit with its sensitivity laid out.

## Reading a fit in words

* `ilm_interpret()` writes out what a model says: each effect on the scale the
  response is measured on, the strength of the evidence, what the diagnostics
  found, and how far to trust the estimates. Methods for `ilm_model`,
  `ilm_dag_model`, `ilm_did` and `ilm_rdd`.
* The prose is **templated, never generated**, so the same fit gives the same
  words and every sentence is testable.
* **Causal language is licensed, not assumed.** A coefficient is an association
  and is called one, unless the object carries a design that identifies an
  effect -- and where a graph licenses it, the interpretation says the licence
  is an assumption the user supplied rather than something the data established.
* Nothing is said about bias that a diagnostic did not measure. Every such
  sentence traces to a check that ran, carries its verdict and names its remedy.
* `ilm_ame()` gives average marginal effects on the response scale by the delta
  method, over the full parameter vector including the covariance parameters,
  since a population-averaged prediction depends on them.
* `ilm_register_insight()` makes `parameters`, `performance` and `report` work
  on an `ilm_model`. `insight` stays in `Suggests`.

## Comparing more than two groups

* `ilm_boot_diff()` compares every pair of levels rather than exactly two, takes
  a formula, and is simultaneous by default. See the 0.0.2.9000 notes.

## Findings behind those changes

* **The adjustment sets and d-separation agree with `dagitty` exactly.** Across
  400 random graphs: 1564 d-separation tests with 0 disagreements, 391 of 391
  identical minimal adjustment sets, and 2532 implied claims every one of which
  dagitty confirms is a d-separation with no redundant member.
* **Conditioning sets built from parents silently drop testable claims.** The
  first version took the parents of each pair and discarded the claim when a
  parent was unobserved. But a minimal separator usually does not contain that
  parent -- most often it is empty -- so genuine claims were being lost. The
  search now runs over subsets of the observed ancestors, where any minimal
  separator must lie, and returns a minimum rather than merely minimal set, so
  each test holds as few things fixed as possible.
* **The parallel-trends check needed the right reference distribution.** It
  compares two groups of *units'* pre-treatment slopes, so units are the
  independent replicates. The large-sample normal rejected 0.060, 0.068 and
  0.050 of the time at 40, 20 and 80 units under trends that really were
  parallel; t on `units - 2` gave 0.048, 0.055 and 0.045. At 800 replicates the
  finished check raised a false alarm 0.054 of the time with null p-values
  uniform by Kolmogorov-Smirnov (p = 0.111), while the ATT itself was unbiased
  (+0.0029) with 0.946 interval coverage.
* **A density check can ask the wrong question and look right.** Comparing
  counts either side of an RD cutoff against 50/50 asks whether the density is
  SYMMETRIC there; the question is whether it is CONTINUOUS. Over 500
  unmanipulated data sets the count split flagged 0.044 of uniform running
  variables but 0.760 of a sloped normal and 1.000 of an exponential, all of
  them perfectly smooth. A local linear density on each side held 0.046 to 0.060
  throughout, and had more power too: with 30% of the units just below the
  cutoff moved above it, the count split caught 0.808 and the local linear fit
  0.984.
* **The marginal effects agree with `marginaleffects`.** On a binomial fit the
  estimates match `avg_slopes()` to 1e-6 and the standard errors to 1e-7, for
  slopes and for a factor contrast, which are different calculations on both
  sides.

## Fixes

* `ilm_did()` read a two-level factor's labels in alphabetical order, so
  `post = factor(x, levels = c("before", "after"))` made "before" the treated
  period and inverted the estimate. A factor's own level order is the author's
  intent and is now what is used.
* `has_rp` and `Drp` are declared in `globalVariables()`. They entered the RTMB
  data list with the Royston-Parmar work and `getAll()` binds them at run time,
  which `codetools` cannot see, so `R CMD check` reported them as undefined.

# illume 0.0.2.9000

Time-to-event, censoring, correlation over irregular time, and a model for the
spread. Everything the diagnostics could only name a remedy for, they can now
also fit.

## Features

* `ilm_model()` fits every supported family through one formula interface, so
  the object, the methods and the diagnostics do not change as the model gets
  more complicated.
* Gaussian models with no random or smooth terms report **exact** t and F
  inference agreeing with `lm()` and `car::Anova()`, rather than the
  large-sample approximations used elsewhere.
* `ilm_anova()` gives Type II and Type III analysis of deviance, testing each
  term jointly across all category dimensions, with Wald, likelihood-ratio and
  parametric bootstrap variants.
* Built-in diagnostics (`ilm_appraise()`, `ilm_rqr_test()`, `ilm_check_ar()`)
  state whether each assumption is consistent with the data and name a remedy
  when it is not.
* Reduced-rank (`rr`) category covariances make larger numbers of outcome
  categories tractable.

## Time structure and cycles

* `ilm_plot_acf()` draws the residual autocorrelation and partial
  autocorrelation by lag, within group, against an envelope built by simulating
  from the fit and refitting. The partial autocorrelation comes from
  Durbin-Levinson applied to the autocorrelations the refits already produced,
  so it costs nothing extra and agrees with `stats::pacf()` to 1e-8.
* The band is drawn at the critical value the p-value is computed against, so a
  point outside the band is exactly a lag the check flags.
* `ilm_check_ar()` now distinguishes a **cycle** from autoregression and names
  the period. The two need opposite remedies and an AR term cannot represent a
  cycle at all.
* `ilm_fourier()` and `ilm_cyclic()` fit that cycle. `ilm_cyclic()` is a cubic
  spline that closes on itself, so the curve and its first two derivatives match
  across the wrap.
* `ilm_check_ar()` and the residual machinery now work for every family.
  Previously the Pearson residual was indexed by `object$J`, which is 2 for each
  univariate family and did not match the one residual series those families
  have.

## Correlation over irregular time

* `ilm_car1()` fits a first-order autoregression whose correlation is `rho`
  raised to the gap between observations, so the spacing need not be regular.
  `ilm_ar1()` is the evenly spaced constructor; neither existed before, and the
  bare `list(idx =, n_group =, Tt =)` the fitter used to take still works.
* `ilm_variogram()` and `ilm_plot_variogram()` bin every within-group pair by
  separation and report the residual correlation in each bin against the
  refit-based envelope. This is what makes a `ilm_car1()` term checkable, and it
  tells apart the three shapes that need three different remedies:
  autoregression, an unabsorbed group effect, and a cycle.
* The latent budget checks no longer fail a gaussian model. With a gaussian
  response and gaussian latents the model is linear-Gaussian and the Laplace
  approximation is exact, so there is no approximation error to run out of
  observations for.

## A flexible survival baseline

* `family = "rp"`, `"rp_odds"` and `"rp_normal"` fit the Royston-Parmar model:
  a restricted cubic spline in log time in place of the straight line the
  parametric families assume, on the hazard, odds or probit scale. The
  coefficients keep their meaning -- a hazard ratio on the hazard scale -- while
  the baseline is free to bend.
* `rp_df` sets how much it can bend, `rp_knots` places the knots directly.
* `ilm_rp_lrt()` tests the spline against the straight-line special case it
  generalises, and names the simpler family when that fits as well.
* `ilm_survival()` and `ilm_plot_survival()` work for it. `predict()` on new
  data does not, and says why: the linear predictor depends on time through the
  baseline, so there is no fitted value for a covariate pattern alone.

## Correlation in space

* `ilm_variogram(coords = )` bins pairs by Euclidean distance rather than by
  separation in time. Time and space are the same statistic over a different
  distance, so it is the same function; a one-dimensional coordinate gives the
  absolute time difference exactly.
* The advice is spatial where the question is. illume fits no spatial
  covariance, so it names what it does have: a tensor-product smooth of the
  coordinates for smooth variation, a random effect for the coarse kind.
* `min_effect`, 0.1 by default, is the smallest departure from the simulated
  null worth a verdict.

## Modelling the spread

* `ilm_model(dispformula = )` models the logarithm of the dispersion:
  `~ group` for a separate spread per level, `~ x` for one that changes with a
  covariate, `~ mu` for a power of the fitted mean. `mu` is a reserved name.
  This is the remedy `ilm_check_variance()` diagnoses, and the check now names
  it.
* Available for every family with a dispersion parameter -- gaussian, negative
  binomial, and the three accelerated failure time families. A Poisson or
  binomial response has none, and asking says so.
* A dispersion model rules out exact t and F inference, for the same reason
  censoring does: both rest on a constant variance.

## Censored responses

* `ilm_censor()` marks observations known only as an interval -- at or below a
  floor, at or above a ceiling -- and `ilm_model(censor = )` fits them as the
  probability of that interval. With no random effects this is the Tobit model.
* `ilm_describe()` names a pile-up at either extreme of a continuous variable,
  and points at `ilm_censor()`, so the problem is visible before the model is
  fitted rather than after.
* Censoring switches off the exact t and F path, because a censored likelihood
  is only asymptotically normal and the `N/(N-p)` correction is an
  ordinary-least-squares result.
* Quantile residuals for censored rows are drawn across the probability of the
  censoring interval, the same randomisation Dunn and Smyth apply to a discrete
  response.
* Simulated replicates are censored by the same limits the data were, so the
  envelope diagnostics are calibrated against the fitted process rather than an
  uncensored one.

## Time to an event

* Three accelerated failure time families -- `"weibull"`, `"lognormal"` and
  `"loglogistic"` -- fitting `log T = X beta + scale * W`, so a coefficient is a
  log time ratio whatever the baseline hazard does.
* `ilm_surv(time, event)` builds the censoring specification in the convention
  `survival::Surv()` uses, which is the opposite of the code stored internally.
* `ilm_survival()` gives the predicted survival curve with intervals, and
  `ilm_plot_survival()` draws it over the Kaplan-Meier estimate with an
  envelope and a verdict. The Kaplan-Meier is computed in-package, and matches
  `survival::survfit()` to 3e-16.
* A family now carries its own distribution function, and the residual
  machinery asks it rather than keeping a third switch over families. Two
  separate switches is how the quantile and Pearson residuals both came to
  assume a multinomial response.

## Comparing more than two groups

* `ilm_boot_diff()` compares every pair of levels, not just two, and takes a
  formula: `ilm_boot_diff(score ~ grp, data = d)` as well as
  `ilm_boot_diff(d, "score", "grp")`. `ref = ` compares every level against one
  instead, and several grouping columns are crossed with `interaction()`.
* The intervals are **simultaneous by default**. Five levels give ten
  comparisons, and ten intervals each nominally 95% do not jointly cover at
  95%; reporting them as though they did is the usual way a pairwise table
  misleads. `adjust = "max_t"` resamples all groups together, standardises each
  comparison by its own bootstrap standard error and takes the `conf` quantile
  of the largest standardised value as one critical value for all of them.
  `"bonferroni"` and `"none"` are there for when they are wanted.
* Each row now carries `p_value` and `p_adj` alongside the interval, read off
  the same bootstrap maximum, so the test and the interval agree.
* `ilm_boot_diff()` is now generic, and its first argument is `x` rather than
  `data`. Positional calls are unaffected; a call naming `data =` for the data
  frame needs the formula form.
* Groups are resampled jointly rather than one comparison at a time -- which is
  what makes a simultaneous statement possible, and changes the random stream,
  so a two-group result at a given `seed` differs numerically from 0.0.1.9000.

## Findings behind those changes

Summary notes; the full tables belong with the methods paper.

* **A `2/sqrt(n)` autocorrelation band is the wrong reference here.** Residuals
  within a group sum to roughly zero, so pairs of them correlate negatively even
  under a correct model: measured at about -0.07 at every lag on a 60-unit,
  10-period gaussian panel. A band centred on zero would have flagged all five
  lags examined on correctly specified data.
* **A cycle was being answered with the wrong remedy.** On a 50 x 36 monthly
  panel carrying an annual cycle and no autoregression, the check reported four
  of five lags outside the envelope and advised an AR(1) term. The residual
  autocorrelation ran +0.45 +0.28 -0.02 -0.27 -0.46 -0.55 -0.43 -0.26 +0.01
  +0.25 +0.45 **+0.50** -- a clean return to a peak at lag 12. The return, not
  the sign change, is what separates a cycle from autoregression: an AR(1)
  fitted without its AR term also goes negative at the longer lags, but does not
  come back.
* **No basis for a cycle is generally best.** At matched degrees of freedom on a
  twelve-phase cycle, Fourier terms and the cyclic spline were within a few AIC
  of each other across four shapes (sinusoid, asymmetric, narrow spike, two
  peaks). On the two shapes that genuinely jumped, `factor(phase)` beat both
  smooth bases by more than 100 AIC. On a long cycle -- day of year with one
  narrow summer feature -- Fourier led below about six degrees of freedom and
  the spline led above it, by at most 23 AIC.
* **Type III main effects depend on the contrast coding.** On `y ~ g * x + z`
  with 600 rows, the `x` row came back at a chi-square of 167 under
  `contr.treatment` and 1086 under `contr.sum` -- same model, same data. Rows
  spanning more than one column (`g`, `g:x`) were invariant, since both codings
  span the same subspace. `ilm_anova()` now warns and names the fix.

* **The log-likelihood was missing a constant that grew with the model.** The
  random-effect priors are written as `0.5 u'Sigma^-1 u + 0.5 log|Sigma|`,
  leaving out the `(1/2)log(2*pi)` each latent scalar contributes; TMB's Laplace
  step then subtracts `(q/2)log(2*pi)` of its own. On a random-intercept fit
  with 50 groups illume reported -595.083 where lme4 and nlme both reported
  -641.030 -- a gap of 45.947 against `25 * log(2*pi) = 45.947` predicted. Every
  extra latent value was worth about 0.92 log-likelihood units for free, so AIC
  and BIC preferred whichever model carried the larger random structure.
* **The latent AR process was in the model but not in `fitted()`.** Residuals
  were therefore taken against a fitted value that omitted it, so every
  residual diagnostic on such a model measured the structure the model had
  already accounted for. Measured: residual correlation at the shortest
  separations went *up* after fitting CAR(1), from 0.23 to 0.49, where it should
  fall toward zero. It now falls inside the envelope.
* **CAR(1) reproduces AR(1) exactly on equal spacing** -- same rho, same
  log-likelihood to 7e-11, from two separately written branches of the
  likelihood -- and matches `nlme`'s `corExp(nugget = TRUE)`, which is the
  marginal form of a latent process plus a residual, to four decimal places in
  log-likelihood and five in the coefficients.
* **Exact-lag matching does not survive irregular times.** On a 60-by-6 panel
  with times drawn from 1..30 it found 53 pairs at lag 1, five at lag 2 and none
  beyond, so five of six lags returned no verdict. Binning by separation uses
  all 900 pairs.

* **`rp_df = 1` reproduces the parametric special case exactly.** Against the
  three accelerated failure time families on the same data, the
  log-likelihoods agree to 1e-10 and the reparameterisation is exact to five
  decimals: the slope on log time is `1 / scale`, and each coefficient is
  `-beta / scale`. Those families are pinned to `survival::survreg()`, so the
  flexible one inherits that validation -- which matters, since neither
  `flexsurv` nor `rstpm2` is available here to check against directly.
* **The hazard ratio survives a baseline the model cannot draw.** On a
  piecewise-constant hazard -- 0.15 before time 3, then 0.9, which no
  parametric family in the package can bend to -- the coefficient came back
  0.711, 0.699 and 0.701 at `rp_df` of 3, 5 and 7, against a truth of 0.700.
* **And the survival-curve check caught the baseline that was too rigid.** At
  `rp_df = 3` it returned `FAIL` with 82% of the curve outside the envelope; at
  5, `OK` with 2%. The largest error in the fitted survivor function fell from
  0.12 to 0.06 over the same change.

* **With thousands of pairs per bin, significance stops meaning anything
  actionable.** On 300 points carrying an exponential field, adding
  `t2(sx, sy)` cut the residual correlation at the shortest distance from 0.242
  to -0.033, dropped AIC from 1165.7 to 1069.4 and narrowed the standard error
  on the covariate from 0.100 to 0.083 -- and four of five bins were still
  flagged, against an envelope 0.02 wide. Hence `min_effect`: the same stance
  the gaussian index takes against normality tests.

* **The dispersion model matches `glmmTMB` and `nlme`.** On a three-group
  design illume, `glmmTMB::glmmTMB(dispformula = ~ g)` and
  `nlme::gls(weights = varIdent())` all returned a log-likelihood of -1993.6974
  and standard errors agreeing to four decimals. On a power-of-the-mean design,
  illume and `nlme::gls(weights = varPower())` agreed on the power to two
  decimals and on the log-likelihood to 1e-3.
* **`~ mu` needs a warm start.** From a cold start the fitted mean is zero for
  every row, `log|mu|` is the log of the numerical floor, and the power
  multiplying it is unidentified there: on data with a true power of 1.0 the fit
  converged falsely with the power at -115818 and the slope collapsed to zero.
  A dispersion model is now started from the ordinary fit.
* **Applying the remedy clears the diagnosis.** On a 900-row design with spreads
  of 0.5 and 2.0 by stratum, `ilm_check_variance()` went from `WARN` with a
  variance ratio of 12.01 to `OK` with 1.06, and AIC fell from 3274.6 to 2653.6.

* **The censored fit matches `survival::survreg()`** to five decimal places in
  the coefficients, 1e-5 in the scale and 1e-6 in the log-likelihood, on both a
  ceiling and a floor. On a 60-unit panel censored at 18% -- which `survreg()`
  cannot fit, having no random effects -- ignoring the ceiling attenuated the
  slope from 1.02 to 0.79 against a truth of 1.0, and the between-unit standard
  deviation from 0.78 to 0.65 against a truth of 0.8.
* **The residual uniformity check does not detect ignored censoring.** Fitting
  the same data without `censor` left the quantile residuals uniform
  (Kolmogorov-Smirnov p = 0.18), because the fitted means vary from row to row
  and spread the pile-up out. The pile-up is visible in the data, not in the
  residuals, which is why `ilm_describe()` now names it.

* **All three AFT families match `survival::survreg()`** to about 1e-7 in the
  coefficients and 1e-10 in the log-likelihood, with standard errors agreeing to
  five decimal places, at 41% censoring. On a study `survreg()` cannot fit --
  80 centres, 44% censored -- a frailty Weibull recovered 0.746 / 0.597 / 0.465
  against truths of 0.7 / 0.6 / 0.5, while dropping the centre term pushed the
  between-centre spread into the scale, which rose from 0.597 to 0.729.
* **A simulated replicate has to be censored the way the study was.** Drawing a
  full event time for every subject and then marking the originally censored
  ones as censored at it is not the same process: measured on a 500-subject
  study, 141 of 190 censored subjects drew a time beyond their own censoring
  time, and the simulated Kaplan-Meier sat above the observed one at every point
  (0.13 against 0.00 in the tail). Every one of the four survival curves shown
  to the first version of the check came back `FAIL`, the correctly specified
  model included. Censoring times are now carried per subject -- known for those
  censored, drawn from the reverse Kaplan-Meier conditioned on exceeding the
  event time for the rest, capped at the end of follow-up -- and the check
  passes correctly specified models and fails the wrong family.
* **A threshold on the gap between the curves would be backwards.** On a
  correctly specified 500-subject model the largest vertical gap to the
  Kaplan-Meier was 0.121, with none of the curve outside the envelope. A Weibull
  fitted to log-logistic data gave a *smaller* gap, 0.093, with 26% outside.
* **The natural AFT residual is not usable when anything is censored.** The
  error on the log-time scale, `(log t - eta) / scale`, is evaluated at the
  censoring time rather than at the event, so it sits systematically low: at 40%
  censoring the mean came to -0.64, -0.84 and -0.56 for the three families. The
  normal score of the quantile residual fills the interval in and stays centred
  and unit-scaled, which is what the diagnostics now use.

* **All-pairs comparisons match `TukeyHSD()` where Tukey is exactly right,
  and hold up where it is not.** On normal, equal-variance, balanced data the
  simultaneous endpoints agree to 0.006 and the adjusted p-values to 0.011.
  Family-wise error over six comparisons at group sizes of 45 to 80 came to
  0.068 with a common variance and 0.072 with variances differing fourfold,
  against a nominal 0.05, where six unadjusted intervals reject 0.302 of the
  time.
* **Under heavy skew the marginal interval fails before the adjustment does.**
  On lognormal data with a spread parameter up to 1.2 at those group sizes, a
  single unadjusted comparison already erred 0.090 of the time, and no interval
  shape moved it: percentile 0.090, BCa 0.089, basic 0.093, normal 0.084. It is
  the bootstrapped mean of a small skewed sample, not the multiplicity, and it
  converges -- at group sizes of 300 to 500 the same design gives 0.049 per
  comparison and 0.062 family-wise. `stat = "median"` is the remedy: on that
  same design, under a null placing the medians together, 0.043 and 0.050.
* **A centred smooth changes what the intercept estimates.** The gaussian
  smooth cell returned an intercept covered 0.880 of the time while every slope
  was nominal, which looked like a defect in the fit and was not. `mgcv`'s
  identifiability constraint centres the basis on the observed data, so the
  intercept is the mean response at the sample average of the smooth rather
  than at its population average, and that average moves from one replicate to
  the next. The sample mean had a standard deviation of 0.0446 across
  replicates against a reported standard error of 0.0766, and
  `sqrt(0.0766^2 + 0.0446^2) = 0.0886` against an observed spread of 0.0875.
  Scored against the moving target, coverage was 0.943. `ilm_model()` now says
  so in its documentation; the slopes are unaffected either way.

## Fixes
* **A reduced-model refit was not the same model.** `ilm_refit_drops()` hands
  its workers a stub of the fitted object rather than the whole thing, and that
  stub carried the family and the weights and nothing else. For a censored
  model the reduced likelihood was therefore computed as if nothing were
  censored -- two likelihoods on different scales, differenced -- and
  `ilm_anova(test = "LRT")` rejected at **100%** under the null. The same held
  for a dispersion model, and a flexible parametric baseline errored outright.
  Everything that makes a model what it is now travels through one function,
  `ilm_refit_like()`, so there is a single place to add the next structure to.
  Measured after the fix, at 200 replicates and a nominal 0.05: 0.045, 0.050,
  0.035 and 0.055 for a Tobit, an accelerated failure time, a dispersion model
  and a flexible baseline.
* `ilm_model()` kept the caller's formula environment. `lme4::nobars()`,
  `mgcv::interpret.gam()` and `stats::reformulate()` each return a formula
  carrying an environment of their own, so a term with a local argument --
  `ilm_fourier(t, 12, K)`, `ns(x, df = d)` -- failed with "object not found"
  when the formula was built inside a function.
* The stored model terms now carry `predvars`. `stats::model.frame()` computes
  them and `stats::terms(formula, data = )` does not, so `predict()` was
  rebuilding `ns()`, `poly()` and `scale()` from whatever rows it was handed:
  `ns()` errored, and `poly()` and `scale()` returned quietly wrong numbers.
* `ilm_describe()` accepts a data frame again. With no `y` it handed the whole
  frame to the categorical branch, which failed with "the condition has length
  > 1" for any frame of more than one column, and returned nonsense rather than
  failing for a frame of exactly one. It now describes every column when they
  are all of a kind, takes several column names, and names
  `ilm_describe_all()` when the kinds are mixed.
* `ilm_plot_model(what = "effect")` draws one curve per interaction partner. It
  previously held the partner at its most common level, drew one of several
  quite different slopes, and did not say so.
* `ilm_rqr()`, `ilm_pearson_ovr()` and `ilm_appraise()` are family-aware
  throughout. The quantile residual errored for gaussian and Poisson models and
  returned plausible-looking but non-uniform values for binomial ones; a test
  that checked only `is.finite()` passed on all of it.

## Validation

Summary of the simulation studies re-run against this version on 2026-09-20.
Full tables and per-cell detail are in `studies/findings/`; that material is
intended for the methods paper rather than for this file.

* **Coverage is nominal for every model type the package fits.** The study grew
  from 11 designs to 23, adding the negative binomial, a random slope, a
  penalised smooth, AR(1) and CAR(1), a Tobit ceiling, all three accelerated
  failure time families with and without a frailty, a flexible parametric
  baseline and a dispersion model. Coverage ran 0.944 to 0.955 against a
  nominal 0.95 at 2000 replicates per cell, with the worst single coefficient
  0.935 in the deliberately under-powered `mn_J5_thin`.
* **Type I error is nominal across 13 designs** -- 0.042 to 0.061 against a
  nominal 0.05, Monte Carlo standard error 0.005, Wald and likelihood-ratio
  agreeing throughout. The cells covering censoring, time to event, a flexible
  baseline and a dispersion model are why this table was expanded: the
  reduced-model refit fault above rejected at 100% under the null, and no
  coverage study of point estimates would ever have seen it.
* **Convergence is now judged on the gradient, not on the optimiser's stopping
  code.** nlminb's "false convergence (8)" means it could not verify a descent
  direction, which is not the same as being away from a stationary point: on
  flexible parametric fits a third of replicates reported it while sitting at a
  gradient of 3.7e-03 and recovering the same coefficients as the replicates
  that reported success. Grading on the code discarded them, 1212 of 2000
  against 1991 on the gradient. The hardest multinomial cell moved the same
  way, 0.332 to 0.416 at 1.5 observations per latent value, with coverage
  unchanged.
* **Agreement with independent implementations holds** at 2.3e-04 or better
  against lme4, glmmTMB and nnet, and the mclogit comparison reproduces: 0.949
  coverage against 0.891 at four observations per cluster, where PQL attenuates
  fixed effects to 73% of their true magnitude.
* **Agreement with brms is unchanged** -- a mean of 0.097 standard errors over
  8 coefficients, a standard error ratio of 0.984, illume at a median 0.59s
  against 45.2s of brms sampling alone.

Coverage in cells with convergence failures is conditional on convergence:
failed fits are excluded, so surviving coverage is optimistic if failure
correlates with extreme estimates.

# illume 0.0.1.9000

## Validation

Summary of the simulation studies run against this version on 2026-09-19. Full
tables and per-cell detail are in `studies/findings/`; that material is intended
for the methods paper rather than for this file.

* **Coverage of Wald confidence intervals is nominal** across 11 designs
  spanning five families, AR(1) correlation and deliberately under-powered
  cells: 0.947 to 0.954 against a nominal 0.95, from 2000 replicates per cell.
* **The latent budget governs convergence, not coverage.** Convergence ran
  33% / 58% / 83% / 98% at 1.5 / 3 / 5 / 10 observations per latent value, while
  coverage among converged fits stayed nominal throughout. The approximation
  declines to return an answer rather than returning an overconfident one.
* **Type I error is nominal** for both Wald and likelihood-ratio tests
  (0.045–0.055 at a nominal 0.05), and null p-values are uniform by
  Kolmogorov–Smirnov test, not merely correct at the 5% threshold.
* **Agreement with independent implementations** is at numerical tolerance:
  1e-15 to 7e-5 against lme4, glmmTMB and nnet, and within 0.10 standard errors
  of brms at roughly 82x the speed.
* **Against `mclogit::mblogit()`**, the only other R package fitting
  random-effects multinomial models, illume holds nominal coverage where PQL
  does not: 0.949 versus 0.891 with four observations per cluster, where PQL
  attenuates fixed effects to 73% of their true magnitude. Note that mclogit
  attains a lower RMSE there, since shrinkage trades bias for variance — a
  defensible trade for prediction but not for inference.

Coverage figures in cells with convergence failures are conditional on
convergence: failed fits are excluded, so surviving coverage is optimistic if
failure correlates with extreme estimates.
