# Package index

## Fitting models

One engine, every family. The formula is lme4’s and the object that
comes back is the same whatever was fitted.

- [`ilm_ar1()`](https://huttoncp.github.io/illume/reference/ilm_ar1.md)
  : A first-order autoregressive structure over evenly spaced time
- [`ilm_car1()`](https://huttoncp.github.io/illume/reference/ilm_car1.md)
  : A continuous-time autoregressive structure over irregular time
- [`ilm_censor()`](https://huttoncp.github.io/illume/reference/ilm_censor.md)
  : Mark censored observations
- [`ilm_cyclic()`](https://huttoncp.github.io/illume/reference/ilm_cyclic.md)
  : Cyclic cubic spline basis for a cycle of known period
- [`ilm_family()`](https://huttoncp.github.io/illume/reference/ilm_family.md)
  : Describes how a response distribution enters the model. The
  random-effect machinery – covariance structures, random slopes,
  smooths, AR(1) – is identical for every family and never sees the
  response; a family supplies only the piece that turns a linear
  predictor into a log-likelihood.
- [`ilm_fit()`](https://huttoncp.github.io/illume/reference/ilm_fit.md)
  : Fit a multinomial mixed model from a design matrix
- [`ilm_fourier()`](https://huttoncp.github.io/illume/reference/ilm_fourier.md)
  : Fourier terms for a cycle of known period
- [`ilm_model()`](https://huttoncp.github.io/illume/reference/ilm_model.md)
  [`ilm_model_formula()`](https://huttoncp.github.io/illume/reference/ilm_model.md)
  : Fit a multinomial linear mixed model
- [`ilm_squeeze()`](https://huttoncp.github.io/illume/reference/ilm_squeeze.md)
  : Move proportions off 0 and 1
- [`ilm_surv()`](https://huttoncp.github.io/illume/reference/ilm_surv.md)
  : Right-censored follow-up, in the survival convention
- [`ilm_thresholds()`](https://huttoncp.github.io/illume/reference/ilm_thresholds.md)
  : Thresholds of a cumulative link fit

## Other designs

Identification strategies and sampling designs that are not a single
regression.

- [`ilm_design()`](https://huttoncp.github.io/illume/reference/ilm_design.md)
  : Describe how a sample was drawn
- [`ilm_did()`](https://huttoncp.github.io/illume/reference/ilm_did.md)
  : Difference in differences
- [`ilm_iv()`](https://huttoncp.github.io/illume/reference/ilm_iv.md) :
  Instrumental-variables regression by two-stage least squares
- [`ilm_iv_ar()`](https://huttoncp.github.io/illume/reference/ilm_iv_ar.md)
  : Anderson-Rubin confidence set for an instrumented coefficient
- [`ilm_mediate()`](https://huttoncp.github.io/illume/reference/ilm_mediate.md)
  : How much of an effect runs through a mediator
- [`ilm_mediate_sens()`](https://huttoncp.github.io/illume/reference/ilm_mediate_sens.md)
  : How strong would unmeasured confounding have to be?
- [`ilm_rdd()`](https://huttoncp.github.io/illume/reference/ilm_rdd.md)
  : Regression discontinuity
- [`ilm_svy_coef()`](https://huttoncp.github.io/illume/reference/ilm_svy_coef.md)
  : Fixed effects with design-based standard errors
- [`ilm_svy_vcov()`](https://huttoncp.github.io/illume/reference/ilm_svy_vcov.md)
  : Linearization variance for a fit from a complex sample

## Diagnostics

Each check reports whether an assumption is consistent with the data,
and names a remedy that exists in this package when it is not.

- [`ilm_appraise()`](https://huttoncp.github.io/illume/reference/ilm_appraise.md)
  : Diagnostic plots for a fitted model
- [`ilm_binned_residuals()`](https://huttoncp.github.io/illume/reference/ilm_binned_residuals.md)
  : Binned residuals for a binary response
- [`ilm_calibration()`](https://huttoncp.github.io/illume/reference/ilm_calibration.md)
  : Calibration of predicted probabilities
- [`ilm_check_ar()`](https://huttoncp.github.io/illume/reference/ilm_check_ar.md)
  : Check for leftover temporal autocorrelation
- [`ilm_check_collinearity()`](https://huttoncp.github.io/illume/reference/ilm_check_collinearity.md)
  : Variance inflation among fixed effects
- [`ilm_check_covariate()`](https://huttoncp.github.io/illume/reference/ilm_check_covariate.md)
  : Test whether residuals shift across a variable
- [`ilm_check_dispersion()`](https://huttoncp.github.io/illume/reference/ilm_check_dispersion.md)
  : Is the response more variable than the model allows?
- [`ilm_check_missing()`](https://huttoncp.github.io/illume/reference/ilm_check_missing.md)
  : What the missing values look like, and whether they matter
- [`ilm_check_omitted()`](https://huttoncp.github.io/illume/reference/ilm_check_omitted.md)
  : Screen variables the model does not use
- [`ilm_check_predictive()`](https://huttoncp.github.io/illume/reference/ilm_check_predictive.md)
  : Compare the observed response with data the model would generate
- [`ilm_check_proportional()`](https://huttoncp.github.io/illume/reference/ilm_check_proportional.md)
  : Does one coefficient per predictor describe every cut?
- [`ilm_check_variance()`](https://huttoncp.github.io/illume/reference/ilm_check_variance.md)
  : Is the residual spread constant?
- [`ilm_check_zeros()`](https://huttoncp.github.io/illume/reference/ilm_check_zeros.md)
  : Are there more zeros than the model expects?
- [`ilm_consistency()`](https://huttoncp.github.io/illume/reference/ilm_consistency.md)
  : Check whether the model can recover itself
- [`ilm_gauss_check()`](https://huttoncp.github.io/illume/reference/ilm_gauss_check.md)
  : How far a variable is from gaussian, and why
- [`ilm_re_mahalanobis()`](https://huttoncp.github.io/illume/reference/ilm_re_mahalanobis.md)
  : Distances of the fitted random effects from zero
- [`ilm_rqr()`](https://huttoncp.github.io/illume/reference/ilm_rqr.md)
  : Randomised quantile residuals for a nominal outcome
- [`ilm_rqr_test()`](https://huttoncp.github.io/illume/reference/ilm_rqr_test.md)
  : Test residuals against a simulated reference
- [`ilm_scores()`](https://huttoncp.github.io/illume/reference/ilm_scores.md)
  : Proper scoring rules for a categorical outcome
- [`ilm_variogram()`](https://huttoncp.github.io/illume/reference/ilm_variogram.md)
  : Residual correlation as a function of separation in time

## ANOVA

Factorial and repeated-measures designs, specified by naming columns
rather than by writing a formula with an error term.

- [`ilm_aov_ez()`](https://huttoncp.github.io/illume/reference/ilm_aov_ez.md)
  : Factorial and repeated-measures ANOVA

## Inference and interpretation

What the model says, on a scale someone can read.

- [`ilm_ame()`](https://huttoncp.github.io/illume/reference/ilm_ame.md)
  : Average marginal effect on the response scale
- [`ilm_anova()`](https://huttoncp.github.io/illume/reference/ilm_anova.md)
  : Analysis of deviance for fixed effects
- [`ilm_boot_ci()`](https://huttoncp.github.io/illume/reference/ilm_boot_ci.md)
  : Bootstrap confidence interval for a statistic
- [`ilm_boot_diff()`](https://huttoncp.github.io/illume/reference/ilm_boot_diff.md)
  : Bootstrap intervals for differences between groups
- [`ilm_coef_table()`](https://huttoncp.github.io/illume/reference/ilm_coef_table.md)
  : Coefficient table with Wald tests
- [`ilm_contrast()`](https://huttoncp.github.io/illume/reference/ilm_contrast.md)
  : Compare estimated marginal means
- [`ilm_denom_df()`](https://huttoncp.github.io/illume/reference/ilm_denom_df.md)
  : Denominator degrees of freedom for one or more contrasts
- [`ilm_effects()`](https://huttoncp.github.io/illume/reference/ilm_effects.md)
  : Effect sizes for the fixed effects
- [`ilm_emmeans()`](https://huttoncp.github.io/illume/reference/ilm_emmeans.md)
  : Estimated marginal means
- [`ilm_interpret()`](https://huttoncp.github.io/illume/reference/ilm_interpret.md)
  : Interpret a fitted model in words
- [`ilm_moderation()`](https://huttoncp.github.io/illume/reference/ilm_moderation.md)
  : Does the effect of X depend on something else?
- [`ilm_pb_lrt()`](https://huttoncp.github.io/illume/reference/ilm_pb_lrt.md)
  : Parametric-bootstrap likelihood-ratio test
- [`ilm_plot_moderation()`](https://huttoncp.github.io/illume/reference/ilm_plot_moderation.md)
  : See a moderation
- [`ilm_robust()`](https://huttoncp.github.io/illume/reference/ilm_robust.md)
  : Fixed effects with sandwich standard errors
- [`ilm_rp_lrt()`](https://huttoncp.github.io/illume/reference/ilm_rp_lrt.md)
  : Was the flexible baseline worth it?
- [`ilm_scenario()`](https://huttoncp.github.io/illume/reference/ilm_scenario.md)
  : What the model says would happen under specified scenarios
- [`ilm_translate()`](https://huttoncp.github.io/illume/reference/ilm_translate.md)
  : Recode a variable against a dictionary held as two vectors
- [`ilm_trends()`](https://huttoncp.github.io/illume/reference/ilm_trends.md)
  : Estimated marginal slopes
- [`ilm_vcov_cluster()`](https://huttoncp.github.io/illume/reference/ilm_vcov_cluster.md)
  : Cluster-robust and heteroskedasticity-robust covariance
- [`ilm_zi_coef()`](https://huttoncp.github.io/illume/reference/ilm_zi_coef.md)
  : Coefficients of the zero part of a zero-inflated or hurdle fit
- [`ilm_zi_prob()`](https://huttoncp.github.io/illume/reference/ilm_zi_prob.md)
  : Fitted excess-zero probability, one value per row

## Design and power

Before the data exist.

- [`ilm_scaffold()`](https://huttoncp.github.io/illume/reference/ilm_scaffold.md)
  : A model built from assumptions instead of data
- [`ilm_power_design()`](https://huttoncp.github.io/illume/reference/ilm_power_design.md)
  : Power for a study that has not been run
- [`ilm_power()`](https://huttoncp.github.io/illume/reference/ilm_power.md)
  : Power for a term, by simulating from a fitted model
- [`ilm_power_n()`](https://huttoncp.github.io/illume/reference/ilm_power_n.md)
  : The sample size a target power implies
- [`plot(`*`<ilm_power>`*`)`](https://huttoncp.github.io/illume/reference/plot.ilm_power.md)
  : Power curve

## Describing data

Descriptive statistics, counts, and the things that are wrong with a
data frame before any model sees it.

- [`ilm_copies()`](https://huttoncp.github.io/illume/reference/ilm_copies.md)
  : Find copied or duplicated rows
- [`ilm_counts()`](https://huttoncp.github.io/illume/reference/ilm_counts.md)
  : Frequency counts of a vector's unique values
- [`ilm_counts_all()`](https://huttoncp.github.io/illume/reference/ilm_counts_all.md)
  : Frequency counts for every column
- [`ilm_counts_tb()`](https://huttoncp.github.io/illume/reference/ilm_counts_tb.md)
  : The most and least frequent values, side by side
- [`ilm_counts_tb_all()`](https://huttoncp.github.io/illume/reference/ilm_counts_tb_all.md)
  : Most and least frequent values for every column
- [`ilm_describe()`](https://huttoncp.github.io/illume/reference/ilm_describe.md)
  : Class-aware description of one variable
- [`ilm_describe_all()`](https://huttoncp.github.io/illume/reference/ilm_describe_all.md)
  : Describe every column of a data frame
- [`ilm_dupes()`](https://huttoncp.github.io/illume/reference/ilm_dupes.md)
  : Duplicated rows only
- [`ilm_frame_issues()`](https://huttoncp.github.io/illume/reference/ilm_frame_issues.md)
  : Problems that belong to pairs of columns
- [`ilm_recode_errors()`](https://huttoncp.github.io/illume/reference/ilm_recode_errors.md)
  : Replace known-bad values with NA or another value
- [`ilm_recode_errors_vec()`](https://huttoncp.github.io/illume/reference/ilm_recode_errors_vec.md)
  : Replace known-bad values in a vector
- [`ilm_wash_df()`](https://huttoncp.github.io/illume/reference/ilm_wash_df.md)
  : Clean up a messy data frame

## Outliers and anomalies

A value extreme for its own column, against a row implausible as a
combination.

- [`ilm_anomaly()`](https://huttoncp.github.io/illume/reference/ilm_anomaly.md)
  : Rows that do not fit the pattern the other rows make
- [`ilm_anomalous()`](https://huttoncp.github.io/illume/reference/ilm_anomalous.md)
  : The rows an anomaly scan flagged
- [`ilm_plot_anomaly()`](https://huttoncp.github.io/illume/reference/ilm_plot_anomaly.md)
  : See an anomaly scan
- [`ilm_outliers()`](https://huttoncp.github.io/illume/reference/ilm_outliers.md)
  : Flag unusual values in a numeric vector
- [`ilm_outliers_all()`](https://huttoncp.github.io/illume/reference/ilm_outliers_all.md)
  : Flag unusual values across a data frame

## Structure: profile, cluster, reduce

ilm_profile() is the front door: it reduces, clusters and describes the
groups in one call, and handles mixed columns. The other two are the
same pipeline taken a step at a time, for when you want the coordinates
or the partition on their own.

- [`ilm_profile()`](https://huttoncp.github.io/illume/reference/ilm_profile.md)
  : Profile a data set: reduce, cluster, and describe the clusters
- [`ilm_var_contrib()`](https://huttoncp.github.io/illume/reference/ilm_var_contrib.md)
  : What each variable contributes to a clustering
- [`ilm_cluster()`](https://huttoncp.github.io/illume/reference/ilm_cluster.md)
  : Cluster observations, choosing the number of clusters
- [`ilm_reduce()`](https://huttoncp.github.io/illume/reference/ilm_reduce.md)
  : Reduce a data frame's variables to a few dimensions
- [`ilm_glrm()`](https://huttoncp.github.io/illume/reference/ilm_glrm.md)
  : A low-rank model with a loss chosen per column
- [`ilm_profile_na()`](https://huttoncp.github.io/illume/reference/ilm_profile_na.md)
  : Profile which values are missing, and for whom
- [`ilm_cluster_na()`](https://huttoncp.github.io/illume/reference/ilm_cluster_na.md)
  : Cluster observations by which values they are missing
- [`ilm_reduce_na()`](https://huttoncp.github.io/illume/reference/ilm_reduce_na.md)
  : Reduce a data frame's missingness pattern to a few dimensions

## Missing data

Diagnosing it, filling it in honestly, and pooling across the
imputations.

- [`ilm_describe_na()`](https://huttoncp.github.io/illume/reference/ilm_describe_na.md)
  : Missingness in one variable
- [`ilm_describe_na_all()`](https://huttoncp.github.io/illume/reference/ilm_describe_na_all.md)
  : Missingness in every variable
- [`ilm_impute()`](https://huttoncp.github.io/illume/reference/ilm_impute.md)
  : Multiple imputation by chained equations
- [`ilm_mi_pool()`](https://huttoncp.github.io/illume/reference/ilm_mi_pool.md)
  : Fit a model across imputations and pool the results
- [`ilm_plot_na()`](https://huttoncp.github.io/illume/reference/ilm_plot_na.md)
  : Missing values in one column, across groups
- [`ilm_plot_na_all()`](https://huttoncp.github.io/illume/reference/ilm_plot_na_all.md)
  : Missing values by column

## Causal models

A DAG, what it implies, and what it licenses you to say.

- [`ilm_adjust_sets()`](https://huttoncp.github.io/illume/reference/ilm_adjust_sets.md)
  : Minimal adjustment sets for an exposure effect
- [`ilm_dag()`](https://huttoncp.github.io/illume/reference/ilm_dag.md)
  : A causal graph
- [`ilm_dag_implied()`](https://huttoncp.github.io/illume/reference/ilm_dag_implied.md)
  : Conditional independencies the graph implies
- [`ilm_dag_model()`](https://huttoncp.github.io/illume/reference/ilm_dag_model.md)
  : Fit the model a causal graph implies
- [`ilm_dag_test()`](https://huttoncp.github.io/illume/reference/ilm_dag_test.md)
  : Test what the graph claims against the data
- [`ilm_dsep()`](https://huttoncp.github.io/illume/reference/ilm_dsep.md)
  : Are two variables d-separated given a conditioning set?

## Plots

Built on tinyplot, named for what they show.

- [`ilm_geom_spec()`](https://huttoncp.github.io/illume/reference/ilm_geom_spec.md)
  : What each geom requires
- [`ilm_pick_geom()`](https://huttoncp.github.io/illume/reference/ilm_pick_geom.md)
  : Which geom would be drawn, and why
- [`ilm_plot()`](https://huttoncp.github.io/illume/reference/ilm_plot.md)
  : Adaptive plot of one or two variables
- [`ilm_plot_acf()`](https://huttoncp.github.io/illume/reference/ilm_plot_acf.md)
  : Autocorrelation and partial-autocorrelation plots with verdicts
- [`ilm_plot_all()`](https://huttoncp.github.io/illume/reference/ilm_plot_all.md)
  : Small multiples of every variable
- [`ilm_plot_bar()`](https://huttoncp.github.io/illume/reference/ilm_plot_bar.md)
  : Bar plot
- [`ilm_plot_boot_diff()`](https://huttoncp.github.io/illume/reference/ilm_plot_boot_diff.md)
  : The bootstrap distribution behind a group difference
- [`ilm_plot_box()`](https://huttoncp.github.io/illume/reference/ilm_plot_box.md)
  : Boxplot
- [`ilm_plot_c()`](https://huttoncp.github.io/illume/reference/ilm_plot_c.md)
  : Combine several plots into one figure
- [`ilm_plot_cluster()`](https://huttoncp.github.io/illume/reference/ilm_plot_cluster.md)
  [`ilm_plot_cluster_na()`](https://huttoncp.github.io/illume/reference/ilm_plot_cluster.md)
  : Map the clusters
- [`ilm_plot_cluster_gap()`](https://huttoncp.github.io/illume/reference/ilm_plot_cluster_gap.md)
  [`ilm_plot_cluster_gap_na()`](https://huttoncp.github.io/illume/reference/ilm_plot_cluster_gap.md)
  : The gap statistic across every k considered
- [`ilm_plot_density()`](https://huttoncp.github.io/illume/reference/ilm_plot_density.md)
  : Density plot
- [`ilm_plot_did()`](https://huttoncp.github.io/illume/reference/ilm_plot_did.md)
  : Event-study plot for a difference in differences
- [`ilm_plot_histogram()`](https://huttoncp.github.io/illume/reference/ilm_plot_histogram.md)
  : Histogram
- [`ilm_plot_line()`](https://huttoncp.github.io/illume/reference/ilm_plot_line.md)
  : Line plot
- [`ilm_plot_missing()`](https://huttoncp.github.io/illume/reference/ilm_plot_missing.md)
  : Proportion missing, by variable
- [`ilm_plot_model()`](https://huttoncp.github.io/illume/reference/ilm_plot_model.md)
  : Diagnostic and summary plots for a fitted model
- [`ilm_plot_profile()`](https://huttoncp.github.io/illume/reference/ilm_plot_profile.md)
  [`ilm_plot_profile_na()`](https://huttoncp.github.io/illume/reference/ilm_plot_profile.md)
  : Map the clusters from a profile
- [`ilm_plot_rdd()`](https://huttoncp.github.io/illume/reference/ilm_plot_rdd.md)
  : Regression discontinuity plot
- [`ilm_plot_reduce()`](https://huttoncp.github.io/illume/reference/ilm_plot_reduce.md)
  [`ilm_plot_reduce_na()`](https://huttoncp.github.io/illume/reference/ilm_plot_reduce.md)
  : Map the observations from a reduction
- [`ilm_plot_reduce_contrib()`](https://huttoncp.github.io/illume/reference/ilm_plot_reduce_contrib.md)
  [`ilm_plot_reduce_contrib_na()`](https://huttoncp.github.io/illume/reference/ilm_plot_reduce_contrib.md)
  : Which variables a dimension is made of
- [`ilm_plot_reduce_scree()`](https://huttoncp.github.io/illume/reference/ilm_plot_reduce_scree.md)
  [`ilm_plot_reduce_scree_na()`](https://huttoncp.github.io/illume/reference/ilm_plot_reduce_scree.md)
  : Scree plot for a reduction
- [`ilm_plot_scatter()`](https://huttoncp.github.io/illume/reference/ilm_plot_scatter.md)
  : Scatter plot
- [`ilm_plot_stat_error()`](https://huttoncp.github.io/illume/reference/ilm_plot_stat_error.md)
  : Group means or medians with an error bar
- [`ilm_plot_survival()`](https://huttoncp.github.io/illume/reference/ilm_plot_survival.md)
  : Survival curve against the Kaplan-Meier estimate
- [`ilm_plot_var()`](https://huttoncp.github.io/illume/reference/ilm_plot_var.md)
  : Plot one or two variables, choosing the geometry
- [`ilm_plot_var_all()`](https://huttoncp.github.io/illume/reference/ilm_plot_var_all.md)
  : Plot every column of a data frame
- [`ilm_plot_var_pairs()`](https://huttoncp.github.io/illume/reference/ilm_plot_var_pairs.md)
  : Pairwise plots
- [`ilm_plot_variogram()`](https://huttoncp.github.io/illume/reference/ilm_plot_variogram.md)
  : Plot a residual variogram
- [`ilm_plot_violin()`](https://huttoncp.github.io/illume/reference/ilm_plot_violin.md)
  : Violin plot

## Simulation

Data with a known structure, and draws from a fitted model.

- [`ilm_fitted()`](https://huttoncp.github.io/illume/reference/ilm_fitted.md)
  : Fitted category probabilities
- [`ilm_sim()`](https://huttoncp.github.io/illume/reference/ilm_sim.md)
  : A simulated mixed-type dataset for testing and examples
- [`ilm_simulate()`](https://huttoncp.github.io/illume/reference/ilm_simulate.md)
  : Simulate new outcomes from a fitted model
- [`ilm_survival()`](https://huttoncp.github.io/illume/reference/ilm_survival.md)
  : Predicted survival curve from a fitted model

## Working with other packages

Registration shims and the methods that let the wider ecosystem dispatch
on an illume fit.

- [`AIC(`*`<ilm_model>`*`)`](https://huttoncp.github.io/illume/reference/AIC.ilm_model.md)
  [`BIC(`*`<ilm_model>`*`)`](https://huttoncp.github.io/illume/reference/AIC.ilm_model.md)
  : Information criteria
- [`Anova.ilm_model()`](https://huttoncp.github.io/illume/reference/Anova.ilm_model.md)
  : car::Anova method
- [`check_model.ilm_model()`](https://huttoncp.github.io/illume/reference/check_model.ilm_model.md)
  : performance::check_model method
- [`coef(`*`<ilm_model>`*`)`](https://huttoncp.github.io/illume/reference/coef.ilm_model.md)
  [`vcov(`*`<ilm_model>`*`)`](https://huttoncp.github.io/illume/reference/coef.ilm_model.md)
  : Coefficients and their covariance matrix
- [`fixef(`*`<ilm_model>`*`)`](https://huttoncp.github.io/illume/reference/fixef.ilm_model.md)
  [`ilm_se_fixef()`](https://huttoncp.github.io/illume/reference/fixef.ilm_model.md)
  : Fixed effects as a matrix
- [`get_coef.ilm_model()`](https://huttoncp.github.io/illume/reference/marginaleffects-methods.md)
  [`set_coef.ilm_model()`](https://huttoncp.github.io/illume/reference/marginaleffects-methods.md)
  [`get_vcov.ilm_model()`](https://huttoncp.github.io/illume/reference/marginaleffects-methods.md)
  [`get_predict.ilm_model()`](https://huttoncp.github.io/illume/reference/marginaleffects-methods.md)
  : marginaleffects interface
- [`ilm_register_insight()`](https://huttoncp.github.io/illume/reference/ilm_register_insight.md)
  : Make illume models readable by the easystats packages
- [`ilm_register_marginaleffects()`](https://huttoncp.github.io/illume/reference/ilm_register_marginaleffects.md)
  : Register the marginaleffects interface
- [`logLik(`*`<ilm_model>`*`)`](https://huttoncp.github.io/illume/reference/logLik.ilm_model.md)
  [`nobs(`*`<ilm_model>`*`)`](https://huttoncp.github.io/illume/reference/logLik.ilm_model.md)
  : Log-likelihood, and number of observations
- [`model_performance.ilm_model()`](https://huttoncp.github.io/illume/reference/model_performance.ilm_model.md)
  : Fit indices for a multinomial mixed model
- [`plot(`*`<ilm_model>`*`)`](https://huttoncp.github.io/illume/reference/plot.ilm_model.md)
  : Plot a fitted model
- [`predict(`*`<ilm_model>`*`)`](https://huttoncp.github.io/illume/reference/predict.ilm_model.md)
  : Predicted category probabilities
- [`print(`*`<ilm_model>`*`)`](https://huttoncp.github.io/illume/reference/print.ilm_model.md)
  : Compact display of a fitted model
- [`summary(`*`<ilm_model>`*`)`](https://huttoncp.github.io/illume/reference/summary.ilm_model.md)
  [`print(`*`<summary.ilm_model>`*`)`](https://huttoncp.github.io/illume/reference/summary.ilm_model.md)
  : Summarise a fitted multinomial mixed model
- [`terms(`*`<ilm_model>`*`)`](https://huttoncp.github.io/illume/reference/ilm_model-accessors.md)
  [`formula(`*`<ilm_model>`*`)`](https://huttoncp.github.io/illume/reference/ilm_model-accessors.md)
  [`model.frame(`*`<ilm_model>`*`)`](https://huttoncp.github.io/illume/reference/ilm_model-accessors.md)
  [`model.matrix(`*`<ilm_model>`*`)`](https://huttoncp.github.io/illume/reference/ilm_model-accessors.md)
  [`getCall(`*`<ilm_model>`*`)`](https://huttoncp.github.io/illume/reference/ilm_model-accessors.md)
  : Standard accessors for a fitted model

## Shared parameters and print methods

Documentation shared across functions, and methods you call by printing
rather than by name.

- [`ilm_progress_arg`](https://huttoncp.github.io/illume/reference/ilm_progress_arg.md)
  : Progress reporting in illume
- [`ilm_selection`](https://huttoncp.github.io/illume/reference/ilm_selection.md)
  : How columns can be chosen
- [`print(`*`<ilm_censor>`*`)`](https://huttoncp.github.io/illume/reference/print.ilm_censor.md)
  : Print a censoring specification
- [`print(`*`<ilm_cor>`*`)`](https://huttoncp.github.io/illume/reference/print.ilm_cor.md)
  : Print a correlation structure
