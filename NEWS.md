# illume 0.0.1.9000

First working version. One engine, `ilm_model()`, fits gaussian, binomial,
Poisson, negative binomial and multinomial models by Laplace approximation via
RTMB, with lme4-style formulas, crossed and nested random intercepts and slopes,
mgcv penalised smooths and AR(1) correlation.

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

## Fixes

* `ilm_model()` kept the caller's formula environment. `lme4::nobars()`,
  `mgcv::interpret.gam()` and `stats::reformulate()` each return a formula
  carrying an environment of their own, so a term with a local argument --
  `ilm_fourier(t, 12, K)`, `ns(x, df = d)` -- failed with "object not found"
  when the formula was built inside a function.
* The stored model terms now carry `predvars`. `stats::model.frame()` computes
  them and `stats::terms(formula, data = )` does not, so `predict()` was
  rebuilding `ns()`, `poly()` and `scale()` from whatever rows it was handed:
  `ns()` errored, and `poly()` and `scale()` returned quietly wrong numbers.
* `ilm_plot_model(what = "effect")` draws one curve per interaction partner. It
  previously held the partner at its most common level, drew one of several
  quite different slopes, and did not say so.
* `ilm_rqr()`, `ilm_pearson_ovr()` and `ilm_appraise()` are family-aware
  throughout. The quantile residual errored for gaussian and Poisson models and
  returned plausible-looking but non-uniform values for binomial ones; a test
  that checked only `is.finite()` passed on all of it.

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
