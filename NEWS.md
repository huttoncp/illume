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
