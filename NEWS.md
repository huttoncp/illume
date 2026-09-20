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
