# Benchmarking and validation

A package that fits models has to earn the numbers it prints. This
vignette is the evidence: what was measured, against what, and where
illume comes off worse.

Nothing here is computed when the vignette builds. The studies behind it
take hours on ten cores, so the results are tabulated from
`studies/findings/` in the source repository, and the scripts that
produced them are in `studies/scripts/`. Each table names the number of
replications it rests on, because a simulation result without one is an
anecdote.

## Why bother

Every defect found in this package so far was found from outside it. Not
by the test suite – by comparing a fit against `lme4`, `nlme`, `survreg`
or `mclogit`, by someone calling a function in a way the tests did not,
or by `R CMD check`. Tests confirm the package still does what it did
yesterday. They cannot tell you that what it did yesterday was right,
because they were written by the same person, on the same understanding,
as the code.

So validation here means an outside reference or a known answer, and the
studies below are organised around four questions:

1.  Does illume agree with implementations that are already trusted?
2.  Do its confidence intervals cover at the rate they claim?
3.  Is its type I error rate correct, and is its power real?
4.  For the one model where illume has a genuine competitor, which is
    better – and in what sense?

## 1. Agreement with established implementations

Where another package fits the same model, the two should agree to
numerical tolerance. This is the strongest available check: a bug that
survives it has to be a bug the reference implementation shares.

160 replications per task, each on a simulated dataset of 1,000
observations; R 4.4.3, illume 0.0.7.9000. `max_diff` is the largest
absolute coefficient difference against any reference fitting the same
model.

| task                     | illume (s) | lme4 (s) | glmmTMB (s) | nnet (s) | max_diff |
|--------------------------|------------|----------|-------------|----------|----------|
| gaussian mixed           | 0.173      | 0.026    | 0.194       | –        | 0.0000   |
| Poisson mixed            | 0.243      | 0.517    | 0.269       | –        | 0.0000   |
| binomial mixed           | 0.360      | 0.568    | 0.499       | –        | 0.0001   |
| multinomial fixed, J = 3 | 0.162      | –        | –           | 0.018    | 0.0011   |
| multinomial fixed, J = 5 | 0.226      | –        | –           | 0.039    | 0.0001   |
| multinomial mixed, J = 3 | 0.575      | –        | –           | –        | –        |
| multinomial mixed, J = 5 | 1.571      | –        | –           | –        | –        |

Largest disagreement with any independent implementation: **1.1e-03**,
on multinomial fixed effects against
[`nnet::multinom`](https://rdrr.io/pkg/nnet/man/multinom.html). The rest
agree to four decimal places or better.

Timing is context, not a contest – these are different optimisers
solving the same problem, and the differences are not large enough to
choose between packages on. `lme4` is faster on the Gaussian case,
illume is faster on the Poisson and binomial ones, and the em-dashes are
the point: the last two rows have no comparison because no other R
package fits a multinomial mixed model by maximum likelihood.

## 2. Do the intervals cover?

A 95% interval should contain the true value 95% of the time. This is
the property inference rests on, and it is cheap to check and easy to
get wrong.

24 model types, up to 2000 replications each. `se_ratio` is the mean
reported standard error divided by the actual standard deviation of the
estimates: 1 means the reported uncertainty is honest.

| model | replications | convergence | coverage | worst coefficient | se_ratio |
|----|----|----|----|----|----|
| Gaussian, fixed | 400 | 1.000 | 0.951 | 0.945 | 1.006 |
| Gaussian, mixed | 600 | 1.000 | 0.950 | 0.944 | 0.995 |
| Gaussian, random slope | 600 | 1.000 | 0.947 | 0.942 | 0.987 |
| Gaussian, smooth | 600 | 1.000 | 0.950 | 0.942 | 0.997 |
| Gaussian, AR(1) | 1200 | 1.000 | 0.945 | 0.943 | 0.984 |
| Gaussian, CAR(1) | 360 | 0.981 | 0.944 | 0.936 | 0.988 |
| binomial, mixed | 600 | 1.000 | 0.953 | 0.951 | 1.002 |
| Poisson, mixed | 600 | 1.000 | 0.949 | 0.945 | 0.996 |
| negative binomial, mixed | 600 | 1.000 | 0.946 | 0.942 | 0.980 |
| zero-inflated Poisson, mixed | 600 | 1.000 | 0.950 | 0.940 | 0.989 |
| hurdle Poisson | 600 | 1.000 | 0.949 | 0.944 | 1.003 |
| beta, fixed | 600 | 1.000 | 0.949 | 0.948 | 0.999 |
| zero-inflated beta | 600 | 1.000 | 0.951 | 0.946 | 0.991 |
| ordinal, fixed | 600 | 1.000 | 0.946 | 0.941 | 0.994 |
| ordinal, mixed | 600 | 1.000 | 0.953 | 0.951 | 1.007 |
| multinomial J = 3, rich | 600 | 0.999 | 0.952 | 0.948 | 0.999 |
| multinomial J = 3, thin | 300 | 0.953 | 0.951 | 0.945 | 0.978 |
| multinomial J = 3, AR(1) | 2000 | 0.768 | 0.954 | 0.949 | 1.009 |
| multinomial J = 5, rich | 600 | 0.844 | 0.949 | 0.943 | 0.982 |
| multinomial J = 5, mid | 360 | 0.611 | 0.951 | 0.939 | 0.980 |
| multinomial J = 5, thin | 240 | 0.416 | 0.948 | 0.935 | 0.954 |
| Tobit, ceiling | 600 | 1.000 | 0.945 | 0.938 | 0.993 |
| AFT Weibull | 600 | 1.000 | 0.944 | 0.940 | 0.976 |
| AFT log-normal | 600 | 1.000 | 0.948 | 0.946 | 0.993 |
| AFT log-logistic | 600 | 1.000 | 0.948 | 0.946 | 0.994 |
| AFT with frailty | 600 | 1.000 | 0.950 | 0.946 | 0.995 |
| Royston-Parmar flexible | 800 | 0.996 | 0.955 | 0.951 | 1.021 |
| dispersion by group | 600 | 1.000 | 0.953 | 0.940 | 1.013 |
| instrumental variables (2SLS) | 800 | 1.000 | 0.946 | 0.941 | 0.993 |
| mediation | 600 | 1.000 | 0.967 | 0.967 | 1.107 |
| survey, stratified | 480 | 1.000 | 0.939 | 0.939 | 0.956 |

Coverage runs from 0.939 to 0.967 across every cell, and the worst
single coefficient anywhere is 0.935. With 600 replications the Monte
Carlo error on a coverage estimate is about 0.009, so most of that
spread is noise.

Two cells are not noise and are worth naming:

**Mediation over-covers, at 0.967 with an `se_ratio` of 1.107.** The
indirect effect’s standard error is about 11% larger than the spread of
the estimates warrants, so its intervals are conservative. That is the
safe direction to be wrong in, but it is still wrong: a real indirect
effect will be declared non-significant more often than 5% of the time.
The cause is the delta method applied to a product of coefficients,
which is known to be conservative in small samples;
`ilm_mediate(boot = TRUE)` is the remedy and is what the mediation
documentation recommends.

**Stratified survey estimation under-covers slightly, at 0.939 with an
`se_ratio` of 0.956.** Reported standard errors are about 4% too small.
The interval is a little narrower than it should be.

### The caveat that matters most

**Coverage in the low-convergence cells is conditional on convergence.**
Failed fits are excluded from the table. The multinomial J = 5 thin cell
converges 41.6% of the time, and its 0.948 coverage describes only those
41.6%. If failure correlates with extreme estimates – and there is every
reason to think it does, since the fits that fail are the ones where the
likelihood is flattest – the surviving coverage is optimistic.

This is not a defect being hidden. It is the reason
[`ilm_model()`](https://huttoncp.github.io/illume/reference/ilm_model.md)
reports convergence status rather than burying it, and the reason the
latent budget diagnostic exists: a multinomial model with 1.5
observations per latent coefficient is being asked a question the data
cannot answer, and illume’s design choice is to fail loudly rather than
return a confident-looking number. See the *Regression models* vignette
for what the budget is and how to read it.

## 3. Size and power

Power is only meaningful if the type I error rate is right, so size is
the primary result. 2000 replications per cell, nominal 0.05.

| model                    | Wald  | LRT   |
|--------------------------|-------|-------|
| Gaussian, fixed          | 0.050 | 0.050 |
| Gaussian, mixed          | 0.054 | 0.054 |
| Gaussian, AR(1)          | 0.055 | 0.054 |
| Gaussian, CAR(1)         | 0.057 | 0.057 |
| binomial, mixed          | 0.048 | 0.050 |
| Poisson, mixed           | 0.053 | 0.053 |
| negative binomial, mixed | 0.054 | 0.054 |
| multinomial J = 3, rich  | 0.055 | 0.056 |
| multinomial J = 3, thin  | 0.045 | 0.055 |
| Tobit, ceiling           | 0.046 | 0.046 |
| AFT Weibull              | 0.042 | 0.043 |
| Royston-Parmar flexible  | 0.048 | 0.049 |
| dispersion by group      | 0.061 | 0.060 |

Monte Carlo error at 2000 replications is 0.005, so the range 0.042 to
0.061 is within about two standard errors of nominal everywhere except
the dispersion-by-group cell, which sits a little high at 0.061.

Power, measured as the rejection rate at the 5% level:

| model                    | d = 0.1 | d = 0.2 | d = 0.3 | d = 0.5 |
|--------------------------|---------|---------|---------|---------|
| Gaussian, fixed          | 0.508   | 0.981   | 1.000   | 1.000   |
| Gaussian, mixed          | 0.679   | 0.997   | 1.000   | 1.000   |
| Gaussian, AR(1)          | 0.887   | 1.000   | 1.000   | 1.000   |
| binomial, mixed          | 0.201   | 0.604   | 0.906   | 1.000   |
| Poisson, mixed           | 0.908   | 1.000   | 1.000   | 1.000   |
| negative binomial, mixed | 0.617   | 0.994   | 1.000   | 1.000   |
| multinomial J = 3, rich  | 0.308   | 0.869   | 0.996   | 1.000   |
| multinomial J = 3, thin  | 0.160   | 0.503   | 0.858   | 0.999   |

These are the numbers
[`ilm_power()`](https://huttoncp.github.io/illume/reference/ilm_power.md)
is checked against. See the *Effect size and power* vignette for the
analytic and simulation routes and how they were validated against
`Superpower`, `simr`, `faux` and `lmerTest`.

## 4. The multinomial mixed model

This is the only model in illume where another R package offers a
genuine alternative, so it gets the closest look.

`mclogit::mblogit()` fits random-effects baseline-category logit models
by penalised quasi-likelihood. PQL is fast and well understood, and is
known to attenuate fixed effects for discrete outcomes – worst when
clusters are small and the random-effect variance is large (Breslow &
Lin 1995; Rodriguez & Goldman 1995). illume uses a Laplace approximation
to the marginal likelihood instead. The question is whether that
theoretical difference shows up in practice and matters.

Comparison is on mclogit’s baseline-category scale, since illume’s
sum-to-zero coefficients are not directly comparable. Attenuation is the
least-squares slope of estimate on truth: 1 means none.

### The textbook regimes

500 replications per regime.

| regime | clusters | per cluster | coverage illume | coverage mclogit | attenuation illume | attenuation mclogit |
|----|----|----|----|----|----|----|
| large | 40 | 25 | 0.953 | 0.949 | 1.005 | 0.960 |
| small | 100 | 4 | 0.949 | 0.891 | **1.004** | **0.731** |

The literature’s prediction holds exactly. With 25 observations per
cluster the two methods agree and both cover. With 4 observations per
cluster mclogit’s coefficients shrink to 73% of their true size and its
95% intervals cover 89.1% of the time – a one-in-nine error rate where
one in twenty was advertised. illume is unaffected on both counts.

### When the data misbehaves

Both regimes above are tidy: balanced clusters, three roughly equal
outcome categories, Gaussian random effects. Real data is none of those,
and a method that only wins on tidy data has not earned much. So a
second study varies the four things most likely to break either
approximation, one at a time and then all together. 400 replications per
regime, 60 clusters of 8 observations, random-effect sd 1.0 except where
stated.

| regime | what changes |
|----|----|
| clean | nothing – the baseline |
| unbalanced | the same 480 rows over wildly uneven clusters, many singletons |
| rare | the rarest outcome category at 3.6%, about 17 of 480 observations |
| flat_re | random-effect sd 0.05, a variance component at the boundary |
| heavy_re | random effects from a two-component mixture, same sd |
| combined | all four at once |

| regime | coverage illume | coverage mclogit | attenuation illume | attenuation mclogit | converged illume | converged mclogit |
|----|----|----|----|----|----|----|
| clean | 0.950 | 0.937 | 1.011 | 0.864 | 1.000 | 1.000 |
| unbalanced | 0.943 | 0.932 | 1.006 | 0.881 | 1.000 | 1.000 |
| rare | 0.931 | 0.853 | 1.462 | 0.816 | 0.973 | 1.000 |
| flat_re | 0.948 | 0.947 | 1.018 | 1.008 | **0.633** | 1.000 |
| heavy_re | 0.936 | 0.902 | 0.982 | 0.846 | 1.000 | 1.000 |
| combined | 0.945 | 0.894 | 1.301 | 0.918 | 0.910 | 1.000 |

illume’s intervals cover better in all six regimes, and the gap widens
where theory says PQL should struggle: 0.931 against 0.853 with a sparse
category, 0.945 against 0.894 with everything combined. mclogit shrinks
fixed effects in five of six regimes and its standard errors never
reflect it. illume is 6 to 13 times faster throughout (median 0.5s
against 3.4 to 6.5s).

Non-Gaussian random effects cost illume less than expected. Both methods
assume Gaussian random effects, so `heavy_re` is misspecification for
both, and illume’s attenuation stays at 0.982 against mclogit’s 0.846.

### Where illume comes off worse

Three results in that table run against the package, and they matter
more than the ones that favour it.

**illume inflates coefficients when a category is sparse.** Attenuation
1.46 in `rare` and 1.30 in `combined` – the opposite direction from
PQL’s shrinkage, but a bias all the same. In `rare`, illume’s mean
absolute bias (0.394) is *worse* than mclogit’s (0.271). The intervals
are wide enough to cover anyway, which is what honest uncertainty buys
you, but a point estimate from a sparse-cell multinomial should not be
read at face value. With 17 observations of a category spread over 60
clusters the likelihood is nearly flat in that direction and the maximum
drifts outward.

**illume’s convergence collapses near a variance boundary.** In
`flat_re`, with a true random-effect sd of 0.05, illume converged on
63.3% of replications against mclogit’s 100%. Its 0.948 coverage there
describes only those 63.3%. Neither method recovers the variance –
illume’s median estimate is 0.099, twice the truth; mclogit’s is 0.004,
driven almost to zero – but only one of them declines to answer.

That last sentence describes illume as it was when the study ran, when a
fit counted only with a positive definite Hessian. It has since reported
the fixed effects when a variance sits at its boundary, holding that
variance at its estimate as lme4 does. On the same 400 datasets, 97.5%
of fits in `flat_re` now give fixed effects with usable standard errors,
and the 128 that did not before cover at 0.944. In the combined regime
the share rises from 91.0% to 99.0%, and the 30 fits recovered there
cover at 0.908 – a little under nominal, which is the price of holding a
variance fixed where there is least to say about it. The other findings
in this section are unchanged.

**mclogit’s RMSE is lower in four of six regimes**, by 3% (`flat_re`) to
14%. illume wins on RMSE only in the two sparse-cell regimes, and
narrowly.

The picture is not that one method dominates. The two fail differently.
PQL is reliable and biased: it always converges, it shrinks, and it does
not say so. The Laplace approximation is honest and fragile: its
intervals mean what they claim, and it refuses more often. Refusing
loudly is the design intent – see the latent budget in the *Regression
models* vignette – but it is still a cost, and a convergence-conditional
table hides it. With a near-boundary variance component or a 3%
category, read illume’s convergence status and not only its
coefficients.

### The caveat that must travel with this result

**mclogit’s RMSE is lower**: 0.336 against illume’s 0.409 in the
small-cluster regime. Shrinkage buys a variance reduction, and if you
are predicting, that trade is defensible and mclogit may well be the
better choice.

It is not defensible for inference, and the reason is specific: the
attenuation bias does not appear in the reported standard errors. A user
reading mclogit’s output in the small-cluster regime sees a coefficient
shrunk by 27% and an interval that does not know it. Nothing in the
printed result says so. That is the distinction this package is built
around – not that illume is more accurate in every sense, but that its
reported uncertainty means what it says.

## 5. Against full Bayes

`brms` fits the same models by Hamiltonian Monte Carlo and is the
closest thing to a gold standard available. Agreement is not expected to
be exact: a maximum likelihood estimate with a Wald interval and a
posterior mean with a credible interval answer different questions, and
the residual gap is prior shrinkage.

3 datasets, 8 coefficients, maximum R-hat 1.008:

- agreement with brms: maximum 0.202 standard errors, mean 0.097
- SE ratio illume/brms: 0.984
- illume median 0.60s; brms sampling only, median 46.2s

A mean gap of a tenth of a standard error is agreement. The SE ratio
near 1 says the Laplace approximation is not understating uncertainty
relative to full posterior sampling on these problems. The 77-fold speed
difference is the price of that sampling, and it is the reason a
frequentist engine is worth having.

## What these studies do not show

Worth stating plainly, because the tables above can look more complete
than they are.

- **Coverage was measured under correct specification.** The
  data-generating model and the fitted model match everywhere except the
  `heavy_re` regime above, which misspecifies the random-effect
  distribution on purpose. Real data is generated by neither, and no
  simulation study of this kind speaks to that.
- **Convergence-conditional results are optimistic** wherever the
  convergence column is below 1.
- **Nothing here validates the exploration half of the package.**
  Clustering, dimension reduction and anomaly detection have their own
  measured claims, reported in the *Profiling* and *Anomaly detection*
  vignettes.
- **Timings are from one machine** (R 4.4.3, Windows, ten cores) and are
  order-of-magnitude statements, not benchmarks.
- **The comparisons are with the versions current at the time of the
  run.** `lme4`, `glmmTMB`, `mclogit` and `brms` all move.

## Reproducing this

The scripts are in `studies/scripts/` in the source repository. Each
takes the number of replications and an output directory, and those that
run in parallel take a number of cores between the two:

``` r
# from studies/scripts/, with the replication counts behind the tables above
Rscript benchmark_freq.R     160    ../runs/0.0.7.9000/bench
Rscript coverage_study.R    2000 10 ../runs/0.0.7.9000/coverage
Rscript power_study.R       2000 10 ../runs/0.0.7.9000/power
Rscript mclogit_compare.R    500 10 ../runs/0.0.7.9000/mclogit
Rscript messy_compare.R      400 10 ../runs/0.0.7.9000/messy
Rscript brms_compare.R         3    ../runs/0.0.7.9000/brms
```

Raw per-replication output is not committed – it is several megabytes
per version and would sit in the history permanently.
`studies/findings/` holds the permanent record, one section per package
version, generated by `summarise_run.R` rather than written by hand.

## See also

- *Regression models* – what
  [`ilm_model()`](https://huttoncp.github.io/illume/reference/ilm_model.md)
  fits, and the latent budget
- *Effect size and power* – the power machinery these numbers validate
- *Missing data* – the imputation comparison, which has its own study
