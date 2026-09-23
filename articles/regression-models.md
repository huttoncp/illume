# Regression models

[`ilm_model()`](https://huttoncp.github.io/illume/reference/ilm_model.md)
fits every model in this package. The formula is lme4’s, the response
distribution is an argument, and the object that comes back is the same
whatever you fitted – so
[`summary()`](https://rdrr.io/r/base/summary.html),
[`ilm_anova()`](https://huttoncp.github.io/illume/reference/ilm_anova.md),
[`predict()`](https://rdrr.io/r/stats/predict.html) and the diagnostics
do not change as the model gets more complicated.

This vignette works through it: a first fit and how to read it, then the
families one at a time, then the structures that go on top of them –
random slopes, smooths, correlation over time and space, non-constant
spread, censoring. The *Introduction* vignette is the shorter
orientation; the *Workflow* vignette is the map of where this stage
sits.

``` r

library(illume)
```

## Why a multinomial mixed model is fitted this way

Fitting a multinomial mixed model has meant choosing which problem to
accept. Penalised quasi-likelihood (`mclogit`) is fast, and in sparse
data it shrinks fixed effects toward zero – 27% attenuation with 100
clusters of 4 – while reporting standard errors that do not know it, so
coverage falls to 0.891 against a nominal 0.95. MCMC (`brms`) is exact
up to Monte Carlo error and takes about 46 seconds of sampling per fit,
which is a long time when you are still deciding what the model is.

illume takes a third route: the Laplace approximation via RTMB, with the
random effects carried as a matrix-normal block. Against brms on the
same data, estimates agree to a mean of 0.10 standard errors with an SE
ratio of 0.984; against mclogit it holds coverage at 0.949 where PQL
reaches 0.891. Median fit time is 0.60 seconds, roughly 75 times faster
than brms’ sampling alone.

One caveat travels with that: mclogit’s RMSE is *lower* (0.336 against
0.409), because shrinkage buys a variance reduction that more than pays
for the bias. For prediction that is a defensible trade. For a
coefficient you intend to interpret it is not, and coverage is the
measure that says so.

On data that misbehave the picture is less tidy. With an outcome
category under 4% of the data, illume *inflates* coefficients rather
than shrinking them, and its mean absolute bias there is worse than
mclogit’s; and in the study, with a random-effect SD near zero, only
63.3% of its fits counted as converged against mclogit’s 100%. (That
study required a positive definite Hessian. illume has since reported
the fixed effects when a variance sits at its boundary, and on the same
datasets 97.5% of fits now give them.) The whole account, including
where illume comes off worse, is in
[`vignette("benchmarking")`](https://huttoncp.github.io/illume/articles/benchmarking.md).

## A first model

We will simulate data where an outcome with three categories depends on
a continuous predictor `x1` and a grouping factor `grp`, measured
repeatedly within subjects.

``` r

n_subj <- 40; per <- 20; N <- n_subj * per
dd <- data.frame(
  subj = factor(rep(seq_len(n_subj), each = per)),
  x1   = rnorm(N),
  grp  = factor(sample(c("a", "b", "c"), N, TRUE))
)

# true structure: an effect of x1 and grp, plus subject-level variation
X  <- model.matrix(~ x1 + grp, dd)
bt <- matrix(c(0.1, -0.2,  0.9, 0.5,  -0.6, 0.4,  0.2, -0.5), 4, 2, byrow = TRUE)
b  <- matrix(rnorm(n_subj * 2), n_subj, 2) %*% diag(c(0.8, 0.6))
eta <- X %*% bt + b[as.integer(dd$subj), ]
P <- exp(eta %*% t(contr.sum(3))); P <- P / rowSums(P)
dd$y <- factor(c("low", "mid", "high")[apply(P, 1, function(p) sample.int(3, 1, prob = p))],
               levels = c("low", "mid", "high"))

table(dd$y)
#> 
#>  low  mid high 
#>  282  206  312
```

The model is written the way you would write it in `lme4`:

``` r

fit <- ilm_model(y ~ x1 + grp + (1 | subj), data = dd,
                 family = "multinomial", verbose = FALSE)
#> Warning: the 'findbars' function has moved to the reformulas package. Please
#> update your imports, or ask an upstream package maintainer to do so.
#> Warning: the 'nobars' function has moved to the reformulas package. Please
#> update your imports, or ask an upstream package maintainer to do so.
fit
#> ilm_model fit: 3 categories (multinomial), 800 obs, 8 fixed + 3 covariance parameters
#>   logLik -549.45 | AIC 1267.9
```

## Reading the output

``` r

summary(fit)
#> Generalized linear mixed model fit by maximum likelihood (Laplace approximation)
#>  Family: multinomial (3 categories: low, mid, high)
#> Formula: y ~ x1 + grp + (1 | subj) 
#> 
#>      AIC      BIC   logLik deviance df.resid
#>   1267.9   1319.5   -623.0   1245.9      789
#> 
#> Random effects:
#>  subj  [us]  40 levels
#>        SD   low
#> low 0.736      
#> mid 0.598 0.173
#> 
#> Number of obs: 800; groups: subj 40
#> 
#> Fixed effects:
#>                 Estimate Std. Error z value Pr(>|z|)    
#> low:(Intercept)  0.29964    0.16332   1.835 0.066554 .  
#> low:x1           1.03963    0.09361  11.106  < 2e-16 ***
#> low:grpb        -0.63624    0.15787  -4.030 5.58e-05 ***
#> low:grpc         0.07219    0.15555   0.464 0.642602    
#> mid:(Intercept) -0.11321    0.15475  -0.732 0.464424    
#> mid:x1           0.48832    0.08796   5.552 2.83e-08 ***
#> mid:grpb         0.50936    0.14912   3.416 0.000636 ***
#> mid:grpc        -0.32462    0.16883  -1.923 0.054505 .  
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> ---
#> Category contrasts are SUM-TO-ZERO: each coefficient is that category's
#> deviation from the across-category average, NOT a contrast against a baseline.
#> 
#> Model checks: all passed.
```

Two things in that output deserve attention.

**The coefficients are sum-to-zero contrasts.** A coefficient such as
`low:x1` is the effect of `x1` on the `low` category *relative to the
average across all categories*, not relative to a baseline category.
This differs from
[`nnet::multinom()`](https://rdrr.io/pkg/nnet/man/multinom.html) and
`brms`, which compare each category against a reference. The summary
prints a reminder, because it is easy to misread if you are used to the
other convention.

**The assumption checks are printed with the coefficients.** That is
unusual for a [`summary()`](https://rdrr.io/r/base/summary.html) method
and it is deliberate, for a reason worth explaining.

## Why the diagnostics are so prominent

This class of model fails *quietly*. A fit can converge, produce
entirely reasonable-looking coefficients with small standard errors and
significance stars, and yet the covariance matrix behind those standard
errors can be unusable – meaning the uncertainty is meaningless even
though nothing in the output looks wrong.

Every fit therefore carries a table of checks:

``` r

head(fit$checks[, c("check", "status", "detail")], 8)
#>                 check status
#> 1     category_counts     OK
#> 2        weights_type     OK
#> 3     re_levels[subj]     OK
#> 4 obs_per_level[subj]     OK
#> 5       latent_budget     OK
#> 6           optimizer     OK
#> 7            gradient     OK
#> 8             hessian     OK
#>                                                                              detail
#> 1                                    smallest category n = 206 of 800 total (J = 3)
#> 2                                     none (all 1) weights; range 1 to 1, total 800
#> 3             40 levels for 3 covariance parameters (13.3 per parameter); us, C = 2
#> 4                                                       20.0 observations per level
#> 5 10.00 observations per latent value (800 observations, 80 latent values: subj 80)
#> 6                                          nlminb code 0 (relative convergence (4))
#> 7                                                         max |gradient| = 1.62e-04
#> 8              positive definite: TRUE; non-finite or non-positive variances: FALSE
```

Each check has one of five statuses. `OK` and `FAIL` are
self-explanatory. `WARN` means something is worth a look. `BOUNDARY`
means a variance sits at zero, where the usual tests do not apply. And
`INCONCLUSIVE` means there is not enough information to judge – which is
reported honestly rather than being passed off as a clean bill of
health.

Anything that is not `OK` also carries a reason and a suggested remedy:

``` r

bad <- fit$checks[fit$checks$status != "OK", ]
if (nrow(bad)) bad[1, c("check", "cause", "suggestion")] else "all checks passed"
#> [1] "all checks passed"
```

### The most useful check

The one to look at first is `latent_budget`. A mixed model introduces
one unobserved value per random effect per group. The method used to fit
the model (the Laplace approximation) is accurate when each of those
values is informed by a reasonable amount of data, and degrades when it
is not.

Importantly this is a *global* constraint. A single term can look
perfectly healthy on its own and the model can still be over-stretched,
because what matters is the total. Per-term checks cannot see that; this
one can.

## Testing fixed effects

[`ilm_anova()`](https://huttoncp.github.io/illume/reference/ilm_anova.md)
gives an analysis of deviance. Each term is tested **jointly across all
category dimensions**, so a predictor with one column is tested on two
degrees of freedom when there are three categories.

``` r

ilm_anova(fit, type = 3)
#> Analysis of Deviance Table (Type III Wald chi-square tests)
#> Response: y   (3 categories, 2 contrast dimensions)
#> Each term is tested jointly across all category dimensions: Df = (columns) x C
#>     Df   Chisq Pr(>Chisq)    
#> x1   2 166.589  < 2.2e-16 ***
#> grp  4  38.953  7.125e-08 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

Type III tests each term with all others present; Type II tests each
term after the terms that do not contain it. They agree when there are
no interactions. Wald tests are the default because they need no
refitting; a likelihood-ratio version is available with `test = "LRT"`
and is more reliable for very strong effects.

For small samples neither is ideal, because both rely on large-sample
approximations. What to do about it depends on the family.

For a **gaussian** mixed model,
[`ilm_denom_df()`](https://huttoncp.github.io/illume/reference/ilm_denom_df.md)
supplies a finite denominator degrees of freedom – Satterthwaite by
default, Kenward-Roger where it applies – which turns the chi-square
into an F:

``` r

ilm_denom_df(fit_gaussian, L)                       # for one contrast
```

Satterthwaite works for every structure this package fits. Kenward-Roger
is only defined where the marginal covariance is linear in the variance
components – random intercepts and slopes with a single residual
variance – and says so rather than returning a number it cannot stand
behind.

For **any other family**, including the multinomial here, neither
approximation is derived and neither is offered.
[`ilm_pb_lrt()`](https://huttoncp.github.io/illume/reference/ilm_pb_lrt.md)
replaces the assumed reference distribution with a simulated one
instead, which needs no such derivation:

``` r

ilm_pb_lrt(fit, "x1", B = 200, ncores = 4)
```

## Maximum likelihood, and when to switch to REML

Variance components estimated by maximum likelihood are biased downward,
because ML does not account for the degrees of freedom spent on the
fixed effects. With many clusters that hardly matters; with few it
carries through into standard errors and degrees of freedom.

`reml = TRUE` corrects it:

``` r

ilm_model(y ~ x1 + grp + (1 | id), data = d, family = "gaussian",
          reml = TRUE)
```

The default is maximum likelihood, and the reason is the order you work
in. A restricted likelihood is the likelihood of contrasts orthogonal to
the design matrix, so changing the fixed effects changes which data it
is the likelihood of. Two REML fits with different fixed effects are not
comparable, and their difference is not a likelihood ratio.

So: **settle the mean structure under the default, then refit with
`reml = TRUE` for the estimates you report.** Once set,
`ilm_anova(test = "LRT")` and
[`ilm_pb_lrt()`](https://huttoncp.github.io/illume/reference/ilm_pb_lrt.md)
refuse rather than quietly comparing things that are not comparable, and
so does
[`ilm_robust()`](https://huttoncp.github.io/illume/reference/ilm_robust.md)
– a sandwich needs per-observation scores, and integrating the
coefficients out leaves none.

[`ilm_dag_model()`](https://huttoncp.github.io/illume/reference/ilm_dag_model.md)
and
[`ilm_aov_ez()`](https://huttoncp.github.io/illume/reference/ilm_aov_ez.md)
default to REML, because in both the fixed effects were fixed before any
data were seen: by the graph in one case and by the design in the other.

## Predictions

By default [`predict()`](https://rdrr.io/r/stats/predict.html) returns a
probability for every observation and category, as
`nnet::multinom(type = "probs")` does.

``` r

head(round(predict(fit), 3))
#>     low   mid  high
#> 1 0.759 0.196 0.045
#> 2 0.639 0.220 0.141
#> 3 0.398 0.279 0.323
#> 4 0.178 0.416 0.406
#> 5 0.786 0.185 0.030
#> 6 0.742 0.241 0.017
```

There is a choice to make here. With `marginal = FALSE` (the default)
the random effects are set to zero, giving probabilities for a *typical*
subject. With `marginal = TRUE` the prediction is averaged over the
distribution of subjects, giving probabilities for the *population*.

``` r

grid <- data.frame(x1 = 0, grp = factor(c("a", "b", "c"), levels = levels(dd$grp)))
round(predict(fit, newdata = grid, marginal = FALSE), 3)
#>     low   mid  high
#> 1 0.439 0.291 0.270
#> 2 0.227 0.473 0.300
#> 3 0.458 0.204 0.338
round(predict(fit, newdata = grid, marginal = TRUE, ndraw = 200), 3)
#>     low   mid  high
#> 1 0.417 0.282 0.301
#> 2 0.235 0.442 0.323
#> 3 0.439 0.206 0.356
```

These differ because averaging and the softmax transform do not commute:
the average of the transformed values is not the transform of the
average. The population-averaged probabilities are pulled toward being
more even across categories. Which you want depends on the question you
are asking.

Standard errors and intervals come from simulation, and the intervals
are percentile intervals, so they always lie between 0 and 1.

``` r

pr <- predict(fit, newdata = grid, se.fit = TRUE, nsim = 100)
round(data.frame(fit = pr$fit[, 1], lower = pr$lower[, 1], upper = pr$upper[, 1]), 3)
#>     fit lower upper
#> 1 0.439 0.325 0.538
#> 2 0.227 0.167 0.327
#> 3 0.458 0.348 0.570
```

## Random slopes, smooths and correlated time

Beyond random intercepts, the formula accepts random slopes and
penalised smooths:

``` r

ilm_model(y ~ x1 + time + (1 + time | subj), data = dd, family = "multinomial")
ilm_model(y ~ x1 + s(xs, k = 10) + (1 | subj), data = dd, family = "multinomial")
ilm_model(y ~ x1 + t2(lon, lat) + (1 | site), data = dd, family = "multinomial")
```

Two rules for smooths: use `t2()` rather than `te()` for tensor
products, and write them unqualified – `s(x)`, not `mgcv::s(x)` –
because mgcv identifies smooths by name.

## Simple models, and exact inference

Nothing requires a random effect.
[`ilm_model()`](https://huttoncp.github.io/illume/reference/ilm_model.md)
fits an ordinary regression, and when the family is gaussian and nothing
is integrated out it reports **exact** inference rather than
large-sample approximations:

``` r

set.seed(1)
sd2 <- data.frame(x = rnorm(200), z = factor(sample(c("a", "b", "c"), 200, TRUE)))
sd2$v <- 2 + 1.5 * sd2$x - 0.8 * (sd2$z == "b") + rnorm(200, 0, 1.2)

f_lm <- ilm_model(v ~ x + z, data = sd2, family = "gaussian", verbose = FALSE)
round(ilm_coef_table(f_lm), 5)
#>             Estimate Std. Error  t value Pr(>|t|)
#> (Intercept)  2.06024    0.14527 14.18256  0.00000
#> x            1.62374    0.09605 16.90496  0.00000
#> zb          -0.89137    0.21705 -4.10674  0.00006
#> zc          -0.12637    0.21410 -0.59024  0.55571
round(summary(lm(v ~ x + z, data = sd2))$coefficients, 5)
#>             Estimate Std. Error  t value Pr(>|t|)
#> (Intercept)  2.06024    0.14527 14.18257  0.00000
#> x            1.62374    0.09605 16.90495  0.00000
#> zb          -0.89137    0.21705 -4.10675  0.00006
#> zc          -0.12637    0.21410 -0.59026  0.55569
```

The two agree to numerical precision, including the residual standard
deviation and the degrees of freedom.
[`ilm_anova()`](https://huttoncp.github.io/illume/reference/ilm_anova.md)
likewise switches from chi-square to F:

``` r

ilm_anova(f_lm, type = 3)
#> Analysis of Deviance Table (Type III F tests)
#> Response: v   (family: gaussian)
#>   Df  F value    Pr(>F)    
#> x  1 285.7775 < 2.2e-16 ***
#> z  2   9.4595 0.0001197 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

This matters because the chi-square test treats the residual variance as
known, which it is not, making it anti-conservative in small samples.
Where an exact reference exists, illume uses it.

It applies only where it is genuinely available. A Poisson or binomial
model without random effects is still a GLM, with no exact small-sample
analogue, so it keeps z and chi-square. So does any model with a random
effect.

## Other families

The same formula syntax, methods and diagnostics apply to every family.
What changes is the response, the link, and whether there is a
dispersion parameter:

``` r

ilm_model(count ~ x + (1 | site), data = dd, family = "poisson")
ilm_model(count ~ x + (1 | site), data = dd, family = "nbinom")   # overdispersed
ilm_model(hit   ~ x + (1 | site), data = dd, family = "binomial")
ilm_model(score ~ x + (1 | site), data = dd, family = "gaussian")
```

[`summary()`](https://rdrr.io/r/base/summary.html) reports a dispersion
parameter where the family has one – the residual standard deviation for
gaussian, the overdispersion parameter for negative binomial – and omits
that section otherwise. Coefficients carry no category prefix outside
the multinomial, and the sum-to-zero note disappears, because there is
only one linear predictor.

A two-category outcome can be stored however your data happen to hold it
– `0`/`1`, a two-level factor, `TRUE`/`FALSE`, or a character column
read straight from a csv. All four give the same fit. For a factor or
character response the *second* level is the one being modelled, exactly
as in [`glm()`](https://rdrr.io/r/stats/glm.html), so
`factor(c("no", "yes"))` models the probability of `"yes"`. If you want
the other direction, relevel the factor. You do not have to remember
which way round it is:
[`summary()`](https://rdrr.io/r/base/summary.html) states the category
it is modelling, and the reference it is measured against.

Two categories belong to `"binomial"` and three or more to
`"multinomial"`; each family says so if it is given the other one’s
data.

If your counts vary more than Poisson allows, `"nbinom"` adds one
parameter `k`, giving variance `mu + mu^2/k`. Smaller `k` means more
overdispersion; large `k` approaches the Poisson.

Zero-inflation, hurdle models and Tweedie are deliberately out of scope.
`glmmTMB` covers those well, and a weaker reimplementation would help
nobody.

## When there are many categories

The number of covariance parameters grows quickly with the number of
categories. With an unstructured covariance and `J` categories you need
`(J-1)J/2` parameters *per grouping factor* – 3 when `J = 3`, but 45
when `J = 10`. Many designs cannot support that.

The remedy is a reduced-rank structure, which describes the covariance
with fewer dimensions:

``` r

ilm_model(y ~ x1 + (1 | subj), data = dd, family = "multinomial",
          re_struct = list(subj = list(type = "rr", rank = 3)))
```

This saves parameters *and* latent values, and the latter is usually the
bigger gain. When a structure is too rich for the data, the checks say
so and name a specific rank to try.

## Checking the fitted model

[`ilm_appraise()`](https://huttoncp.github.io/illume/reference/ilm_appraise.md)
draws diagnostic panels built for a categorical outcome rather than
borrowed from tools designed for continuous responses:

``` r

ilm_appraise(fit)
```

Two further checks are worth knowing about, because broad residual
summaries have very little power to detect the problems they address.
[`ilm_check_omitted()`](https://huttoncp.github.io/illume/reference/ilm_check_omitted.md)
tests variables you did *not* include, and
[`ilm_check_ar()`](https://huttoncp.github.io/illume/reference/ilm_check_ar.md)
tests for leftover correlation over time:

``` r

ilm_check_omitted(fit, dd)
ilm_check_ar(fit, time = dd$time, group = dd$subj)
```

And
[`ilm_consistency()`](https://huttoncp.github.io/illume/reference/ilm_consistency.md)
asks the most direct question available: can the model recover itself
from data it generated? If not, the fitting method is not reliable for
this model and this amount of data.

``` r

ilm_consistency(fit, B = 50, ncores = 4)
```

## Structure over time

Repeated measurements on the same unit are rarely independent, and there
are two quite different reasons why – which matter because they need
opposite remedies.

[`ilm_plot_acf()`](https://huttoncp.github.io/illume/reference/ilm_plot_acf.md)
draws the residual autocorrelation and partial autocorrelation by lag,
within group, with a verdict on each lag:

``` r

ilm_plot_acf(fit, time = dd$time, group = dd$subj, maxlag = 14)
```

Three things make this not
[`stats::acf()`](https://rdrr.io/r/stats/acf.html):

- Pairs are matched on exact time differences **within group**. A panel
  dataset stacked long is many short series, and pairing across the
  joins manufactures correlation that is not there.
- The band is not `2/sqrt(n)`. Residuals within a group sum to roughly
  zero, so pairs of them correlate *negatively* even when the model is
  exactly right – measured, around -0.07 at every lag on a 60-unit
  panel. A band centred on zero calls correct models wrong. The band
  here comes from simulating from the fit and refitting.
- The band is drawn at the critical value the p-value is computed
  against, so a point outside the band is exactly a lag the check flags.
  The picture and the table cannot disagree.

### When the times are not evenly spaced

[`ilm_ar1()`](https://huttoncp.github.io/illume/reference/ilm_ar1.md)
needs a common step, because “one step” has to mean something. When
visits fall where they fall, use
[`ilm_car1()`](https://huttoncp.github.io/illume/reference/ilm_car1.md):
the correlation between two observations is `rho` raised to the gap
between them, which is `nlme`’s `corCAR1()` and reduces to AR(1) when
every gap is one.

``` r

ilm_model(y ~ x + (1 | subj), data = dd, family = "gaussian",
          ar = ilm_car1(dd$day, dd$subj))
```

Both constructors report their **latent budget** before you fit, because
that is the constraint that bites. These are latent-variable structures:
one latent value exists per distinct unit-and-time, so genuinely
continuous times give one latent per observation. For a gaussian
response that is fine — the model is linear-Gaussian, the Laplace
approximation is exact, and coverage measured over 400 replicates at
exactly one observation per latent came back at 0.955, 0.937 and 0.949
against a nominal 0.95. For a categorical response it is not, and the
fit says so. Rounding `time` to a coarser grid, so observations share a
latent, is the usual fix.

### Checking a correlation structure over irregular time

[`ilm_check_ar()`](https://huttoncp.github.io/illume/reference/ilm_check_ar.md)
matches pairs at exact lags, which stops working once the times are
irregular: on a sixty-unit panel with six observations each drawn from
thirty possible times, it found 53 pairs at lag 1, five at lag 2 and
none beyond.
[`ilm_variogram()`](https://huttoncp.github.io/illume/reference/ilm_variogram.md)
bins every within-group pair by how far apart it is instead, and uses
all 900.

``` r

ilm_variogram(fit, time = dd$day, group = dd$subj)
```

It distinguishes three shapes, because they need three different
remedies:

| what you see | what it is | what to do |
|----|----|----|
| decays from the shortest separation | autoregression | [`ilm_car1()`](https://huttoncp.github.io/illume/reference/ilm_car1.md) |
| the same level at every separation | a group effect | add or widen a random effect |
| alternates in sign | a cycle | [`ilm_fourier()`](https://huttoncp.github.io/illume/reference/ilm_fourier.md), [`ilm_cyclic()`](https://huttoncp.github.io/illume/reference/ilm_cyclic.md) |

The envelope again comes from refitting simulated data, and again is not
centred on zero: residuals within a group sum to roughly zero, so pairs
of them correlate negatively even under a correct model, and more so the
further apart they are.

### Correlation in space

The same function takes coordinates instead of a time, because
separation in space and separation in time are the same statistic over a
different distance:

``` r

ilm_variogram(fit, coords = dd[, c("easting", "northing")])
```

The remedy it names is different, and deliberately so. illume fits no
spatial covariance — no Matern field, no exponential one — so the advice
does not pretend otherwise. What it points at is putting the structure
in the **mean**: a tensor-product smooth of the coordinates, `t2(x, y)`,
absorbs smooth spatial variation, and a random effect for a spatial
grouping absorbs the coarse kind.

Measured on 300 points carrying an exponential field, residual
correlation ran 0.242 at the shortest distance down through 0.071 to
−0.072 at the longest. Adding `t2(sx, sy)` cut the first of those to
−0.033, dropped AIC from 1165.7 to 1069.4, and narrowed the standard
error on the covariate of interest from 0.100 to 0.083.

### A departure too small to act on

Each bin rests on thousands of pairs, so the envelope narrows until any
imperfection clears it. In the example above the residual correlation
after the smooth was about 0.03 against an envelope 0.02 wide — a
seventh of what it had been, and nothing anyone would act on — and four
of five bins were still flagged.

`min_effect` is the floor, 0.1 by default: a bin outside the envelope by
less than that is reported as `OK`. It is the same stance the gaussian
index takes against normality tests. With the floor in place the field
is still caught before the smooth is added and the diagnosis is clean
after it; with `min_effect = 0` both are flagged and the verdict stops
distinguishing them.

`type = "semivariance"` reports the same information in the form
[`nlme::Variogram()`](https://rdrr.io/pkg/nlme/man/Variogram.html) uses,
since for standardised residuals the semivariance at separation `d` is
`1 - r(d)`.

### A cycle is not autocorrelation

Both fill the autocorrelation plot with flagged lags. Only one of them
is fixed by an AR term.

The partial autocorrelation is what separates them, and so is the
*shape* of the ordinary one. Autoregression decays through zero and
stays there. A cycle comes back up to a positive peak at its period – a
spike at lag 12 in monthly data, lag 7 in daily data.
[`ilm_check_ar()`](https://huttoncp.github.io/illume/reference/ilm_check_ar.md)
reads that difference and names the period rather than reaching for
AR(1), which cannot represent a cycle at all.

The remedy is ordinary fixed effects; there is no new likelihood
involved, only the right columns:

``` r

head(round(ilm_fourier(1:12, period = 12, K = 2), 3), 3)
#>       sin1  cos1  sin2 cos2
#> [1,] 0.500 0.866 0.866  0.5
#> [2,] 0.866 0.500 0.866 -0.5
#> [3,] 1.000 0.000 0.000 -1.0
```

``` r

# a smooth annual cycle in monthly data
ilm_model(y ~ x + ilm_fourier(month, 12, K = 2) + (1 | subj), data = dd,
          family = "gaussian")

# local flexibility instead, for a long cycle such as day-of-year
ilm_model(y ~ x + ilm_cyclic(doy, 365.25, df = 8) + (1 | subj), data = dd,
          family = "gaussian")
```

[`ilm_cyclic()`](https://huttoncp.github.io/illume/reference/ilm_cyclic.md)
is a cubic spline that closes on itself: the value and its first two
derivatives match across the wrap, so December runs into January without
a step. An ordinary `bs()` or `ns()` on a phase variable does not know
the two ends of the cycle are the same point, and leaves a jump there.

Which to use is less consequential than it looks. Measured at matched
degrees of freedom on a twelve-phase cycle, Fourier terms and the cyclic
spline came within a few AIC of each other across four shapes. Where the
shape genuinely *jumped*, a plain `factor(month)` beat both by more than
a hundred, because no smooth basis buys a step cheaply. The spline pulls
ahead only on long cycles, and only above about six degrees of freedom,
by at most 23 AIC. Pick on the shape you expect, not on the basis.

## When the spread is not constant

[`ilm_check_variance()`](https://huttoncp.github.io/illume/reference/ilm_check_variance.md)
reports two patterns, and both have the same remedy: model the logarithm
of the dispersion instead of treating it as one number.

``` r

# a separate spread per group
ilm_model(y ~ x + (1 | id), data = dd, family = "gaussian",
          dispformula = ~ site)

# spread that grows with the fitted mean; `mu` is the one reserved name
ilm_model(y ~ x + (1 | id), data = dd, family = "gaussian",
          dispformula = ~ mu)
```

These are `nlme`’s `varIdent()` and `varPower()`, and `glmmTMB`’s
`dispformula`, written one way. On a three-group design illume,
`glmmTMB` and `gls` returned the same log-likelihood to five decimals
and the same standard errors; on a power-of-the-mean design illume and
`gls` agreed on the power to two decimals.

It works for any family that has a dispersion to model — the negative
binomial’s overdispersion and an accelerated failure time model’s scale
as well as a gaussian standard deviation. A Poisson or binomial response
has no dispersion parameter at all, its spread being fixed by its mean,
and asking for one says so.

Two things follow. Exact t and F inference is switched off, because it
rests on a constant variance and the variance is now itself estimated —
`gls()` reports t here by convention, but this package reserves t for
where it is exact. And the fit is warm started from the ordinary one,
which matters most for `~ mu`: from a cold start the fitted mean is zero
for every row, `log|mu|` is the log of the numerical floor, and the
power multiplying it can go anywhere. Measured on data with a true power
of 1.0, a cold start converged falsely with the power at -115818 and the
slope collapsed to zero.

## Floors, ceilings and detection limits

A reading of `0` from an assay whose detection limit is 0 does not mean
zero. It means “somewhere at or below the limit”, and fitting it as a
zero pulls the mean down and the variance in. Same for a scale that tops
out, a score capped at 100, a measurement beyond an instrument’s range.

[`ilm_describe()`](https://huttoncp.github.io/illume/reference/ilm_describe.md)
names the pile-up before you model it:

``` r

set.seed(1)
y <- 0.5 + rnorm(1500, 0, 1.4)
d0 <- data.frame(unbounded = y, capped = pmin(y, 1.2))
ilm_describe(d0)[, c("variable", "gauss", "gauss_note")]
#>    variable gauss
#> 1 unbounded     1
#> 2    capped     0
#>                                                                                 gauss_note
#> 1                                                                                         
#> 2 31% of values sit exactly at the maximum (1.2): a ceiling, see ilm_censor(); left-skewed
```

and
[`ilm_censor()`](https://huttoncp.github.io/illume/reference/ilm_censor.md)
is the remedy. A censored observation contributes the probability of the
interval it is known to lie in rather than a density at a value it was
never observed at:

``` r

ilm_model(y ~ x + (1 | subj), data = dd, family = "gaussian",
          censor = ilm_censor(dd$y, upper = 1.2))
```

With no random effects this is the Tobit model, and it agrees with
[`survival::survreg()`](https://rdrr.io/pkg/survival/man/survreg.html)
to five decimal places in the coefficients and to 1e-6 in the
log-likelihood. What illume adds is that the random effects still work.
Measured on a 60-unit panel censored at 18%, ignoring the ceiling
attenuated the slope from 1.02 to 0.79 against a truth of 1.0, and the
between-unit standard deviation from 0.78 to 0.65 against a truth of
0.8.

Two things follow automatically. Exact t and F inference is switched
off, because a censored likelihood mixes densities with tail
probabilities and is only asymptotically normal — `survreg()` reports
Wald for the same reason. And the quantile residual of a censored row is
drawn across the probability of its interval, the same randomisation
Dunn and Smyth use for a discrete response; without it every censored
row lands on one quantile and the QQ plot shows a spike at the limit
that has nothing to do with fit.

## Time to an event

The same interval-probability mechanism gives parametric survival
models. Three families are available, and all three are the same model
with a different error distribution:

``` math
\log T = X\beta + \text{scale} \times W
```

`"weibull"`, `"lognormal"` and `"loglogistic"`, for `W` extreme-value,
normal and logistic respectively.

``` r

ilm_model(time ~ treatment + age + (1 | centre), data = dd, family = "weibull",
          censor = ilm_surv(dd$time, dd$event))
```

[`ilm_surv()`](https://huttoncp.github.io/illume/reference/ilm_surv.md)
takes `event` in the convention
[`survival::Surv()`](https://rdrr.io/pkg/survival/man/Surv.html) uses —
`1` when the event happened, `0` when the subject was still event-free
at the end of follow-up. It exists because that is the opposite of the
code stored internally, and getting it backwards fits the model to the
wrong subjects without looking wrong.

A coefficient is a **log time ratio**: `+0.5` means survival times are
`exp(0.5)` times longer. That reading holds whatever the baseline hazard
does, which is the appeal of an accelerated failure time model over a
proportional hazards one.

All three agree with
[`survival::survreg()`](https://rdrr.io/pkg/survival/man/survreg.html)
to about seven decimal places in the coefficients and ten in the
log-likelihood, with the same standard errors. What illume adds is the
random effect. On an 80-centre study with 44% censoring, a frailty
Weibull recovered a slope of 0.746, a scale of 0.597 and a
between-centre standard deviation of 0.465 against truths of 0.7, 0.6
and 0.5; fitting the same data without the centre term pushed the
between-centre spread into the scale, which came back at 0.729.

### Survival curves, and whether the family was right

``` r

ilm_plot_survival(fit, dd$time, dd$event)
```

The fitted curve is drawn over the Kaplan-Meier estimate, which assumes
nothing about the shape of the baseline. Where the two part company, the
parametric assumption is what is wrong, and the footnote says to try
another family and compare by AIC.

The grey band is not a confidence interval for the curve. It is the
spread of the Kaplan-Meier across datasets simulated from the fit and
refitted, which is the reference every check in this package uses. It
matters here more than most: on a correctly specified 500-subject model
the largest vertical gap between the two curves was 0.12, which looks
alarming and is not — the Kaplan-Meier of 500 subjects wobbles by that
much in the tail, and the verdict came back `OK` with none of the curve
outside the band. The same measure on a Weibull fitted to log-logistic
data gave a *smaller* gap, 0.09, and a verdict of `FAIL` with 26% of the
curve outside. A threshold on the gap would have got both backwards.

[`ilm_survival()`](https://huttoncp.github.io/illume/reference/ilm_survival.md)
returns the curve itself, with intervals, for plotting elsewhere.

Residuals work as they do everywhere else. The quantile residual of a
censored subject is drawn across the probability of the interval they
are known to lie in, so it stays uniform under a correct model —
measured across all three families at 40% censoring, Kolmogorov-Smirnov
p-values ran from 0.37 to 0.90.

## Interactions

Interactions are ordinary terms: `x * z` and `g * x` fit, their columns
are grouped into one multi-degree-of-freedom row in
[`ilm_anova()`](https://huttoncp.github.io/illume/reference/ilm_anova.md),
and
[`ilm_check_collinearity()`](https://huttoncp.github.io/illume/reference/ilm_check_collinearity.md)
reports them as single terms.

Two things are worth knowing.

**Type III depends on how factors are coded.** With R’s default
treatment contrasts, the main-effect row for a variable involved in an
interaction tests its effect *at the reference level* of the factor it
interacts with, not averaged over it. On a 600-row model the `x` row
came back at a chi-square of 167 under `contr.treatment` and 1086 under
`contr.sum` – same model, same data. Only the second is what a Type III
main effect is supposed to mean, so
[`ilm_anova()`](https://huttoncp.github.io/illume/reference/ilm_anova.md)
warns and names the fix:

``` r

ilm_model(y ~ g * x + (1 | subj), data = dd, family = "gaussian",
          contrasts = list(g = "contr.sum"))
```

`type = 2` sidesteps the question entirely, at the cost of testing a
different (and for main effects, usually more useful) hypothesis.

**An effect plot cannot show one curve.** When a term is in an
interaction it has no single effect, so `ilm_plot_model(fit, "effect")`
draws one curve per level of the interacting factor – or, for a numeric
partner, one per quartile – and says so underneath rather than quietly
showing you the modal level:

``` r

ilm_plot_model(fit, "effect", term = "x")
```

## From association to effect

Everything above estimates associations. Three designs turn one into a
causal claim, each naming an assumption and letting it be partly
checked:

``` r

# the regression a causal graph implies, with the graph itself tested
g <- ilm_dag("dag { x [exposure] ; y [outcome] ; z -> x -> y ; z -> y }")
ilm_dag_model(g, d)

# a treatment that switches on for some units at some time
ilm_did(panel, "y", "unit", "time", treated = "treated", post = "post")

# a treatment assigned by crossing a threshold
ilm_rdd(d, "score", "running", cutoff = 0)
```

The rule that makes these safe to automate: the **design** fixes the
mean structure and nothing searches over it. Only the error structure is
adjusted when a diagnostic asks.
[`vignette("causal-models")`](https://huttoncp.github.io/illume/articles/causal-models.md)
covers all three.

## Reading a fit back in words

[`ilm_interpret()`](https://huttoncp.github.io/illume/reference/ilm_interpret.md)
writes out what a model says – each effect on the scale the response is
measured on, how strong the evidence is, what the diagnostics found, and
how far to trust the estimates:

``` r

ilm_interpret(fit)
```

The prose is templated rather than generated, so the same fit gives the
same words. It says nothing about bias that a diagnostic did not
measure, and it uses causal language only when the object carries a
design that licenses it – otherwise every effect is “associated with”.

[`ilm_ame()`](https://huttoncp.github.io/illume/reference/ilm_ame.md)
gives the piece that makes this readable for a model with a link
function: the average marginal effect on the response scale. An odds
ratio of 2 can be a 3-point change in probability or a 20-point one
depending on where the data sit, so the change in probability is what
gets reported.

``` r

ilm_ame(fit)
```

## Ordered outcomes

An ordered response treated as unordered throws the ordering away and
pays for it in parameters: a multinomial spends `J - 1` coefficients per
predictor where a cumulative link model spends one. Treated as a number
it invents a spacing the categories do not have – “agree” is not twice
“neutral”.

``` r

fit <- ilm_model(rating ~ dose + age, data = dd, family = "ordinal")
ilm_thresholds(fit)
```

The model is `P(Y <= j) = F(threshold_j - eta)`, so the linear predictor
is **subtracted** from each cut and a positive coefficient pushes
probability towards the **higher** categories. Reading it the other way
reverses every conclusion, which is why
[`summary()`](https://rdrr.io/r/base/summary.html) prints the convention
beside the thresholds. `"ordinal_probit"` and `"ordinal_cloglog"` change
the link.

One coefficient describing every cut is the proportional-odds
assumption, and it is an assumption:

``` r

ilm_check_proportional(fit)
```

If a term acts differently at different cuts, the remedy is
`family = "multinomial"`, which spends more coefficients and assumes
nothing about the ordering. It nests this model, so
[`AIC()`](https://rdrr.io/r/stats/AIC.html) says whether the extra
parameters earn their place – and the ordering is real information, so
giving all of it up to accommodate one predictor may cost more than it
buys.

## Counts with too many zeros

More zeros than a count model expects is the one misspecification a
residual plot cannot show, because a zero is a perfectly ordinary value
of a count.

``` r

ilm_check_zeros(fit)
fit <- ilm_model(visits ~ age, data = dd, family = "poisson", ziformula = ~ insured)
ilm_zi_coef(fit)
```

`zi_type` picks between two models that are **not** the same model.
`"inflated"` is a mixture: some units are structural zeros and the count
part can produce zeros of its own, so a zero has two possible origins
and the probability is the structural share only. `"hurdle"` is two
processes: whether the response clears zero, and how far past it goes,
the second fitted to a count that cannot be zero. Which applies is a
question about the subject – whether there is a subpopulation that was
never at risk – rather than one the fit can settle, though
[`AIC()`](https://rdrr.io/r/stats/AIC.html) does separate them.

## Proportions

A proportion that is not a count of anything – percent cover, share of
time, a score already scaled to the unit interval – has no denominator
for a binomial fit, and a gaussian one puts mass outside `[0, 1]` and
assumes a spread that does not shrink at the ends.

``` r

fit <- ilm_model(cover ~ elevation, data = dd, family = "beta")
```

The parameterisation is a mean and a **precision**, so a larger
dispersion parameter means *less* spread – the opposite of every other
one here. `dispformula` lets the precision depend on a covariate.

Values sitting exactly at 0 or 1 have no likelihood, because the density
has no mass there, and the error says so along with the two honest
responses: they are a separate process, in which case model them as one
with `ziformula`, or they are rounding, in which case
[`ilm_squeeze()`](https://huttoncp.github.io/illume/reference/ilm_squeeze.md)
applies the Smithson-Verkuilen shift and reports how far it moved
things.

## When a regressor is not exogenous

A treatment people chose rather than were assigned, a price set in
response to demand, an exposure measured with error. If an instrument is
available – something that moves the regressor without acting on the
outcome any other way –
[`ilm_iv()`](https://huttoncp.github.io/illume/reference/ilm_iv.md) uses
it.

``` r

fit <- ilm_iv(earnings ~ schooling + age | distance + age, data = dd)
```

Everything to the right of the bar is the **full** list of exogenous
variables, instruments and covariates alike. Leaving a covariate out
does not mean it is fine, it means it is endogenous with no instrument,
and that is an error rather than a silent one.

Three diagnostics print unasked: the first-stage F, against both the
familiar threshold of 10 and the roughly 104 a conventional 5% t-test
actually needs; Durbin-Wu-Hausman, asking whether the regressor was
endogenous at all, since if it was not then least squares was unbiased
and far more precise; and Sargan’s J when there is more than one
instrument. When the instrument is weak,
[`ilm_iv_ar()`](https://huttoncp.github.io/illume/reference/ilm_iv_ar.md)
gives a confidence set that stays valid however weak it is – and can
come back unbounded, which is the honest answer a Wald interval hides.

## A complex sample

Sampling weights are not replicate counts. Handed to a model that treats
weights as frequencies they give standard errors too small by at least
`sqrt(n / sum(w))`, because the likelihood believes it saw `sum(w)`
observations.

``` r

des <- ilm_design(dd, weights = ~ wt, ids = ~ psu, strata = ~ stratum)
fit <- ilm_model(y ~ x, data = dd, family = "binomial", design = des)
ilm_svy_coef(fit)
```

The coefficients are fine in a fixed-effects fit – a weighted likelihood
is design-consistent for the population parameter – which is what makes
the mistake easy to miss: everything looks right except the uncertainty,
and the uncertainty looks better than right. The degrees of freedom are
the number of sampling units less the number of strata, not the number
of rows.

## Mediation

How much of an effect runs **through** something else.

``` r

fm <- ilm_model(stress ~ treat + age, data = dd, family = "gaussian")
fy <- ilm_model(outcome ~ treat * stress + age, data = dd, family = "gaussian")
md <- ilm_mediate(fm, fy, treat = "treat", mediator = "stress")
ilm_mediate_sens(md)
```

The Baron-Kenny product of coefficients is the special case where the
outcome model is linear and has no treatment-by-mediator interaction.
With an interaction the indirect effect differs by arm and one product
cannot be two numbers; with a binary outcome the product is not on the
scale of the effect at all.
[`ilm_mediate()`](https://huttoncp.github.io/illume/reference/ilm_mediate.md)
estimates the counterfactual definitions, and reproduces the product
exactly where the product is right.

The decomposition needs no unmeasured confounding of the **mediator**
and the outcome, which the data cannot check because the mediator was
not randomised.
[`ilm_mediate_sens()`](https://huttoncp.github.io/illume/reference/ilm_mediate_sens.md)
asks the answerable question instead: how strong would such confounding
have to be before the indirect effect went away, and how that compares
with the predictors you did measure.

## Cluster-robust standard errors

When the mean structure is credible but the variance structure is not,
and you would rather not say what the right one is:

``` r

ilm_robust(fit, ~ school)
```

The defaults are CR2 with a t reference on Bell-McCaffrey degrees of
freedom, because that is what covers: with a cluster-level predictor and
10 clusters, the usual CR0-with-a-normal-reference covers a nominal 95%
interval 80% of the time and this covers 95%. Those degrees of freedom
are far below the number of clusters – a treatment assigned to three
schools of thirty does not have twenty-nine behind it – and reporting
them is the warning.
