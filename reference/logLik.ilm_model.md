# Log-likelihood, and number of observations

Returns the **Laplace-approximate** marginal log-likelihood – the
likelihood with the random effects integrated out, using the
approximation described in
[`ilm_fit()`](https://huttoncp.github.io/illume/reference/ilm_fit.md).
The `df` attribute counts fixed effects *and* covariance parameters,
which is the convention
[`AIC.ilm_model()`](https://huttoncp.github.io/illume/reference/AIC.ilm_model.md)
and
[`BIC.ilm_model()`](https://huttoncp.github.io/illume/reference/AIC.ilm_model.md)
use.

Returns the **Laplace-approximate** marginal log-likelihood – the
likelihood with the random effects integrated out, using the
approximation described in
[`ilm_fit()`](https://huttoncp.github.io/illume/reference/ilm_fit.md).
The `df` attribute counts fixed effects *and* covariance parameters,
which is the convention
[`AIC.ilm_model()`](https://huttoncp.github.io/illume/reference/AIC.ilm_model.md)
and
[`BIC.ilm_model()`](https://huttoncp.github.io/illume/reference/AIC.ilm_model.md)
use.

## Usage

``` r
# S3 method for class 'ilm_model'
logLik(object, ...)

# S3 method for class 'ilm_model'
nobs(object, ...)
```

## Arguments

- object:

  A fitted `"ilm_model"` object.

- ...:

  Unused.

## Value

An object of class `"logLik"`.

An object of class `"logLik"`.

## Details

Because it is an approximation, differences between models are most
trustworthy when the random structure is held fixed, so that the
approximation error largely cancels.

The value is on the same scale as
[`lme4::lmer()`](https://rdrr.io/pkg/lme4/man/lmer.html) and
[`nlme::lme()`](https://rdrr.io/pkg/nlme/man/lme.html), and agrees with
both to numerical tolerance for a gaussian model they can also fit. It
therefore includes the normalising constants, so
[`AIC()`](https://rdrr.io/r/stats/AIC.html) and
[`BIC()`](https://rdrr.io/r/stats/AIC.html) may be compared across
models with *different* random structures as well as the same one.

One constant is left out: for a binomial response with trials – a
proportion weighted by its number of trials – the binomial coefficients,
`sum(lchoose(trials, successes))`, which
[`glm()`](https://rdrr.io/r/stats/glm.html), `lme4` and `glmmTMB`
include. It is a constant of the data, so it changes no comparison
between models of the same data, but the log-likelihoods of the two
differ by exactly that much.

Because it is an approximation, differences between models are most
trustworthy when the random structure is held fixed, so that the
approximation error largely cancels.

## See also

[`AIC.ilm_model()`](https://huttoncp.github.io/illume/reference/AIC.ilm_model.md),
[`BIC.ilm_model()`](https://huttoncp.github.io/illume/reference/AIC.ilm_model.md),
[`ilm_anova()`](https://huttoncp.github.io/illume/reference/ilm_anova.md).

[`AIC.ilm_model()`](https://huttoncp.github.io/illume/reference/AIC.ilm_model.md),
[`BIC.ilm_model()`](https://huttoncp.github.io/illume/reference/AIC.ilm_model.md),
[`ilm_anova()`](https://huttoncp.github.io/illume/reference/ilm_anova.md).
