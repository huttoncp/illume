# Information criteria

AIC and BIC computed from the Laplace-approximate marginal
log-likelihood, with degrees of freedom counting both fixed effects and
covariance parameters.

AIC and BIC computed from the Laplace-approximate marginal
log-likelihood, with degrees of freedom counting both fixed effects and
covariance parameters.

## Usage

``` r
# S3 method for class 'ilm_model'
AIC(object, ..., k = 2)

# S3 method for class 'ilm_model'
BIC(object, ..., n = c("obs", "groups"))
```

## Arguments

- object:

  A fitted `"ilm_model"` object.

- ...:

  Unused.

- k:

  Numeric. Penalty per parameter; 2 gives AIC.

- n:

  `"obs"` or `"groups"`; [`BIC()`](https://rdrr.io/r/stats/AIC.html)
  only.

## Value

A numeric scalar.

A numeric scalar.

## Three caveats worth knowing

**The likelihood is approximate.** Differences between models are most
trustworthy when the random structure is the same in both, so that the
approximation error largely cancels.

**Marginal AIC favours simpler random structures.** It is well suited to
comparing fixed-effect specifications; comparisons of *random*
structures are biased toward the simpler model. Conditional AIC, which
uses effective degrees of freedom, is the right tool there and is not
implemented.

**BIC's sample size is ambiguous in a mixed model.** BIC penalises by
`log(n)`, but what counts as `n`? The number of observations, or the
number of independent groups? For covariance parameters the effective
sample size is closer to the group count. The default `n = "obs"`
matches lme4 and glmmTMB so values are comparable with those packages;
`n = "groups"` gives the alternative. They can select different models,
and that is a genuine disagreement rather than a bug.

A variance sitting at zero also makes the degrees of freedom an
overcount, since that parameter is not freely estimated, making both
criteria conservative.

**The likelihood is approximate.** Differences between models are most
trustworthy when the random structure is the same in both, so that the
approximation error largely cancels.

**Marginal AIC favours simpler random structures.** It is well suited to
comparing fixed-effect specifications; comparisons of *random*
structures are biased toward the simpler model. Conditional AIC, which
uses effective degrees of freedom, is the right tool there and is not
implemented.

**BIC's sample size is ambiguous in a mixed model.** BIC penalises by
`log(n)`, but what counts as `n`? The number of observations, or the
number of independent groups? For covariance parameters the effective
sample size is closer to the group count. The default `n = "obs"`
matches lme4 and glmmTMB so values are comparable with those packages;
`n = "groups"` gives the alternative. They can select different models,
and that is a genuine disagreement rather than a bug.

A variance sitting at zero also makes the degrees of freedom an
overcount, since that parameter is not freely estimated, making both
criteria conservative.

## From a REML fit

A restricted likelihood is the likelihood of contrasts orthogonal to the
fixed-effect design, so changing the fixed effects changes which data it
is the likelihood of. AIC and BIC from REML fits therefore compare only
models with the same fixed effects – different random or dispersion
structures, say. To compare fixed effects, refit each model with
`reml = FALSE`. Both functions say so in a message when given a REML
fit; they cannot refuse, since a single criterion does not know what it
will be compared with. A gaussian model with a dispersion formula and no
random terms is fitted by REML by default (see
[`ilm_model()`](https://huttoncp.github.io/illume/reference/ilm_model.md)).

## References

Vaida, F., & Blanchard, S. (2005). Conditional Akaike information for
mixed-effects models. *Biometrika*, 92(2), 351–370.

Delattre, M., Lavielle, M., & Poursat, M.-A. (2014). A note on BIC in
mixed-effects models. *Electronic Journal of Statistics*, 8, 456–475.

Vaida, F., & Blanchard, S. (2005). Conditional Akaike information for
mixed-effects models. *Biometrika*, 92(2), 351–370.

Delattre, M., Lavielle, M., & Poursat, M.-A. (2014). A note on BIC in
mixed-effects models. *Electronic Journal of Statistics*, 8, 456–475.
