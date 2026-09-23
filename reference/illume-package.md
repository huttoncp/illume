# illume: A Unified Engine for Exploration and Frequentist Inference

A frequentist regression workflow aimed at **inference** rather than
prediction, together with the exploratory steps that precede it. One
engine fits gaussian, binomial, Poisson, negative binomial, multinomial,
beta, ordinal, zero-inflated, hurdle and survival models, so every model
returns the same object and shares the same methods and diagnostics.

## What problem this solves

Fitting a range of models in R usually means a range of packages, each
with its own object, its own
[`anova()`](https://rdrr.io/r/stats/anova.html) semantics, its own
[`predict()`](https://rdrr.io/r/stats/predict.html) arguments and its
own idea of what a diagnostic plot is. `illume` fits them with one
engine, so the interface does not change as the model gets more
complicated – from `y ~ x` through to a multinomial mixed model with
smooths and correlated errors.

Where an exact small-sample reference exists it is used: a gaussian
model with no random effects reports t and F tests agreeing with
[`stats::lm()`](https://rdrr.io/r/stats/lm.html), rather than the
large-sample approximations the Laplace machinery would otherwise give.

The **multinomial** family closes a specific gap.
[`nnet::multinom()`](https://rdrr.io/pkg/nnet/man/multinom.html) handles
unordered categorical outcomes but has no random effects;
[`lme4::glmer()`](https://rdrr.io/pkg/lme4/man/glmer.html) and `glmmTMB`
have rich random-effect machinery but no multinomial family; `brms` fits
the full model but is far slower and answers a different inferential
question.

## Families

`"gaussian"`, `"binomial"`, `"poisson"`, `"nbinom"`, `"beta"`,
`"multinomial"`, the ordinal families, and the survival families
`"weibull"`, `"lognormal"`, `"loglogistic"` and the flexible `"rp"`; see
[`ilm_family()`](https://huttoncp.github.io/illume/reference/ilm_family.md).
A count or a proportion can have a zero part, inflated or hurdle
(`ziformula`). By default the family is read off the response, and the
fit says which it chose and why. Tweedie is out of scope, because
`glmmTMB` covers it well.

## How the model is written

With `J` outcome categories the model has `C = J - 1` free linear
predictors. `illume` uses **sum-to-zero** coding across categories, so
each coefficient is that category's deviation from the average across
categories, *not* a contrast against a baseline category. This differs
from [`nnet::multinom()`](https://rdrr.io/pkg/nnet/man/multinom.html)
and `brms`, which use baseline-category coding, and matters when reading
output: see
[`ilm_coef_table()`](https://huttoncp.github.io/illume/reference/ilm_coef_table.md).

## Why the diagnostics are prominent

This model class fails *quietly*. A fit can return entirely plausible
coefficients while its covariance matrix is unusable, so the standard
errors are meaningless even though nothing looks wrong. Every fit
therefore carries a table of checks (see `fit$checks`), each with a
status, a plain-language reason, and a suggested remedy.
[`summary.ilm_model()`](https://huttoncp.github.io/illume/reference/summary.ilm_model.md)
prints them.

## Two spellings

Every exported function whose name starts `ilm_` also answers to `iml_`,
an easy transposition to type:
[`iml_model()`](https://huttoncp.github.io/illume/reference/ilm_model.md)
is
[`ilm_model()`](https://huttoncp.github.io/illume/reference/ilm_model.md)
itself, not a wrapper, and opens the same help page.

## Main entry points

- [`ilm_model()`](https://huttoncp.github.io/illume/reference/ilm_model.md):

  fit a model from a formula, for any family

- [`ilm_family()`](https://huttoncp.github.io/illume/reference/ilm_family.md):

  what each response distribution assumes

- [`summary.ilm_model()`](https://huttoncp.github.io/illume/reference/summary.ilm_model.md):

  coefficient table plus assumption checks

- [`ilm_anova()`](https://huttoncp.github.io/illume/reference/ilm_anova.md):

  Type II/III analysis of deviance for fixed effects

- [`predict.ilm_model()`](https://huttoncp.github.io/illume/reference/predict.ilm_model.md):

  category probabilities, conditional or population-averaged

- [`ilm_appraise()`](https://huttoncp.github.io/illume/reference/ilm_appraise.md):

  residual diagnostic panels

- [`ilm_pb_lrt()`](https://huttoncp.github.io/illume/reference/ilm_pb_lrt.md):

  parametric-bootstrap likelihood-ratio test

## References

Agresti, A. (2013). *Categorical Data Analysis*, 3rd ed. Wiley. (Chapter
8 covers baseline-category logit models.)

Kristensen, K., Nielsen, A., Berg, C. W., Skaug, H., & Bell, B. M.
(2016). TMB: Automatic differentiation and Laplace approximation.
*Journal of Statistical Software*, 70(5), 1–21.
[doi:10.18637/jss.v070.i05](https://doi.org/10.18637/jss.v070.i05)

Bates, D., Maechler, M., Bolker, B., & Walker, S. (2015). Fitting linear
mixed-effects models using lme4. *Journal of Statistical Software*,
67(1), 1–48.
[doi:10.18637/jss.v067.i01](https://doi.org/10.18637/jss.v067.i01)

Brooks, M. E., Kristensen, K., van Benthem, K. J., Magnusson, A., Berg,
C. W., Nielsen, A., Skaug, H. J., Maechler, M., & Bolker, B. M. (2017).
glmmTMB balances speed and flexibility among packages for zero-inflated
generalized linear mixed modeling. *The R Journal*, 9(2), 378–400.

Wood, S. N. (2017). *Generalized Additive Models: An Introduction with
R*, 2nd ed. Chapman & Hall/CRC.

## See also

Useful links:

- <https://github.com/huttoncp/illume>

- Report bugs at <https://github.com/huttoncp/illume/issues>

## Author

**Maintainer**: Craig Hutton <craig.hutton@gmail.com>
