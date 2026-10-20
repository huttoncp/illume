# A fitted model's response distribution, as functions

The density, distribution function, quantile function, random generator
and mean of a fit's response, in illume's own parameterisation, taking
linear predictors and internal parameters – a draw from
[`ilm_draws()`](https://huttoncp.github.io/illume/reference/ilm_draws.md),
say – rather than means.

## Usage

``` r
ilm_dist(object)
```

## Arguments

- object:

  A fitted `"ilm_model"`.

## Value

An object of class `"ilm_dist"`: a list with the functions
`d(y, lp, par, log = FALSE, status = NULL, size = 1)`,
`p(q, lp, par, lower.tail = TRUE, size = 1)` (`NULL` for a multinomial
response), `q(prob, lp, par, size = 1)` (likewise),
`r(lp, par, size = 1)` and `mean(lp, par, size = 1)`, and `family`,
`discrete`, `zero` (the kind of zero part, or `NULL`), `categories`,
`lp_names` and `par_names`.

## Details

**Inputs.** `lp` is a named list of linear predictors on the link scale:
`mu`, the mean's (an `n`-by-`C` matrix for a multinomial response);
`zi`, the zero part's, for a zero-inflated or hurdle model; and `disp`,
the dispersion model's, for a model with `dispformula` – its data
columns only, since the dispersion's power of the mean is added here
from `par["mu_pow"]`. `par` holds the family's own parameters, named as
the blocks of
[`ilm_draws()`](https://huttoncp.github.io/illume/reference/ilm_draws.md):
`logdisp` (the log dispersion), `zeta_raw` (an ordered response's
thresholds, as fitted: a first value and log increments) and `mu_pow`.
Each entry of `lp` and `par` may have one value or one per element, so a
vector of rows by draws is scored in one call.

**Censoring.** `status` is `-1` for a value known only to be at or below
`y`, `1` for one at or above it, and `0` for an observed one, as
[`ilm_censor()`](https://huttoncp.github.io/illume/reference/ilm_censor.md)
codes them; a censored row's density is the probability of its interval.

**A binomial response** is a proportion, with `size` its number of
trials. The density then includes the binomial coefficient, which the
fit's likelihood leaves out (see
[`logLik.ilm_model()`](https://huttoncp.github.io/illume/reference/logLik.ilm_model.md)).

**Categories.** For an ordered or multinomial response `y` is a
category: its position among `categories`, or its label.
[`mean()`](https://rdrr.io/r/base/mean.html) gives each category's
probability. A multinomial response has no ordering, so no distribution
or quantile function.

The flexible parametric survival families are not covered: their linear
predictor is a function of time itself; see
[`ilm_survival()`](https://huttoncp.github.io/illume/reference/ilm_survival.md).

## See also

[`ilm_draws()`](https://huttoncp.github.io/illume/reference/ilm_draws.md),
[`ilm_matrices()`](https://huttoncp.github.io/illume/reference/ilm_matrices.md).

## Examples

``` r
set.seed(1)
d <- data.frame(x = rnorm(200))
d$y <- rpois(200, exp(0.5 + 0.4 * d$x))
fit <- ilm_model(y ~ x, data = d, family = "poisson", verbose = FALSE)
dd <- ilm_dist(fit)
lp <- list(mu = as.numeric(fit$X %*% fit$beta))
## the fit's log-likelihood, row by row
sum(dd$d(d$y, lp, par = NULL, log = TRUE))
#> [1] -323.5795
logLik(fit)
#> 'log Lik.' -323.5795 (df=2)
```
