# Fit a multinomial mixed model from a design matrix

The computational engine. Most users should call
[`ilm_model()`](https://huttoncp.github.io/illume/reference/ilm_model.md)
with a formula instead; this is the entry point when you already have a
design matrix, and it is what the formula interface calls internally.

## Usage

``` r
ilm_fit(
  X,
  y,
  J = NULL,
  re_list,
  re_struct = NULL,
  ar = NULL,
  ylevels = NULL,
  weights = NULL,
  family = "multinomial",
  verbose = TRUE,
  restarts = 3L,
  joint = FALSE,
  censor = NULL,
  Zd = NULL,
  disp_mu = FALSE,
  rp = NULL,
  Zzi = NULL,
  zi_type = c("inflated", "hurdle"),
  reml = FALSE,
  boundary = c("hold", "avoid")
)
```

## Arguments

- X:

  Numeric design matrix for the fixed effects, with `N` rows.

- y:

  Integer vector of length `N` giving the observed category, coded `1`
  to `J`.

- J:

  Integer. Number of outcome categories, for the multinomial family
  only; `NULL` for every other family.

- re_list:

  Named list of random terms. Each element is a grouping vector (random
  intercept), a `list(group =, Z =)` (random slopes), or a
  `list(basis =)` (a smooth, from `ilm_smooth()`).

- re_struct:

  Optional named list of category covariance structures, named by term
  as `re_list` is. A term it leaves out gets `"us"`, so only the terms
  that differ need naming. Set `d_cor = FALSE` within an element to drop
  an intercept-slope correlation.

- ar:

  Optional AR(1) specification: `list(idx =, n_group =, Tt =)`.

- ylevels:

  Optional character vector of category labels, used in output.

- weights:

  Optional numeric vector of **frequency** weights: the number of
  replicate observations each row stands for, exactly as in a binomial
  [`glm()`](https://rdrr.io/r/stats/glm.html) fitted to grouped data. A
  weighted fit is then numerically identical to expanding each row,
  standard errors included.

  Valid only when the pooled observations share every random effect in
  the model – you may aggregate within subject and covariate pattern,
  but not across subjects. These are **not** sampling weights; see the
  `weights_type` check and the note below.

- family:

  Response distribution: a name, or the object returned by
  [`ilm_family()`](https://huttoncp.github.io/illume/reference/ilm_family.md).
  See
  [`ilm_family()`](https://huttoncp.github.io/illume/reference/ilm_family.md)
  for what each one assumes.

- verbose:

  Logical. Print the checks while fitting.

- restarts:

  Integer. Number of optimiser restarts from the previous solution,
  which helps on difficult surfaces.

- joint:

  Logical. Also compute the joint precision over fixed and random
  parameters. Needed by
  [`predict.ilm_model()`](https://huttoncp.github.io/illume/reference/predict.ilm_model.md)
  to propagate uncertainty in penalised smooth coefficients;
  [`ilm_model()`](https://huttoncp.github.io/illume/reference/ilm_model.md)
  switches it on automatically when the model contains smooths.

- censor:

  Optional censoring specification from
  [`ilm_censor()`](https://huttoncp.github.io/illume/reference/ilm_censor.md),
  marking observations known only as an interval – at or below a floor,
  at or above a ceiling. Supported for the gaussian family, where it
  gives a Tobit model.

- Zd:

  Optional design matrix for the dispersion model, one row per
  observation. Its columns become a linear predictor for the logarithm
  of the dispersion, so its intercept replaces the single dispersion
  parameter. Built by
  [`ilm_model()`](https://huttoncp.github.io/illume/reference/ilm_model.md)
  from `dispformula`.

- disp_mu:

  Logical. Add a term in the logarithm of the fitted mean to that
  predictor, making the dispersion a power of the mean.

- rp:

  Internal. The flexible parametric baseline: knots, the derivative
  design and which columns of `X` hold the spline. Built by
  [`ilm_model()`](https://huttoncp.github.io/illume/reference/ilm_model.md)
  from `rp_df`.

- Zzi:

  Optional design matrix for the zero part of a count model, one row per
  observation, modelling the logit of an excess-zero probability. Built
  by
  [`ilm_model()`](https://huttoncp.github.io/illume/reference/ilm_model.md)
  from `ziformula`.

- zi_type:

  `"inflated"` for a mixture, where a zero may have come either from the
  zero process or from the count itself, or `"hurdle"`, where every zero
  comes from the zero process and the positives come from a count that
  cannot be zero. Ignored when `Zzi` is `NULL`.

- reml:

  Logical. Integrate the fixed effects out along with the random ones,
  giving restricted maximum likelihood. Under a flat prior that integral
  IS the restricted likelihood, and for a linear-gaussian model the
  Laplace approximation to it is exact – so this is REML rather than an
  approximation to it. Gaussian responses only; elsewhere the integral
  is still well defined but has none of REML's properties, so it is
  refused. See
  [`ilm_model()`](https://huttoncp.github.io/illume/reference/ilm_model.md)
  for when to switch it on.

- boundary:

  `"hold"` (maximum likelihood; a covariance that reaches the edge of
  its range is held there) or `"avoid"` (the boundary-avoiding penalty
  of Chung et al. 2013, 2015). See
  [`ilm_model()`](https://huttoncp.github.io/illume/reference/ilm_model.md)
  for what each implies.

## Value

An object of class `"ilm_model"`: a list whose most useful elements are
`checks` (the diagnostic table), `Sigma` (fitted category covariances),
`beta` (fixed effects as a p-by-C matrix), `opt` (optimiser output),
`sdr` ([`TMB::sdreport()`](https://rdrr.io/pkg/TMB/man/sdreport.html)
output) and `ok` (whether every check passed).

## How the model is parameterised

With `J` outcome categories there are `C = J - 1` free linear
predictors. Coefficients use **sum-to-zero** coding, so each one is that
category's deviation from the average across categories, not a contrast
against a baseline. This differs from
[`nnet::multinom()`](https://rdrr.io/pkg/nnet/man/multinom.html) and
`brms`.

Random effects are held as a matrix rather than a long vector, with a
*matrix-normal* prior. Writing it this way means the covariance over
groups and the covariance over categories never have to be combined into
one large Kronecker product, which keeps both memory and computation
manageable and lets AR(1), spatial and i.i.d. structures share a single
code path.

## Choosing a category covariance structure

Each random term gets its own structure through `re_struct`:

- `"us"`:

  unstructured: every variance and correlation free. Costs `C(C+1)/2`
  parameters, which grows quickly – 45 at `J = 10`.

- `"diag"`:

  diagonal: categories uncorrelated. Costs `C`.

- `"rr"` with `rank = r`:

  reduced rank: `Sigma = Lambda Lambda'` with `Lambda` of size C-by-r.
  Costs fewer parameters **and** fewer latent values, which is usually
  the better trade when a term is stretched.

If a structure is too rich for the data the pre-fit checks will say so
and name a specific rank to try.

## Models with no random effects

`re_list` may be empty. The latent parameter block is then omitted
entirely rather than declared with length zero, and for a gaussian
response the fit reports exact t and F inference instead of the
large-sample approximations. See
[`ilm_model()`](https://huttoncp.github.io/illume/reference/ilm_model.md).

## Fitting method

The random effects are integrated out by the **Laplace approximation**,
implemented through `RTMB`. Exact quadrature is not an option here: with
crossed random effects the integral does not factorise, so its dimension
is the total number of latent values, often in the thousands. Laplace is
an approximation, and its accuracy depends on having enough data per
latent value – see the `latent_budget` entry in `fit$checks`, and
[`ilm_consistency()`](https://huttoncp.github.io/illume/reference/ilm_consistency.md)
to test it directly on your own fit.

## A warning about sampling weights

Passing inverse-probability or survey weights here is not supported and
will mislead you. In a *mixed* model the weight multiplies the
conditional likelihood inside the Laplace integral, so each cluster
appears to carry more information than it does; both the standard errors
and the coefficients are affected. Proper handling needs scaled
level-specific weights and a design-based variance estimator, neither of
which is implemented. The `weights_type` check flags weights that look
like sampling weights.

## References

Kristensen, K., Nielsen, A., Berg, C. W., Skaug, H., & Bell, B. M.
(2016). TMB: Automatic differentiation and Laplace approximation.
*Journal of Statistical Software*, 70(5), 1–21.

Agresti, A. (2013). *Categorical Data Analysis*, 3rd ed. Wiley.

Pfeffermann, D., Skinner, C. J., Holmes, D. J., Goldstein, H., &
Rasbash, J. (1998). Weighting for unequal selection probabilities in
multilevel models. *Journal of the Royal Statistical Society, Series B*,
60(1), 23–40.

Rabe-Hesketh, S., & Skrondal, A. (2006). Multilevel modelling of complex
survey data. *Journal of the Royal Statistical Society, Series A*,
169(4), 805–827.

## See also

[`ilm_model()`](https://huttoncp.github.io/illume/reference/ilm_model.md)
for the formula interface,
[`summary.ilm_model()`](https://huttoncp.github.io/illume/reference/summary.ilm_model.md),
[`ilm_consistency()`](https://huttoncp.github.io/illume/reference/ilm_consistency.md).
