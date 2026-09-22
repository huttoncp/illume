# Describes how a response distribution enters the model. The random-effect machinery – covariance structures, random slopes, smooths, AR(1) – is identical for every family and never sees the response; a family supplies only the piece that turns a linear predictor into a log-likelihood.

Describes how a response distribution enters the model. The
random-effect machinery – covariance structures, random slopes, smooths,
AR(1) – is identical for every family and never sees the response; a
family supplies only the piece that turns a linear predictor into a
log-likelihood.

## Usage

``` r
ilm_family(
  family = c("gaussian", "binomial", "poisson", "nbinom", "beta", "multinomial",
    "ordinal", "ordinal_probit", "ordinal_cloglog", "weibull", "lognormal",
    "loglogistic", "rp", "rp_odds", "rp_normal")
)
```

## Arguments

- family:

  Character: one of `"gaussian"`, `"binomial"`, `"poisson"`, `"nbinom"`,
  `"multinomial"`.

## Value

A list describing the family, with elements `name`, `link`, `n_disp`
(number of dispersion parameters), `disp_names`, `C_of()` (number of
linear predictor dimensions), `nll()`, `linkinv()` and `sim()`.

## Details

Supported families:

- `"gaussian"`:

  continuous response, identity link, one dispersion parameter (the
  residual standard deviation).

- `"binomial"`:

  two-category response, logit link. The response may be 0/1, a
  two-level factor, `TRUE`/`FALSE`, or a character column – all four are
  treated identically. For a factor or character response the **second**
  level is the one being modelled, as in
  [`stats::glm()`](https://rdrr.io/r/stats/glm.html), so
  `factor(c("no", "yes"))` models the probability of `"yes"`. A
  proportion with `weights` giving the number of trials also works.

- `"poisson"`:

  non-negative counts, log link, no dispersion parameter.

- `"nbinom"`:

  counts with more variability than Poisson allows, log link, one
  dispersion parameter. Variance is `mu + mu^2/k`, so smaller `k` means
  more overdispersion; as `k` grows it approaches the Poisson.

- `"multinomial"`:

  three or more unordered categories. The only family needing more than
  one linear predictor: with `J` categories it uses `J - 1` dimensions
  with sum-to-zero coding.

Zero-inflation, hurdle models, Tweedie and other specialised families
are deliberately out of scope. `glmmTMB` covers those well and there is
nothing to gain from a weaker reimplementation.

## On numerical stability

The binomial log-likelihood uses `logspace_add()` rather than
`log(1 + exp(eta))`. They are mathematically identical, but the naive
form overflows to infinity once `eta` exceeds about 709, which an
optimiser can easily reach while exploring. The stable form is exact
throughout.

## References

McCullagh, P., & Nelder, J. A. (1989). *Generalized Linear Models*, 2nd
ed. Chapman & Hall.

Hilbe, J. M. (2011). *Negative Binomial Regression*, 2nd ed. Cambridge
University Press.

## Examples

``` r
ilm_family("poisson")$n_disp      # Poisson has no dispersion parameter
#> [1] 0
ilm_family("nbinom")$disp_names
#> [1] "log_k"
```
