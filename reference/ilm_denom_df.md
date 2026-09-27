# Denominator degrees of freedom for one or more contrasts

Turns a Wald statistic into an F or t test by working out how much the
estimated standard error is itself worth trusting. Without this a mixed
model reports a chi-square, which treats the variance components as
known and is anti-conservative when there are few clusters.

## Usage

``` r
ilm_denom_df(
  object,
  L,
  method = c("auto", "satterthwaite", "kenward-roger", "residual", "asymptotic"),
  h = 1e-05
)
```

## Arguments

- object:

  A fitted `"ilm_model"` object.

- L:

  A contrast vector of length `p`, or a `q`-by-`p` contrast matrix.

- method:

  One of `"auto"`, `"satterthwaite"`, `"kenward-roger"`, `"residual"`,
  `"asymptotic"`. `"auto"` gives the exact residual df when nothing was
  integrated out, Satterthwaite for a gaussian mixed model, and `Inf`
  otherwise.

- h:

  Step size for the numerical derivatives.

## Value

A list with `df`, the `method` actually used, and for Kenward-Roger the
adjusted covariance `V` and `scale`, the factor by which the Wald F
formed with `V` is multiplied before it is referred to F(q, df).

## Which method

`"satterthwaite"` matches a scaled chi-square to the estimated variance
of the contrast. It is the default because it costs first derivatives
only, it applies to every structure this package fits, and it is what
`lmerTest` reports.

`"kenward-roger"` also inflates the covariance for the uncertainty in
the variance components, which matters most when clusters are few, and
gives an F test its own denominator df and scale factor (Kenward and
Roger 1997). It is computed as `pbkrtest` computes it, and agrees with
it. It needs a REML fit (`reml = TRUE`) of a gaussian model whose random
terms are intercepts and slopes – correlated or not, nested or crossed –
with one residual variance and no weights, and it forms the covariance
of all the observations, so it stops above 4000 rows. Asked for
elsewhere it stops and says why. For one contrast its scale is exactly
1: the t statistic is the estimate over its standard error from the
adjusted covariance `V`, on the df returned.

`"residual"` is `N - p`, which is exact with nothing integrated out and
optimistic otherwise. `"asymptotic"` returns `Inf`, recovering the
chi-square.

## What it is not for

Both approximations are derived for LINEAR mixed models. For a
non-gaussian family the small-sample problem is a different one and
neither answers it;
[`ilm_pb_lrt()`](https://huttoncp.github.io/illume/reference/ilm_pb_lrt.md)
simulates the null rather than approximating its reference, and is the
remedy there.

## References

Kenward, M. G., & Roger, J. H. (1997). Small sample inference for fixed
effects from restricted maximum likelihood. *Biometrics*, 53(3),
983–997.

Halekoh, U., & Hojsgaard, S. (2014). A Kenward-Roger approximation and
parametric bootstrap methods for tests in linear mixed models: the R
package pbkrtest. *Journal of Statistical Software*, 59(9), 1–30.

## See also

[`ilm_anova()`](https://huttoncp.github.io/illume/reference/ilm_anova.md),
[`ilm_emmeans()`](https://huttoncp.github.io/illume/reference/ilm_emmeans.md),
[`ilm_pb_lrt()`](https://huttoncp.github.io/illume/reference/ilm_pb_lrt.md).

## Examples

``` r
set.seed(1)
d <- ilm_sim(n_id = 25)
f <- ilm_model(score ~ income + (1 | id), data = d, family = "gaussian",
               verbose = FALSE)
L <- c(0, 1)                      # the slope on income
ilm_denom_df(f, L)
#> $df
#> [1] 282.5906
#> 
#> $method
#> [1] "satterthwaite"
#> 
#> $V
#> NULL
#> 
```
