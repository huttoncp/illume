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
adjusted covariance `V`.

## Which method

`"satterthwaite"` matches a scaled chi-square to the estimated variance
of the contrast. It is the default because it costs first derivatives
only, it applies to every structure this package fits, and it is what
`lmerTest` reports.

`"kenward-roger"` also inflates the covariance for the uncertainty in
the variance components, which matters most when clusters are few. It
costs second derivatives, and it is only defined where the marginal
covariance is linear in the variance components – random intercepts and
slopes with a single residual variance. Asked for elsewhere it stops and
says why.

`"residual"` is `N - p`, which is exact with nothing integrated out and
optimistic otherwise. `"asymptotic"` returns `Inf`, recovering the
chi-square.

## What it is not for

Both approximations are derived for LINEAR mixed models. For a
non-gaussian family the small-sample problem is a different one and
neither answers it;
[`ilm_pb_lrt()`](https://craig-hutton.github.io/illume/reference/ilm_pb_lrt.md)
simulates the null rather than approximating its reference, and is the
remedy there.

## See also

[`ilm_anova()`](https://craig-hutton.github.io/illume/reference/ilm_anova.md),
[`ilm_emmeans()`](https://craig-hutton.github.io/illume/reference/ilm_emmeans.md),
[`ilm_pb_lrt()`](https://craig-hutton.github.io/illume/reference/ilm_pb_lrt.md).

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
