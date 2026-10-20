# Coefficients of the zero part of a zero-inflated or hurdle fit

The second linear predictor a `ziformula` adds, on the logit scale: a
positive coefficient means more excess zeros. What "excess" means
depends on which model was fitted – under `zi_type = "inflated"` these
are the structural zeros only, over and above the ones the count part
produces by itself, while under `"hurdle"` they are every zero there is.

## Usage

``` r
ilm_zi_coef(object, level = 0.95)
```

## Arguments

- object:

  A fitted
  [`ilm_model()`](https://huttoncp.github.io/illume/reference/ilm_model.md)
  with a `ziformula`.

- level:

  Confidence level for the interval.

## Value

A data frame of `term`, `estimate`, `se`, `lower`, `upper`, `z` and `p`,
with an `odds_ratio` column since the scale is a logit.

## See also

[`ilm_zi_prob()`](https://huttoncp.github.io/illume/reference/ilm_zi_prob.md)
for the fitted probabilities themselves,
[`ilm_check_zeros()`](https://huttoncp.github.io/illume/reference/ilm_check_zeros.md)
for whether a zero part is called for at all.

## Examples

``` r
set.seed(1); n <- 400
d <- data.frame(x = rnorm(n), z = rnorm(n))
lam <- exp(0.6 + 0.4 * d$x)
d$y <- ifelse(runif(n) < plogis(-0.5 + 0.8 * d$z), 0, rpois(n, lam))
f <- ilm_model(y ~ x, data = d, family = "poisson", ziformula = ~ z,
               verbose = FALSE)
ilm_zi_coef(f)
#>          term   estimate        se      lower       upper odds_ratio         z
#> 1 (Intercept) -0.3965041 0.1685590 -0.7268736 -0.06613458  0.6726675 -2.352317
#> 2           z  1.0212491 0.1789195  0.6705733  1.37192501  2.7766611  5.707868
#>              p
#> 1 1.865688e-02
#> 2 1.144001e-08
```
