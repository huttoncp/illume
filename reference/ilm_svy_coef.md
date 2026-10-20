# Fixed effects with design-based standard errors

Fixed effects with design-based standard errors

## Usage

``` r
ilm_svy_coef(
  object,
  design = NULL,
  level = 0.95,
  lonely = c("adjust", "certainty", "fail")
)
```

## Arguments

- object:

  A fitted
  [`ilm_model()`](https://huttoncp.github.io/illume/reference/ilm_model.md)
  from a complex sample.

- design:

  An
  [`ilm_design()`](https://huttoncp.github.io/illume/reference/ilm_design.md),
  or `NULL` to use the fit's own.

- level:

  Confidence level.

- lonely:

  See
  [`ilm_svy_vcov()`](https://huttoncp.github.io/illume/reference/ilm_svy_vcov.md).

## Value

A data frame of `term`, `estimate`, `se`, `se_model`, `t`, `p`, `lower`
and `upper`, on the design degrees of freedom.

## See also

[`ilm_design()`](https://huttoncp.github.io/illume/reference/ilm_design.md),
[`ilm_svy_vcov()`](https://huttoncp.github.io/illume/reference/ilm_svy_vcov.md).

## Examples

``` r
set.seed(1); n <- 600
d <- data.frame(x = rnorm(n), psu = rep(1:60, each = 10),
                st = rep(1:2, each = 300))
d$w <- ifelse(d$st == 1, 40, 10)
d$y <- 1 + 0.5 * d$x + rep(rnorm(60, 0, 0.8), each = 10) + rnorm(n)
des <- ilm_design(d, weights = ~ w, ids = ~ psu, strata = ~ st)
fit <- ilm_model(y ~ x, data = d, family = "gaussian", design = des,
                 verbose = FALSE)
ilm_svy_coef(fit)
#> Fixed effects with design-based (linearization) standard errors
#>   60 primary sampling unit(s) in 2 stratum(a); 58 degrees of freedom
#>         term estimate      se se_model df     t        p  lower  upper
#>  (Intercept)  0.97200 0.12970  0.01111 58 7.496 4.29e-10 0.7125 1.2320
#>            x  0.39787 0.07624  0.01134 58 5.218 2.53e-06 0.2453 0.5505
#> 
#>   design se / model-based se: 6.72 to 11.67
#>   The model-based column is there to be ignored: a weighted likelihood
#>   believes it saw sum(w) observations, so its errors are too small by
#>   at least sqrt(n / sum(w)), and by more once the sample is clustered.
```
