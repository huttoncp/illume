# Thresholds of a cumulative link fit

The cut points on the latent scale, with their standard errors. There
are `J - 1` of them for `J` categories and they are increasing by
construction: the fit estimates the first and the logarithms of the
gaps, so an ordering violation is not something the optimiser can reach.

## Usage

``` r
ilm_thresholds(object, level = 0.95)
```

## Arguments

- object:

  A fitted
  [`ilm_model()`](https://craig-hutton.github.io/illume/reference/ilm_model.md)
  with an ordinal family.

- level:

  Confidence level.

## Value

A data frame of `cut`, `estimate`, `se`, `lower` and `upper`.

## Details

Read them against the linear predictor, which is **subtracted**:
`P(Y <= j) = F(threshold_j - eta)`. A positive coefficient therefore
pushes probability towards the higher categories.

## See also

[`ilm_check_proportional()`](https://craig-hutton.github.io/illume/reference/ilm_check_proportional.md)
for the assumption that makes one coefficient per predictor enough.

## Examples

``` r
set.seed(1); n <- 400
d <- data.frame(x = rnorm(n))
z <- 0.8 * d$x + rlogis(n)
d$y <- factor(cut(z, c(-Inf, -1, 1, Inf), labels = c("lo", "mid", "hi")),
              ordered = TRUE)
f <- ilm_model(y ~ x, data = d, family = "ordinal", verbose = FALSE)
ilm_thresholds(f)
#>      cut   estimate        se      lower      upper
#> 1 lo|mid -0.7810335 0.1112162 -0.9990132 -0.5630538
#> 2 mid|hi  0.9322873 0.1138474  0.7091505  1.1554241
```
