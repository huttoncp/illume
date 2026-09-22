# Compare estimated marginal means

Differences between the means from
[`ilm_emmeans()`](https://huttoncp.github.io/illume/reference/ilm_emmeans.md),
with intervals and a multiplicity adjustment.

## Usage

``` r
ilm_contrast(
  object,
  method = c("pairwise", "trt.vs.ctrl", "poly"),
  ref = NULL,
  adjust = c("max_t", "bonferroni", "none"),
  level = 0.95,
  nsim = 20000L,
  seed = 1L
)
```

## Arguments

- object:

  An
  [`ilm_emmeans()`](https://huttoncp.github.io/illume/reference/ilm_emmeans.md)
  result.

- method:

  `"pairwise"`, `"trt.vs.ctrl"` or `"poly"`. The last assumes the levels
  are ordered and equally spaced.

- ref:

  Reference level for `"trt.vs.ctrl"`.

- adjust:

  `"max_t"`, `"bonferroni"` or `"none"`.

- level:

  Confidence level for the family.

- nsim:

  Draws used to find the studentized-maximum critical value.

- seed:

  Random seed for that simulation.

## Value

A data frame with `contrast`, `estimate`, `se`, `lower`, `upper`,
`p_value`, `p_adj` and `adjust`.

## The adjustment

Comparing every pair of five groups is ten tests, and ten intervals each
nominally 95% do not jointly cover at 95%. The default is the
single-step studentized maximum: the contrasts' joint covariance is
known exactly here, `C L V L' C'`, so the reference distribution is
simulated from it directly rather than bootstrapped. This is the same
construction
[`ilm_boot_diff()`](https://huttoncp.github.io/illume/reference/ilm_boot_diff.md)
uses on raw data, where it reproduced
[`TukeyHSD()`](https://rdrr.io/r/stats/TukeyHSD.html) to 0.006 on the
design Tukey is exact for.

`"bonferroni"` is the conservative fallback and `"none"` is there for
comparisons chosen in advance. With a single contrast there is nothing
to adjust and the column reads `"none"` whatever was asked.

## See also

[`ilm_emmeans()`](https://huttoncp.github.io/illume/reference/ilm_emmeans.md),
[`ilm_boot_diff()`](https://huttoncp.github.io/illume/reference/ilm_boot_diff.md)
for the same comparison made without a model.

## Examples

``` r
set.seed(1); n <- 200
d <- data.frame(g = factor(sample(c("a", "b", "c"), n, TRUE)), x = rnorm(n))
d$y <- 1 + 0.5 * (d$g == "b") + 0.2 * d$x + rnorm(n)
fit <- ilm_model(y ~ g + x, data = d, family = "gaussian", verbose = FALSE)
ilm_contrast(ilm_emmeans(fit, "g"))
#> <ilm_contrast> 3 comparison(s), adjust = max_t 
#> 
#>  contrast  estimate     se   lower  upper  p_adj
#>     b - a  0.009823 0.1764 -0.4069 0.4265 0.9984
#>     c - a -0.226155 0.1780 -0.6467 0.1943 0.4142
#>     c - b -0.235978 0.1768 -0.6538 0.1818 0.3800
#> 
#>   Intervals hold jointly at 95% across all 3 comparisons.
```
