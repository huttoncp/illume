# Fixed effects with sandwich standard errors

The same coefficients, re-tested against a covariance that does not
assume the variance structure is right. See
[`ilm_vcov_cluster()`](https://huttoncp.github.io/illume/reference/ilm_vcov_cluster.md)
for what the corrections do.

## Usage

``` r
ilm_robust(
  object,
  cluster = NULL,
  type = c("CR2", "CR1", "CR0"),
  df = c("bm", "G-1", "normal"),
  level = 0.95
)
```

## Arguments

- object:

  A fitted
  [`ilm_model()`](https://huttoncp.github.io/illume/reference/ilm_model.md)
  with no random effects.

- cluster:

  Clustering, as in
  [`ilm_vcov_cluster()`](https://huttoncp.github.io/illume/reference/ilm_vcov_cluster.md).
  `NULL` gives heteroskedasticity-robust standard errors.

- type:

  `"CR2"`, `"CR1"` or `"CR0"`.

- df:

  `"bm"` for Bell-McCaffrey (the default), `"G-1"`, `"normal"`, or a
  single number to use for every coefficient.

- level:

  Confidence level.

## Value

A data frame of `term`, `estimate`, `se`, `df`, `t`, `p`, `lower`,
`upper`, with the model-based standard error alongside for comparison.

## The degrees of freedom are the point

A sandwich is consistent in the number of clusters, so with few of them
the normal reference is badly optimistic. The default here is a t
distribution on the Bell-McCaffrey degrees of freedom, computed per
coefficient from the design rather than fixed at `G - 1`. Those can be
far smaller than `G - 1` when one cluster dominates a predictor – a
treatment assigned to three schools out of thirty does not have
twenty-nine degrees of freedom behind it, and the number this reports is
the warning.

## See also

[`ilm_vcov_cluster()`](https://huttoncp.github.io/illume/reference/ilm_vcov_cluster.md).

## Examples

``` r
set.seed(1); G <- 30; n <- 12
d <- data.frame(g = factor(rep(seq_len(G), each = n)), x = rnorm(G * n))
d$y <- 1 + 0.4 * d$x + rep(rnorm(G, 0, 1.5), each = n) + rnorm(G * n)
f <- ilm_model(y ~ x, data = d, family = "gaussian", verbose = FALSE)
ilm_robust(f, ~ g)
#> Fixed effects with cluster-robust standard errors (CR2)
#>   30 clusters
#>         term estimate      se se_model   df     t        p  lower  upper
#>  (Intercept)  0.66577 0.27230  0.09531 30.0 2.445 0.020600 0.1097 1.2220
#>            x  0.35726 0.08816  0.09855 26.9 4.052 0.000387 0.1763 0.5382
#> 
#>   se / model-based se: 0.89 to 2.86
#> 
#>   Under 40 clusters, the correction and the reference distribution
#>   both matter. CR2 with Bell-McCaffrey degrees of freedom is the
#>   least bad combination and is what was used here.
```
