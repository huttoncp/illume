# How far a variable is from gaussian, and why

Reports an agreement index on 0-1 and, when it is low, the reason.

## Usage

``` r
ilm_gauss_check(x, min_n = 20L, cap = 0.12)
```

## Arguments

- x:

  A numeric vector.

- min_n:

  Below this many observations the index is not computed.

- cap:

  The excess distance at which agreement reaches 0. The default 0.12 is
  set so clear departures (lognormal 0.110, t(3) 0.105, bimodal 0.108,
  Poisson 0.148) land at the bottom of the scale.

## Value

A list with `gauss` (0-1), `ks_d` (the raw distance) and `gauss_note`
(empty when agreement is high).

## Details

This deliberately does not report a normality test p-value. Any test of
exact normality rejects everything once `n` is large, so the p-value
answers a question nobody asked. What matters is how far from normal a
variable is, which is an effect size.

The measure is the Kolmogorov distance between the data and the
best-fitting normal: the largest amount, on the cumulative probability
scale, by which a normal model misstates the data. Sampling noise alone
produces a distance of about `0.6/sqrt(n)`, so the index charges only
the excess beyond what normal data of the same size would give. It is an
agreement index, not a probability.

## References

Lilliefors, H. W. (1967). On the Kolmogorov-Smirnov test for normality
with mean and variance unknown. Journal of the American Statistical
Association, 62(318), 399-402.

## Examples

``` r
set.seed(1)
ilm_gauss_check(rnorm(500))
#> $gauss
#> [1] 1
#> 
#> $ks_d
#> [1] 0.037
#> 
#> $gauss_note
#> [1] ""
#> 
ilm_gauss_check(rlnorm(500))
#> $gauss
#> [1] 0
#> 
#> $ks_d
#> [1] 0.2189
#> 
#> $gauss_note
#> [1] "bounded at zero; right-skewed"
#> 
ilm_gauss_check(c(rnorm(250), rnorm(250, 5)))
#> $gauss
#> [1] 0.135
#> 
#> $ks_d
#> [1] 0.1441
#> 
#> $gauss_note
#> [1] "multimodal (check for subgroups); light-tailed / flat"
#> 
```
