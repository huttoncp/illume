# Cluster-robust and heteroskedasticity-robust covariance

A sandwich covariance for the fixed effects: the model's own covariance
on the outside, and what the data actually did on the inside. Use it
when the mean structure is credible but the variance structure is not –
spread that changes with a predictor, or observations that are
correlated within a school, firm, household or subject.

## Usage

``` r
ilm_vcov_cluster(object, cluster = NULL, type = c("CR2", "CR1", "CR0"))
```

## Arguments

- object:

  A fitted
  [`ilm_model()`](https://craig-hutton.github.io/illume/reference/ilm_model.md)
  with no random effects.

- cluster:

  A one-sided formula naming variables in the model frame (`~ school`,
  or `~ school + year` for their interaction), or a vector with one
  entry per row. `NULL` treats every row as its own cluster.

- type:

  `"CR2"`, `"CR1"` or `"CR0"`; see above.

## Value

A `p x p` covariance matrix, with the number of clusters and the type in
its attributes.

## Which correction, and why it matters

A sandwich is consistent in the number of **clusters**, not the number
of rows, and with few clusters the raw version is biased downwards:

- `"CR0"` is the raw sandwich, with no correction.

- `"CR1"` multiplies it by `G/(G-1) * (N-1)/(N-p)`, the constant Stata
  applies. It is better than nothing and does not depend on the design.

- `"CR2"` (the default) adjusts each cluster's residuals by that
  cluster's own leverage, after Bell and McCaffrey. It is the one that
  works when the clusters are unbalanced, which is when the problem is
  worst.

With `cluster = NULL` each row is its own cluster and these become the
familiar heteroskedasticity-robust HC0, HC1 and HC2.

## How few is too few

There is no threshold that makes the problem go away, but below roughly
40 clusters the correction and the reference distribution both start to
matter a great deal, and below about 15 no adjustment reliably rescues
the interval – the estimator is being asked to learn a covariance from a
dozen numbers.
[`ilm_robust()`](https://craig-hutton.github.io/illume/reference/ilm_robust.md)
reports the count and says so. The honest alternative there is a design
with more clusters, or a model that says what the correlation is rather
than working around it.

## References

Bell, R. M. and McCaffrey, D. F. (2002). Bias reduction in standard
errors for linear regression with multi-stage samples. *Survey
Methodology* 28, 169-181.

Cameron, A. C. and Miller, D. L. (2015). A practitioner's guide to
cluster-robust inference. *Journal of Human Resources* 50, 317-372.

## See also

[`ilm_robust()`](https://craig-hutton.github.io/illume/reference/ilm_robust.md)
for a coefficient table built on it.

## Examples

``` r
set.seed(1); G <- 40; n <- 10
d <- data.frame(g = factor(rep(seq_len(G), each = n)), x = rnorm(G * n))
d$y <- 1 + 0.5 * d$x + rep(rnorm(G, 0, 1.5), each = n) + rnorm(G * n)
f <- ilm_model(y ~ x, data = d, family = "gaussian", verbose = FALSE)
sqrt(diag(vcov(f)))                          # model-based, too small here
#> (Intercept)           x 
#>  0.08739040  0.09017861 
sqrt(diag(ilm_vcov_cluster(f, ~ g)))         # cluster-robust
#> (Intercept)           x 
#>  0.22044028  0.08571531 
```
