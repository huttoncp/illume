# Expectation of a function of a normal variable

`E[h(eta + sd * Z)]` for `Z` standard normal, by the Gauss-Hermite
quadrature `predict(groups = "population")` uses for a population
average – the same nodes and the same rule for how many – so code built
on a fit averages over a latent spread exactly as the fit's own
predictions do.

## Usage

``` r
ilm_normal_expect(h, eta, sd, n = NULL)
```

## Arguments

- h:

  A vectorised function of one argument.

- eta, sd:

  The mean and the standard deviation of the normal, recycled to a
  common length.

- n:

  Optional number of nodes, in place of the rule.

## Value

The expectation for each element: a vector shaped as `h(eta)`.

## Details

The number of nodes grows with the spread: 40 times the smallest whole
number at least `sd^2`, and no more than 800, with one rule for the
whole call, set by the largest `sd`. The rules are computed once and
cached. Measured against
[`stats::integrate()`](https://rdrr.io/r/stats/integrate.html), through
an inverse link that is better than 1e-12 with a logit up to `sd = 6`,
and better than 1e-9 with a complementary log-log up to `sd = 5`. For
any other integrand, check the accuracy it needs; `n` sets the number of
nodes directly.

## See also

[`predict.ilm_model()`](https://huttoncp.github.io/illume/reference/predict.ilm_model.md),
whose `groups = "population"` uses the same rule.

## Examples

``` r
## a population-averaged probability through a logit
ilm_normal_expect(stats::plogis, eta = 0.5, sd = 1.2)
#> [1] 0.5960873
## the probability of a count of 3 in a new group, under a log link
ilm_normal_expect(function(e) stats::dpois(3, exp(e)), eta = 1, sd = 0.6)
#> [1] 0.1559136
```
