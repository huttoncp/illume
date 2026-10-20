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
cached.

The rule is centred where `h(eta + sd * z)` times the normal density
peaks, when that is more than two standard deviations from zero. For a
bounded `h`, such as an inverse link, it rarely is, and the rule is the
plain one. For an `h` that grows, the mass moves away: under a log link
it sits at `z = sd`, where a rule centred at zero has almost no nodes.

Accuracy, measured against a fine grid over `eta` from -12 to 8:

- better than 1e-12 through a logit up to `sd = 6`;

- better than 5e-9 through a complementary log-log up to `sd = 5`, and
  5e-8 at 6;

- better than 1e-13 for `exp`, whose expectation is `exp(eta + sd^2/2)`,
  up to `sd = 20`.

Past that, an integrand growing as fast as `exp` overflows double
precision at the nodes, and the result is `Inf`, never `NaN`. For any
other integrand, check the accuracy it needs; `n` sets the number of
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
