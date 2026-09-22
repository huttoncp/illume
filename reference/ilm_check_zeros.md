# Are there more zeros than the model expects?

More zeros than a count model expects is the one misspecification that
cannot be read off a residual plot, because a zero is a perfectly
ordinary value of a count. This simulates from the fit and asks how
often it produces as many zeros as were seen.

## Usage

``` r
ilm_check_zeros(object, B = 200L, seed = 1L)
```

## Arguments

- object:

  A fitted `"ilm_model"` object.

- B:

  Simulated datasets. A simulated p-value cannot fall below
  `1 / (B + 1)`, so B must exceed 100 for the strongest verdict to be
  reachable at all.

- seed:

  Random seed.

## Value

Invisibly, a list with the observed and expected zero counts, a p-value
and a `status`.

## Details

The remedy is [ilm_model(ziformula =
)](https://huttoncp.github.io/illume/reference/ilm_model.md), which adds
a second linear predictor for the probability of an excess zero. Whether
that should be a mixture (`zi_type = "inflated"`) or two processes
(`"hurdle"`) is a question about what the zeros mean rather than one
this check can answer.

## Examples

``` r
set.seed(1)
d <- ilm_sim()
f <- ilm_model(downtime ~ income + (1 | id), data = d, family = "poisson",
               verbose = FALSE)
#> Warning: NA/NaN function evaluation
#> Warning: NA/NaN function evaluation
#> Warning: NA/NaN function evaluation
#> Warning: NA/NaN function evaluation
#> Warning: NA/NaN function evaluation
#> Warning: NA/NaN function evaluation
ilm_check_zeros(f, B = 200)
#> ilm_check_zeros: observed 576 zeros against about 237 expected; refit with ilm_model(ziformula = ~ 1), or ziformula = ~ x if the excess depends on a predictor. Use zi_type = "hurdle" when every zero comes from one process and "inflated" when some units were never at risk
```
