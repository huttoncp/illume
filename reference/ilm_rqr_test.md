# Test residuals against a simulated reference

Takes any single-number summary of the residuals and compares it to the
distribution that summary takes when the model is true, obtained by
simulating from the fit and refitting each replicate.

## Usage

``` r
ilm_rqr_test(
  object,
  B = 30L,
  ncores = 1L,
  seed = 1L,
  stat = function(u) mean(u),
  verbose = TRUE
)
```

## Arguments

- object:

  A fitted `"ilm_model"` object.

- B:

  Integer. Simulated datasets.

- ncores:

  Integer. Worker processes.

- seed:

  Integer. Random seed.

- stat:

  Function taking the residual vector and returning one number.

- verbose:

  Logical. Print the result.

## Value

Invisibly, a list with `status`, the observed statistic, the simulated
null, and a z score and p-value.

## Details

The refitting matters: it builds in the fact that in-sample residuals
are not exactly uniform (see
[`ilm_rqr()`](https://craig-hutton.github.io/illume/reference/ilm_rqr.md)),
so a correctly specified model is not flagged for a discrepancy that
fitting itself created.

## What this can and cannot detect

Broad summaries of the residuals have essentially **no power** against
an omitted covariate. This is structural rather than a matter of tuning:
the reference distribution is "this model refitting data it generated",
and a misspecified model reproduces its own behaviour faithfully, so the
misspecification cancels from both sides of the comparison.

A statistic aimed at specific structure does have power. Use
[`ilm_check_covariate()`](https://craig-hutton.github.io/illume/reference/ilm_check_covariate.md)
or
[`ilm_check_omitted()`](https://craig-hutton.github.io/illume/reference/ilm_check_omitted.md)
to point the test at particular variables, including ones the model does
not contain.
