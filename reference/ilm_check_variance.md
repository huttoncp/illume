# Is the residual spread constant?

The assumption illume had no check for. Heteroscedasticity does not
usually bias the coefficients much, but it does bias their standard
errors, which is what an inference package is for.

## Usage

``` r
ilm_check_variance(
  object,
  by = NULL,
  B = 100L,
  seed = 1L,
  ncores = 1L,
  plot = TRUE,
  verbose = TRUE,
  progress = NULL
)
```

## Arguments

- object:

  A fitted `"ilm_model"` object.

- by:

  Optional name of a column in the model frame, or a vector, whose
  levels might have different variances.

- B:

  Simulated datasets. A simulated p-value cannot fall below
  `1 / (B + 1)`, so B must exceed 100 for the strongest verdict to be
  reachable.

- seed:

  Random seed.

- ncores:

  Worker processes for the refits.

- plot:

  Draw the scale-location panel.

- verbose:

  Print the verdict.

- progress:

  Show a progress bar. Defaults to
  [`interactive()`](https://rdrr.io/r/base/interactive.html), so a bar
  appears when someone is watching and nothing is written in a script or
  a knitted document. See
  [illumex::ilm_progress_arg](https://huttoncp.github.io/illumex/reference/ilm_progress_arg.html).

## Value

Invisibly, a list with the observed statistics, their simulated nulls,
p-values, a `status`, a suggested remedy, and `by`: the column the
groups were taken from, or `NA` when `by` was a vector or not given.
Pass the list to
[`ilm_remedies()`](https://huttoncp.github.io/illume/reference/ilm_remedies.md)
to have the remedy written out.

## Details

Two patterns are tested, because they call for different fixes: spread
that trends with the fitted value, and spread that differs between the
levels of a grouping variable. Both are measured on randomised quantile
residuals, which are standard normal when the model and its variance are
both correct, and both are compared against a null simulated from the
fitted model rather than against an assumed distribution.

## See also

[`ilm_rqr_test()`](https://huttoncp.github.io/illume/reference/ilm_rqr_test.md)
for other residual statistics,
[`ilm_check_dispersion()`](https://huttoncp.github.io/illume/reference/ilm_check_dispersion.md)
when the whole response is over-dispersed rather than unevenly
dispersed.

## Examples

``` r
set.seed(1)
d <- ilm_sim(n_id = 25)
f <- ilm_model(score ~ income + (1 | id), data = d, family = "gaussian",
               verbose = FALSE)
ilm_check_variance(f, B = 20, plot = FALSE)
#> Warning: B = 20 puts the smallest achievable p-value at 0.048, so a FAIL verdict is unreachable. Use B >= 100.
#> spread vs fitted: rho = 0.043, p = 0.4118
#> OK
```
