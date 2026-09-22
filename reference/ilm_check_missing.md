# What the missing values look like, and whether they matter

Reports how much is missing, which variables go missing together, and
what missingness is related to – then says what that implies for the
analysis.

## Usage

``` r
ilm_check_missing(
  data,
  y = NULL,
  covariates = NULL,
  min_effect = 0.1,
  alpha = 0.05,
  adjust = "holm",
  verbose = TRUE
)
```

## Arguments

- data:

  A data frame.

- y:

  Optional outcome column, named as a string. Naming it is what
  separates "complete cases are fine" from "complete cases are biased",
  so it is worth naming. A formula is also accepted and is usually the
  clearer way to write it: `outcome ~ x + z` sets the outcome and takes
  the right-hand side as `covariates`.

- covariates:

  Variables to hold fixed when asking whether missingness depends on the
  outcome. Default is every other column; pass the model's predictors
  when they are a subset.

- min_effect:

  Smallest association worth reporting, as a correlation or Cramer's V.
  With several thousand rows an association of 0.02 is significant and
  means nothing.

- alpha:

  Level for the adjusted p-values.

- adjust:

  Multiplicity adjustment across pairs, passed to
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html).

- verbose:

  Narrate the findings.

## Value

An object of class `"ilm_missing"`: `variables`, `patterns`,
`associations`, `monotone` and a `verdict`.

## What can and cannot be established

**MCAR** – missingness unrelated to anything – is testable, and is
tested here by asking whether each variable's missingness is associated
with the observed values of the others.

**MAR versus MNAR is not testable.** The data that would separate them
are precisely the data that are missing. No amount of pattern analysis
settles it, and this function does not pretend otherwise: it reports
what missingness is associated with among the things you *can* see, and
says plainly that a dependence on the unseen values themselves cannot be
ruled out.

## Why "not MCAR" does not mean "impute"

For a regression of `y` on covariates, dropping incomplete rows is
unbiased whenever missingness is independent of `y` **given** the
covariates – even when it depends on the covariates strongly. That is a
far weaker condition than MCAR, and it is why naming `y` changes the
advice:

- missingness related to the covariates but not to `y`: complete cases
  are unbiased. Imputation can recover precision but is not needed for
  correctness, and it adds assumptions.

- missingness related to `y` itself: complete cases are biased, and
  [`ilm_impute()`](https://craig-hutton.github.io/illume/reference/ilm_impute.md)
  is the remedy.

## See also

[`ilm_impute()`](https://craig-hutton.github.io/illume/reference/ilm_impute.md),
[`ilm_describe_na_all()`](https://craig-hutton.github.io/illume/reference/ilm_describe_na_all.md),
[`ilm_plot_missing()`](https://craig-hutton.github.io/illume/reference/ilm_plot_missing.md).

## Examples

``` r
d <- ilm_sim()
ilm_check_missing(d, y = "score", verbose = FALSE)
#> <ilm_missing> 792 of 900 rows complete (88.0%)
#> 
#>   missing by variable
#>     lab_value               108   12.0%
#> 
#>   pattern: monotone (2 distinct)
#>   verdict: MCAR_NOT_REJECTED 
#>     consistent with MCAR; complete cases are unbiased 
#> 
#>   does missingness depend on `score` once the others are held fixed?
#>     lab_value        partial r = 0.029, p = 0.44
#> 
#>   MNAR cannot be ruled out by any test; that needs a sensitivity analysis.
```
