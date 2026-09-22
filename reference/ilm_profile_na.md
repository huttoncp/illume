# Profile which values are missing, and for whom

The missingness counterpart to
[`ilm_profile()`](https://craig-hutton.github.io/illume/reference/ilm_profile.md):
runs
[`ilm_reduce_na()`](https://craig-hutton.github.io/illume/reference/ilm_reduce_na.md)
then
[`ilm_cluster_na()`](https://craig-hutton.github.io/illume/reference/ilm_cluster_na.md),
and describes each cluster of rows by which columns' *missingness* sets
it apart – "cluster 2 is characterised by dim 1 (missing: income,
missing: age)" rather than by those columns' values.

## Usage

``` r
ilm_profile_na(
  data,
  cols = NULL,
  ndim = 5,
  ...,
  vtest_threshold = 1.96,
  top_n_vars = 2
)
```

## Arguments

- data:

  A data frame.

- cols:

  Columns to use. A character vector of names, a regular expression, a
  predicate function such as `is.numeric`, or `NULL` for all of them –
  see
  [ilm_selection](https://craig-hutton.github.io/illume/reference/ilm_selection.md).

- ndim:

  Number of dimensions to keep.

- ...:

  Passed to
  [`ilm_cluster_na()`](https://craig-hutton.github.io/illume/reference/ilm_cluster_na.md).

- vtest_threshold:

  Smallest `|vtest|` for a dimension to count.

- top_n_vars:

  How many top-loading indicators to name per dimension.

## Value

An object of class `"ilm_profile_na"`, which is also an `"ilm_profile"`.

## Details

This is how a structured gap shows itself: a block of variables that go
missing together points at a shared cause, such as a section of a form
everyone in one group skipped, which is a different problem from values
going missing one at a time.
[`ilm_check_missing()`](https://craig-hutton.github.io/illume/reference/ilm_check_missing.md)
then says whether any of it threatens the model you intend to fit.

## See also

[`ilm_check_missing()`](https://craig-hutton.github.io/illume/reference/ilm_check_missing.md),
[`ilm_impute()`](https://craig-hutton.github.io/illume/reference/ilm_impute.md),
[`ilm_profile()`](https://craig-hutton.github.io/illume/reference/ilm_profile.md).

## Examples

``` r
p <- ilm_profile_na(airquality, k_max = 4, B = 25, seed = 1)
#> ilm_reduce_na(): dropping column(s) whose missingness never varies (always or never missing): Wind, Temp, Month, Day
cat(p$summary, sep = "\n")
#> Cluster 1 (n = 35, 22.9% of the data, stable) is characterised by dim 2 (missing: Ozone, missing: Solar.R).
#> Cluster 2 (n = 2, 1.3% of the data, stable) is characterised by dim 1 (missing: Ozone, missing: Solar.R). It is a small cluster, 1.3% of observations: possibly a real minority pattern, possibly a data problem, but worth looking at either way.
#> Cluster 3 (n = 111, 72.5% of the data, stable) is characterised by dim 1 (missing: Ozone, missing: Solar.R).
#> Cluster 4 (n = 5, 3.3% of the data, stable) is characterised by dim 2 (missing: Ozone, missing: Solar.R). It is a small cluster, 3.3% of observations: possibly a real minority pattern, possibly a data problem, but worth looking at either way.
```
