# Profile a data set: reduce, cluster, and describe the clusters

Runs
[`ilm_reduce()`](https://huttoncp.github.io/illume/reference/ilm_reduce.md),
then
[`ilm_cluster()`](https://huttoncp.github.io/illume/reference/ilm_cluster.md)
on the dimensions it produces, then says what distinguishes each cluster
– which dimensions it sits unusually far along, and which original
variables those dimensions are made of. The result is a sentence per
cluster rather than a column of numbers.

## Usage

``` r
ilm_profile(
  data,
  cols = NULL,
  ndim = 5,
  method = c("pcamix", "glrm"),
  ...,
  vtest_threshold = 1.96,
  top_n_vars = 2,
  var_contrib = TRUE,
  var_contrib_B = 199L
)
```

## Arguments

- data:

  A data frame.

- cols:

  Columns to use. A character vector of names, a regular expression, a
  predicate function such as `is.numeric`, or `NULL` for all of them –
  see
  [ilm_selection](https://huttoncp.github.io/illume/reference/ilm_selection.md).

- ndim:

  Number of dimensions to keep.

- method:

  `"pcamix"` (the default) for PCA, MCA or FAMD depending on the column
  types, in closed form. `"glrm"` fits a generalized low rank model
  instead, which uses a loss appropriate to each column's type rather
  than squared error on one-hot indicators, and reconstructs a category
  as a category. It costs an iterative fit, and on all-numeric data the
  two are the same model – see
  [`ilm_glrm()`](https://huttoncp.github.io/illume/reference/ilm_glrm.md)
  for when it is worth that.

- ...:

  Passed to
  [`ilm_cluster()`](https://huttoncp.github.io/illume/reference/ilm_cluster.md),
  for instance `k`, `k_max`, `method`, `B` or `seed`.

- vtest_threshold:

  Smallest `|vtest|` for a dimension to count toward a cluster's
  description.

- top_n_vars:

  How many top-loading variables to name per dimension.

- var_contrib:

  Run
  [`ilm_var_contrib()`](https://huttoncp.github.io/illume/reference/ilm_var_contrib.md)
  on the result and report it. Nothing in this pipeline selects
  variables, and an irrelevant one is not neutral: on three known
  clusters with two informative columns, adding two pure-noise factors
  took recovery from 0.301 to 0.006. So the check runs here rather than
  waiting to be asked for, and a clustering that turns out to be one
  variable's levels under another name raises a warning – the cluster
  descriptions would otherwise be read at face value, and they would all
  be true and all about that variable. Costs roughly 70% on top of the
  clustering at 500 rows and five columns.

- var_contrib_B:

  Permutations for that check.

## Value

An object of class `"ilm_profile"`: `reduce`, `cluster`,
`characterization` (one row per dimension that characterises a cluster),
`var_contrib` (or `NULL`) and `summary`, one sentence per cluster.

## How a cluster gets characterised

By a v-test: for cluster `c` and dimension `d`,
`(mean_c - mean_overall) / (sd_overall / sqrt(n_c))`. It measures how
surprising the cluster's average position on that dimension would be if
membership had nothing to do with it.

It is **not a p-value**. The clusters were found from the very
coordinates being tested, so the usual sampling argument does not apply
and the 1.96 default is a threshold that behaves sensibly rather than a
5% test. Treat a characterisation as a description of the partition you
have, not as evidence that the partition is real – the stability and
silhouette figures in
[`ilm_cluster()`](https://huttoncp.github.io/illume/reference/ilm_cluster.md)
are what speak to that.

For numeric variables the reported direction combines the v-test's sign
with the variable's own correlation on that dimension, since a dimension
can run either way. Categorical variables are named without a direction.

## See also

[`ilm_reduce()`](https://huttoncp.github.io/illume/reference/ilm_reduce.md),
[`ilm_cluster()`](https://huttoncp.github.io/illume/reference/ilm_cluster.md),
[`ilm_plot_profile()`](https://huttoncp.github.io/illume/reference/ilm_plot_profile.md),
[`ilm_profile_na()`](https://huttoncp.github.io/illume/reference/ilm_profile_na.md).

## Examples

``` r
p <- ilm_profile(mtcars, k_max = 5, B = 25, seed = 1)
cat(p$summary, sep = "\n")
#> Cluster 1 (n = 5, 15.6% of the data, stable) is characterised by dim 2 (low qsec, high gear).
#> Cluster 2 (n = 7, 21.9% of the data, stable) is characterised by dim 3 (high carb, high vs); dim 2 (high qsec, low gear).
#> Cluster 3 (n = 12, 37.5% of the data, stable) is characterised by dim 1 (high cyl, high disp).
#> Cluster 4 (n = 8, 25.0% of the data, stable) is characterised by dim 1 (low cyl, low disp).
```
