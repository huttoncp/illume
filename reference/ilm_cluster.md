# Cluster observations, choosing the number of clusters

Groups the rows of an
[`ilm_reduce()`](https://huttoncp.github.io/illume/reference/ilm_reduce.md)
result, or any numeric coordinates, by k-means or hierarchical
clustering. When `k` is not given it is searched for over `1..k_max` by
the gap statistic. Every cluster gets a stability score from bootstrap
resampling, every observation gets a silhouette width, and clusters
holding less than `small_cluster_frac` of the data are flagged.

## Usage

``` r
ilm_cluster(
  x,
  k = NULL,
  k_max = 10,
  method = c("kmeans", "hclust"),
  dist_method = "euclidean",
  hclust_method = "ward.D2",
  nstart = 25,
  B = 100,
  gap_method = "firstSEmax",
  small_cluster_frac = 0.05,
  ambiguous_threshold = 0.1,
  seed = NULL,
  progress = NULL
)
```

## Arguments

- x:

  An
  [`ilm_reduce()`](https://huttoncp.github.io/illume/reference/ilm_reduce.md)
  result, or a numeric matrix or data frame of coordinates with one row
  per observation.

- k:

  Number of clusters. Omit to choose it by the gap statistic.

- k_max:

  Largest `k` to consider in that search. Ignored when `k` is given.

- method:

  `"kmeans"` or `"hclust"`.

- dist_method:

  Distance for `method = "hclust"`. Euclidean by default, which is right
  here because the clustering runs on continuous dimension-reduction
  coordinates rather than raw mixed-type columns.

- hclust_method:

  Linkage for `method = "hclust"`; `"ward.D2"` by default.

- nstart:

  Random starts for k-means, guarding against a poor local optimum.

- B:

  Bootstrap resamples, used both for the gap statistic and for the
  stability score.

- gap_method:

  Rule for reading `k` off the gap curve; see
  [`cluster::maxSE()`](https://rdrr.io/pkg/cluster/man/clusGap.html).

- small_cluster_frac:

  Clusters holding less than this fraction of the data are flagged.

- ambiguous_threshold:

  Silhouette width below which an observation is flagged individually.

- seed:

  Random seed.

- progress:

  Show a progress bar. Defaults to
  [`interactive()`](https://rdrr.io/r/base/interactive.html), so a bar
  appears when someone is watching and nothing is written in a script or
  a knitted document. See
  [ilm_progress_arg](https://huttoncp.github.io/illume/reference/ilm_progress_arg.md).

## Value

An object of class `"ilm_cluster"`: `method`, `k`, `gap` (the search, or
`NULL` when `k` was given), `clusters` (one row per cluster, with
`size`, `pct`, `jaccard`, `stability` labelled `"stable"` above 0.75,
`"moderate"` from 0.5 to 0.75 and `"unstable"` below, following Hennig,
plus `mean_silhouette` and `anomalous`), `ind_cluster` (one row per
observation, with `silhouette`, `is_small_cluster`, `is_ambiguous` and
`is_anomalous`), and `coords`.

## Two different things called anomalous

A **small cluster** is a property of the cluster: a handful of rows that
sit apart from everything else. It may be a real minority pattern or it
may be a data problem, and either way it is worth looking at.

An **ambiguous observation** is a property of the row, taken from its
silhouette width – its average distance to its own cluster against its
average distance to the nearest other one. Near 1 means clearly placed,
near 0 means sitting on a boundary, negative means it is closer to
another cluster than to its own. This is independent of cluster size: a
point can be ambiguous inside a large and otherwise stable cluster,
which a size-based flag alone would never show. Kaufman and Rousseeuw
treat widths below 0.25 as no substantial structure; the default
threshold here is stricter at 0.1.

## References

Hennig, C. (2007). Cluster-wise assessment of cluster stability.
Computational Statistics and Data Analysis 52(1).

Tibshirani, R., Walther, G. and Hastie, T. (2001). Estimating the number
of clusters in a data set via the gap statistic. JRSS B 63(2).

## See also

[`ilm_reduce()`](https://huttoncp.github.io/illume/reference/ilm_reduce.md),
[`ilm_profile()`](https://huttoncp.github.io/illume/reference/ilm_profile.md),
[`ilm_plot_cluster()`](https://huttoncp.github.io/illume/reference/ilm_plot_cluster.md).

## Examples

``` r
cl <- ilm_cluster(ilm_reduce(mtcars), k_max = 5, B = 25, seed = 1)
cl
#> <ilm_cluster> method = kmeans, k = 4 (chosen by gap statistic)
#> 
#>   cluster   size    pct  jaccard  stability     sil
#>         1      5   15.6    0.834     stable   0.281
#>         2      7   21.9    0.944     stable   0.363
#>         3     12   37.5    1.000     stable   0.625
#>         4      8   25.0    0.881     stable   0.487
```
