# What each variable contributes to a clustering

Reports, for every variable, how strongly it separates the clusters that
were found. A variable near the bottom is contributing distance without
contributing structure, and dropping it usually improves the result
rather than losing information.

## Usage

``` r
ilm_var_contrib(x, data, B = 199L, seed = 1L)
```

## Arguments

- x:

  An
  [`ilm_profile()`](https://huttoncp.github.io/illume/reference/ilm_profile.md)
  result, or an
  [`ilm_cluster()`](https://huttoncp.github.io/illume/reference/ilm_cluster.md)
  result that carries its reduction.

- data:

  The data frame the clustering came from.

- B:

  Permutations for the reference.

- seed:

  Random seed.

## Value

A data frame, one row per variable, with the separation statistic, the
permutation mean, a p-value and a verdict. Sorted strongest first.

## Details

The statistic is the **between-cluster share of variance** – eta squared
for a numeric variable, Cramer's V for a categorical one – so both are
on a 0 to 1 scale and directly comparable. A permutation reference is
computed by shuffling the cluster labels, which preserves each
variable's own distribution and destroys only its relationship with the
clustering, and the variable is called uninformative when it does not
beat that.

Neither FAMD nor any distance used here selects variables, and an
uninformative one is not harmless: on three known clusters with two
informative columns, adding two pure-noise factors took recovery from
0.301 to 0.006. Checking is cheap and the remedy is `cols`.

## What this cannot tell you

The clustering was fitted to these variables, so a variable the
clustering **used** will separate the clusters it helped make, whether
or not it means anything. On three real clusters plus one pure-noise
numeric column, the noise column scored 0.441 against the informative
columns' 0.576 and 0.550 – visibly weaker, and still "better than
chance", because k-means had split on it.

So this ranks contribution; it does not prove relevance. Two readings
are safe and are what it is for: a variable the clustering **ignored**
(no better than a shuffled label) is adding distance and nothing else,
and a variable that separates the clusters almost perfectly while
nothing else does has not explained the clustering – it *is* the
clustering, and the partition is its levels under another name. That
second case is reported separately because from inside a clustering it
looks like the best possible result.

Settling relevance properly means refitting without each variable and
seeing whether the partition survives, which costs one clustering per
variable and is not done here.

## See also

[`ilm_profile()`](https://huttoncp.github.io/illume/reference/ilm_profile.md),
[`ilm_cluster()`](https://huttoncp.github.io/illume/reference/ilm_cluster.md),
[`ilm_reduce()`](https://huttoncp.github.io/illume/reference/ilm_reduce.md).

## Examples

``` r
# \donttest{
p <- ilm_profile(mtcars, k = 3)
ilm_var_contrib(p, mtcars)
#> Variable contribution to 3 clusters (199 permutations)
#>  variable    type separation permuted     p p_adj                verdict
#>        vs numeric      0.881    0.070 0.005 0.005 carries the clustering
#>       cyl numeric      0.833    0.057 0.005 0.005 carries the clustering
#>      disp numeric      0.716    0.065 0.005 0.005 carries the clustering
#>      gear numeric      0.637    0.066 0.005 0.005 carries the clustering
#>        hp numeric      0.594    0.073 0.005 0.005 carries the clustering
#>      qsec numeric      0.572    0.062 0.005 0.005 carries the clustering
#>       mpg numeric      0.559    0.059 0.005 0.005 carries the clustering
#>      carb numeric      0.553    0.067 0.005 0.005 carries the clustering
#>        wt numeric      0.531    0.059 0.005 0.005 carries the clustering
#>        am numeric      0.516    0.064 0.005 0.005 carries the clustering
#>      drat numeric      0.492    0.065 0.005 0.005 carries the clustering
#> 
#>   Every variable separates the clusters better than chance.
#>   Read the ranking, not the verdict: the clustering was fitted to
#>   these variables, so one it USED separates the clusters it helped
#>   make whether or not it means anything. A pure-noise column scored
#>   0.441 here against 0.576 for a real one.
# }
```
