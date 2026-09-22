# Test what the graph claims against the data

Every missing edge in a DAG is a testable claim: those two variables
should be independent given the conditioning set. This fits that claim
and reports whether the data agree. It is the step that turns a DAG from
an assumption into something the data can argue with.

## Usage

``` r
ilm_dag_test(
  g,
  data,
  min_effect = 0.1,
  alpha = 0.05,
  adjust = "holm",
  verbose = TRUE
)
```

## Arguments

- g:

  An
  [`ilm_dag()`](https://huttoncp.github.io/illume/reference/ilm_dag.md).

- data:

  A data frame holding the graph's observed variables.

- min_effect:

  Smallest partial correlation worth calling a contradiction.

- alpha:

  Level for the adjusted p-values.

- adjust:

  Multiplicity adjustment across claims, passed to
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html). A DAG
  can imply dozens of claims and testing them all at 0.05 guarantees
  false alarms, so the default is Holm.

- verbose:

  Print progress.

## Value

A data frame with one row per claim: `x`, `y`, `given`, `n`, `estimate`
(partial correlation), `p_value`, `p_adj`, `verdict` and `note`. The
overall verdict is attached as the `verdict` attribute.

## Details

A failed claim does not say which of the two is at fault. The graph may
be missing an arrow, or the data may not measure what the graph's node
names suppose – a mismeasured covariate, a selected sample, a variable
that means something different from what it is called. Both are worth
knowing before any effect estimate is taken seriously.

Each claim is tested by a likelihood-ratio test of the exposure term in
a generalised linear model for one of the pair given the other and the
conditioning set, with the family taken from the response's type. That
is a screen rather than a full model: it assumes the conditional mean is
linear on the link scale, so a curved dependence can pass. It will not
test a claim between two variables that are both multi-level factors,
and says so in the `note` column rather than quietly dropping it.

**Two thresholds, not one.** With a few thousand rows, a correlation of
0.03 is significant and means nothing, so a claim counts as contradicted
only when it is both statistically significant after adjustment and
larger than `min_effect` in partial correlation. This is the stance
[`ilm_variogram()`](https://huttoncp.github.io/illume/reference/ilm_variogram.md)
takes for the same reason.

## See also

[`ilm_dag_implied()`](https://huttoncp.github.io/illume/reference/ilm_dag_implied.md)
for the claims themselves.

## Examples

``` r
set.seed(1)
n <- 400
z <- rnorm(n); x <- 0.6 * z + rnorm(n); y <- 0.5 * x + 0.4 * z + rnorm(n)
d <- data.frame(x = x, y = y, z = z, w = rnorm(n))
## w is unconnected, so every claim involving it should hold
g <- ilm_dag("dag { x [exposure] ; y [outcome] ; z -> x -> y ; z -> y }")
ilm_dag_test(g, cbind(d["x"], d["y"], d["z"]), verbose = FALSE)
#> [1] x        y        n        estimate p_value  p_adj    verdict  note    
#> [9] given   
#> <0 rows> (or 0-length row.names)
```
