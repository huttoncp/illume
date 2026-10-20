# Minimal adjustment sets for an exposure effect

Every set of observed variables that, adjusted for, identifies the total
effect of `exposure` on `outcome`, with no redundant member. An empty
result is a finding rather than a failure: it says the effect is not
identifiable from the variables available, which is the most useful
thing a DAG can tell you before any model is fitted.

## Usage

``` r
ilm_adjust_sets(
  g,
  exposure = NULL,
  outcome = NULL,
  observed = NULL,
  max_size = NULL
)
```

## Arguments

- g:

  An
  [`ilm_dag()`](https://huttoncp.github.io/illume/reference/ilm_dag.md).

- exposure, outcome:

  Variable names; taken from the graph if it declares them.

- observed:

  Variables available for adjustment. Defaults to every node the graph
  does not mark unobserved. Pass `names(data)` to restrict to what was
  actually measured.

- max_size:

  Largest set to look for. The default searches every subset of the
  candidate pool when that is affordable and warns when it is not.

## Value

A list of character vectors, one per minimal set, shortest first.
`character(0)` as the only element means the empty set suffices – no
adjustment is needed. A zero-length list means no admissible set exists.

## Details

Several minimal sets are common and worth having. Adjusting for each in
turn should give the same exposure effect if the graph is right, so the
spread across them is a sensitivity analysis that costs only compute –
[`ilm_dag_model()`](https://huttoncp.github.io/illume/reference/ilm_dag_model.md)
runs it automatically.

This is the back-door criterion, which is sufficient but not quite
complete: the full adjustment criterion admits a small number of further
sets containing descendants of the exposure that lie off every causal
path. Where they differ, the sets returned here are a subset of the
valid ones, never a superset, so a set reported here is always
admissible.

## References

Pearl, J. (1995). Causal diagrams for empirical research. Biometrika
82(4).

## Examples

``` r
g <- ilm_dag("dag {
  smoking [exposure]
  cancer  [outcome]
  smoking -> cancer
  genotype -> smoking
  genotype -> cancer
}")
ilm_adjust_sets(g)
#> [[1]]
#> [1] "genotype"
#> 
#> attr(,"exposure")
#> [1] "smoking"
#> attr(,"outcome")
#> [1] "cancer"
#> attr(,"pool")
#> [1] "genotype"
#> attr(,"unobserved")
#> character(0)
## and with the confounder unmeasured, nothing identifies the effect
ilm_adjust_sets(g, observed = c("smoking", "cancer"))
#> list()
#> attr(,"exposure")
#> [1] "smoking"
#> attr(,"outcome")
#> [1] "cancer"
#> attr(,"pool")
#> character(0)
#> attr(,"unobserved")
#> [1] "genotype"
```
