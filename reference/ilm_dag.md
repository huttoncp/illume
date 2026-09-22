# A causal graph

Holds a directed acyclic graph, which variables are unobserved, and
which are the exposure and the outcome. This is the object
[`ilm_adjust_sets()`](https://huttoncp.github.io/illume/reference/ilm_adjust_sets.md)
and
[`ilm_dag_test()`](https://huttoncp.github.io/illume/reference/ilm_dag_test.md)
work from.

## Usage

``` r
ilm_dag(x, exposure = NULL, outcome = NULL, latent = NULL)
```

## Arguments

- x:

  A dagitty-style string, a data frame of edges, or a `dagitty` object.

- exposure, outcome:

  Variable names. Optional if the string or object already declares
  them; required otherwise before an adjustment set can be found.

- latent:

  Names of unobserved variables, added to any the input declares.

## Value

An object of class `"ilm_dag"`.

## Details

Three input forms, all giving the same object:

- a dagitty-style string, with or without the `dag { }` wrapper;

- a data frame of edges with columns `from` and `to`;

- a `dagitty` object, if that package is installed.

A bidirected edge `X <-> Y` states that something unobserved causes
both, so it is stored as an unobserved node with an arrow into each.
That keeps one representation – a DAG, possibly with unobserved nodes –
rather than two.

## See also

[`ilm_adjust_sets()`](https://huttoncp.github.io/illume/reference/ilm_adjust_sets.md)
for what the graph licenses,
[`ilm_dag_test()`](https://huttoncp.github.io/illume/reference/ilm_dag_test.md)
for whether the data agree with it.

## Examples

``` r
g <- ilm_dag("dag {
  smoking [exposure]
  cancer  [outcome]
  smoking -> tar -> cancer
  genotype -> smoking
  genotype -> cancer
}")
g
#> <ilm_dag> 4 variables, 4 edges
#>   exposure: smoking
#>   outcome:  cancer
#>   edges:
#>     smoking -> tar
#>     tar -> cancer
#>     genotype -> smoking
#>     genotype -> cancer
```
