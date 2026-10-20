# Conditional independencies the graph implies

One testable claim per missing edge: two variables with no arrow between
them are independent given some set of the others. If the graph is right
the data should agree; where they do not, either the graph is wrong or
the data are not what they are taken to be.

## Usage

``` r
ilm_dag_implied(g, observed = NULL)
```

## Arguments

- g:

  An
  [`ilm_dag()`](https://huttoncp.github.io/illume/reference/ilm_dag.md).

- observed:

  Restrict to claims every variable of which was measured.

## Value

A data frame with `x`, `y` and a list column `given`.

## Details

The conditioning set is minimal – no member can be dropped without the
claim failing. Conditioning on the parents of both variables always
works and is the usual textbook basis, but it is often far larger than
it needs to be, and every unnecessary covariate costs power in the test
and invites a modelling error of its own. Where the two variables are
already independent with nothing held fixed, the set is empty.

## See also

[`ilm_dag_test()`](https://huttoncp.github.io/illume/reference/ilm_dag_test.md),
which tests these against data.

## Examples

``` r
g <- ilm_dag("dag { z -> x -> y ; z -> w }")
ilm_dag_implied(g)
#>   x y given
#> 1 z y     x
#> 2 x w     z
#> 3 y w     x
```
