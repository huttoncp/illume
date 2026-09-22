# Are two variables d-separated given a conditioning set?

The graphical criterion behind every claim a DAG makes about the data.
Tested by the moralisation route: take the ancestral subgraph of
everything mentioned, marry the parents, drop the directions, remove the
conditioning set, and ask whether a path remains.

## Usage

``` r
ilm_dsep(g, x, y, z = character())
```

## Arguments

- g:

  An
  [`ilm_dag()`](https://huttoncp.github.io/illume/reference/ilm_dag.md).

- x, y:

  Variable names.

- z:

  Conditioning set; may be empty.

## Value

`TRUE` if `x` and `y` are d-separated given `z`.

## Examples

``` r
g <- ilm_dag("dag { z -> x -> y ; z -> y }")
ilm_dsep(g, "x", "y", character())   # FALSE, x causes y
#> [1] FALSE
ilm_dsep(g, "x", "y", "z")           # still FALSE, the direct edge remains
#> [1] FALSE
```
