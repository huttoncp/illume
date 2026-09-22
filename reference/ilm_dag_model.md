# Fit the model a causal graph implies

Takes a DAG and a data frame and does what an analyst would: checks the
graph against the data, finds what has to be adjusted for, picks a
response distribution, looks for grouping structure, fits, runs the
diagnostics, and fixes the error structure if an assumption fails. It
narrates each step.

## Usage

``` r
ilm_dag_model(
  dag,
  data,
  exposure = NULL,
  outcome = NULL,
  family = NULL,
  cluster = NULL,
  auto_error = TRUE,
  test_dag = TRUE,
  max_sets = 8L,
  verbose = TRUE,
  reml = TRUE,
  ...
)
```

## Arguments

- dag:

  An
  [`ilm_dag()`](https://huttoncp.github.io/illume/reference/ilm_dag.md),
  or anything
  [`ilm_dag()`](https://huttoncp.github.io/illume/reference/ilm_dag.md)
  accepts.

- data:

  A data frame.

- exposure, outcome:

  Variable names; taken from the graph if declared.

- family:

  Response distribution. Inferred from the outcome when `NULL`, and the
  inference is stated rather than assumed.

- cluster:

  Grouping variable(s) for random intercepts. `NULL` looks for them
  among variables the graph does not mention; `character(0)` fits none.

- auto_error:

  Adjust the error structure when a diagnostic asks.

- test_dag:

  Test the graph's implied conditional independencies first.

- max_sets:

  Stop after this many adjustment sets.

- verbose:

  Narrate each step.

- reml:

  Estimate the variance components by restricted maximum likelihood.
  Defaults to `TRUE` here, unlike
  [`ilm_model()`](https://huttoncp.github.io/illume/reference/ilm_model.md),
  because the graph fixed the adjustment set before any data were looked
  at: the fixed effects are not being selected, so the one thing REML
  forbids – comparing likelihoods across different fixed structures –
  never arises, and its unbiased variance components are simply better.
  Ignored for non-gaussian responses, where restricted likelihood has no
  meaning.

- ...:

  Passed to
  [`ilm_model()`](https://huttoncp.github.io/illume/reference/ilm_model.md).

## Value

An object of class `"ilm_dag_model"`: the graph, the sets, the fits, an
`effects` table with one row per set, the DAG test, and `steps`.

## What is fixed and what is searched

**The graph fixes the mean structure.** The exposure and its adjustment
set are determined by the causal question, and nothing here searches
over them. That is the point of supplying a DAG: it commits to the
specification before the data are looked at, so the intervals reported
at the end mean what they say. Searching over which covariates to
include and then reporting the winner's p-values is post-selection
inference, and no amount of care elsewhere repairs it.

**Only the error structure is adjusted**, and only when a diagnostic
asks for it – a model for the dispersion when the spread is not
constant, a random effect for grouping the graph does not speak to.
These change what the standard errors are, not which effect is being
estimated. Every step taken is recorded in `$steps` and printed as it
happens.

## Several adjustment sets

If more than one minimal set is admissible, every one is fitted and
reported. They estimate the same quantity, so if the graph is right they
should agree; the spread across them is a sensitivity analysis that
comes free. A set that disagrees markedly is worth more attention than
any single point estimate.

## When the effect is not identified

If no set of measured variables suffices, that is the result, and this
function returns it rather than fitting something that cannot answer the
question. It is the most useful thing a DAG can tell you, and it can
only be said before the modelling, not after.

## See also

[`ilm_adjust_sets()`](https://huttoncp.github.io/illume/reference/ilm_adjust_sets.md),
[`ilm_dag_test()`](https://huttoncp.github.io/illume/reference/ilm_dag_test.md),
[`ilm_model()`](https://huttoncp.github.io/illume/reference/ilm_model.md).

## Examples

``` r
set.seed(1); n <- 300
z <- rnorm(n); x <- 0.5 * z + rnorm(n); y <- 0.4 * x + 0.6 * z + rnorm(n)
d <- data.frame(x = x, y = y, z = z)
g <- ilm_dag("dag { x [exposure] ; y [outcome] ; z -> x -> y ; z -> y }")
ilm_dag_model(g, d, verbose = FALSE)
#> <ilm_dag_model> x -> y 
#>   family: gaussian
#>   graph vs data: UNTESTED (0 claims tested) 
#> 
#>   effect of x by adjustment set:
#>     [1] z                            x                0.4415  ( 0.3225,  0.5605)  p = <1e-04
```
