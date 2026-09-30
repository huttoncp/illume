# The rows a model used

How many rows the model was given, how many it used, how many rows with
missing values it dropped, and which columns those were missing in.

## Usage

``` r
ilm_rows_used(object, ...)

# S3 method for class 'ilm_model'
ilm_rows_used(object, ...)

# S3 method for class 'ilm_dag_model'
ilm_rows_used(object, ...)
```

## Arguments

- object:

  A fitted `"ilm_model"`, or an `"ilm_dag_model"`.

- ...:

  Unused.

## Value

For a model, a list of class `"ilm_rows_used"`: `n_input` (the rows of
`data`; `NA` when `data` was not a data frame), `n_used`, `n_dropped`,
`dropped_by` (a named integer vector: for each column the model reads,
the dropped rows missing in it) and `n_zero_weight`. For an
`"ilm_dag_model"`, a data frame with one row per adjustment set, and the
per-column counts as a list column.

## Details

[`ilm_model()`](https://huttoncp.github.io/illume/reference/ilm_model.md)
drops a row missing in any column the model reads (the default
`na.action = na.omit`). Nothing else drops a row: there is no `subset`
argument, and a row with a frequency weight of zero stays in the fit –
it counts as used, contributing nothing, and is counted as
`n_zero_weight`.

Which rows went is stats' own record, `stats::na.action(fit)`: their
positions in `data`, with the row names. The rows used are the others.

A row missing in two columns is dropped once and counted under each, so
`dropped_by` can sum to more than `n_dropped`.

An
[`ilm_dag_model()`](https://huttoncp.github.io/illume/reference/ilm_dag_model.md)
fits one model per adjustment set, and sets adjust for different
columns, so each can lose different rows: its method gives one row per
set.

## See also

[`stats::na.action()`](https://rdrr.io/r/stats/na.action.html),
[`ilm_model()`](https://huttoncp.github.io/illume/reference/ilm_model.md).

## Examples

``` r
set.seed(1)
d <- data.frame(x = rnorm(50), z = rnorm(50))
d$y <- 1 + d$x + rnorm(50)
d$x[1:4] <- NA; d$z[3:6] <- NA
fit <- ilm_model(y ~ x + z, data = d, family = "gaussian", verbose = FALSE)
ilm_rows_used(fit)
#> Rows: 44 of 50 analysed; 6 dropped for missing values (x 4, z 4; a row missing in several columns counts in each)
stats::na.action(fit)
#> 1 2 3 4 5 6 
#> 1 2 3 4 5 6 
#> attr(,"class")
#> [1] "omit"
```
