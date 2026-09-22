# Plot every column of a data frame

[`ilm_plot_var()`](https://craig-hutton.github.io/illume/reference/ilm_plot_var.md)
once per column, arranged in a grid.

## Usage

``` r
ilm_plot_var_all(
  data,
  var2 = NULL,
  by = NULL,
  cols = NULL,
  nrow = NULL,
  ncol = NULL,
  verbose = FALSE,
  ...
)
```

## Arguments

- data:

  A data frame.

- var2:

  Optional name of a second column.

- by:

  Optional grouping column.

- cols:

  Columns to use. A character vector of names, a regular expression, a
  predicate function such as `is.numeric`, or `NULL` for all of them –
  see
  [ilm_selection](https://craig-hutton.github.io/illume/reference/ilm_selection.md).

- nrow, ncol:

  Panel grid. Default is as square as it goes.

- verbose:

  Say which plot was chosen and why.

- ...:

  Passed to the underlying plot function.

## Value

`NULL`, invisibly.

## See also

[`ilm_plot_all()`](https://craig-hutton.github.io/illume/reference/ilm_plot_all.md),
[`ilm_plot_var_pairs()`](https://craig-hutton.github.io/illume/reference/ilm_plot_var_pairs.md).

## Examples

``` r
ilm_plot_var_all(mtcars, cols = c("mpg", "cyl", "wt", "gear"))
```
