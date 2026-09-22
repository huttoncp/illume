# Plot one or two variables, choosing the geometry

Picks from the columns' types: a continuous variable alone gets a
density, a categorical one a bar chart, two continuous ones a scatter,
one of each a boxplot, two categorical ones a grouped bar chart.

## Usage

``` r
ilm_plot_var(data, var1, var2 = NULL, by = NULL, verbose = FALSE, ...)
```

## Arguments

- data:

  A data frame.

- var1:

  Name of the primary column.

- var2:

  Optional name of a second column.

- by:

  Optional grouping column.

- verbose:

  Say which plot was chosen and why.

- ...:

  Passed to the underlying plot function.

## Value

`NULL`, invisibly.

## Details

[`ilm_plot()`](https://craig-hutton.github.io/illume/reference/ilm_plot.md)
does the same job through a `geom = "auto"` argument. This is the same
choice made by naming variables rather than a geometry, and it says
which it picked when asked.

## See also

[`ilm_plot()`](https://craig-hutton.github.io/illume/reference/ilm_plot.md),
[`ilm_plot_var_all()`](https://craig-hutton.github.io/illume/reference/ilm_plot_var_all.md).

## Examples

``` r
ilm_plot_var(mtcars, "mpg")

ilm_plot_var(mtcars, "mpg", "cyl", verbose = TRUE)
#> both continuous -- scatter
```
