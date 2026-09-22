# Combine several plots into one figure

Sets up a panel grid and evaluates each plot-producing expression into
it.

## Usage

``` r
ilm_plot_c(..., nrow = NULL, ncol = NULL)
```

## Arguments

- ...:

  Plot-producing expressions, for instance
  `ilm_plot_histogram(mtcars, "mpg")`.

- nrow, ncol:

  Panel grid. Default is as square as it goes.

## Value

`NULL`, invisibly.

## Details

It takes **expressions**, not plot objects, and that is base graphics
rather than a shortcut: a base plot draws immediately and leaves nothing
behind to combine, so the drawing has to happen inside the grid.

## See also

[`ilm_plot_var_all()`](https://huttoncp.github.io/illume/reference/ilm_plot_var_all.md).

## Examples

``` r
ilm_plot_c(
  ilm_plot_histogram(mtcars, "mpg"),
  ilm_plot_box(mtcars, "mpg", x = "cyl"),
  nrow = 1
)
```
