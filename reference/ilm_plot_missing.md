# Proportion missing, by variable

Proportion missing, by variable

## Usage

``` r
ilm_plot_missing(data, ...)
```

## Arguments

- data:

  A data frame.

- ...:

  Passed to
  [`graphics::barplot()`](https://rdrr.io/r/graphics/barplot.html).

## Value

Invisibly, the named vector of proportions.

## See also

[`ilm_describe_na_all()`](https://craig-hutton.github.io/illume/reference/ilm_describe_na_all.md)
for the same information as a table.

## Examples

``` r
ilm_plot_missing(ilm_sim())
```
