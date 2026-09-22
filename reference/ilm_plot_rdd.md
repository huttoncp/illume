# Regression discontinuity plot

Binned means of the outcome against the running variable, with the
fitted line either side of the cutoff. The canonical figure: it shows
the jump, and it shows whether the jump is the only thing going on.

## Usage

``` r
ilm_plot_rdd(
  x,
  nbins = 20L,
  main = "Regression discontinuity",
  xlab = NULL,
  ylab = NULL,
  ...
)
```

## Arguments

- x:

  An
  [`ilm_rdd()`](https://craig-hutton.github.io/illume/reference/ilm_rdd.md)
  result.

- nbins:

  Bins on each side.

- main, xlab, ylab:

  Labels.

- ...:

  Passed to
  [`graphics::plot()`](https://rdrr.io/r/graphics/plot.default.html).

## Value

`x`, invisibly.

## See also

[`ilm_rdd()`](https://craig-hutton.github.io/illume/reference/ilm_rdd.md).

## Examples

``` r
set.seed(1); n <- 2000
r <- runif(n, -1, 1)
y <- 0.5 * r + 0.8 * (r >= 0) + rnorm(n, 0, 0.5)
fit <- ilm_rdd(data.frame(r = r, y = y), "y", "r", cutoff = 0,
               verbose = FALSE)
ilm_plot_rdd(fit)
```
