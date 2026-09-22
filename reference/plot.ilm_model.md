# Plot a fitted model

A shorthand for
[`ilm_plot_model()`](https://craig-hutton.github.io/illume/reference/ilm_plot_model.md),
which is the fuller interface.

## Usage

``` r
# S3 method for class 'ilm_model'
plot(x, what = "coef", ...)
```

## Arguments

- x:

  A fitted `"ilm_model"` object.

- what:

  Which plot to draw; see
  [`ilm_plot_model()`](https://craig-hutton.github.io/illume/reference/ilm_plot_model.md).

- ...:

  Passed to
  [`ilm_plot_model()`](https://craig-hutton.github.io/illume/reference/ilm_plot_model.md).

## Value

Invisibly, the data behind the plot.

## Examples

``` r
set.seed(1)
d <- ilm_sim()
f <- ilm_model(score ~ income + (1 | id), data = d, family = "gaussian",
               verbose = FALSE)
plot(f)
```
