# Diagnostic and summary plots for a fitted model

Diagnostic and summary plots for a fitted model

## Usage

``` r
ilm_plot_model(
  model,
  what = "coef",
  term = NULL,
  conf = 0.95,
  colour = "black",
  fill = "grey70",
  alpha = NULL,
  size = 1,
  main = NULL,
  ...
)
```

## Arguments

- model:

  A fitted `"ilm_model"` object.

- what:

  `"coef"` for fixed effects with intervals, `"effect"` for the partial
  effect of one predictor, `"random"` for a caterpillar plot of the
  estimated group effects, `"residual"` to hand off to
  [`ilm_appraise()`](https://huttoncp.github.io/illume/reference/ilm_appraise.md),
  or `"acf"` to hand off to
  [`ilm_plot_acf()`](https://huttoncp.github.io/illume/reference/ilm_plot_acf.md),
  which also needs `time` and `group`.

- term:

  Which predictor or grouping factor to show, for `"effect"` and
  `"random"`.

- conf:

  Confidence level for intervals.

- colour, fill, alpha, size:

  Appearance.

- main:

  Plot title.

- ...:

  Passed to the underlying plot.

## Value

Invisibly, the data behind the plot.

## See also

[`ilm_appraise()`](https://huttoncp.github.io/illume/reference/ilm_appraise.md)
for the full residual panel.

## Examples

``` r
set.seed(1)
d <- ilm_sim()
f <- ilm_model(score ~ income + grp + (1 | id), data = d,
               family = "gaussian", verbose = FALSE)
ilm_plot_model(f, "coef")

ilm_plot_model(f, "random")
```
