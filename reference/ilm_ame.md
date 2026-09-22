# Average marginal effect on the response scale

What a coefficient means where the response lives. For a gaussian model
that is the coefficient itself; for anything with a link it is not, and
the difference matters: an odds ratio of 2 can be a 3-point change in
probability or a 20-point one depending on where the data sit.

## Usage

``` r
ilm_ame(object, terms = NULL, eps = 1e-04)
```

## Arguments

- object:

  An
  [`ilm_model()`](https://craig-hutton.github.io/illume/reference/ilm_model.md).

- terms:

  Which predictors. Default is every fixed-effect term.

- eps:

  Relative step for the numerical derivatives.

## Value

A data frame with `term`, `level`, `estimate`, `se`, `lower`, `upper`,
and `kind` (`"slope"` or `"contrast"`).

## Details

For a numeric predictor this is the average of the derivative of the
fitted mean with respect to it, taken over the rows actually observed.
For a factor it is the average change in fitted mean from moving every
row to that level from the reference, which is a contrast rather than a
derivative.

Standard errors come from the delta method: the effect is differentiated
numerically with respect to the parameters and combined with their
covariance. The parameter vector includes the covariance parameters,
because a population-averaged prediction depends on them and pretending
otherwise would understate the uncertainty.

## See also

[`ilm_interpret()`](https://craig-hutton.github.io/illume/reference/ilm_interpret.md),
and the `marginaleffects` package, which does this and a great deal more
once
[`ilm_register_marginaleffects()`](https://craig-hutton.github.io/illume/reference/ilm_register_marginaleffects.md)
is called.

## Examples

``` r
set.seed(1); n <- 400
d <- data.frame(x = rnorm(n), g = factor(sample(c("a", "b"), n, TRUE)))
d$y <- rbinom(n, 1, plogis(0.4 * d$x + 0.6 * (d$g == "b")))
fit <- ilm_model(y ~ x + g, data = d, family = "binomial", verbose = FALSE)
#> Warning: the ‘findbars’ function has moved to the reformulas package. Please update your imports, or ask an upstream package maintainer to do so.
#> Warning: the ‘nobars’ function has moved to the reformulas package. Please update your imports, or ask an upstream package maintainer to do so.
ilm_ame(fit)
#>   term level     kind   estimate         se      lower     upper
#> 1    x  <NA>    slope 0.09084394 0.02377348 0.04424877 0.1374391
#> 2    g     b contrast 0.18714683 0.04793616 0.09319367 0.2811000
```
