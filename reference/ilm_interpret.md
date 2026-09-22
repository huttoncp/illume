# Interpret a fitted model in words

Writes out what a fit says: each effect on the scale the response is
measured on, how strong the evidence is, what the diagnostics found, and
what that implies for taking the estimates at face value. Intended for
reporting to people who will not read a coefficient table.

## Usage

``` r
ilm_interpret(object, ...)

# S3 method for class 'ilm_model'
ilm_interpret(object, causal = NULL, ame = TRUE, digits = 3, ...)

# S3 method for class 'ilm_dag_model'
ilm_interpret(object, causal = NULL, ame = TRUE, digits = 3, ...)

# S3 method for class 'ilm_did'
ilm_interpret(object, causal = NULL, ame = FALSE, digits = 3, ...)

# S3 method for class 'ilm_rdd'
ilm_interpret(object, causal = NULL, ame = FALSE, digits = 3, ...)
```

## Arguments

- object:

  An
  [`ilm_model()`](https://huttoncp.github.io/illume/reference/ilm_model.md),
  [`ilm_dag_model()`](https://huttoncp.github.io/illume/reference/ilm_dag_model.md),
  [`ilm_did()`](https://huttoncp.github.io/illume/reference/ilm_did.md)
  or
  [`ilm_rdd()`](https://huttoncp.github.io/illume/reference/ilm_rdd.md).

- ...:

  Unused.

- causal:

  Force causal or associational language. `NULL` decides from the
  design, which is what you want.

- ame:

  Report average marginal effects on the response scale. Costs a
  delta-method calculation; worth it for anything with a link function.

- digits:

  Rounding.

## Value

An object of class `"ilm_interpretation"`: a list of sections, which
[`print()`](https://rdrr.io/r/base/print.html) renders as wrapped text.

## Causal language is licensed, not assumed

A regression coefficient is an association. This says "associated with"
unless the object carries a design that identifies an effect – an
[`ilm_dag_model()`](https://huttoncp.github.io/illume/reference/ilm_dag_model.md)
with a valid adjustment set, an
[`ilm_did()`](https://huttoncp.github.io/illume/reference/ilm_did.md),
an [`ilm_rdd()`](https://huttoncp.github.io/illume/reference/ilm_rdd.md)
– in which case it says so and names what licenses it. That is the
single most common error in reporting a model, and the one place an
automatic interpreter could do real damage, so it is handled explicitly.

## What it will not do

It says nothing about bias that a diagnostic did not measure. Every such
sentence traces to a check that was actually run, and names the remedy
that check names. Where no diagnostic has spoken, it is silent rather
than reassuring: a clean report from this function means the checks that
ran passed, not that the model is right.

The prose is templated. The same fit gives the same words every time.

## See also

[`ilm_ame()`](https://huttoncp.github.io/illume/reference/ilm_ame.md),
[`ilm_appraise()`](https://huttoncp.github.io/illume/reference/ilm_appraise.md).

## Examples

``` r
set.seed(1); n <- 300
d <- data.frame(x = rnorm(n), g = factor(sample(c("a", "b"), n, TRUE)))
d$y <- 0.5 * d$x + 0.8 * (d$g == "b") + rnorm(n)
fit <- ilm_model(y ~ x + g, data = d, family = "gaussian", verbose = FALSE)
ilm_interpret(fit, ame = FALSE)
#> INTERPRETATION
#> ==============
#> 
#> The model
#>   A gaussian model of y, fitted to 300 observations.
#> 
#> What it says
#>   x: very strong evidence that a higher x is associated with a higher value
#>   of y (estimate 0.450, 95% interval 0.321 to 0.579, p = <1e-04).
#> 
#>   gb: very strong evidence that being b rather than a is associated with a
#>   higher value of y (estimate 0.892, 95% interval 0.643 to 1.141, p =
#>   <1e-04).
#> 
#> What the checks found
#>   All 5 fitting checks passed: the optimiser converged, the gradient is at
#>   zero and the information matrix is usable. These say the fit is sound,
#>   not that the model is right -- for that, run ilm_appraise().
#> 
#> How far to trust it
#>   These are associations. Nothing here rules out a common cause of a
#>   predictor and the response, so an effect could differ in size or sign
#>   from what is reported. To say more you need a design that identifies one:
#>   ilm_dag_model() with an adjustment set, ilm_did(), ilm_rdd().
#> 
#>   This model has no random or smooth terms, so its t and F tests are exact
#>   rather than large-sample approximations.
#> 
```
