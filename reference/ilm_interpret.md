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

# S3 method for class 'ilm_power'
ilm_interpret(
  object,
  causal = NULL,
  ame = FALSE,
  digits = 2,
  target = 0.8,
  ...
)

# S3 method for class 'ilm_contrast'
ilm_interpret(object, causal = NULL, ame = FALSE, digits = 3, ...)

# S3 method for class 'ilm_profile'
ilm_interpret(object, causal = NULL, ame = FALSE, digits = 3, ...)
```

## Arguments

- object:

  An
  [`ilm_model()`](https://huttoncp.github.io/illume/reference/ilm_model.md),
  [`ilm_dag_model()`](https://huttoncp.github.io/illume/reference/ilm_dag_model.md),
  [`ilm_did()`](https://huttoncp.github.io/illume/reference/ilm_did.md)
  or
  [`ilm_rdd()`](https://huttoncp.github.io/illume/reference/ilm_rdd.md);
  or the result of
  [`ilm_power()`](https://huttoncp.github.io/illume/reference/ilm_power.md)
  or
  [`ilm_power_design()`](https://huttoncp.github.io/illume/reference/ilm_power_design.md),
  which is written up as a power analysis, of
  [`ilm_contrast()`](https://huttoncp.github.io/illume/reference/ilm_contrast.md),
  written up as the comparisons it makes, or of
  [`ilm_profile()`](https://huttoncp.github.io/illumex/reference/ilm_profile.html)
  or
  [`ilm_profile_na()`](https://huttoncp.github.io/illumex/reference/ilm_profile_na.html),
  written up as the profile of each cluster.

- ...:

  Unused.

- causal:

  Force causal or associational language. `NULL` decides from the
  design, which is what you want.

- ame:

  Give the effects intervals on the response scale, by the delta method.
  With `FALSE` the predictions are still reported, without them, which
  is quicker.

- digits:

  Rounding.

- target:

  For a power analysis, the power a study is to reach.

## Value

An object of class `"ilm_interpretation"`: a list of sections, which
[`print()`](https://rdrr.io/r/base/print.html) renders as wrapped text.

## How an effect is said

In the response's own units, over a change in the predictor a reader can
picture: across the middle half of a number (its quartiles), level by
level for a factor. "Across the middle half of age, predicted income is
49.4 at 29 against 64.7 at 54: 15.2 higher (95% interval 12.3 to 18.1)",
where the coefficient alone, 0.61, says nothing until the reader knows
what a unit of age is and how much age varies. The predictions are
averaged over the rows the model was fitted to, as
[`ilm_ame()`](https://huttoncp.github.io/illume/reference/ilm_ame.md)'s
effects are, with any random effects at zero: in a mixed model they are
a typical group's, and with a nonlinear link the text says so, since
averaging over the groups instead – `predict(groups = "population")` –
gives different numbers. A probability is given as one, and a difference
of two in percentage points, with the odds ratio after it for those who
want it. The evidence is the term's joint test, so a factor with several
levels, or a predictor in a multinomial model, gets one verdict rather
than one per coefficient. A term that is not a plain variable – an
interaction, a spline – is described by its coefficients.

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
#>   x: very strong evidence (p < 0.001) that x is associated with y. Across
#>   the middle half of x, predicted y is 0.0787 at -0.59 against 0.644 at
#>   0.667: 0.565 higher. That is 0.45 per unit of x.
#> 
#>   g: very strong evidence (p < 0.001) that g is associated with y.
#>   Predicted y is -0.0988 for 'a' and 0.793 for 'b': 0.892 higher for 'b'.
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
