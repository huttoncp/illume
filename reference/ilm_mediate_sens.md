# How strong would unmeasured confounding have to be?

The one assumption a mediation analysis cannot check is that nothing
unmeasured causes both the mediator and the outcome. This does not test
it, because nothing can. It asks the answerable question instead: given
an unmeasured confounder of stated strength, what would the indirect
effect have been, and how strong would it have to be for that effect to
vanish.

## Usage

``` r
ilm_mediate_sens(
  object,
  lambda_m = seq(0, 1.5, length.out = 16L),
  lambda_y = seq(0, 1.5, length.out = 16L)
)
```

## Arguments

- object:

  An
  [`ilm_mediate()`](https://craig-hutton.github.io/illume/reference/ilm_mediate.md)
  result on linear models.

- lambda_m, lambda_y:

  Grids of confounder strengths.

## Value

An object of class `"ilm_mediate_sens"`: the grid of corrected indirect
effects, with the zero-crossing threshold attached.

## Reading the result

`lambda_m` and `lambda_y` are the confounder's coefficients in the
mediator and outcome models, with the confounder standardised, so they
are on the scale of those models' own coefficients: a `lambda_y` of 0.5
means an unmeasured variable with the same pull on the outcome as a
predictor whose coefficient is 0.5.

The `threshold` is the product `lambda_m * lambda_y` at which the
corrected indirect effect reaches zero. Compare it against the products
of the coefficients you DID measure: if a confounder half as strong as
your strongest covariate would overturn the finding, the finding is
fragile, and if it would take one twice as strong as anything observed,
it is not.

## References

VanderWeele, T. J. (2015). *Explanation in Causal Inference*. Oxford
University Press, chapter 3.

## See also

[`ilm_mediate()`](https://craig-hutton.github.io/illume/reference/ilm_mediate.md).

## Examples

``` r
set.seed(1); n <- 400
d <- data.frame(x = rbinom(n, 1, 0.5))
d$m <- 0.3 + 0.7 * d$x + rnorm(n)
d$y <- 1 + 0.4 * d$x + 0.6 * d$m + rnorm(n)
fm <- ilm_model(m ~ x, data = d, family = "gaussian", verbose = FALSE)
fy <- ilm_model(y ~ x + m, data = d, family = "gaussian", verbose = FALSE)
md <- ilm_mediate(fm, fy, treat = "x", mediator = "m", sims = 200)
ilm_mediate_sens(md)
#> Sensitivity of the indirect effect to unmeasured mediator-outcome confounding
#>   observed ACME: 0.4873
#>   it reaches zero when lambda_m * lambda_y = 0.6626
#>   so a confounder equally strong in both models would need
#>   lambda = 0.814
#> 
#>   For scale, the outcome model's own coefficients run from 0.222 to 0.608.
#>   If a confounder that strong is plausible, the indirect effect is not
#>   established; if it is larger than anything you measured, it is.
```
