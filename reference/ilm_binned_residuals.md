# Binned residuals for a binary response

A raw residual from a 0/1 outcome takes one of two values and tells you
almost nothing. Averaging residuals within bins of fitted probability
makes the pattern visible: points outside the band mark regions of the
predictor space where the model is systematically wrong.

## Usage

``` r
ilm_binned_residuals(object, nbins = NULL, plot = TRUE)
```

## Arguments

- object:

  A fitted `"ilm_model"` object with a binomial family.

- nbins:

  Number of bins. The default follows the usual `sqrt(n)` rule.

- plot:

  Draw the plot.

## Value

Invisibly, a data frame of bin midpoints, mean residuals, bounds and a
`status`.

## References

Gelman, A. and Hill, J. (2007). Data Analysis Using Regression and
Multilevel/Hierarchical Models. Cambridge University Press, ch. 5.

## Examples

``` r
set.seed(1)
d <- ilm_sim()
f <- ilm_model(flag ~ income + (1 | id), data = d, family = "binomial",
               verbose = FALSE)
ilm_binned_residuals(f)
```
