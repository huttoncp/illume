# Fitted excess-zero probability, one value per row

Fitted excess-zero probability, one value per row

## Usage

``` r
ilm_zi_prob(object, newdata = NULL)
```

## Arguments

- object:

  A fitted
  [`ilm_model()`](https://craig-hutton.github.io/illume/reference/ilm_model.md)
  with a `ziformula`.

- newdata:

  Optional data to predict for; the model frame by default.

## Value

Numeric vector of probabilities between 0 and 1.

## See also

[`ilm_zi_coef()`](https://craig-hutton.github.io/illume/reference/ilm_zi_coef.md).

## Examples

``` r
set.seed(1); n <- 300
d <- data.frame(x = rnorm(n))
d$y <- ifelse(runif(n) < 0.3, 0, rpois(n, exp(0.5 + 0.3 * d$x)))
f <- ilm_model(y ~ x, data = d, family = "poisson", ziformula = ~ 1,
               verbose = FALSE)
summary(ilm_zi_prob(f))
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>  0.2808  0.2808  0.2808  0.2808  0.2808  0.2808 
```
