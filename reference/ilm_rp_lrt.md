# Was the flexible baseline worth it?

Compares a flexible parametric fit against the straight-line special
case it generalises – Weibull on the hazard scale, log-logistic on the
odds scale, lognormal on the normal scale – by refitting with
`rp_df = 1`.

## Usage

``` r
ilm_rp_lrt(object, verbose = TRUE)
```

## Arguments

- object:

  A fitted `"ilm_model"` with a flexible parametric baseline.

- verbose:

  Print the verdict.

## Value

Invisibly, a list with both log-likelihoods, the statistic, its degrees
of freedom, a p-value and a `status`.

## Details

The comparison is a likelihood-ratio test on `rp_df - 1` degrees of
freedom. It is a test of the baseline SHAPE, not of any covariate: a
large statistic says the straight line was wrong, which is the thing the
flexible model exists to fix. Because the knots are placed from the data
rather than fixed in advance, treat the p-value as indicative.

## See also

[`ilm_plot_survival()`](https://huttoncp.github.io/illume/reference/ilm_plot_survival.md),
which checks the same thing against the Kaplan-Meier estimate rather
than against the special case.

## Examples

``` r
set.seed(1)
d <- data.frame(x = rnorm(300))
d$time <- exp(1.5 + 0.8 * d$x + 0.6 * log(rexp(300)))
d$event <- 1L
f <- ilm_model(time ~ x, data = d, family = "rp", rp_df = 3,
               censor = ilm_surv(d$time, d$event), verbose = FALSE)
ilm_rp_lrt(f)
#> 
#> flexible baseline against the weibull special case
#>   logLik -679.4763 vs -680.3160, chi-square 1.68 on 2 df, p = 0.432
#>   a straight line fits as well: family = "weibull" is simpler and says the same
```
