# Predicted survival curve from a fitted model

The probability of surviving past each time, for one or more covariate
patterns, with confidence intervals.

## Usage

``` r
ilm_survival(object, newdata = NULL, times = NULL, conf = 0.95)
```

## Arguments

- object:

  A fitted `"ilm_model"` with an accelerated failure time family.

- newdata:

  Covariate patterns, one row each. With none, the curve is drawn at the
  median of each numeric predictor and the commonest level of each
  factor.

- times:

  Times at which to evaluate. With none, a grid spanning the observed
  follow-up.

- conf:

  Confidence level.

## Value

A data frame with `row`, `time`, `surv`, `lower` and `upper`.

## Details

All three accelerated failure time families have a closed-form survivor
function in `z = (log t - eta) / scale`, so the curve and its interval
come straight from the fitted linear predictor. Intervals are computed
on the complementary log-log scale, `log(-log S)`, which is linear in
`eta` and so keeps the interval inside `(0, 1)` without truncation.

Uncertainty in the scale parameter is not included: the interval
reflects uncertainty in the coefficients only. It is therefore slightly
narrow, most noticeably in the tail beyond the last observed event.

## See also

[`ilm_plot_survival()`](https://huttoncp.github.io/illume/reference/ilm_plot_survival.md),
which draws it against the Kaplan-Meier estimate,
[`ilm_surv()`](https://huttoncp.github.io/illume/reference/ilm_surv.md).

## Examples

``` r
set.seed(1)
d <- data.frame(x = rnorm(200))
tt <- exp(1.5 + 0.8 * d$x + 0.7 * log(rexp(200)))
ct <- rexp(200, rate = 1 / (2 * median(tt)))
d$time <- pmin(tt, ct); d$event <- as.integer(tt <= ct)
f <- ilm_model(time ~ x, data = d, family = "weibull",
               censor = ilm_surv(d$time, d$event), verbose = FALSE)
head(ilm_survival(f, newdata = data.frame(x = c(-1, 1))))
#>   row       time      surv     lower     upper
#> 1   1 0.02367112 0.9969355 0.9959874 0.9976598
#> 2   1 0.23229399 0.9385701 0.9203042 0.9527572
#> 3   1 0.44091686 0.8621862 0.8234503 0.8929776
#> 4   1 0.64953972 0.7804836 0.7227639 0.8276276
#> 5   1 0.85816259 0.6986832 0.6251829 0.7605540
#> 6   1 1.06678545 0.6197232 0.5342923 0.6940207
```
