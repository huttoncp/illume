# Does one coefficient per predictor describe every cut?

A cumulative link model says a predictor shifts all `J - 1` cut points
by the same amount. Under a logit link that is the proportional-odds
assumption, and it is what makes the model economical – and wrong, when
it does not hold, in a way that is invisible in the coefficient table.

## Usage

``` r
ilm_check_proportional(
  object,
  B = 199L,
  alpha = 0.05,
  seed = 1L,
  progress = NULL
)
```

## Arguments

- object:

  A fitted
  [`ilm_model()`](https://huttoncp.github.io/illume/reference/ilm_model.md)
  with an ordinal family.

- B:

  Simulated datasets for the reference. A simulated p-value cannot fall
  below `1 / (B + 1)`.

- alpha:

  Level for the per-term verdicts.

- seed:

  Random seed.

- progress:

  Show a progress bar; see
  [ilm_progress_arg](https://huttoncp.github.io/illume/reference/ilm_progress_arg.md).

## Value

Invisibly, a data frame with one row per fixed-effect term giving the
spread statistic, a simulated p-value, the largest gap between cuts and
a `status`, with the per-cut coefficients in a `coefs` attribute.

## Details

This is the Brant idea, done by refitting rather than by a score test.
Each cut becomes a binary outcome (`y > j`), a separate model is fitted
to each, and the statistic is how far the resulting coefficients spread
apart relative to their standard errors.

## Why the reference is simulated

Those separate fits share their data, so their estimates are correlated,
and the spread statistic is not chi-squared however tempting the shape.
Treating it as one was measured here at a **31.5% false-alarm rate** for
a binary predictor whose effect was perfectly proportional – with two
cuts the statistic is `(b1 - b2)^2` over an average variance, where the
honest denominator is `v1 + v2 - 2c`, and the covariance `c` that is
missing depends on the predictor.

The reference is therefore simulated from the fitted proportional-odds
model, which carries that correlation because it arises from the same
refitting. Over 400 replicates per design, with the assumption holding:

                       flagged at 0.05        p-values under the null
      n = 400   x 0.032, g 0.040          p<.05  0.042 / 0.056
      n = 1500  x 0.052, g 0.072          p<.25  0.255 / 0.244
                                          p<.50  0.510 / 0.484

and with `x` acting on the first cut only, the offending term is flagged
in 400 of 400 samples at both sizes, while the term that does act
proportionally is flagged 1.5% and 2.8% of the time – it localises the
violation rather than condemning the model. The 7.2% at the largest size
is a little above nominal and is reported here rather than rounded away.

## What to do when it fails

The remedy is `family = "multinomial"`, which spends `J - 1`
coefficients per predictor and assumes nothing about the ordering. It
nests this model, so [`AIC()`](https://rdrr.io/r/stats/AIC.html) will
say whether the extra parameters earn their place. When only one or two
terms offend, that is the trade to weigh: the ordering is real
information, and giving all of it up to accommodate one predictor may
cost more than it buys.

## References

Brant, R. (1990). Assessing proportionality in the proportional odds
model for ordinal logistic regression. *Biometrics* 46, 1171-1178.

## See also

[`ilm_thresholds()`](https://huttoncp.github.io/illume/reference/ilm_thresholds.md),
[`ilm_model()`](https://huttoncp.github.io/illume/reference/ilm_model.md).

## Examples

``` r
set.seed(1); n <- 500
d <- data.frame(x = rnorm(n))
z <- 0.8 * d$x + rlogis(n)
d$y <- factor(cut(z, c(-Inf, -1, 1, Inf), labels = c("lo", "mid", "hi")),
              ordered = TRUE)
f <- ilm_model(y ~ x, data = d, family = "ordinal", verbose = FALSE)
ilm_check_proportional(f, B = 99)
#> Warning: B = 99 puts the smallest achievable p-value at 0.01, so a FAIL verdict is unreachable (it needs p < 0.01) and even a flagrant violation can only read WARN. Use B >= 100.
#> Proportional-odds check
#>   reference simulated from the fitted model, B = 99
#> 
#>  term stat    p max_gap status
#>     x 0.02 0.79   0.025     OK
```
