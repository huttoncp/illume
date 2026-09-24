# Diagnostic plots for a fitted model

Six panels, chosen for the family: randomised quantile residuals against
the normal; residuals against the fitted value (for a multinomial, the
predicted probability of the category observed); calibration of the
predicted probabilities with a simulated band for a binomial or
multinomial outcome, and scale-location otherwise; the observed response
against simulations from the fit (category frequencies for a
multinomial); a random-effects distance plot; and the model check
verdicts.

## Usage

``` r
ilm_appraise(object, nbins = 10L, B = 200L, seed = 1L, ...)
```

## Arguments

- object:

  A fitted `"ilm_model"` object.

- nbins:

  Integer. Bins for the calibration panel.

- B:

  Integer. Simulated datasets for the bands.

- seed:

  Integer. Random seed.

- ...:

  Unused.

## Value

Invisibly, the residuals, calibration data and random-effect distances.

## Details

The residuals are
[`ilm_rqr()`](https://huttoncp.github.io/illume/reference/ilm_rqr.md)'s,
which are uniform under a correct model in every family – including a
nominal outcome, where the usual residual plots have no clear meaning.

[`performance::check_model()`](https://easystats.github.io/performance/reference/check_model.html)
draws these panels too.

## See also

[`ilm_check_ar()`](https://huttoncp.github.io/illume/reference/ilm_check_ar.md)
for temporal correlation,
[`ilm_check_omitted()`](https://huttoncp.github.io/illume/reference/ilm_check_omitted.md)
for omitted variables.
