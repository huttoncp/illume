# Diagnostic plots for a fitted model

Six panels: a quantile-residual normal plot, residuals against fitted
probability, per-category calibration with a simulated band, observed
against simulated category frequencies, a random-effects distance plot,
and the model check verdicts.

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

These are built for a nominal categorical outcome rather than adapted
from tools designed for continuous responses, because the usual residual
plots have no clear meaning here.
[`performance::check_model()`](https://easystats.github.io/performance/reference/check_model.html)
routes to this function.

## See also

[`ilm_check_ar()`](https://huttoncp.github.io/illume/reference/ilm_check_ar.md)
for temporal correlation,
[`ilm_check_omitted()`](https://huttoncp.github.io/illume/reference/ilm_check_omitted.md)
for omitted variables.
