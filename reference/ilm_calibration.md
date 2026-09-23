# Calibration of predicted probabilities

Asks a direct question: among observations the model gave roughly a 30%
chance of category A, did about 30% turn out to be category A?
Predictions that pass this are said to be well calibrated, and for a
categorical outcome this is often more informative than any residual
plot.

## Usage

``` r
ilm_calibration(object, nbins = 10L, B = 200L, seed = 1L)
```

## Arguments

- object:

  A fitted `"ilm_model"` object.

- nbins:

  Integer. Number of bins of predicted probability.

- B:

  Integer. Simulated datasets for the band.

- seed:

  Integer. Random seed.

## Value

A list with one element per category, each holding mean predicted and
observed proportions per bin plus the simulated interval.

## Details

The comparison band is simulated from the model itself rather than taken
from a formula, so it reflects how much scatter is expected at this
sample size.

For a binomial, multinomial or ordinal fit – every family whose
prediction is a probability for each category. An ordinal fit is checked
category by category, from the probabilities its thresholds imply.
