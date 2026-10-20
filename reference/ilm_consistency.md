# Check whether the model can recover itself

Simulates datasets from the fitted model, refits each one, and asks
whether the re-estimates centre on the original values. A systematic
offset indicates that the Laplace approximation is biased for this model
and this amount of data.

## Usage

``` r
ilm_consistency(fit, B = 50L, seed = 1L, ncores = 1L, verbose = TRUE)
```

## Arguments

- fit:

  A fitted `"ilm_model"` object.

- B:

  Integer. Simulated datasets. Detectable bias scales roughly as
  `1/sqrt(B)`, so `B = 30` resolves offsets of about 0.37 standard
  deviations.

- seed:

  Integer. Random seed.

- ncores:

  Integer. Worker processes. Every replicate is simulated up front on
  the main process, so results do not depend on this.

- verbose:

  Logical. Print the table.

## Value

Invisibly, a list with the comparison table, the raw draws, the number
of usable refits and the refit rate.

## How to read the result

The **refit rate comes first and matters most**. If many replicates fail
to refit, the model cannot reliably recover itself from data it
generated, and the replicates that did converge are a
success-conditioned sample – their apparent lack of bias understates the
problem rather than excusing it. A low refit rate is a finding, not a
technical hiccup.

Parameters at the **zero boundary** are reported separately from bias. A
variance or loading fitted at essentially zero cannot recentre on
itself: re-estimates are bounded below by zero, so their mean must
exceed it. That is arithmetic, not approximation error, and the usual
symmetric test does not apply. Such a parameter is flagged `"BOUNDARY"`,
which tells you the component is unsupported by the data rather than
that the fitting method failed.

With many parameters tested, some will exceed a two-sigma threshold by
chance; the printed output states how many to expect.

## References

Joe, H. (2008). Accuracy of Laplace approximation for discrete response
mixed models. *Computational Statistics & Data Analysis*, 52(12),
5066–5074.

Self, S. G., & Liang, K.-Y. (1987). Asymptotic properties of maximum
likelihood estimators and likelihood ratio tests under nonstandard
conditions. *Journal of the American Statistical Association*, 82(398),
605–610. (On why boundary parameters need separate treatment.)
