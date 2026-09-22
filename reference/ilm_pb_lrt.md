# Parametric-bootstrap likelihood-ratio test

Tests one fixed-effect term by comparing the observed likelihood-ratio
statistic to its **simulated** null distribution, instead of to a
chi-square approximation.

## Usage

``` r
ilm_pb_lrt(
  object,
  term,
  B = 200L,
  ncores = 1L,
  seed = 1L,
  restarts = 1L,
  verbose = TRUE
)
```

## Arguments

- object:

  A fitted `"ilm_model"` object, fitted through the formula interface.

- term:

  Character name or integer index of the term to test.

- B:

  Integer. Bootstrap replicates. 200 or more is advisable; each costs
  two model fits.

- ncores:

  Integer. Worker processes.

- seed:

  Integer. Random seed.

- restarts:

  Integer. Optimiser restarts within replicates.

- verbose:

  Logical. Print the result.

## Value

Invisibly, a list with `LR`, `df`, `p_boot`, `p_chisq`, the simulated
null distribution, and calibration diagnostics.

## Why bother

The usual p-value for a likelihood-ratio test assumes the statistic
follows a chi-square distribution, which is a large-sample result.
Linear mixed models have corrections for this (Kenward-Roger,
Satterthwaite) but **no such correction exists for a multinomial GLMM**,
so the chi-square reference is all the asymptotic route offers.

How wrong can it be? On a design with 22 groups, the test that should
reject a true null 5% of the time rejected it 8.5% of the time. The
bootstrap replaces the assumed reference with the statistic's actual
behaviour under this model and this design.

## How it works

Fit the full model and the model without the term. Simulate many
datasets **from the reduced model**, because that is what the null
hypothesis says generated the data. Refit both models to each simulated
dataset to build the null distribution, then see where the observed
statistic falls.

## Reading the output

With `B` replicates no p-value can fall below `1/(B+1)`, so the smallest
reportable value at `B = 200` is about 0.005. The function says when a
result has hit that floor, rather than implying more precision than the
simulation can deliver. It also reports the implied type-I error of the
asymptotic test, which tells you directly whether the chi-square
shortcut would have been adequate for your data.

## References

Halekoh, U., & Hojsgaard, S. (2014). A Kenward-Roger approximation and
parametric bootstrap methods for tests in linear mixed models: the R
package pbkrtest. *Journal of Statistical Software*, 59(9), 1–30.

## See also

[`ilm_anova()`](https://huttoncp.github.io/illume/reference/ilm_anova.md).
