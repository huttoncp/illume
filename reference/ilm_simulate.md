# Simulate new outcomes from a fitted model

Draws fresh random effects from their estimated distributions, forms the
linear predictor, and samples a new response for each observation from
the model's family – a category for a multinomial or ordinal outcome.

## Usage

``` r
ilm_simulate(fit, nsim = 1L, seed = NULL)
```

## Arguments

- fit:

  A fitted `"ilm_model"` object.

- nsim:

  Integer. Number of simulated datasets.

- seed:

  Integer. Random seed, for reproducibility.

## Value

An integer matrix with one column per simulated dataset, each holding
category codes.

## Why this is central

Four separate tools depend on it: population-averaged predictions,
parametric-bootstrap tests, the simulated envelopes in every residual
diagnostic, and
[`ilm_consistency()`](https://huttoncp.github.io/illume/reference/ilm_consistency.md).
Simulation is how this package builds reference distributions, rather
than relying on theoretical ones that may not hold.

[`TMB::checkConsistency()`](https://rdrr.io/pkg/TMB/man/checkConsistency.html)
cannot be used in its place. It relies on the objective function being
able to simulate its own data, but this model's likelihood contains no
distribution calls that TMB can invert – the response enters as a fixed
matrix and the random-effect densities are written out by hand – so
there is nothing in the computation graph for TMB to redraw.

## See also

[`ilm_consistency()`](https://huttoncp.github.io/illume/reference/ilm_consistency.md),
[`ilm_pb_lrt()`](https://huttoncp.github.io/illume/reference/ilm_pb_lrt.md),
[`ilm_appraise()`](https://huttoncp.github.io/illume/reference/ilm_appraise.md).
