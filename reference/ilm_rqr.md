# Randomised quantile residuals for a nominal outcome

Turns categorical outcomes into residuals that should be uniformly
distributed when the model is right.

## Usage

``` r
ilm_rqr(object, conditional = TRUE, seed = 1L)
```

## Arguments

- object:

  A fitted `"ilm_model"` object.

- conditional:

  Logical. Evaluate at the fitted random effects.

- seed:

  Integer. Random seed; the transform uses randomisation.

## Value

A numeric vector on the unit interval.

## The obstacle, and the way around it

Quantile residuals need a cumulative distribution function, and a CDF
needs an ordering. Nominal categories have none – "red", "green", "blue"
cannot be put in order – which is why `DHARMa` cannot be applied
directly to this model.

The resolution is that you do not need an ordering of *categories*, only
of some number, and the model supplies one. For each observation the
**log score** `-log(p)` at the observed category is a real-valued
quantity whose distribution under the model is completely known: it
takes the value `-log(p_j)` with probability `p_j`. The standard
randomised transform then applies, and the fitted probabilities
themselves provide the ordering.

## Do not test these against a uniform distribution

Verified by simulation: when the probabilities are known exactly the
residuals are uniform, as intended. But **in sample** they are not,
because the random effects were fitted to the same data, which pulls the
observed log score down. A Kolmogorov-Smirnov test against uniform
therefore rejects a correctly specified model. Use
[`ilm_rqr_test()`](https://huttoncp.github.io/illume/reference/ilm_rqr_test.md),
which builds its reference by refitting simulated data and so accounts
for this automatically.

## References

Dunn, P. K., & Smyth, G. K. (1996). Randomized quantile residuals.
*Journal of Computational and Graphical Statistics*, 5(3), 236–244.

Hartig, F. (2022). DHARMa: Residual diagnostics for hierarchical
regression models. R package. (The simulation-based approach this
borrows from.)

## See also

[`ilm_rqr_test()`](https://huttoncp.github.io/illume/reference/ilm_rqr_test.md),
[`ilm_appraise()`](https://huttoncp.github.io/illume/reference/ilm_appraise.md).
