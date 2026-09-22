# Regression discontinuity

Estimates the jump in an outcome at the point where treatment switches
on, by fitting a straight line on each side of the cutoff within a
bandwidth and taking the difference at the cutoff.

## Usage

``` r
ilm_rdd(
  data,
  y,
  running,
  cutoff = 0,
  treatment = NULL,
  h = NULL,
  kernel = "triangular",
  poly = 1L,
  covariates = NULL,
  family = NULL,
  bw_range = c(0.5, 0.75, 1, 1.5, 2),
  placebo = NULL,
  verbose = TRUE,
  ...
)
```

## Arguments

- data:

  A data frame.

- y:

  Outcome column.

- running:

  Running (forcing) variable.

- cutoff:

  Value of `running` at which treatment switches on.

- treatment:

  Optional column of actual treatment received, used to check whether
  assignment is sharp.

- h:

  Bandwidth. `NULL` uses a rule-of-thumb starting point.

- kernel:

  `"triangular"` (default), `"uniform"` or `"epanechnikov"`.

- poly:

  Polynomial order on each side; 1 is local linear.

- covariates:

  Further columns for the mean structure.

- family:

  Response distribution; inferred from `y` when `NULL`.

- bw_range:

  Multipliers of `h` for the sensitivity sweep.

- placebo:

  Placebo cutoffs, as quantiles of the running variable either side of
  the real cutoff. `NULL` picks a few.

- verbose:

  Narrate each step.

- ...:

  Passed to
  [`ilm_model()`](https://huttoncp.github.io/illume/reference/ilm_model.md).

## Value

An object of class `"ilm_rdd"`: `jump`, `fit`, `bandwidth`, `density`,
`balance`, `placebo` and the settings used.

## What this reports, and what it does not

A local linear fit with a triangular kernel, the bandwidth sensitivity
laid out, and the design checks below. It does **not** implement the
bias-corrected robust intervals of Calonico, Cattaneo and Titiunik; the
`rdrobust` package does that and does it well. The interval here is the
ordinary one for a weighted local fit, which understates uncertainty
when the bandwidth is large enough for curvature to matter – which is
exactly what `$bandwidth` is there to show you.

The default is local **linear**. A high-order global polynomial produces
estimates driven by points far from the cutoff and is a well-documented
way to manufacture a discontinuity (Gelman and Imbens 2019). `poly` is
available and going above 2 warns.

## The design checks

A regression discontinuity is only as good as the claim that crossing
the cutoff is the *only* thing that changes there.

- **Density** – if people can move themselves across the cutoff, the
  ones just above are not comparable to the ones just below. The test is
  whether the density is *discontinuous* at the cutoff, not whether it
  is symmetric about it: a local linear density is fitted on each side
  and compared there. Over 500 unmanipulated data sets it flagged 0.046
  to 0.060 of them across uniform, normal and exponential running
  variables, and caught 0.984 of the cases where 30% of the units just
  below the cutoff had been moved above it.

- **Covariate balance** – anything measured before treatment should not
  jump at the cutoff. If it does, something other than treatment changes
  there.

- **Placebo cutoffs** – the same estimate at points where nothing
  happened. These should be null; a comparable jump elsewhere means the
  method is finding jumps in noise.

- **Bandwidth sensitivity** – the estimate across a range of bandwidths.
  An effect that appears only in a narrow window is not an effect.

## Fuzzy assignment

This estimates a **sharp** design, where crossing the cutoff determines
treatment. If a treatment column is supplied and assignment turns out
fuzzy – the probability of treatment jumps but not from 0 to 1 – this
says so and declines to report the sharp estimate as if it were the
effect of treatment. A fuzzy design needs an instrumental-variables
estimator, which illume does not have.

## References

Gelman, A. and Imbens, G. (2019). Why high-order polynomials should not
be used in regression discontinuity designs. Journal of Business and
Economic Statistics 37(3).

Calonico, S., Cattaneo, M. D. and Titiunik, R. (2014). Robust
nonparametric confidence intervals for regression-discontinuity designs.
Econometrica 82(6).

## See also

[`ilm_plot_rdd()`](https://huttoncp.github.io/illume/reference/ilm_plot_rdd.md),
[`ilm_did()`](https://huttoncp.github.io/illume/reference/ilm_did.md),
[`ilm_model()`](https://huttoncp.github.io/illume/reference/ilm_model.md).

## Examples

``` r
set.seed(1); n <- 2000
r <- runif(n, -1, 1)
y <- 0.5 * r + 0.8 * (r >= 0) + rnorm(n, 0, 0.5)
ilm_rdd(data.frame(r = r, y = y), "y", "r", cutoff = 0, verbose = FALSE)
#> <ilm_rdd> y at r = 0 
#>   local linear fit, triangular kernel, h = 0.2353 (258 below, 202 above)
#>   family: gaussian 
#> 
#>   jump   0.6949  ( 0.4685,  0.9213)  p = <1e-04
#>   across bandwidths 0.5x-2x: 0.6711 to 0.7936
#> 
#>   design checks
#>     density at the cutoff: OK (258 below, 202 above)
#>     placebo cutoffs:      0 of 4 show a jump
#> 
#>   Interval is the ordinary one for a weighted local fit. For
#>   bias-corrected robust intervals see the rdrobust package.
```
