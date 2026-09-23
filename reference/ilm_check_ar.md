# Check for leftover temporal autocorrelation

Tests whether residuals from observations close together in time within
the same group are more similar than the model predicts. If they are,
the model is missing temporal structure.

## Usage

``` r
ilm_check_ar(
  object,
  time,
  group,
  maxlag = 8L,
  B = 30L,
  ncores = 1L,
  seed = 1L,
  verbose = TRUE,
  plot = FALSE,
  progress = NULL
)
```

## Arguments

- object:

  A fitted `"ilm_model"` object.

- time:

  Integer time index, one per observation. Must be evenly spaced.

- group:

  Grouping variable, one per observation.

- maxlag:

  Integer. Largest lag to examine.

- B:

  Integer. Simulated datasets. No p-value can fall below `1/(B+1)`, so
  `B` of about 100 is needed before the strongest verdict becomes
  reachable.

- ncores:

  Integer. Worker processes.

- seed:

  Integer. Random seed.

- verbose:

  Logical. Print the table.

- plot:

  Logical. Draw the autocorrelation and its band.

- progress:

  Show a progress bar. Defaults to
  [`interactive()`](https://rdrr.io/r/base/interactive.html), so a bar
  appears when someone is watching and nothing is written in a script or
  a knitted document. See
  [illumex::ilm_progress_arg](https://rdrr.io/pkg/illumex/man/ilm_progress_arg.html).

## Value

Invisibly, a list with the per-lag `table`, the matching `pacf` table,
the observed values and per-component z scores, the simulated null and
the number of replicates that refitted. The whole thing can be handed
straight to
[`ilm_plot_acf()`](https://huttoncp.github.io/illume/reference/ilm_plot_acf.md)
to avoid refitting.

## Which residual, and why it matters

This uses signed one-against-the-rest Pearson residuals rather than the
log-score quantile residuals used elsewhere. The reason is worth
understanding: the log score depends only on the probability assigned to
the category that actually occurred, so it registers the *size* of a
discrepancy but not its **direction**. Autocorrelation is directional –
consecutive observations being pushed the same way – so a
direction-blind residual cannot see it. Measured side by side, the
log-score residual had essentially no power here while the signed
Pearson residual separated the cases cleanly.

## What it cannot tell you

It detects autocorrelation the model does not **account for**, so it
fires on a model fitted without an AR term to data that have one. It
cannot tell you that an AR term you *did* fit is poorly estimated,
because simulating from that fit reproduces its own autocorrelation and
the discrepancy cancels out.

## The comparison band

Built by refitting simulated data rather than using the usual
`plus or minus 2/sqrt(n)` lines, which assume independent observations
and would over-flag here. Because the statistic takes a maximum across
residual components, its reference is calibrated by leave-one-out to
avoid an artificial bias toward significance.

## References

Rue, H., & Held, L. (2005). *Gaussian Markov Random Fields: Theory and
Applications*. Chapman & Hall/CRC. (Background on latent autoregressive
structures.)

## See also

[`ilm_plot_acf()`](https://huttoncp.github.io/illume/reference/ilm_plot_acf.md)
for the picture,
[`ilm_appraise()`](https://huttoncp.github.io/illume/reference/ilm_appraise.md),
[`ilm_check_omitted()`](https://huttoncp.github.io/illume/reference/ilm_check_omitted.md).

## Examples

``` r
set.seed(1)
d <- ilm_sim(n_id = 20, n_period = 8)
f <- ilm_model(score ~ income + (1 | id), data = d, family = "gaussian",
               verbose = FALSE)
ilm_check_ar(f, time = as.integer(factor(d$date)), group = d$id,
             maxlag = 3, B = 12)
#> 
#> residual autocorrelation vs simulated envelope (12 refits)
#>  lag component estimate null_mean      lo     hi  maxz      p status n_pairs
#>    1  residual  -0.1768   -0.1278 -0.3871 0.1315 -0.57 0.4615     OK     140
#>    2  residual  -0.1317   -0.0836 -0.2908 0.1236 -0.44 0.9231     OK     120
#>    3  residual  -0.0654   -0.1369 -0.3296 0.0558  0.83 0.5385     OK     100
#> 
#> >> residual autocorrelation is consistent with the fitted model.
```
