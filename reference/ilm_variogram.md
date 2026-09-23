# Residual correlation as a function of separation in time

Bins every pair of observations within a group by how far apart they
are, and reports the residual correlation in each bin against an
envelope built by simulating from the fitted model and refitting. This
is the check for a
[`ilm_car1()`](https://huttoncp.github.io/illume/reference/ilm_car1.md)
term, and the one to use whenever observation times are irregular.

## Usage

``` r
ilm_variogram(
  object,
  time,
  group,
  coords = NULL,
  breaks = 8L,
  B = 100L,
  type = "correlation",
  min_effect = 0.1,
  ncores = 1L,
  seed = 1L,
  max_pairs = 2e+05,
  plot = TRUE,
  verbose = TRUE,
  progress = NULL
)
```

## Arguments

- object:

  A fitted `"ilm_model"` object.

- time:

  Time, one value per observation. Give this or `coords`.

- group:

  Grouping variable, one value per observation. Pairs are only formed
  within a group. Optional for a spatial variogram, where the default
  treats every observation as comparable with every other.

- coords:

  Spatial coordinates, one to three columns, for a variogram over
  distance rather than over time.

- breaks:

  Number of separation bins, or explicit break points. Bins are chosen
  to hold roughly equal numbers of pairs, so a bin far out is not
  estimated from a handful of them.

- B:

  Simulated datasets behind the envelope. A p-value cannot fall below
  `1/(B+1)`, so about 100 is needed before `FAIL` is reachable.

- type:

  `"correlation"` or `"semivariance"`.

- min_effect:

  Smallest departure from the simulated null worth a verdict. A bin
  outside the envelope by less than this is reported as `OK`. Each bin
  rests on thousands of pairs, so the envelope narrows until any
  imperfection clears it, and significance stops being the same thing as
  something to act on. Set to `0` to flag on the envelope alone.

- ncores:

  Worker processes for the refits.

- seed:

  Random seed.

- max_pairs:

  Cap on the number of within-group pairs used. A group of `m`
  observations contributes `m(m-1)/2` of them, so a few long series can
  run to millions; beyond the cap a random subset is taken, which
  changes the precision but not what is being estimated.

- plot:

  Draw the variogram.

- verbose:

  Print the table.

- progress:

  Show a progress bar. Defaults to
  [`interactive()`](https://rdrr.io/r/base/interactive.html), so a bar
  appears when someone is watching and nothing is written in a script or
  a knitted document. See
  [ilm_progress_arg](https://huttoncp.github.io/illume/reference/ilm_progress_arg.md).

## Value

Invisibly, a list with the per-bin `table`, the simulated null, and the
number of replicates that refitted.

## Why not [`ilm_check_ar()`](https://huttoncp.github.io/illume/reference/ilm_check_ar.md)

That function matches pairs at exact lags, which is right for evenly
spaced data and close to useless without it. On a 60-unit panel with six
observations each, drawn from thirty possible times, exact matching
found 53 pairs at lag 1, five at lag 2 and none beyond – five of six
lags returned no verdict. Binning by distance uses all 900 pairs.

## Reading it

Correlation that starts high and decays toward zero as separation grows
is what an autoregressive process looks like, and is what
[`ilm_car1()`](https://huttoncp.github.io/illume/reference/ilm_car1.md)
fits. Correlation that is flat and positive at every separation is a
group effect the random intercept has not absorbed. Correlation that
rises again at some separation is a cycle; see
[`ilm_fourier()`](https://huttoncp.github.io/illume/reference/ilm_fourier.md).

The envelope is not centred on zero, and should not be. Residuals within
a group sum to roughly zero, so pairs of them correlate negatively even
when the model is exactly right, and the further apart they are the more
strongly this bites.

## Relation to the semivariogram

For standardised residuals the semivariance at separation `d` is
`1 - r(d)`, so `type = "semivariance"` plots the same information the
other way up, in the form
[`nlme::Variogram()`](https://rdrr.io/pkg/nlme/man/Variogram.html) uses.

## References

Pinheiro, J. C., & Bates, D. M. (2000). *Mixed-Effects Models in S and
S-PLUS*. Springer. (Chapter 5 covers the residual variogram.)

## See also

[`ilm_car1()`](https://huttoncp.github.io/illume/reference/ilm_car1.md)
for the remedy,
[`ilm_check_ar()`](https://huttoncp.github.io/illume/reference/ilm_check_ar.md)
for evenly spaced data,
[`ilm_plot_acf()`](https://huttoncp.github.io/illume/reference/ilm_plot_acf.md).

## Examples

``` r
set.seed(1)
d <- data.frame(id = factor(rep(1:20, each = 5)),
                t = as.vector(replicate(20, sort(sample(1:20, 5)))),
                x = rnorm(100))
d$y <- d$x + rnorm(100)
f <- ilm_model(y ~ x + (1 | id), data = d, family = "gaussian",
               verbose = FALSE)
#> ilm_model(): the random-effect covariance of `id` sits at the edge of its range -- a variance of zero or a correlation of +/-1 -- where the data cannot resolve it. The fixed effects and their standard errors are still usable; summary() says what else is. If the term belongs in the model, boundary = "avoid" keeps it inside its range with a small penalty: it is then assumed nonzero rather than estimated at zero, so do not test whether it is; its variance comes out larger, and for a binary or categorical outcome the fixed effects a little further from zero -- markedly so when a category is rare.
ilm_variogram(f, d$t, d$id, breaks = 4, B = 15, plot = FALSE)
#> Warning: B = 15 puts the smallest achievable p-value at 0.062, so a FAIL verdict is unreachable. Use B >= 100.
#> 
#> residual correlation by separation (15 refits, 200 pairs)
#>  bin separation component estimate null_mean      lo     hi  maxz      p status
#>    1      2.096  residual  -0.2914   -0.0394 -0.2110 0.1322 -2.89 0.0625     OK
#>    2      4.938  residual  -0.1124   -0.0152 -0.2037 0.1734 -0.94 0.5000     OK
#>    3      8.029  residual   0.1819   -0.0389 -0.4175 0.3396  1.49 0.1875     OK
#>    4     12.200  residual  -0.1688   -0.0531 -0.2898 0.1836 -0.92 0.5000     OK
#>  n_pairs
#>       52
#>       64
#>       34
#>       50
#> 
#> >> correlation falls away with separation as the fitted model says it should.
```
