# Bootstrap intervals for differences between groups

Reports each difference with its own interval, rather than a set of
intervals to compare by eye: judging a difference by whether separate
intervals overlap is conservative and lossy.

## Usage

``` r
ilm_boot_diff(x, ...)

# S3 method for class 'formula'
ilm_boot_diff(x, data = NULL, ...)

# S3 method for class 'data.frame'
ilm_boot_diff(
  x,
  y = NULL,
  group = NULL,
  stat = "mean",
  R = 2000L,
  conf = 0.95,
  ci_type = "percentile",
  adjust = "max_t",
  ref = NULL,
  seed = NULL,
  progress = NULL,
  ...
)
```

## Arguments

- x:

  A data frame, or a formula `y ~ g` with `data` supplied.

- ...:

  Passed between methods.

- data:

  The data frame, when `x` is a formula.

- y:

  Name of the numeric column. Ignored when `x` is a formula.

- group:

  Name of the grouping column, or several to cross with
  [`interaction()`](https://rdrr.io/r/base/interaction.html). Ignored
  when `x` is a formula.

- stat:

  `"mean"`, `"median"`, `"sd"`, `"var"`, or a function taking a numeric
  vector.

- R:

  Bootstrap replicates.

- conf:

  Confidence level. Under `"max_t"` and `"bonferroni"` this is the level
  of the whole family of comparisons, not of one interval.

- ci_type:

  `"percentile"`, `"bca"`, `"normal"` or `"basic"`. Applies to
  `"bonferroni"` and `"none"`; `"max_t"` has its own construction.

- adjust:

  `"max_t"`, `"bonferroni"` or `"none"`; see Details.

- ref:

  Optional level to compare every other level against, giving `J - 1`
  comparisons instead of `J * (J - 1) / 2`.

- seed:

  Random seed.

- progress:

  Show a progress bar. Defaults to
  [`interactive()`](https://rdrr.io/r/base/interactive.html), so a bar
  appears when someone is watching and nothing is written in a script or
  a knitted document. See
  [ilm_progress_arg](https://huttoncp.github.io/illume/reference/ilm_progress_arg.md).

## Value

One row per comparison, with `from`, `to`, `observed`, `lower`, `upper`,
`p_value`, `p_adj`, `excludes_zero` and the settings used. Differences
are `to` minus `from`.

## Details

With more than two groups every pair is compared, and the intervals are
**simultaneous by default**. Ten comparisons each at a nominal 95% do
not jointly cover at 95%, and reporting them as though they did is the
usual way a pairwise table misleads. `adjust` chooses how that is
handled:

- `"max_t"` (default) resamples all groups together, standardises each
  comparison by its own bootstrap standard error, and takes the `conf`
  quantile of the largest absolute standardised value across comparisons
  as one critical value for all of them. This is the single-step
  studentized maximum. It assumes neither normality nor a common
  variance, which is what separates it from Tukey's range test. Its
  intervals are symmetric about the observed difference, so `ci_type`
  does not apply to them.

- `"bonferroni"` builds each interval at level `1 - (1 - conf) / m` in
  the requested `ci_type`. More conservative than `"max_t"`, but it
  keeps the shape of a percentile or BCa interval, which matters on
  skewed data.

- `"none"` builds each interval at `conf`. Correct for one comparison,
  optimistic for several; use it when the comparisons were chosen in
  advance.

With a single comparison there is nothing to adjust, and the `adjust`
column of the result reads `"none"` whatever was asked for.

`p_value` and `p_adj` are built the same way the interval is, so the two
cannot contradict each other: under `"max_t"` both come off the
studentized maximum, and otherwise both come off the bootstrap
distribution directly. They are read off the same `R` replicates from
opposite directions – a quantile and a tail proportion – so a comparison
sitting within one replicate of the critical value can still land either
side. Raise `R` if a result is that close to the line. With `ci_type` of
`"basic"`, `"normal"` or `"bca"` the correspondence is close rather than
exact, since those reshape the same replicates while the p-value does
not.

## What this is calibrated for

Against [`stats::TukeyHSD()`](https://rdrr.io/r/stats/TukeyHSD.html) on
normal, equal-variance, balanced data – the case Tukey is exactly right
for – the endpoints agree to 0.006 and the adjusted p-values to 0.011.
Family-wise error over six comparisons at group sizes of 45 to 80 came
to 0.068 on normal data with a common variance and 0.072 with variances
differing fourfold, against a nominal 0.05.

Heavy skew at small n is the case to know about. On lognormal data with
a spread parameter up to 1.2 and those same group sizes, a SINGLE
unadjusted comparison already erred 0.090 of the time, and no `ci_type`
moved it (percentile 0.090, BCa 0.089, basic 0.093, normal 0.084). That
is the bootstrapped mean of a heavily skewed small sample, not the
multiplicity adjustment, which is still doing its work: those six
comparisons reject 0.302 of the time unadjusted against 0.126 under
`"max_t"`. It converges – at group sizes of 300 to 500 the same design
gives 0.049 per comparison and 0.062 family-wise.

`stat = "median"` is the remedy when the data look like that. On the
same lognormal design, under a null placing the medians rather than the
means together, it gave 0.043 per comparison and 0.050 family-wise.

## References

Efron, B. and Tibshirani, R. J. (1993). An Introduction to the
Bootstrap. Chapman and Hall.

Westfall, P. H. and Young, S. S. (1993). Resampling-Based Multiple
Testing. Wiley.

## See also

[`ilm_boot_ci()`](https://huttoncp.github.io/illume/reference/ilm_boot_ci.md)
for a single group.

## Examples

``` r
d <- ilm_sim()
## all pairs, simultaneous by default
ilm_boot_diff(d, "score", "grp", R = 300, seed = 1)
#>   stat group  from    to   observed      lower     upper conf   R ci_type
#> 1 mean   grp alpha  beta -0.7764371 -1.8206914 0.2678173 0.95 300   max_t
#> 2 mean   grp alpha gamma -0.3748473 -1.6604311 0.9107366 0.95 300   max_t
#> 3 mean   grp alpha delta  3.3021040  0.2588339 6.3453740 0.95 300   max_t
#> 4 mean   grp  beta gamma  0.4015898 -0.8919237 1.6951033 0.95 300   max_t
#> 5 mean   grp  beta delta  4.0785410  1.0519031 7.1051789 0.95 300   max_t
#> 6 mean   grp gamma delta  3.6769512  0.4048276 6.9490749 0.95 300   max_t
#>   adjust n_comparisons n_from n_to   p_value       p_adj p_superiority
#> 1  max_t             6    404  329 0.0800000 0.266666667        0.0400
#> 2  max_t             6    404  164 0.4900000 0.920000000        0.2600
#> 3  max_t             6    404    3 0.0000000 0.026666667        1.0000
#> 4  max_t             6    329  164 0.4666667 0.890000000        0.7667
#> 5  max_t             6    329    3 0.0000000 0.006666667        1.0000
#> 6  max_t             6    164    3 0.0000000 0.016666667        1.0000
#>   excludes_zero
#> 1         FALSE
#> 2         FALSE
#> 3          TRUE
#> 4         FALSE
#> 5          TRUE
#> 6          TRUE
## the same thing through a formula
ilm_boot_diff(score ~ grp, data = d, R = 300, seed = 1)
#>   stat group  from    to   observed      lower     upper conf   R ci_type
#> 1 mean   grp alpha  beta -0.7764371 -1.8206914 0.2678173 0.95 300   max_t
#> 2 mean   grp alpha gamma -0.3748473 -1.6604311 0.9107366 0.95 300   max_t
#> 3 mean   grp alpha delta  3.3021040  0.2588339 6.3453740 0.95 300   max_t
#> 4 mean   grp  beta gamma  0.4015898 -0.8919237 1.6951033 0.95 300   max_t
#> 5 mean   grp  beta delta  4.0785410  1.0519031 7.1051789 0.95 300   max_t
#> 6 mean   grp gamma delta  3.6769512  0.4048276 6.9490749 0.95 300   max_t
#>   adjust n_comparisons n_from n_to   p_value       p_adj p_superiority
#> 1  max_t             6    404  329 0.0800000 0.266666667        0.0400
#> 2  max_t             6    404  164 0.4900000 0.920000000        0.2600
#> 3  max_t             6    404    3 0.0000000 0.026666667        1.0000
#> 4  max_t             6    329  164 0.4666667 0.890000000        0.7667
#> 5  max_t             6    329    3 0.0000000 0.006666667        1.0000
#> 6  max_t             6    164    3 0.0000000 0.016666667        1.0000
#>   excludes_zero
#> 1         FALSE
#> 2         FALSE
#> 3          TRUE
#> 4         FALSE
#> 5          TRUE
#> 6          TRUE
## every level against one reference, on medians
ilm_boot_diff(score ~ grp, data = d, stat = "median", ref = "alpha",
              R = 300, seed = 1)
#>     stat group  from    to observed      lower     upper conf   R ci_type
#> 1 median   grp alpha  beta   -0.760 -2.3104532 0.7904532 0.95 300   max_t
#> 2 median   grp alpha gamma   -0.755 -2.1822794 0.6722794 0.95 300   max_t
#> 3 median   grp alpha delta    3.570 -0.5570813 7.6970813 0.95 300   max_t
#>   adjust n_comparisons n_from n_to   p_value      p_adj p_superiority
#> 1  max_t             3    404  329 0.3300000 0.70666667        0.2067
#> 2  max_t             3    404  164 0.2466667 0.62666667        0.1367
#> 3  max_t             3    404    3 0.0100000 0.09333333        0.9900
#>   excludes_zero
#> 1         FALSE
#> 2         FALSE
#> 3         FALSE
```
