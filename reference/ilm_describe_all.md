# Describe every column of a data frame

Applies
[`ilm_describe()`](https://craig-hutton.github.io/illume/reference/ilm_describe.md)
to each column and returns one table per class, so columns with
different summaries are not forced into a common shape.

## Usage

``` r
ilm_describe_all(
  data,
  by = NULL,
  digits = 3,
  gauss = c("index", "ks_d", "both", "none"),
  probs = c(0, 0.5, 1),
  skew = FALSE,
  kurt = FALSE,
  dispersion = TRUE,
  class = "all",
  rare_n = 5L,
  cap = 0.12,
  min_n = 20L
)
```

## Arguments

- data:

  A data frame, or a vector when `y` is `NULL`.

- by:

  Optional character vector of grouping columns.

- digits:

  Rounding for numeric columns, quantiles included.

- gauss:

  One of `"index"` (the 0-1 agreement index), `"ks_d"` (the raw
  Kolmogorov distance behind it), `"both"`, or `"none"`.

- probs:

  Quantiles to report, as proportions. `c(0, 0.5, 1)` gives `p0`, `p50`
  and `p100`; any vector works and the column names follow it.

- skew, kurt:

  Add moment-based skewness and excess kurtosis.

- dispersion:

  Add the variance-to-mean ratio for non-negative integer variables.
  Meaningful only if you intend to model the variable as a count.

- class:

  `"all"`, or one or more of `"numeric"`, `"categorical"`, `"logical"`,
  `"time"`.

- rare_n:

  Levels with fewer observations than this count as rare.

- cap, min_n:

  Control the `gauss` index; see
  [`ilm_gauss_check()`](https://craig-hutton.github.io/illume/reference/ilm_gauss_check.md).

## Value

A named list of data frames, one per class present, plus `constant` when
any column has a single value; a single data frame if only one section
results.

## Details

Constant columns are split into their own `constant` section rather than
padding every other table with rows of `NA`: a variable with one value
tells you nothing a statistic can show, and it will break a model
matrix.

## See also

[`ilm_frame_issues()`](https://craig-hutton.github.io/illume/reference/ilm_frame_issues.md)
for problems belonging to pairs of columns.

## Examples

``` r
d <- ilm_sim()
ilm_describe_all(d)
#> $numeric
#>    variable obs   n  na         sum      mean        sd      se      p0
#> 1        id 900 900   0    34200.00    38.000    21.661   0.722    1.00
#> 2     score 900 900   0    44709.09    49.677     5.880   0.196   29.99
#> 3    income 900 900   0 33361609.93 37068.455 22902.677 763.423 6356.18
#> 4    visits 900 900   0     2779.00     3.088     1.830   0.061    0.00
#> 5    claims 900 900   0     3773.00     4.192     4.657   0.155    0.00
#> 6  downtime 900 900   0     1339.00     1.488     2.386   0.080    0.00
#> 7 lab_value 900 792 108     5840.86     7.375     1.153   0.041    4.04
#>         p50      p100 p_zero dispersion gauss
#> 1    38.000     75.00  0.000     12.347 0.719
#> 2    49.835     65.34  0.000         NA 1.000
#> 3 31137.270 182886.22  0.000         NA 0.305
#> 4     3.000     10.00  0.053      1.085 0.042
#> 5     3.000     27.00  0.177      5.172 0.000
#> 6     0.000     11.00  0.640      3.826 0.000
#> 7     7.330     11.66  0.000         NA 1.000
#>                                                                                 gauss_note
#> 1                                                     bounded at zero; light-tailed / flat
#> 2                                                                                         
#> 3                                                               right-skewed; heavy-tailed
#> 4                                           discrete (11 distinct values); bounded at zero
#> 5 18% of values sit exactly at the minimum (0): a floor, see ilm_censor(); bounded at zero
#> 6                                           discrete (12 distinct values); bounded at zero
#> 7                                                                                         
#> 
#> $time
#>   variable obs   n na n_unique      start        end span_days spacing regular
#> 1     date 900 900  0       12 2024-01-01 2024-12-01       335 monthly    TRUE
#>   n_gaps n_dup_times   tz
#> 1      0         888 <NA>
#>                                                       note
#> 1 repeated dates: 12 unique across 900 rows (long format?)
#> 
#> $categorical
#>   variable obs   n na n_empty n_unique ordered p_max n_rare n_unused
#> 1      grp 900 900  0       0        4   FALSE 0.449      1        1
#> 2     site 900 900  0      46        6   FALSE 0.308      0        0
#>   case_variants                      counts_tb
#> 1             0 alpha_404, beta_329, gamma_164
#> 2             1 North_277, South_265, East_217
#>                                               note
#> 1 1 unused levels; 1 rare levels (separation risk)
#> 2 46 empty strings (not NA); 1 case/space variants
#> 
#> $logical
#>    variable obs   n na n_TRUE n_FALSE p_TRUE                            note
#> 1      flag 900 900  0    459     441  0.510                                
#> 2 consented 900 900  0    898       2  0.998 near-constant (separation risk)
#> 
#> $constant
#>   variable       class obs   n na value
#> 1   cohort categorical 900 900  0  2024
#> 
ilm_describe_all(d, class = "numeric", skew = TRUE)
#> ilm_describe_all: 1 constant column(s) set aside (cohort); use class = "all" to see them
#>    variable obs   n  na         sum      mean        sd      se      p0
#> 1        id 900 900   0    34200.00    38.000    21.661   0.722    1.00
#> 2     score 900 900   0    44709.09    49.677     5.880   0.196   29.99
#> 3    income 900 900   0 33361609.93 37068.455 22902.677 763.423 6356.18
#> 4    visits 900 900   0     2779.00     3.088     1.830   0.061    0.00
#> 5    claims 900 900   0     3773.00     4.192     4.657   0.155    0.00
#> 6  downtime 900 900   0     1339.00     1.488     2.386   0.080    0.00
#> 7 lab_value 900 792 108     5840.86     7.375     1.153   0.041    4.04
#>         p50      p100   skew p_zero dispersion gauss
#> 1    38.000     75.00  0.000  0.000     12.347 0.719
#> 2    49.835     65.34 -0.077  0.000         NA 1.000
#> 3 31137.270 182886.22  1.707  0.000         NA 0.305
#> 4     3.000     10.00  0.561  0.053      1.085 0.042
#> 5     3.000     27.00  1.946  0.177      5.172 0.000
#> 6     0.000     11.00  1.553  0.640      3.826 0.000
#> 7     7.330     11.66 -0.038  0.000         NA 1.000
#>                                                                                 gauss_note
#> 1                                                     bounded at zero; light-tailed / flat
#> 2                                                                                         
#> 3                                                               right-skewed; heavy-tailed
#> 4                                           discrete (11 distinct values); bounded at zero
#> 5 18% of values sit exactly at the minimum (0): a floor, see ilm_censor(); bounded at zero
#> 6                                           discrete (12 distinct values); bounded at zero
#> 7                                                                                         
```
