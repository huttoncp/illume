# Rows that do not fit the pattern the other rows make

Finds observations that are implausible as a **combination** of values,
even when no single value is extreme.
[`ilm_outliers()`](https://huttoncp.github.io/illume/reference/ilm_outliers.md)
asks whether a number is far out in its own column; this asks whether a
row is far from the structure the columns share.

## Usage

``` r
ilm_anomaly(
  data,
  cols = NULL,
  method = c("reconstruction", "iforest"),
  rank = NULL,
  trim = 0.25,
  ntrees = 500L,
  B = 39L,
  alpha = 0.05,
  seed = 1L,
  keep_data = TRUE,
  progress = NULL
)
```

## Arguments

- data:

  A data frame or matrix.

- cols:

  Columns to use; see
  [ilm_selection](https://huttoncp.github.io/illume/reference/ilm_selection.md).
  Under the default method, numeric columns only – the reconstruction is
  a projection, and a category has no residual along a direction.
  Anything else is dropped with a note. `method = "iforest"` uses every
  column.

- method:

  `"reconstruction"` (the default) or `"iforest"`. The reconstruction is
  unbeaten at what it is for – against planted anomalies it scores 1.000
  for an extreme value and 0.999 for a jointly-implausible numeric
  combination – and it is the only one here with a calibrated null, so
  it reports FDR-adjusted p-values. It is also blind to the categories:
  0.771 when a category contradicts the numbers, and 0.505, a coin toss,
  for a category pairing that never otherwise occurs. An isolation
  forest splits on factors directly and reaches 0.964 and 0.886 on those
  two, giving up little elsewhere (0.975, 0.950) – but it returns a
  score with no null behind it, so `p` and `p_adj` are `NA` and `alpha`
  becomes the share of rows you are calling anomalous rather than an
  error rate being controlled. Needs the isotree package.

- rank:

  Number of directions; `NULL` uses parallel analysis.

- trim:

  Share of the worst-fitting rows held out of the fit. `0` fits every
  row, which lets the anomalies define the structure they are scored
  against.

- ntrees:

  Trees in the isolation forest. Ignored by the default method.

- B:

  Simulated null datasets for the reference. A Benjamini-Hochberg
  adjusted p-value cannot fall much below `1 / B`, so a single anomaly
  among many rows needs `B` comfortably above `1 / alpha`.

- alpha:

  Flagging level, after a Benjamini-Hochberg adjustment across rows –
  one row in twenty at 0.05 would be 50 rows in a thousand, which is a
  list nobody reads.

- seed:

  Random seed.

- keep_data:

  Keep the scanned data on the result, so that
  [`ilm_anomalous()`](https://huttoncp.github.io/illume/reference/ilm_anomalous.md),
  [`ilm_profile()`](https://huttoncp.github.io/illume/reference/ilm_profile.md),
  [`ilm_cluster()`](https://huttoncp.github.io/illume/reference/ilm_cluster.md)
  and the describe functions can take the object directly. Set `FALSE`
  if the frame is large and you only want the scores.

- progress:

  Show a progress bar; see
  [ilm_progress_arg](https://huttoncp.github.io/illume/reference/ilm_progress_arg.md).

## Value

A data frame with one row per observation: `row`, `score`, `p`, `p_adj`,
`flag`, and `driver`, the column contributing most to the score. The
rank and the residual matrix are attributes.

## Details

A set of correlated columns puts most of its variation in a few
directions. Those directions are estimated, each row is projected onto
them, and the score is what is left over.

## How many directions

By parallel analysis, and deliberately not by the cross-validation
[`ilm_impute()`](https://huttoncp.github.io/illume/reference/ilm_impute.md)
uses on the same decomposition, because the two answer different
questions. Imputation wants the rank that best predicts a missing cell.
This wants the number of directions that are real shared structure, so
that what is left over is residual rather than signal it failed to fit.

The difference is large. On a rank-2 structure in eight columns,
cross-validation chose 6 or 7 – on clean data as well as contaminated –
and those extra components span the directions the anomalies depart
along. Detection fell from 0.975 to 0.560. Parallel analysis chose 2
every time.

## Why the fit is trimmed

The anomalies sit in the same data the directions are estimated from, so
they pull the directions towards themselves and are then reconstructed
well. `trim` excludes the worst-fitting rows from the fit and refits, so
a row is scored against a structure it did not help define. Every row is
still scored, the trimmed ones included.

It earns its place only where the anomalies share a direction, which is
where they can form a component between them. Detection with `trim = 0`
against the default 0.25, over 12 datasets each:

      anomalies along random directions      along ONE shared direction
       2%   0.925  vs  0.917                  2%   0.950  vs  0.942
       5%   0.967  vs  0.977                  5%   0.797  vs  0.887
      10%   0.970  vs  0.978                 10%   0.587  vs  0.670

Holding rows out in folds instead, which was the first design here, does
not work and is not offered: measured at four contamination levels it
matched in-sample scoring to three decimals, because most of the
anomalies remain in every training fold.

## What it cannot do

When a large share of rows depart along the **same** direction they are
not anomalies, they are a subpopulation, and a rank-k fit of the whole
data legitimately includes their direction. Detection degrades
accordingly: with anomalies sharing one direction, 0.89 of them were
found at 5% contamination, 0.67 at 10% and 0.41 at 20%. That is the
method reaching its limit rather than failing quietly, and
[`ilm_cluster()`](https://huttoncp.github.io/illume/reference/ilm_cluster.md)
is the tool for a second group, since finding one is what it is for.

## What the score is compared with

Not a chi-squared distribution. Datasets with the same structure and no
anomalies are simulated and scored the same way. `B` can be modest
because each simulated dataset contributes `n` null scores rather than
one.

Measured on 100 matrices of 400 by 8 with a rank-2 structure and **no
anomalies at all**: 11 rows of 40,000 were flagged, a rate of 0.00028,
with 95 of the 100 datasets producing none. With anomalies present, the
rows that were not planted were flagged at 0.0015 to 0.0030 across every
design tried above.

`flag` uses `p_adj` and is what those numbers describe. The raw `p` runs
a little hot – about 0.064 of clean rows fall below 0.05 rather than
0.05 of them – so read it as a ranking rather than as a test.

## References

Hawkins, D. M. (1974). The detection of errors in multivariate data
using principal components. *Journal of the American Statistical
Association* 69, 340-344.

Horn, J. L. (1965). A rationale and test for the number of factors in
factor analysis. *Psychometrika* 30, 179-185.

## See also

[`ilm_outliers()`](https://huttoncp.github.io/illume/reference/ilm_outliers.md)
for the one-column-at-a-time question,
[`ilm_cluster()`](https://huttoncp.github.io/illume/reference/ilm_cluster.md)
when the unusual rows turn out to be a group,
[`ilm_impute()`](https://huttoncp.github.io/illume/reference/ilm_impute.md)
which fits the same decomposition to fill values in.

## Examples

``` r
set.seed(1)
n <- 300
f <- rnorm(n)
d <- data.frame(a = f + rnorm(n, 0, .3), b = 2 * f + rnorm(n, 0, .3),
                c = -f + rnorm(n, 0, .3))
## a row that is ordinary in every column but breaks the pattern
d[1, ] <- c(1.2, -2.4, 1.2)
head(ilm_anomaly(d), 3)
#> Multivariate anomalies: 1 of 3 rows flagged at 0.05
#>   rank 1 over 3 columns, fitted on the best 75% of rows, 11700 null scores
#>  row score        p  p_adj  flag driver
#>    1 3.746 8.55e-05 0.0256  TRUE      a
#>  195 0.870 2.05e-03 0.3080 FALSE      a
#>   75 0.756 4.19e-03 0.3440 FALSE      c
#> 
#>   `driver` is the column contributing most to each row's score.
#>   A flagged row is a combination the other rows do not make. It is
#>   not necessarily an error, and deleting it because a method said so
#>   is how real effects get removed.
```
