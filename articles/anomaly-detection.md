# Anomaly detection

Two functions ask different questions, and the difference is the whole
point.

``` r

library(illume)
```

## One column, or the combination

[`ilm_outliers()`](https://huttoncp.github.io/illume/reference/ilm_outliers.md)
asks whether a **value** is extreme for its own column.

``` r

ilm_outliers(d$height)
ilm_outliers_all(d)            # every column at once
```

[`ilm_anomaly()`](https://huttoncp.github.io/illume/reference/ilm_anomaly.md)
asks whether a **row** is a plausible combination.

``` r

ilm_anomaly(d)
```

Someone 150 cm tall is unremarkable. Someone weighing 110 kg is
unremarkable. Someone who is both is not, and nothing in either column’s
distribution says so. On simulated data with ten rows pushed off the
correlation structure but inside every column’s range, the
column-at-a-time scan caught **none** of them.

The method: a set of correlated columns puts most of its variation in a
few directions. Fit those directions, project each row onto them, and
score how much of the row is left over. A row that respects the
correlations needs only a few numbers to describe it; one that does not
leaves a large remainder.

``` r

r <- ilm_anomaly(d)
head(r)
#>  row score        p  p_adj  flag driver
#>    1 3.746 8.55e-05 0.0256  TRUE      a
```

`driver` is the column contributing most to that row’s score. An anomaly
nobody can explain is not actionable, so it is reported alongside the
score.

## Three things that had to be measured

This looks like a short function and is not. Two of the three decisions
below overturned the first version of it.

### The rank is not the imputation rank

[`ilm_impute()`](https://huttoncp.github.io/illume/reference/ilm_impute.md)
chooses the number of directions by cross-validating held-out cells,
which is exactly right for filling a value in and wrong here. On a
rank-2 structure in eight columns it chose **6 or 7** – on clean data as
well as contaminated – and those extra components span the very
directions the anomalies depart along. Detection fell from 0.975 to
0.560.

[`ilm_anomaly()`](https://huttoncp.github.io/illume/reference/ilm_anomaly.md)
uses parallel analysis instead: permute each column independently, which
destroys everything the columns share while leaving each one’s own
distribution alone, and keep a component only when it beats what the
permuted data produces. That chose 2 on every dataset tried,
contaminated or not.

The two functions answer different questions – what best predicts a
missing cell, against how many directions are real shared structure – so
they use different criteria, deliberately.

### The fit is trimmed

The anomalies sit in the same data the directions are estimated from, so
they pull those directions towards themselves and are then reconstructed
well. The method hides what it is looking for.

`trim = 0.25` fits on the best-fitting three quarters of rows and
refits, so a row is scored against a structure it did not help define.
It earns its place only where the anomalies share a direction, which is
where they can form a component between them:

      anomalies along random directions      along ONE shared direction
       2%   0.925  vs  0.917                  2%   0.950  vs  0.942
       5%   0.967  vs  0.977                  5%   0.797  vs  0.887
      10%   0.970  vs  0.978                 10%   0.587  vs  0.670

Holding rows out in folds instead – which was the first design – does
nothing, and is not offered. At four contamination levels it matched
in-sample scoring to three decimals, because four fifths of the
anomalies remain in every training fold.

### The reference is simulated, and calibrated

The score is not chi-squared: the noise scale is estimated, the rank
came from the same data, and the residual is taken against an estimated
subspace. Datasets with the same structure and no anomalies are
simulated and scored through the identical pipeline.

Getting that simulation right mattered twice. The noise has to be **per
column** – standardising gives every column unit variance, so one that
loads weakly on the shared directions keeps proportionally more of its
variance off them – and its level has to be matched to the observed
median rather than read off shrunken in-sample residuals. A first
version got both wrong and flagged **38.9%** of the rows of data
containing no anomalies at all.

As it stands: across 100 clean datasets of 400 rows, 11 rows in 40,000
were flagged, with 95 of the 100 producing none; with anomalies present,
98% of the planted rows were found at 2% contamination and 99.7% at 5%.

`B` defaults to 39 because a Benjamini-Hochberg adjusted p-value cannot
fall much below `1 / B`. At `B = 19` a lone anomaly among 300 rows could
not be flagged however extreme it was, and the function says so rather
than returning an empty list.

## Seeing the scan

A ranked table cannot say the one thing that decides what to do next:
whether the flagged rows **stand apart** from the rest, or are only the
top few percent of a smooth continuum that the threshold happened to
cut. Five genuine outliers and the upper tail of a heavy-tailed column
can produce the same list.
[`ilm_plot_anomaly()`](https://huttoncp.github.io/illume/reference/ilm_plot_anomaly.md)
shows the difference, four ways:

``` r

ilm_plot_anomaly(r)                    # scores against rank, and the reference
ilm_plot_anomaly(r, "drivers")         # which column drives the flags
ilm_plot_anomaly(r, "map")             # where the flagged rows sit
ilm_plot_anomaly(r, "row", row = 17)   # one row, column by column
```

**Scores.** Every row’s score against its rank, with a line at the
number flagged and the band the scan simulated: where 95% of datasets
with the same structure and no anomalies put the score at each rank.
Genuine outliers are a few points well above the band, after which the
curve drops back into it. A continuum has no break, and the plot says
so. The choice then is to say plainly that you are reporting a fixed
share of the rows, or to treat the tail as structure and look for it
with
[`ilm_profile()`](https://huttoncp.github.io/illume/reference/ilm_profile.md).

The band **restarts at the line**, and the step is deliberate. Setting
the flagged rows aside moves every other row up that many ranks, so
against the band as simulated a perfectly clean remainder scores too
high for a long stretch. Over 60 scans with 2% planted anomalies, 0.69
of the ranks after the line sat above that band – not far short of the
0.82 to 1.00 of noise that genuinely has heavy tails, which is the
confusion the plot exists to prevent. Right of the line the band is the
reference for the unflagged rows on their own.

The verdict is a reading of that band, not a test. With planted
anomalies it said the flagged rows stand clear in 89 and 97 scans of 100
and never said there was no break; with t-distributed noise it said
there was no break in between half and two thirds of the scans that
flagged anything. Nothing is read off the band when nothing is flagged,
because it inherits the reference’s slightly light upper tail – the same
thing
[`?ilm_anomaly`](https://huttoncp.github.io/illume/reference/ilm_anomaly.md)
reports as the raw p-value running hot.

**Drivers.** How many flagged rows each column drives, beside what the
drivers of the other rows would predict. One column driving nearly every
flag is usually a problem *in* that column – a unit, a sentinel code, a
misplaced decimal – rather than in how the columns combine.

**Map.** The rows on the first two dimensions of
[`ilm_reduce()`](https://huttoncp.github.io/illume/reference/ilm_reduce.md),
the flagged ones in red and sized by score. A reconstruction anomaly is
odd *off* the main dimensions, so it need not sit at the edge of this
map: the view shows whether the flagged rows are alike, not whether they
are extreme.

**Row.** One row’s z-scores beside its residuals, each residual divided
by its column’s typical residual so that 3 means unusual on both. A
small \|z\| beside a large residual is a row that is odd only as a
combination – the case the reconstruction exists for, and the one a
column-at-a-time scan cannot see.

An isolation forest has no null behind its score, so its `"scores"` view
has no band and labels its line as the share you chose to call
anomalous; its `"row"` view needs a fitted structure and is not
available.

## What it cannot do

When a large share of rows depart along the **same** direction, they are
not anomalies – they are a subpopulation, and a rank-`k` fit of the
whole data legitimately includes their direction. Detection degrades
accordingly: 0.89 at 5% contamination, 0.67 at 10%, 0.41 at 20%.

That is the method reaching its limit rather than failing quietly, and
the tool for a second group is
[`ilm_cluster()`](https://huttoncp.github.io/illume/reference/ilm_cluster.md),
which is looking for exactly that.

``` r

ilm_profile(d)                 # is it a group rather than a scattering?
```

## A flagged row is not an error

It is a combination the other rows do not make. It may be a data-entry
mistake, or a real and interesting case, or the most important
observation in the dataset. Deleting it because a method said so is how
real effects get removed, and the print says so every time.

## See also

[`?ilm_anomaly`](https://huttoncp.github.io/illume/reference/ilm_anomaly.md),
[`?ilm_plot_anomaly`](https://huttoncp.github.io/illume/reference/ilm_plot_anomaly.md),
[`?ilm_outliers`](https://huttoncp.github.io/illume/reference/ilm_outliers.md),
[`vignette("profiling")`](https://huttoncp.github.io/illume/articles/profiling.md)
for the dimension reduction underneath, and
[`vignette("workflow")`](https://huttoncp.github.io/illume/articles/workflow.md)
for where this sits.
