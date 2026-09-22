# Flag unusual values in a numeric vector

Scores every value by how far it sits from the centre, in the units of
whichever rule is chosen, and flags those past a threshold.

## Usage

``` r
ilm_outliers(
  y,
  method = c("iqr", "mad", "zscore"),
  threshold = NULL,
  na.rm = TRUE
)
```

## Arguments

- y:

  A numeric vector.

- method:

  `"iqr"`, `"mad"` or `"zscore"`.

- threshold:

  Flagging threshold. Defaults per rule: 1.5 for `"iqr"` (the boxplot
  convention), 3.5 for `"mad"` (Iglewicz and Hoaglin), 3 for `"zscore"`.

- na.rm:

  Compute the reference statistics with missing values removed. The
  result always has one row per element of `y`, with `NA` for `score`
  and `is_outlier` wherever `y` is `NA`.

## Value

A data frame with `value`, `method`, `threshold`, `score` (distance from
the centre in the rule's own units) and `is_outlier`.

## The three rules

- `"iqr"` – distance beyond the nearer quartile in interquartile ranges.
  The classic Tukey fence, and the same rule a boxplot's whiskers draw.

- `"mad"` – distance from the median in median absolute deviations. More
  robust than `"zscore"` on skewed or heavy-tailed data, because neither
  the median nor the MAD is itself dragged around by the values being
  detected, where the mean and standard deviation both are.

- `"zscore"` – distance from the mean in standard deviations. Included
  because it is expected, but it is the weakest of the three for exactly
  that reason: a large outlier inflates the standard deviation and so
  hides itself.

## What a flag is not

It is not a verdict that the value is wrong. A flagged value may be a
recording error, a member of a different population, or an ordinary draw
from a heavy tail – and the third is common. Dropping flagged rows
because they are flagged changes the estimand and biases whatever is
fitted next; if the tail is real, the remedy is a model that expects it,
which for illume means a different family or [ilm_model(dispformula =
)](https://craig-hutton.github.io/illume/reference/ilm_model.md).

## References

Iglewicz, B. and Hoaglin, D. C. (1993). How to Detect and Handle
Outliers. ASQC Quality Press.

## See also

[`ilm_outliers_all()`](https://craig-hutton.github.io/illume/reference/ilm_outliers_all.md)
for a whole data frame,
[`ilm_describe()`](https://craig-hutton.github.io/illume/reference/ilm_describe.md)
for the distribution a flag should be read against.

## Examples

``` r
table(ilm_outliers(mtcars$hp)$is_outlier)
#> 
#> FALSE  TRUE 
#>    31     1 
head(ilm_outliers(mtcars$hp, method = "mad"))
#>   value method threshold score is_outlier
#> 1   110    mad       3.5 0.169      FALSE
#> 2   110    mad       3.5 0.169      FALSE
#> 3    93    mad       3.5 0.389      FALSE
#> 4   110    mad       3.5 0.169      FALSE
#> 5   175    mad       3.5 0.674      FALSE
#> 6   105    mad       3.5 0.233      FALSE
```
