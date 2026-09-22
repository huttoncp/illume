# How columns can be chosen

Wherever illume takes a `cols` or `by` argument, it accepts any of four
things. All are ordinary values, so any of them can be held in a
variable and passed along – which a bare-name interface cannot.

## Details

- `NULL`, meaning every eligible column.

- A **character vector** of names: `c("age", "income")`. Every name must
  exist, and a typo is an error naming the columns that do exist.

- A **regular expression**, as a single string that is not itself a
  column name: `"^score_"` takes every column whose name starts with
  `score_`. Matching nothing is an error rather than a silent empty
  selection.

- A **predicate function** applied to each column: `is.numeric`, or
  `function(v) is.numeric(v) && !anyNA(v)`.

The one ambiguity is a single string, which could be a name or a
pattern. A name wins: `cols = "score"` selects the column called `score`
even if `scores_2024` also exists. Anything that is not a column name is
read as a pattern.

## Examples

``` r
d <- data.frame(id = 1:5, score_a = rnorm(5), score_b = rnorm(5),
                label = letters[1:5])
ilm_outliers_all(d, cols = "^score_")
#> [1] row_id     variable   value      score      is_outlier
#> <0 rows> (or 0-length row.names)
ilm_outliers_all(d, cols = is.numeric)
#> [1] row_id     variable   value      score      is_outlier
#> <0 rows> (or 0-length row.names)
```
