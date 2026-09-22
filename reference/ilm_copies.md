# Find copied or duplicated rows

Find copied or duplicated rows

## Usage

``` r
ilm_copies(
  data,
  ...,
  filter = c("all", "dupes", "first", "last", "unique"),
  na_last = TRUE,
  sort_by = TRUE
)
```

## Arguments

- data:

  A data frame.

- ...:

  Columns defining a duplicate. All columns if none are given.

- filter:

  `"all"` appends `copy_number` and `n_copies`; `"dupes"` keeps only
  repeated rows; `"first"`, `"last"` and `"unique"` filter rows.

- na_last:

  Sort missing values last.

- sort_by:

  Sort the result by the key columns.

## Value

A data frame; the columns depend on `filter`.

## See also

[`ilm_dupes()`](https://craig-hutton.github.io/illume/reference/ilm_dupes.md)
for the common case.

## Examples

``` r
d <- data.frame(a = c(1, 1, 2), b = c("x", "x", "y"))
ilm_copies(d)
#>   a b copy_number n_copies
#> 1 1 x           1        2
#> 2 1 x           2        2
#> 3 2 y           1        1
ilm_copies(d, filter = "unique")
#>   a b
#> 1 2 y
```
