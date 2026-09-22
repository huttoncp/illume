# Clean up a messy data frame

Removes empty rows and columns, standardises names to snake_case, and
converts text columns that are really numeric or logical.

## Usage

``` r
ilm_wash_df(
  data,
  clean_names = TRUE,
  retype = TRUE,
  drop_empty = TRUE,
  column_to_rownames = FALSE,
  names_col = NULL
)
```

## Arguments

- data:

  A data frame.

- clean_names:

  Standardise column names to snake_case.

- retype:

  Convert character columns that are really numeric or logical.

- drop_empty:

  Drop rows and columns that are entirely missing or blank.

- column_to_rownames:

  Use a column's values as row names.

- names_col:

  The column to use when `column_to_rownames = TRUE`.

## Value

A cleaned data frame.

## Details

Conversion happens only when every non-missing value converts cleanly,
so a single stray `"n/a"` cannot silently turn a column into `NA`s.

## Examples

``` r
m <- data.frame("Col One" = c("1", "2", ""), someFlag = c("TRUE", "FALSE", ""),
                check.names = FALSE)
ilm_wash_df(m)
#>   col_one some_flag
#> 1       1      TRUE
#> 2       2     FALSE
```
