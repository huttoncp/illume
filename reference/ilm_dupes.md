# Duplicated rows only

A shorthand for `ilm_copies(data, ..., filter = "dupes")`.

## Usage

``` r
ilm_dupes(data, ..., na_last = TRUE)
```

## Arguments

- data:

  A data frame.

- ...:

  Columns defining a duplicate. All columns if none are given.

- na_last:

  Sort missing values last.

## Value

A data frame of repeated rows with `n_copies` appended.

## Examples

``` r
ilm_dupes(data.frame(a = c(1, 1, 2), b = c("x", "x", "y")))
#>   a b n_copies
#> 1 1 x        2
#> 2 1 x        2
```
