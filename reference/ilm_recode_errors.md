# Replace known-bad values with NA or another value

Sentinels such as 999, -1 or `"unknown"` mean missing in one column and
can be a legitimate measurement in another, which is why the replacement
can be restricted to particular rows and columns.

## Usage

``` r
ilm_recode_errors(
  data,
  errors,
  replacement = NA,
  rows = NULL,
  cols = NULL,
  ind = NULL
)
```

## Arguments

- data:

  A vector, data frame or matrix.

- errors:

  Values to recode.

- replacement:

  What to put in their place. `NA` by default.

- rows, cols:

  Restrict the replacement (data frame or matrix input).

- ind:

  Restrict the replacement (vector input).

## Value

An object of the same shape as `data`.

## Examples

``` r
ilm_recode_errors(c(1, 2, 999, -1), errors = c(999, -1))
#> [1]  1  2 NA NA
d <- data.frame(x = c(1, 999, 3), y = c(999, 2, 3))
ilm_recode_errors(d, errors = 999, cols = "x")
#>    x   y
#> 1  1 999
#> 2 NA   2
#> 3  3   3
```
