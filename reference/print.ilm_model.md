# Compact display of a fitted model

A one-line description plus the log-likelihood and AIC. Flags
`[CHECKS FAILED]` when any assumption check did not pass, so a failure
is visible even from a bare
[`print()`](https://rdrr.io/r/base/print.html).

## Usage

``` r
# S3 method for class 'ilm_model'
print(x, ...)
```

## Arguments

- x:

  A fitted `"ilm_model"` object.

- ...:

  Unused.

## Value

`x`, invisibly.
