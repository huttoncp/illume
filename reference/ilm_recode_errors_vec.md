# Replace known-bad values in a vector

The vector path of
[`ilm_recode_errors()`](https://huttoncp.github.io/illume/reference/ilm_recode_errors.md),
exposed for use inside other pipelines. Factor levels that are recoded
away are dropped.

## Usage

``` r
ilm_recode_errors_vec(x, errors, replacement = NA)
```

## Arguments

- x:

  A vector.

- errors:

  Values to recode.

- replacement:

  What to put in their place.

## Value

A vector of the same length as `x`.

## Examples

``` r
ilm_recode_errors_vec(factor(c("a", "b", "unknown")), errors = "unknown")
#> [1] a    b    <NA>
#> Levels: a b
```
