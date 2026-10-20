# Fixed effects as a matrix

The same numbers as
[`coef.ilm_model()`](https://huttoncp.github.io/illume/reference/coef.ilm_model.md)
but shaped as predictors by categories, which is usually easier to read
for a multinomial model. `ilm_se_fixef()` returns matching standard
errors.

The same numbers as
[`coef.ilm_model()`](https://huttoncp.github.io/illume/reference/coef.ilm_model.md)
but shaped as predictors by categories, which is usually easier to read
for a multinomial model. `ilm_se_fixef()` returns matching standard
errors.

## Usage

``` r
# S3 method for class 'ilm_model'
fixef(object, ...)

ilm_se_fixef(object)
```

## Arguments

- object:

  A fitted `"ilm_model"` object.

- ...:

  Unused.

## Value

A numeric matrix with one row per design column and one column per
category dimension.

A numeric matrix with one row per design column and one column per
category dimension.
