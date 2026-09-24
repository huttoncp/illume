# car::Anova method

Lets `car::Anova(fit, type = 3)` dispatch to
[`ilm_anova()`](https://huttoncp.github.io/illume/reference/ilm_anova.md),
so it uses the correct joint blocking across categories instead of
`car`'s default, which would test only one category.

## Usage

``` r
# S3 method for class 'ilm_model'
Anova(mod, type = c("II", "III", 2, 3), test.statistic = "Chisq", ...)
```

## Arguments

- mod:

  A fitted `"ilm_model"` object.

- type:

  `"II"`, `"III"`, `2` or `3`.

- test.statistic:

  Ignored; present for compatibility.

- ...:

  Passed to
  [`ilm_anova()`](https://huttoncp.github.io/illume/reference/ilm_anova.md).

## Value

An `"anova"` data frame.

## Details

Registered with car's generic when car is loaded. It used to be exported
as an ordinary function, which R's method lookup does not reach, so
[`car::Anova()`](https://rdrr.io/pkg/car/man/Anova.html) fell back to
its own tests: on a multinomial fit with a numeric predictor and a
three-level factor, 1 and 2 degrees of freedom where the joint tests
have 2 and 4.
