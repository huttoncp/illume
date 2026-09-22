# car::Anova method

Lets `car::Anova(fit, type = 3)` dispatch to
[`ilm_anova()`](https://craig-hutton.github.io/illume/reference/ilm_anova.md),
so it uses the correct joint blocking across categories instead of
`car`'s default, which would test only one category.

## Usage

``` r
Anova.ilm_model(
  mod,
  type = c("II", "III", 2, 3),
  test.statistic = "Chisq",
  ...
)
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
  [`ilm_anova()`](https://craig-hutton.github.io/illume/reference/ilm_anova.md).

## Value

An `"anova"` data frame.
