# Summarise a fitted multinomial mixed model

Prints the model, fit statistics, random-effect covariances, a
coefficient table with Wald tests, and – unusually for a
[`summary()`](https://rdrr.io/r/base/summary.html) method – the
**assumption checks**.

Prints the model, fit statistics, random-effect covariances, a
coefficient table with Wald tests, and – unusually for a
[`summary()`](https://rdrr.io/r/base/summary.html) method – the
**assumption checks**.

## Usage

``` r
# S3 method for class 'ilm_model'
summary(object, ...)

# S3 method for class 'summary.ilm_model'
print(x, digits = 4, max_corr_dim = 6L, ...)
```

## Arguments

- object, x:

  A fitted `"ilm_model"` object.

- ...:

  Unused.

- digits:

  Integer. Significant digits in the coefficient table.

- max_corr_dim:

  Integer. Largest `C` for which correlation matrices are printed in
  full.

## Value

[`summary()`](https://rdrr.io/r/base/summary.html) returns an object of
class `"summary.ilm_model"`; its print method returns it invisibly.

[`summary()`](https://rdrr.io/r/base/summary.html) returns an object of
class `"summary.ilm_model"`; its print method returns it invisibly.

## Details

The checks appear here on purpose. This model class fails quietly: the
coefficient table can look completely ordinary, significance stars and
all, while the covariance matrix behind it is unusable and every
standard error is meaningless. Putting the verdicts behind a separate
function would mean the people most likely to be misled are the least
likely to look. When any check fails,
[`summary()`](https://rdrr.io/r/base/summary.html) says so directly
beneath the coefficients.

The printed note about sum-to-zero contrasts is also deliberate: a
coefficient here is a deviation from the across-category average, not a
contrast against a baseline category, and readers used to
[`nnet::multinom()`](https://rdrr.io/pkg/nnet/man/multinom.html) will
otherwise assume the latter.

Correlation matrices are suppressed above `max_corr_dim` categories,
where they become too large to read; use `fit$Sigma` for the full
matrices.

The checks appear here on purpose. This model class fails quietly: the
coefficient table can look completely ordinary, significance stars and
all, while the covariance matrix behind it is unusable and every
standard error is meaningless. Putting the verdicts behind a separate
function would mean the people most likely to be misled are the least
likely to look. When any check fails,
[`summary()`](https://rdrr.io/r/base/summary.html) says so directly
beneath the coefficients.

The printed note about sum-to-zero contrasts is also deliberate: a
coefficient here is a deviation from the across-category average, not a
contrast against a baseline category, and readers used to
[`nnet::multinom()`](https://rdrr.io/pkg/nnet/man/multinom.html) will
otherwise assume the latter.

Correlation matrices are suppressed above `max_corr_dim` categories,
where they become too large to read; use `fit$Sigma` for the full
matrices.

## See also

[`ilm_coef_table()`](https://huttoncp.github.io/illume/reference/ilm_coef_table.md),
[`ilm_appraise()`](https://huttoncp.github.io/illume/reference/ilm_appraise.md).

[`ilm_coef_table()`](https://huttoncp.github.io/illume/reference/ilm_coef_table.md),
[`ilm_appraise()`](https://huttoncp.github.io/illume/reference/ilm_appraise.md).
