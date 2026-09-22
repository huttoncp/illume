# Mark censored observations

Describes which observations are known only as an interval: values at or
below a floor, at or above a ceiling, or both.

## Usage

``` r
ilm_censor(y, lower = NA, upper = NA, status = NULL)
```

## Arguments

- y:

  The response.

- lower:

  Floor. Values at or below it are left-censored.

- upper:

  Ceiling. Values at or above it are right-censored.

- status:

  Alternatively, codes given directly: `-1` left-censored, `0` observed,
  `1` right-censored.

## Value

An integer vector of codes, classed `"ilm_censor"`, carrying the limits
it was built from.

## Details

Pass the result to `ilm_model(censor = )`. With `lower` or `upper`
given, the codes are derived from the response, and the limits are
carried along so that the simulation-based diagnostics can censor their
replicates the same way.

`status` is for the case where the limits vary by observation and only
the outcome is known – right-censored follow-up in a survival study,
most commonly. Diagnostics then hold the censoring pattern fixed across
replicates rather than re-drawing it, which is stated where it matters.

## See also

[`ilm_model()`](https://huttoncp.github.io/illume/reference/ilm_model.md),
[`ilm_describe()`](https://huttoncp.github.io/illume/reference/ilm_describe.md),
whose `p_zero` column is often the first sign of a floor.

## Examples

``` r
y <- c(0, 0, 1.4, 2.9, 5, 5)
ilm_censor(y, lower = 0, upper = 5)
#> censoring: 2 observed, 2 left, 2 right (of 6)
#>   limits: lower 0, upper 5
#>   67% censored; the fit rests mostly on interval probabilities
table(ilm_censor(y, lower = 0))
#> 
#> -1  0 
#>  2  4 
```
