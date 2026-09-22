# Power for a study that has not been run

Builds a scaffold from your assumptions and runs
[`ilm_power()`](https://huttoncp.github.io/illume/reference/ilm_power.md)
on it, in one call. `n_unit` counts participants when the formula has a
grouping bar and rows when it does not, so the numbers you give are the
numbers you would write in a protocol.

## Usage

``` r
ilm_power_design(
  formula,
  design,
  n_unit,
  family = "gaussian",
  coefs = NULL,
  cells = NULL,
  sd = NULL,
  re_sd = NULL,
  re_cor = NULL,
  icc = NULL,
  within = NULL,
  contrasts = NULL,
  term = NULL,
  effect = NULL,
  sims = 200L,
  alpha = 0.05,
  seed = 1L,
  progress = NULL,
  verbose = TRUE
)
```

## Arguments

- formula:

  A model formula, with random-effect bars if the design has repeated
  measures.

- design:

  Named list describing each variable in the design. A character vector
  becomes a factor with those levels; a numeric vector of length two or
  more becomes numeric values to cross; a function is called with the
  number of values needed and must return that many, which is the escape
  hatch for any distribution you like.

- n_unit:

  Integer vector of unit counts to trace power across.

- family:

  Family name or an
  [`ilm_family()`](https://huttoncp.github.io/illume/reference/ilm_family.md)
  object.

- coefs:

  Named numeric vector of coefficients on the link scale.

- cells:

  Expected cell means on the response scale: either a named numeric
  vector whose names are the factor levels joined by `"."`, or a data
  frame with one column per design factor and a column `mean`.

- sd:

  **Residual** standard deviation, for families that have one – not the
  standard deviation of the outcome. The random effects sit on top of
  it, so with `icc` the outcome's own spread is `sd / sqrt(1 - icc)`:
  `sd = 6, icc = 0.5` is data whose standard deviation is 8.49. The
  print method reports both. Superpower and faux take the **total**
  standard deviation, and so does a paper you might read one off, so
  converting one of those is `sd = sd_total * sqrt(1 - icc)`.

- re_sd:

  Random-effect standard deviations, named by grouping factor. For a
  term with a random slope, give a vector of standard deviations in the
  order the bar lists them; see `re_cor`.

- re_cor:

  Correlation between the random effects within a term, named by
  grouping factor. A single number for a two-column term, or a
  correlation matrix. Defaults to zero.

- icc:

  Intraclass correlation, as an alternative to `re_sd` for a random
  intercept: the random-effect variance becomes `icc / (1 - icc)` times
  the residual variance.

- within:

  Character vector of design variables that vary within a unit.

- contrasts:

  Passed to
  [`ilm_model()`](https://huttoncp.github.io/illume/reference/ilm_model.md);
  also used when solving `cells`, so the two always agree.

- term:

  The term to test, named as the variable or as the coefficient; see
  [`ilm_power()`](https://huttoncp.github.io/illume/reference/ilm_power.md).

- effect:

  Values for that coefficient on the link scale. Defaults to whatever
  the assumptions imply for it.

- sims, alpha, seed, progress:

  Passed to
  [`ilm_power()`](https://huttoncp.github.io/illume/reference/ilm_power.md).

- verbose:

  Logical. Report what was built.

## Value

An
[`ilm_power()`](https://huttoncp.github.io/illume/reference/ilm_power.md)
result, with an extra `n_unit` column.

## Details

For anything beyond a power curve – seeing what the assumptions imply,
simulating a data set, checking marginal means – build the scaffold with
[`ilm_scaffold()`](https://huttoncp.github.io/illume/reference/ilm_scaffold.md)
and use it directly.

## See also

[`ilm_scaffold()`](https://huttoncp.github.io/illume/reference/ilm_scaffold.md),
[`ilm_power()`](https://huttoncp.github.io/illume/reference/ilm_power.md),
[`ilm_power_n()`](https://huttoncp.github.io/illume/reference/ilm_power_n.md).

## Examples

``` r
# \donttest{
ilm_power_design(y ~ arm, design = list(arm = c("control", "treatment")),
                 n_unit = c(60, 120, 200),
                 cells = c(control = 12, treatment = 14.5), sd = 4,
                 term = "arm", sims = 50)
#> <ilm_scaffold>  parameters ASSUMED, not estimated
#>   formula    : y ~ arm
#>   family     : gaussian
#>   size       : 200 rows
#>   between    : arm
#> 
#>   assumed coefficients
#>  (Intercept) armtreatment 
#>         12.0          2.5 
#> 
#>   residual sd : 4
#> 
#>   The standard errors on this object come from ONE realisation of
#>   the design at this size. For power, use ilm_power() or
#>   ilm_power_design(), which refit many simulated studies.
#> Simulated power for armtreatment (gaussian family, alpha = 0.05)
#>   50 replicates per cell; the fitted study had 200 rows
#>  n_unit   n effect power mc_lower mc_upper power_converged converged
#>      60  60    2.5  0.62    0.482    0.741            0.62         1
#>     120 120    2.5  0.94    0.838    0.979            0.94         1
#>     200 200    2.5  1.00    0.929    1.000            1.00         1
#> 
#>   mc_lower/mc_upper is a Wilson interval on the power ESTIMATE: at
#>   0.80 from 50 replicates the standard error is 0.057, so a
#>   sample size read off this curve is a range. More replicates narrow it;
#>   nothing else does.
# }
```
