# A model built from assumptions instead of data

Builds a design grid from a study specification, attaches the parameters
you assume, and returns a fitted-shaped `"ilm_model"`. Use it to plan a
study: pass it to
[`ilm_power()`](https://huttoncp.github.io/illume/reference/ilm_power.md)
for a power curve, to
[`ilm_simulate()`](https://huttoncp.github.io/illume/reference/ilm_simulate.md)
for data a pipeline can be rehearsed on, or to
[`ilm_emmeans()`](https://huttoncp.github.io/illume/reference/ilm_emmeans.md)
and
[`ilm_trends()`](https://huttoncp.github.io/illume/reference/ilm_trends.md)
to see what your assumptions actually imply before you rely on them.

## Usage

``` r
ilm_scaffold(
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
  seed = 1L,
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

  Integer. Independent units: participants if the formula has a grouping
  bar, otherwise rows.

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

- seed:

  Integer seed for the grid's own realisation.

- verbose:

  Logical. Report what was built.

## Value

An `"ilm_model"` that also carries class `"ilm_scaffold"`.

## Saying what you assume

Give **either** `coefs` or `cells`, not both.

`coefs` is a named vector of regression coefficients, on the link scale,
named as the model matrix names them (`"armtreatment"`, not `"arm"`). It
is exact and works for every model, including continuous predictors and
interactions with them.

`cells` is the expected mean in each cell of the design, on the
**response** scale, which is how most people actually hold an
assumption: "controls average 12, the treated group averages 14.5".
illume solves back to coefficients. If the cell means you give cannot be
produced by the formula you gave – crossed means under an additive
formula, say – it says so and names the term that is missing, rather
than quietly fitting the closest thing it can.

## Saying how big the study is

`n_unit` counts **independent units**: participants when the formula has
a grouping bar, rows when it does not. A variable named in `within` is
crossed inside each unit; everything else is allocated across units as
evenly as the numbers allow. A variable that appears in a bar's
left-hand side – the `time` of `(1 + time | id)` – is within by
construction and does not need naming.

## What is NOT here

The standard errors on a scaffold come from one realisation of the
design at its own size. They are not a property of your assumptions, and
reading power off them would be reading one coin flip. That is what
[`ilm_power()`](https://huttoncp.github.io/illume/reference/ilm_power.md)
is for: it refits many simulated studies and counts.

## See also

[`ilm_power_design()`](https://huttoncp.github.io/illume/reference/ilm_power_design.md)
for the common case in one call;
[`ilm_power()`](https://huttoncp.github.io/illume/reference/ilm_power.md),
[`ilm_simulate()`](https://huttoncp.github.io/illume/reference/ilm_simulate.md),
[`ilm_emmeans()`](https://huttoncp.github.io/illume/reference/ilm_emmeans.md),
[`ilm_trends()`](https://huttoncp.github.io/illume/reference/ilm_trends.md).

## Examples

``` r
# \donttest{
## a parallel-arm trial, stated as cell means
s <- ilm_scaffold(y ~ arm, design = list(arm = c("control", "treatment")),
                  n_unit = 120, cells = c(control = 12, treatment = 14.5),
                  sd = 4, verbose = FALSE)
coef(s)
#>  (Intercept) armtreatment 
#>         12.0          2.5 

## a treatment-by-time design with repeated measures
s2 <- ilm_scaffold(y ~ arm * time + (1 | id),
                   design = list(arm = c("control", "treatment"),
                                 time = c("pre", "post")),
                   within = "time", n_unit = 60,
                   cells = c(control.pre = 12, control.post = 12.2,
                             treatment.pre = 12, treatment.post = 14.5),
                   sd = 4, icc = 0.5, verbose = FALSE)
ilm_emmeans(s2, c("arm", "time"))
#> <ilm_emm> marginal means of arm x time (link scale, equal weights)
#> 
#>        arm time estimate    se  lower upper
#>    control  pre     12.0 1.064  9.916 14.08
#>  treatment  pre     12.0 1.064  9.916 14.08
#>    control post     12.2 1.064 10.116 14.28
#>  treatment post     14.5 1.064 12.416 16.58
#>   Compare them with ilm_contrast().
# }
```
