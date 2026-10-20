# Random effects of a fitted model

The predicted random effects – the conditional modes – labelled by the
groups they belong to: one row per level of each grouping factor and
each of its coefficients, one per basis function of a smooth, and one
per cell of a correlation over time.

## Usage

``` r
ilm_ranef(object, ...)

# S3 method for class 'ilm_model'
ranef(object, ...)

# S3 method for class 'ilm_ranef'
as.data.frame(x, row.names = NULL, optional = FALSE, lme4 = FALSE, ...)
```

## Arguments

- object:

  A fitted `"ilm_model"`.

- ...:

  Unused.

- x:

  An `"ilm_ranef"` object.

- row.names, optional:

  Unused; for the generic.

- lme4:

  Logical. Give `lme4`'s long-format names – `grpvar`, `term`, `grp`,
  `condval`, `condsd` – in place of illume's own.

## Value

A data frame of class `"ilm_ranef"`, with columns

- `type`:

  `"re"`, `"smooth"`, or the correlation over time: `"ar1"`, `"car1"`,
  `"rw1"`.

- `term`:

  the term's name in the fit, as in `object$Sigma`.

- `factor`:

  the grouping variable.

- `level`:

  the group, a factor with the fitted levels; a smooth's basis functions
  are `b1`, `b2`, ...

- `dim`:

  the coefficient – `"(Intercept)"`, a slope's variable – prefixed by
  the category for a multinomial outcome.

- `time`, `cell`:

  for a correlation over time, the cell's time and its row of
  [`ilm_cells()`](https://huttoncp.github.io/illume/reference/ilm_cells.md);
  `NA` otherwise.

- `row`:

  the position in `object$obj$env$par`.

- `mode`, `sd`:

  the conditional mode and its standard deviation.

## Details

`sd` is each mode's conditional standard deviation with the variance
components at their estimates, from the inner Hessian of the Laplace
approximation: the quantity `lme4` calls `condVar`. With `reml = TRUE`
the fixed effects are integrated out alongside the random ones, so there
it also carries their uncertainty. A fit read back from disk no longer
has the compiled objective, and gives `NA`.

`row` is each value's position in the fit's full parameter vector,
`object$obj$env$par`, which is also the order of its joint precision.
Under REML the fixed effects sit in that vector's random block too; they
are not random effects, and are not listed.

For a reduced-rank term the rows are its latent factor scores;
multiplied by the loadings, `object$Lambda`, they give the effect on
each category. A random walk's first cell in each group is held at zero
rather than estimated, so it has no `row`, and a mode and `sd` of zero.

With `nlme` or `lme4` attached, `ranef(fit)` gives the same.

## See also

[`ilm_varcorr()`](https://huttoncp.github.io/illume/reference/ilm_varcorr.md)
for the variance components,
[`ilm_cells()`](https://huttoncp.github.io/illume/reference/ilm_cells.md).

## Examples

``` r
set.seed(1)
d <- data.frame(id = factor(rep(sprintf("s%02d", 1:10), each = 6)),
                x = rnorm(60))
d$y <- 1 + 0.5 * d$x + rnorm(10)[d$id] + rnorm(60)
fit <- ilm_model(y ~ x + (1 | id), data = d, family = "gaussian",
                 verbose = FALSE)
head(ilm_ranef(fit))
#> <ilm_ranef> 6 predicted random effects (conditional modes)
#>   re: id: 6 
#> 
#>  type term factor level         dim row     mode     sd
#>    re   id     id   s01 (Intercept)   4  1.53000 0.3115
#>    re   id     id   s02 (Intercept)   5 -0.76030 0.3115
#>    re   id     id   s03 (Intercept)   6  0.40380 0.3115
#>    re   id     id   s04 (Intercept)   7  0.05781 0.3115
#>    re   id     id   s05 (Intercept)   8 -1.32600 0.3115
#>    re   id     id   s06 (Intercept)   9 -0.29910 0.3115
```
