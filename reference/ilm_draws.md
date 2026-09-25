# Joint draws of every parameter of a fitted model

Draws the fixed effects, the variance and dispersion parameters, the
random effects and the cells of a correlation over time together, from
the Gaussian approximation to their joint distribution that the fit's
standard errors come from.

## Usage

``` r
ilm_draws(
  object,
  nsim = 1000L,
  seed = NULL,
  given = c("none", "theta"),
  blocks = NULL,
  natural = TRUE
)
```

## Arguments

- object:

  A fitted `"ilm_model"`.

- nsim:

  Number of draws.

- seed:

  Optional random seed.

- given:

  `"none"` to draw everything, or `"theta"` to hold the variance
  parameters at their estimates.

- blocks:

  Optional character vector of the blocks to return, named as in
  `names(object$obj$env$par)`: `"beta"`, `"theta"`, `"bvec"`,
  `"logdisp"`, `"B_ar"`, `"lchol_ar"`, `"rho_raw"` and so on.

- natural:

  Logical. Also transform each draw's variance parameters to the natural
  scale, as
  [`ilm_varcorr()`](https://huttoncp.github.io/illume/reference/ilm_varcorr.md)
  reports them.

## Value

An object of class `"ilm_draws"`, a list with

- `draws`:

  a matrix, one row per parameter and one column per draw, on the fit's
  internal scale (log standard deviations and the like), with rows named
  by block.

- `mode`:

  the centre: the estimates and conditional modes.

- `map`:

  a data frame saying what each row is: `row`, `block`, `term` (a
  coefficient's or parameter's name, or a random term's), and for a
  random effect its `level` (the group's label), `dim` and, for a
  correlation over time, `cell` (the row of
  [`ilm_cells()`](https://huttoncp.github.io/illume/reference/ilm_cells.md)).

- `natural`:

  with `natural = TRUE`, each draw's variance parameters on their
  natural scale: `re`, one array per random term of its covariance by
  draw; `latent`, a correlation over time's covariance by draw with its
  `rho` (and `range`) or, for a random walk, `var_per_time`; and
  `dispersion`, one row per dispersion parameter. At the estimate these
  are exactly
  [`ilm_varcorr()`](https://huttoncp.github.io/illume/reference/ilm_varcorr.md)'s.

- `held`:

  the number of directions held at a boundary, and the terms they belong
  to.

- `given`:

  as asked.

## Details

**The distribution.** The draws are centred on the estimates and the
conditional modes, with covariance the inverse of the joint precision –
the Laplace approximation's curvature in every parameter at once. For
the fixed parameters that is exactly `vcov(fit, full = TRUE)`; for the
random effects it adds the uncertainty the fixed parameters bring to the
conditional SDs of
[`ilm_ranef()`](https://huttoncp.github.io/illume/reference/ilm_ranef.md).

**At a boundary.** When a random term's covariance sits at the edge of
its range, the fit holds the direction the data cannot resolve at its
estimate (see `fit$hessian_held`), and so do the draws: exactly, by
conditioning the joint distribution on those directions. `held` says how
many.

**`given = "theta"`** holds the variance parameters of the random terms
and of a correlation over time at their estimates, and draws everything
else from its distribution given them – the approach of the `merTools`
package. The dispersion and the family's other parameters are still
drawn.

**`blocks`** returns only the named blocks of the parameter vector – say
`c("beta", "bvec")` – which saves memory; the draw is joint either way.

A fit made without `joint = TRUE` has no joint precision stored, and it
is formed here from the fit's compiled objective; a fit read back from
disk no longer has that, and has to be refitted.

## See also

[`ilm_ranef()`](https://huttoncp.github.io/illume/reference/ilm_ranef.md),
whose `row` indexes these draws' full layout;
[`ilm_varcorr()`](https://huttoncp.github.io/illume/reference/ilm_varcorr.md);
[`ilm_cells()`](https://huttoncp.github.io/illume/reference/ilm_cells.md).

## Examples

``` r
set.seed(1)
d <- data.frame(id = factor(rep(1:12, each = 6)), x = rnorm(72))
d$y <- 1 + 0.5 * d$x + rnorm(12)[d$id] + rnorm(72)
fit <- ilm_model(y ~ x + (1 | id), data = d, family = "gaussian",
                 verbose = FALSE)
dr <- ilm_draws(fit, nsim = 200, seed = 1)
dim(dr$draws)
#> [1]  16 200
head(dr$map)
#>   row block        term level         dim cell
#> 1   1  beta (Intercept)  <NA>        <NA>   NA
#> 2   2  beta           x  <NA>        <NA>   NA
#> 3   3 theta   id:L[1,1]  <NA>        <NA>   NA
#> 4   4  bvec          id     1 (Intercept)   NA
#> 5   5  bvec          id     2 (Intercept)   NA
#> 6   6  bvec          id     3 (Intercept)   NA
```
