# A low-rank model with a loss chosen per column

The same idea as principal components – describe every row with a
handful of numbers – but with a loss appropriate to each column's type
rather than squared error everywhere. A binary column gets logistic
loss, a categorical one multinomial across its levels, a count Poisson,
and a numeric one squared error, which is where this reduces to ordinary
PCA.

## Usage

``` r
ilm_glrm(
  data,
  cols = NULL,
  rank = 2L,
  loss = NULL,
  lambda = NULL,
  weights = NULL,
  maxit = 300L,
  tol = 1e-07,
  seed = 1L,
  progress = NULL
)
```

## Arguments

- data:

  A data frame.

- cols:

  Columns to use; see
  [ilm_selection](https://craig-hutton.github.io/illume/reference/ilm_selection.md).

- rank:

  Number of dimensions.

- loss:

  Optional named character vector overriding the automatic choice per
  column: `"quadratic"`, `"logistic"`, `"multinomial"` or `"poisson"`.
  An **ordered** factor defaults to `"multinomial"`, which ignores the
  ordering but never invents a spacing; `"quadratic"` on its integer
  codes uses the ordering and does invent one, and is available if that
  is the trade you want.

- lambda:

  Ridge penalty on the factors, or `NULL` (the default) to choose it by
  holding out observed cells. It is not a detail: a row with few
  observed cells has as many scores as observations and fits them
  exactly, so without enough penalty the reconstruction of everything
  else in that row goes wherever the algebra sends it. On a 400 by 4
  matrix with 20% missing, held-out error ran 0.82 at `lambda = 0.1`,
  0.68 at 2, and 1.16 at 25 – against 0.78 for an iterative SVD and 1.18
  for column means. The best value depends on the size, the width and
  the missingness, which is why it is chosen rather than fixed.

- weights:

  Optional non-negative row weights. A row counted twice contributes
  twice to the loss, which is how
  [`ilm_impute()`](https://craig-hutton.github.io/illume/reference/ilm_impute.md)
  draws a bootstrap replicate without resampling the rows themselves and
  losing the ones it has to reconstruct.

- maxit:

  Maximum alternating passes.

- tol:

  Relative change in the objective at which to stop.

- seed:

  Random seed for the starting point.

- progress:

  Show a progress bar; see
  [ilm_progress_arg](https://craig-hutton.github.io/illume/reference/ilm_progress_arg.md).

## Value

An object of class `"ilm_glrm"` with `scores` (one row per observation),
`archetypes`, the per-column `loss`, the `objective` trace and a
`fitted` data frame reconstructing each column on its own scale.

## When this is worth the iterative fit

On all-numeric data it is not: with quadratic loss everywhere it *is*
PCA, and
[`ilm_reduce()`](https://craig-hutton.github.io/illume/reference/ilm_reduce.md)
gets there in closed form. It earns its cost when categorical columns
matter, because FAMD reaches them by one-hot encoding and applying
squared loss to the indicators, which is a Gaussian approximation to
something that is not Gaussian, and can reconstruct a category as -0.3.

## Missing values are not a special case

A cell that was not observed contributes nothing to the loss, so the fit
uses whatever is there and the reconstruction fills the rest in. That is
why the same machinery serves
[`ilm_impute()`](https://craig-hutton.github.io/illume/reference/ilm_impute.md).

## References

Udell, M., Horn, C., Zadeh, R. and Boyd, S. (2016). Generalized low rank
models. *Foundations and Trends in Machine Learning* 9, 1-118.

## See also

[`ilm_reduce()`](https://craig-hutton.github.io/illume/reference/ilm_reduce.md)
for the closed-form FAMD route and [ilm_reduce(method =
"glrm")](https://craig-hutton.github.io/illume/reference/ilm_reduce.md)
to use this one inside it.

## Examples

``` r
set.seed(1); n <- 200
f <- rnorm(n)
d <- data.frame(a = f + rnorm(n, 0, .4), b = -f + rnorm(n, 0, .4),
                g = factor(ifelse(f + rnorm(n, 0, .5) > 0, "hi", "lo")))
g <- ilm_glrm(d, rank = 1)
head(g$scores)
#>         dim1
#> 1 -0.3688959
#> 2  0.1534611
#> 3 -0.2614024
#> 4  0.5985325
#> 5  0.1067899
#> 6 -0.2188636
```
