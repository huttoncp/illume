# The design matrices of a fitted model, for new rows

Everything needed to assemble a prediction for new rows from the fit's
parameters or from draws of them
([`ilm_draws()`](https://huttoncp.github.io/illume/reference/ilm_draws.md)):
the fixed design, each random term's design and the group each row
belongs to, each smooth's penalised basis, the zero part's and the
dispersion model's designs, and each row's place among the cells of a
correlation over time.

## Usage

``` r
ilm_matrices(object, newdata, time = NULL, group = NULL)
```

## Arguments

- object:

  A fitted `"ilm_model"`.

- newdata:

  A data frame of new rows, with the columns the model uses.

- time, group:

  For a correlation over time built from vectors, the new rows' times
  and groups. Not needed when it was built by name.

## Value

A list with

- `X`:

  the fixed design, columns as `colnames(object$X)`.

- `re`:

  one element per grouping term: `Z`, its design (the intercept and
  slope columns, named by dimension as
  [`ilm_ranef()`](https://huttoncp.github.io/illume/reference/ilm_ranef.md)
  and
  [`ilm_draws()`](https://huttoncp.github.io/illume/reference/ilm_draws.md)'
  map name them, `"(Intercept)"` for a random intercept alone); `level`,
  each row's group label; `group`, its position among the fitted levels,
  `NA` for a new group; `new_group`; and `factor`, the grouping
  variable.

- `smooth`:

  one penalised basis per smooth term.

- `zi`, `disp`:

  the zero part's and the dispersion model's designs, when the model has
  them. Their columns are named as the coefficients they multiply are,
  in `coef(object, full = TRUE)` and in
  [`ilm_draws()`](https://huttoncp.github.io/illume/reference/ilm_draws.md)'
  map: `"zi:(Intercept)"`, `"disp:x"`.

- `ar`:

  with a correlation over time, a data frame with a row per new row:
  `group`, `time`, `cell`, `prev_cell`, `next_cell`, `dt_prev`,
  `dt_next` and `new_group`.

## Details

**Groups.** A row's group is matched to the fitted levels by label; a
group the fit has not seen has no code and is marked `new_group`.

**Over time.** `ar` places each row among its group's cells, as numbered
by
[`ilm_cells()`](https://huttoncp.github.io/illume/reference/ilm_cells.md):
`cell` if the row's time is one of them, and the nearest cells before
and after, `prev_cell` and `next_cell`, with the time from each,
`dt_prev` and `dt_next`, in the units of the fit's time. A time past a
group's last cell has a `prev_cell` and no `next_cell` – the forecast
case. An AR(1) term's grid is fixed, so a time that falls between its
steps is an error rather than rounded to one. The time and the group are
read from the columns the term was built from, `ilm_rw1(~ time | group)`
and the like; for a term built from vectors, pass them as `time` and
`group`.

## See also

[`ilm_draws()`](https://huttoncp.github.io/illume/reference/ilm_draws.md),
[`ilm_cells()`](https://huttoncp.github.io/illume/reference/ilm_cells.md),
[`ilm_ranef()`](https://huttoncp.github.io/illume/reference/ilm_ranef.md).

## Examples

``` r
set.seed(1)
d <- data.frame(id = factor(rep(c("a", "b", "c"), each = 10)), t = rep(1:10, 3))
d$y <- rnorm(30) + rep(cumsum(rnorm(10, 0, 0.5)), 3)
fit <- ilm_model(y ~ 1, data = d, family = "gaussian",
                 ar = ilm_rw1(~ t | id), verbose = FALSE)
#> ilm_model(): the random-effect covariance of `ar` sits at the edge of its range -- a variance of zero or a correlation of +/-1 -- where the data cannot resolve it. The fixed effects and their standard errors are still usable; summary() says what else is. If the term belongs in the model, boundary = "avoid" keeps it inside its range with a small penalty: it is then assumed nonzero rather than estimated at zero, so do not test whether it is; its variance comes out larger, and for a binary or categorical outcome the fixed effects a little further from zero -- markedly so when a category is rare.
nd <- data.frame(t = c(5, 12, 3), id = c("a", "a", "z"))
ilm_matrices(fit, nd)$ar
#>   group time cell prev_cell next_cell dt_prev dt_next new_group
#> 1     a    5    5         4         6       1       1     FALSE
#> 2     a   12   NA        10        NA       2      NA     FALSE
#> 3     z    3   NA        NA        NA      NA      NA      TRUE
```
