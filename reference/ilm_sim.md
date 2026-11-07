# Simulated datasets from a model

Whole datasets drawn from a fitted model or an
[`ilm_scaffold()`](https://huttoncp.github.io/illume/reference/ilm_scaffold.md):
each row's own columns, with a simulated response in place of the
observed one. A shortcut for
[`ilm_simulate()`](https://huttoncp.github.io/illume/reference/ilm_simulate.md),
whose matrix of responses would otherwise have to be bound back to the
rows by hand.

## Usage

``` r
ilm_sim(object, nsim = 1L, seed = NULL, n = NULL, effect = NULL, term = NULL)
```

## Arguments

- object:

  A fitted
  [`ilm_model()`](https://huttoncp.github.io/illume/reference/ilm_model.md),
  or an
  [`ilm_scaffold()`](https://huttoncp.github.io/illume/reference/ilm_scaffold.md).

- nsim:

  Integer. Number of datasets.

- seed:

  Integer. Random seed, for reproducibility. With `NULL` the draws come
  from the current random stream, as
  [`ilm_simulate()`](https://huttoncp.github.io/illume/reference/ilm_simulate.md)'s
  do.

- n:

  The number of rows in each simulated study, as in
  [`ilm_power()`](https://huttoncp.github.io/illume/reference/ilm_power.md);
  `NULL` keeps the rows the model was fitted to.

- effect:

  The tested coefficient's value on the link scale, or a multiple of a
  term's coefficients, as in
  [`ilm_power()`](https://huttoncp.github.io/illume/reference/ilm_power.md);
  `NULL` keeps the fitted values.

- term:

  The coefficient or term `effect` sets, named as in
  [`ilm_power()`](https://huttoncp.github.io/illume/reference/ilm_power.md).

## Value

A data frame when `nsim = 1`, and otherwise a list of `nsim` data
frames.

## The rows

The columns are the ones the model reads: the response, the variables
its terms are built from, its grouping factors and weights, and the
variables of an offset, a dispersion model and a zero part. The rows are
the ones it was fitted to, so a row it dropped for a missing value is
not there. For a scaffold they are its design grid.

With neither `n` nor `effect`, the responses are
[`ilm_simulate()`](https://huttoncp.github.io/illume/reference/ilm_simulate.md)'s,
draw for draw with the same `seed`. The rows and their group labels are
the fitted ones, and each dataset draws new random effects for those
groups, so it is a new sample of groups under the old labels.

## A study of another size or effect

`n` draws a study of about that many rows, as
[`ilm_power()`](https://huttoncp.github.io/illume/reference/ilm_power.md)
draws one. A fitted model's rows are resampled, whole clusters at a time
when it has a grouping factor, and a scaffold's design is laid out again
at that size. A cluster drawn twice is two clusters: the unit is
renumbered by the copy drawn, 1, 2, and so on, a grouping factor nested
in it is prefixed with its copy (`"2.b"`), and one crossed with it, such
as items every participant sees, keeps its labels. A row with a
frequency weight is drawn as that many rows, each with weight 1.

`effect` draws the data with one coefficient, or one term, changed, as
in
[`ilm_power()`](https://huttoncp.github.io/illume/reference/ilm_power.md):
for a single coefficient, its value on the link scale; for a term of
several coefficients, a multiple of them, so `effect = 0.5` halves it.
`term` names it, and defaults to the first coefficient that is not an
intercept.

## The response

The simulated response goes in the response's own column, of the same
type: a factor stays a factor with the same levels, a logical stays
logical, and a whole-number count stays integer. The response has to be
a column the formula names, not a calculation on one such as `log(y)`;
for that,
[`ilm_simulate()`](https://huttoncp.github.io/illume/reference/ilm_simulate.md)
returns the response on the scale the model fits.

A censored response is censored as
[`ilm_simulate()`](https://huttoncp.github.io/illume/reference/ilm_simulate.md)
censors it, and a column beside it says what
[`ilm_model()`](https://huttoncp.github.io/illume/reference/ilm_model.md)'s
`censor` needs to fit it again:

- a survival time gets `.event`, `1` where the event was drawn and `0`
  where the time was cut at the subject's censoring time, as
  [`ilm_surv()`](https://huttoncp.github.io/illume/reference/ilm_surv.md)
  takes it: `ilm_surv(d$time, d$.event)`;

- codes given to
  [`ilm_censor()`](https://huttoncp.github.io/illume/reference/ilm_censor.md)
  as `status` are held fixed by the draw, and come back as `.censored`:
  `ilm_censor(d$y, status = d$.censored)`;

- a floor or ceiling needs no column, because the drawn values are held
  at it and
  [`ilm_censor()`](https://huttoncp.github.io/illume/reference/ilm_censor.md)
  finds them again: `ilm_censor(d$y, lower = 0)`.

## See also

[`ilm_simulate()`](https://huttoncp.github.io/illume/reference/ilm_simulate.md)
for the responses alone,
[`ilm_power()`](https://huttoncp.github.io/illume/reference/ilm_power.md)
to count detections over many simulated studies,
[`ilm_refit()`](https://huttoncp.github.io/illume/reference/ilm_refit.md)
to refit the same model, and
[`illumex::ilm_example_data()`](https://huttoncp.github.io/illumex/reference/ilm_example_data.html)
for example data of a fixed design.

## Examples

``` r
set.seed(1)
d <- data.frame(id = factor(rep(1:12, each = 6)), x = rnorm(72))
d$y <- 1 + 0.5 * d$x + rnorm(12)[d$id] + rnorm(72)
fit <- ilm_model(y ~ x + (1 | id), data = d, family = "gaussian",
                 verbose = FALSE)
head(ilm_sim(fit, seed = 2))
#>            y          x id
#> 1 -0.7448145 -0.6264538  1
#> 2 -0.9491267  0.1836433  1
#> 3  0.9619524 -0.8356286  1
#> 4 -1.4237852  1.5952808  1
#> 5  0.6902534  0.3295078  1
#> 6 -0.4705379 -0.8204684  1

## rehearse the analysis on a study of 120 rows, with half the slope
d2 <- ilm_sim(fit, n = 120, effect = 0.25, term = "x", seed = 3)
coef(ilm_model(y ~ x + (1 | id), data = d2, family = "gaussian",
               verbose = FALSE))
#> (Intercept)           x 
#>   0.5277337   0.2235750 

## from a scaffold: three draws of a planned two-arm trial
s <- ilm_scaffold(y ~ arm, design = list(arm = c("control", "treatment")),
                  n_unit = 60, cells = c(control = 12, treatment = 14.5),
                  sd = 4, verbose = FALSE)
sims <- ilm_sim(s, nsim = 3, seed = 4)
sapply(sims, function(d) diff(tapply(d$y, d$arm, mean)))
#> treatment treatment treatment 
#>  1.330933  2.165573  1.838666 
```
