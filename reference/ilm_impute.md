# Multiple imputation by chained equations

Fills in missing values `m` times, each fill drawn rather than fitted,
so that the analysis afterwards can price in not having observed them.
Each incomplete variable is regressed on the others with
[`ilm_model()`](https://huttoncp.github.io/illume/reference/ilm_model.md),
using the family its own type calls for, and the missing entries are
drawn from the resulting predictive distribution.

## Usage

``` r
ilm_impute(
  data,
  m = 20L,
  predictors = NULL,
  exclude = NULL,
  maxit = 5L,
  method = c("auto", "fcs", "lowrank", "glrm"),
  ncp = NULL,
  single = FALSE,
  seed = NULL,
  verbose = TRUE,
  progress = NULL
)
```

## Arguments

- data:

  A data frame.

- m:

  Number of imputations. 20 or more is cheap insurance; the classic
  advice of 5 dates from when it was not.

- predictors:

  Columns to use as predictors in the imputation models. Default is
  every column except those being imputed in that step.

- exclude:

  Columns never to impute and never to use, such as an identifier.

- maxit:

  Cycles through the variables per imputation.

- method:

  `"auto"` uses chained equations and falls back to the `"glrm"` fits a
  generalized low rank model instead, which uses a loss suited to each
  column's type and so can impute CATEGORICAL columns, which `"lowrank"`
  leaves alone – see
  [`illumex::ilm_glrm()`](https://huttoncp.github.io/illumex/reference/ilm_glrm.html).
  low-rank route when they cannot be fitted; `"fcs"` and `"lowrank"`
  force one. See the section below for what each costs.

- ncp:

  Rank for the low-rank route. `NULL` chooses it by cross-validation
  over held-out observed cells.

- single:

  Return a single completed data frame instead. Warns.

- seed:

  Random seed.

- verbose:

  Narrate progress.

- progress:

  Show a progress bar. Defaults to
  [`interactive()`](https://rdrr.io/r/base/interactive.html), so a bar
  appears when someone is watching and nothing is written in a script or
  a knitted document. See
  [illumex::ilm_progress_arg](https://huttoncp.github.io/illumex/reference/ilm_progress_arg.html).

## Value

An object of class `"ilm_mids"` holding the `m` completed data sets, or
a data frame when `single = TRUE`.

## Why several, and why drawn

Filling in one set of values and analysing it as though those values had
been observed treats a guess as data. The estimates may be fine; the
standard errors are not, because nothing in them knows that part of the
data set was invented. Intervals come out too narrow by an amount that
grows with the proportion missing.

Two things fix that. Each imputed value is **drawn** from the predictive
distribution – the coefficients are drawn from their sampling
distribution and the response's own randomness is added on top – so the
fills differ from one imputation to the next in the way the unknown
values really could. And the `m` completed data sets are analysed
separately and pooled by Rubin's rules in
[`ilm_mi_pool()`](https://huttoncp.github.io/illume/reference/ilm_mi_pool.md),
where the spread of the estimates across imputations becomes part of the
reported uncertainty.

`single = TRUE` returns one completed data set. It exists because it is
occasionally what is wanted – a plot, a rough look – and it warns,
because inference computed from it will be overconfident.

## What it assumes

That the data are missing at random given the variables supplied: the
chance a value is missing may depend on what is observed, but not on the
missing value itself once the observed variables are accounted for. That
assumption is **not testable**
([`illumex::ilm_check_missing()`](https://huttoncp.github.io/illumex/reference/ilm_check_missing.html)
explains why), and imputation does not make it true. Including variables
that predict both the missingness and the missing value makes it more
plausible.

Imputing an outcome from nothing but the predictors of the model you
intend to fit adds no information at all – the imputations lie on the
very surface you are about to estimate. Auxiliary variables are what
make imputing an outcome worthwhile.

## Which method, and what it costs

`method = "auto"` uses chained equations and falls back to the low-rank
route only when the per-variable regressions cannot be fitted – when
there are more predictors than rows on which a variable was observed.
That ordering is not a preference, it is what the ablation says. Hiding
known cells and scoring the imputations against them, with a noise floor
of 0.500:

      design                          mean-fill    fcs   lowrank
      n=200 p=8  rank 3                   2.289  1.177     1.963
      n=200 p=8  full rank                2.739  2.622     3.581
      n=400 p=12 rank 4, 30% missing      2.010  1.054     1.683
      n=60  p=80 rank 3  (p > n)          1.920      -     0.887

Where chained equations can be fitted it is clearly better: the
per-variable regressions use everything, where a rank-k fit discards
whatever falls outside those k components. Where `p > n` it cannot be
fitted at all, and the low-rank route gets within twice the noise floor.

The middle row is the warning. With no low-rank structure to find, the
reconstruction imposes one, and does **worse than filling in column
means**. `ilm_impute()` warns when cross-validation picks the largest
rank it was offered, which is the signature of that case.

Reconstruction accuracy is not the whole story either. It says how close
the filled values are, not whether inference computed afterwards is
calibrated – single imputation scores well on the first and fails the
second. The coverage table below is the one that matters for inference.

## What this is calibrated for

200 replicates, 400 rows, `y = 0.5x + 0.3z + e`, with `x` made missing
three ways. Coverage of a nominal 95% interval for the coefficient on
`x`, Monte Carlo error 0.015:

                           full   complete   multiple    single
                           data      cases  imputation  imputation
      MCAR, 30% missing
        bias            -0.0028    -0.0012     -0.0035    -0.0069
        coverage          0.955      0.960       0.970      0.890
      MAR on a covariate, 40% missing
        bias            -0.0028     0.0001     -0.0017     0.0002
        coverage          0.955      0.980       0.970      0.835
      MAR on the OUTCOME, 41% missing
        bias            -0.0028    -0.0989     -0.0090    -0.0058
        coverage          0.955      0.615       0.955      0.790

Three things to read off it. Complete cases are **unbiased** when
missingness depends on a covariate, even at 40% missing – which is why
[`illumex::ilm_check_missing()`](https://huttoncp.github.io/illumex/reference/ilm_check_missing.html)
distinguishes that case and tells you not to bother imputing. Complete
cases **fail badly** when missingness depends on the outcome: a bias of
-0.099 is a fifth of the effect, and coverage collapses to 0.615.
Multiple imputation repairs exactly that case, 0.955.

And single imputation is the cautionary column. Its point estimates are
respectable throughout – the bias is no worse than multiple imputation's
– but its coverage runs 0.790 to 0.890, because nothing in its standard
errors knows that part of the data was invented. That is the whole
reason `m` defaults to more than one.

Pooled coverage sits slightly high, 0.970 where 0.950 is nominal, in the
two cells where complete cases were already valid. Rubin's rules are
known to be mildly conservative, and erring wide is the right direction
for a method whose job is to stop intervals being too narrow.

## References

Rubin, D. B. (1987). Multiple Imputation for Nonresponse in Surveys.
Wiley.

van Buuren, S. and Groothuis-Oudshoorn, K. (2011). mice: Multivariate
Imputation by Chained Equations in R. Journal of Statistical Software
45(3).

Josse, J. and Husson, F. (2016). missMDA: A Package for Handling Missing
Values in Multivariate Data Analysis. Journal of Statistical Software
70(1). The low-rank route follows their multiple imputation PCA: a
regularised iterative fit, bootstrapped so the imputations differ, then
pooled.

## See also

[`ilm_mi_pool()`](https://huttoncp.github.io/illume/reference/ilm_mi_pool.md)
to analyse them,
[`illumex::ilm_check_missing()`](https://huttoncp.github.io/illumex/reference/ilm_check_missing.html)
to decide whether you need to.

## Examples

``` r
set.seed(1); n <- 200
d <- data.frame(x = rnorm(n), z = rnorm(n))
d$y <- 0.5 * d$x + 0.3 * d$z + rnorm(n)
d$x[sample(n, 40)] <- NA
imp <- ilm_impute(d, m = 5, seed = 1, verbose = FALSE)
imp
#> <ilm_mids> 5 imputation(s) of 200 rows
#>   imputed:
#>     x                  gaussian     40 value(s)
#> 
#>   Analyse with ilm_mi_pool(), which combines the 5 fits by Rubin's rules.
#>   Using one of these on its own would understate the uncertainty.
```
