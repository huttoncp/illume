# Analysis of deviance for fixed effects

Tests each fixed-effect term jointly: all of its coefficients at once,
and for a multinomial outcome **across all category dimensions**.

## Usage

``` r
ilm_anova(
  object,
  type = 2,
  test = c("Wald", "LRT"),
  ncores = 1L,
  restarts = 2L,
  recode = TRUE,
  df = "auto",
  statistic = c("F", "Chisq")
)
```

## Arguments

- object:

  A fitted `"ilm_model"` object, fitted through the formula interface.

- type:

  `2` or `"II"` (the default), or `3` or `"III"`. Type II tests each
  term against everything not containing it, which does not depend on
  how the factors are coded. Type III tests each term against every
  other, which does – see `recode`.

- test:

  `"Wald"` or `"LRT"`. A likelihood-ratio test compares
  maximum-likelihood fits, so a REML fit is refitted by maximum
  likelihood for it, with a message; the fit passed in is unchanged. For
  a gaussian mixed or dispersion model the Wald F test needs no refit
  (see the section on unequal spreads).

- ncores:

  Integer. Worker processes for the refits.

- restarts:

  Integer. Optimiser restarts in refits.

- recode:

  When `type = 3` and a term inside an interaction is coded in a way
  that makes its main-effect row something other than a Type III test,
  refit with
  [`stats::contr.sum()`](https://rdrr.io/r/stats/contrast.html) for
  those factors and say so. `FALSE` tests the model exactly as coded and
  warns instead. The fit passed in is never modified either way.

- df:

  For the Wald tests of a gaussian mixed model, or of a gaussian model
  with a dispersion formula and no random terms, the denominator degrees
  of freedom of each F: `"auto"` (Satterthwaite's, as `lmerTest`
  reports), `"satterthwaite"`, `"kenward-roger"` (their df and scaled F,
  on their adjusted covariance; REML fits) or `"asymptotic"`. Other fits
  are unaffected. For a term of several columns, Satterthwaite's df
  combines the columns' own df by Fai and Cornelius (1996), which is
  defined only when each has more than 2 df; otherwise the term's
  `DenDF` and p-value are `NA` and the table's heading says why.
  Kenward-Roger's df, or testing fewer contrasts at once, may give one.

- statistic:

  `"F"` (the default) or `"Chisq"`. For the Wald tests of a gaussian
  mixed model, or of a gaussian dispersion model with no random terms,
  `"F"` gives an F test on the denominator df above and `"Chisq"` the
  large-sample chi-square the table gave before; with few groups the
  chi-square is too liberal. Every other fit keeps its own reference:
  exact F with nothing integrated out, chi-square otherwise.

## Value

An `"anova"` data frame with one row per fixed-effect term. The columns
are `Df`, `Chisq` and `Pr(>Chisq)` in general; `Df`, `F value` and
`Pr(>F)` when the model admits exact inference – a gaussian model with
no random or smooth terms, where the residual variance is estimated
rather than assumed known (see
[`ilm_model()`](https://huttoncp.github.io/illume/reference/ilm_model.md));
and `NumDF`, `DenDF`, `F value` and `Pr(>F)` for the Wald tests of a
gaussian mixed model, or of a gaussian dispersion model with no random
terms, with `statistic = "F"`. A model with no terms to test, such as an
intercept-only model, returns a table with zero rows rather than an
error: "there is nothing to test" is an answer, not a failure.

## Why this exists rather than car

:Anova: [`car::Anova()`](https://rdrr.io/pkg/car/man/Anova.html) reads
the `assign` attribute of the model matrix, which has one entry per
design column. A multinomial fit has `C` coefficients per column, so
`car` would match a short vector against a long one and silently test
only the first category, reporting `df = 1` where the correct joint test
has `df = C`. It returns a plausible table that answers the wrong
question. `illume` registers its own method so that
[`car::Anova()`](https://rdrr.io/pkg/car/man/Anova.html) dispatches here
and gets the right answer.

## Type II versus Type III

Type III tests each term with every other term in the model. Type II
tests each term after all terms that do **not** contain it, respecting
marginality. They agree when there are no interactions.

With interactions present, Type II is usually preferred for main
effects, and it has the practical advantage of not depending on how
factors are coded. Type III does depend on the coding and is only well
defined with orthogonal contrasts such as
[`stats::contr.sum()`](https://rdrr.io/r/stats/contrast.html).

## Unequal spreads, and Welch's test

A gaussian model with a spread per group, `dispformula = ~ group`, and
no random terms gets F tests on Satterthwaite's denominator df by
default, combined over a term's columns by Fai and Cornelius (1996).
Fitted by REML, which is that model's default, the F test of the group
term has Welch's (1951) numerator exactly. With two groups it is Welch's
t test: the same statistic, df and p-value as
`t.test(var.equal = FALSE)`. With three or more, Welch divides the
statistic by \\1 + 2(k - 2)\Lambda/(k^2 - 1)\\ and uses \\(k^2 -
1)/(3\Lambda)\\ denominator df, where \\\Lambda = \sum_i (1 -
w_i/W)^2/(n_i - 1)\\, \\w_i = n_i/s_i^2\\ and \\W = \sum_i w_i\\; the df
here come instead from Satterthwaite's approximation, which applies to
any design, not only one factor. So with three or more groups the F here
is somewhat larger than Welch's and its df differ, and
[`stats::oneway.test()`](https://rdrr.io/r/stats/oneway.test.html) gives
Welch's own.

In a simulation over 14 designs, under equal means at a nominal 5%, this
F rejected a median 5.8% of the time (4% to 10%), against Welch's 5.1%
in the one-way designs. With normal errors, Welch held 5% better with
four or six groups of six to eight observations (6.1% against 6.7%, and
5.2% against 6.6%); a one-degree-of-freedom term beside other terms, and
a spread modelled on a covariate in 20 observations, came out near 7% to
8%. Fitted by maximum likelihood the F is larger again, because the
group variances are divided by \\n_i\\ rather than \\n_i - 1\\, and too
liberal: a median 8.2%, up to 14.5%. `statistic = "Chisq"` gives the
large-sample chi-square, which rejected a median 11.5% on a
maximum-likelihood fit, up to 26%.

## Wald versus likelihood ratio

Wald tests are fast and need no refitting for Type III. Likelihood-ratio
tests refit a reduced model for each term, which is slower but avoids
the Hauck-Donner effect, where a Wald statistic can shrink for very
strong effects. If a term is important to your conclusions, confirm it
with the LRT; for small samples use
[`ilm_pb_lrt()`](https://huttoncp.github.io/illume/reference/ilm_pb_lrt.md)
instead, since both rely on large-sample approximations.

## References

Fox, J., & Weisberg, S. (2019). *An R Companion to Applied Regression*,
3rd ed. Sage. (Chapter 5 explains Type II and Type III tests.)

Fai, A. H.-T., & Cornelius, P. L. (1996). Approximate F-tests of
multiple degree of freedom hypotheses in generalized least squares
analyses of unbalanced split-plot experiments. *Journal of Statistical
Computation and Simulation*, 54(4), 363–378.

Hauck, W. W., & Donner, A. (1977). Wald's test as applied to hypotheses
in logit analysis. *Journal of the American Statistical Association*,
72(360), 851–853.

Satterthwaite, F. E. (1946). An approximate distribution of estimates of
variance components. *Biometrics Bulletin*, 2(6), 110–114.

Welch, B. L. (1951). On the comparison of several mean values: An
alternative approach. *Biometrika*, 38(3/4), 330–336.

## See also

[`ilm_pb_lrt()`](https://huttoncp.github.io/illume/reference/ilm_pb_lrt.md),
[`ilm_coef_table()`](https://huttoncp.github.io/illume/reference/ilm_coef_table.md).
