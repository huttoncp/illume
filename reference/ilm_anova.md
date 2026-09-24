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
  recode = TRUE
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

  `"Wald"` or `"LRT"`.

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

## Value

An `"anova"` data frame with one row per fixed-effect term. The columns
are `Df`, `Chisq` and `Pr(>Chisq)` in general, or `Df`, `F value` and
`Pr(>F)` when the model admits exact inference – a gaussian model with
no random or smooth terms, where the residual variance is estimated
rather than assumed known (see
[`ilm_model()`](https://huttoncp.github.io/illume/reference/ilm_model.md)).
A model with no terms to test, such as an intercept-only model, returns
a table with zero rows rather than an error: "there is nothing to test"
is an answer, not a failure.

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

Hauck, W. W., & Donner, A. (1977). Wald's test as applied to hypotheses
in logit analysis. *Journal of the American Statistical Association*,
72(360), 851–853.

## See also

[`ilm_pb_lrt()`](https://huttoncp.github.io/illume/reference/ilm_pb_lrt.md),
[`ilm_coef_table()`](https://huttoncp.github.io/illume/reference/ilm_coef_table.md).
