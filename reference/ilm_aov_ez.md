# Factorial and repeated-measures ANOVA

Specify the design by naming columns – the participant, the outcome,
which factors vary between participants and which vary within them –
rather than by writing a formula with an error term in it. Reports F on
a denominator degrees of freedom, mean squared error, generalized eta
squared, and a sphericity correction where one is needed.

## Usage

``` r
ilm_aov_ez(
  id,
  dv,
  data,
  between = NULL,
  within = NULL,
  covariate = NULL,
  observed = NULL,
  engine = c("aov", "mixed"),
  reml = TRUE,
  fun_aggregate = NULL,
  posthoc = "auto",
  alpha = 0.05,
  correction = c("auto", "GG", "HF", "none"),
  verbose = TRUE
)
```

## Arguments

- id:

  Column identifying the participant.

- dv:

  Column holding the numeric outcome.

- data:

  A long data frame: one row per participant per occasion.

- between:

  Factors that vary BETWEEN participants, as a character vector.

- within:

  Factors that vary WITHIN participants.

- covariate:

  Numeric covariates. Must be constant within participant.

- observed:

  Variables measured rather than manipulated, affecting the `ges`
  denominator. Covariates are included automatically.

- engine:

  `"aov"` (the default) computes the omnibus table classically and fits
  a mixed model alongside. `"mixed"` makes the mixed model the primary
  analysis, translating this specification into the corresponding
  [`ilm_model()`](https://huttoncp.github.io/illume/reference/ilm_model.md)
  formula.

- reml:

  Passed to
  [`ilm_model()`](https://huttoncp.github.io/illume/reference/ilm_model.md)
  for the companion fit. `TRUE` here, because the design fixed the fixed
  effects before any data were seen.

- fun_aggregate:

  Applied when a participant has several rows per cell. Defaults to the
  mean, with a message saying it happened.

- posthoc:

  Follow a significant effect with pairwise comparisons. `"auto"` does
  it for significant effects involving a factor with three or more
  levels; `FALSE` never does; a character vector names effects.

- alpha:

  Threshold for deciding what `posthoc` follows up, and for flagging
  sphericity.

- correction:

  Sphericity correction applied to the reported p-value: `"auto"`
  (Greenhouse-Geisser when Mauchly's test rejects), `"GG"`, `"HF"`, or
  `"none"`.

- verbose:

  Narrate the checks.

## Value

An object of class `"ilm_aov_ez"`: `anova` (the table), `sphericity`,
`fit` (the companion
[`ilm_model()`](https://huttoncp.github.io/illume/reference/ilm_model.md)),
`posthoc`, `notes`, and the design.

## What it reports, and why not just a mixed model

On a balanced complete design a mixed model and a classical
repeated-measures ANOVA answer the same question. They report it
differently, and one of them reports things the other cannot:
[`ilm_model()`](https://huttoncp.github.io/illume/reference/ilm_model.md)
gives a Wald chi-square, because once random effects are integrated out
there is no exact residual degrees of freedom to divide by. An F, a mean
squared error, a generalized eta squared and a Greenhouse-Geisser
correction do not fall out of that fit. So the omnibus table here is
computed classically, and a mixed model is fitted alongside for the
things the classical route cannot do – marginal means, contrasts, simple
slopes, diagnostics, and an analysis that keeps participants with
incomplete data rather than dropping them.

## When the classical answer is the wrong one

Three situations make it the wrong tool, and each is checked and
reported with the alternative named rather than left to be noticed:

- **Sphericity fails.** The within-participant variances and covariances
  are not what the F test assumes. Greenhouse-Geisser and Huynh-Feldt
  corrections are reported, and so is the other option: an unstructured
  mixed model does not assume sphericity, so it needs no correction at
  all.

- **Participants have missing cells.** Classical repeated measures has
  to drop them completely, because a participant missing one occasion
  contributes to no within-participant contrast. How many were dropped
  is reported. A mixed model uses them.

- **A covariate varies within participant.** Refused, because a single
  coefficient for a time-varying covariate blends two effects that can
  point in opposite directions – how participants who score higher on
  average differ, and what happens when a participant scores higher than
  usual.

## Generalized eta squared

`ges` is comparable across designs in a way that partial eta squared is
not, which is why it is the default here and in `afex`. Its denominator
includes the sums of squares of any **measured** rather than manipulated
variable. Covariates are measured by definition, so they are added to
`observed` automatically; name any measured factors there too.

## See also

[`ilm_model()`](https://huttoncp.github.io/illume/reference/ilm_model.md)
for the general engine,
[`ilm_emmeans()`](https://huttoncp.github.io/illume/reference/ilm_emmeans.md)
and
[`ilm_contrast()`](https://huttoncp.github.io/illume/reference/ilm_contrast.md)
for follow-ups,
[`ilm_trends()`](https://huttoncp.github.io/illume/reference/ilm_trends.md)
when the interaction is with a continuous predictor.

## Examples

``` r
set.seed(1)
n <- 24
d <- expand.grid(id = factor(seq_len(n)), time = factor(1:3))
d$grp <- factor(rep(c("ctl", "trt"), each = n / 2))[as.integer(d$id)]
d$score <- 10 + 2 * (d$grp == "trt") * as.integer(d$time) +
  rnorm(nrow(d)) + rep(rnorm(n), 3)
a <- ilm_aov_ez("id", "score", d, between = "grp", within = "time",
                verbose = FALSE)
a
#> Analysis of Variance
#>   score ~ grp (between) + time (within)
#>   24 participants
#> 
#>    Effect    df   MSE       F   ges      p    
#>       grp 1, 22 2.743 159.353 0.813 <1e-04 ***
#>      time 2, 44 0.915  29.761 0.351 <1e-04 ***
#>  grp:time 2, 44 0.915  32.088 0.369 <1e-04 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Sphericity (Mauchly):
#>    Effect     W     p    GG HF
#>      time 0.997 0.968 0.997  1
#>  grp:time 0.997 0.968 0.997  1
#> 
#> Pairwise comparisons:
#> 
#> -- time --
#> <ilm_contrast> 3 comparison(s), adjust = max_t 
#> 
#>  contrast estimate     se  lower upper   p_adj
#>     2 - 1   0.8775 0.2762 0.2310 1.524 0.00385
#>     3 - 1   2.1202 0.2762 1.4736 2.767 0.00000
#>     3 - 2   1.2427 0.2762 0.5961 1.889 0.00000
#> 
#>   Intervals hold jointly at 95% across all 3 comparisons.
#> 
#> -- grp:time --
#>   ~ grp within each time ~
#>     [time = 1]
#> <ilm_contrast> 1 comparison(s), adjust = none 
#> 
#>   contrast estimate    se lower upper     p_adj
#>  trt - ctl     2.62 0.504 1.632 3.608 2.022e-07
#> 
#>     [time = 2]
#> <ilm_contrast> 1 comparison(s), adjust = none 
#> 
#>   contrast estimate    se lower upper     p_adj
#>  trt - ctl    5.133 0.504 4.145 6.121 2.353e-24
#> 
#>     [time = 3]
#> <ilm_contrast> 1 comparison(s), adjust = none 
#> 
#>   contrast estimate    se lower upper     p_adj
#>  trt - ctl     7.03 0.504 6.042 8.018 3.264e-44
#> 
#>   ~ time within each grp ~
#>     [grp = ctl]
#> <ilm_contrast> 3 comparison(s), adjust = max_t 
#> 
#>  contrast estimate     se   lower  upper  p_adj
#>     2 - 1 -0.37908 0.3906 -1.2934 0.5353 0.5988
#>     3 - 1 -0.08503 0.3906 -0.9994 0.8293 0.9748
#>     3 - 2  0.29406 0.3906 -0.6203 1.2084 0.7349
#> 
#>   Intervals hold jointly at 95% across all 3 comparisons.
#> 
#>     [grp = trt]
#> <ilm_contrast> 3 comparison(s), adjust = max_t 
#> 
#>  contrast estimate     se lower upper p_adj
#>     2 - 1    2.134 0.3906 1.220 3.048     0
#>     3 - 1    4.325 0.3906 3.411 5.240     0
#>     3 - 2    2.191 0.3906 1.277 3.106     0
#> 
#>   Intervals hold jointly at 95% across all 3 comparisons.
#> 
#> 
#> Notes:
#>   * pairwise comparisons computed for: time, grp:time. These come from the    mixed model, so they use all the data and a single pooled error; they    are not the classical stratum-specific error terms.
```
