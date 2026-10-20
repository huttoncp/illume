# Moderation: does the effect hold for everyone?

A treatment effect is an average. The question underneath it — did it
work the same way for everyone — is a different one, and asking it badly
is one of the most reliable ways to publish something that will not
replicate.

This vignette is about asking it well.

``` r

library(illume)
```

## The trap, with a number on it

Suppose you fit a model, find an effect, and then wonder whether it
differs by site, by age, by baseline severity. You test each one and
report the strongest.

That is not three tests. It is **a search**, and the p-value from a
search is not the p-value from a test. Measured directly: six candidate
moderators, a model with no moderation in it whatsoever, 150
replications.

      procedure                                   rejects at a nominal 5%
      pick the strongest and report it                     0.293
      the same, Holm-adjusted                              0.027
      the same, Benjamini-Hochberg                         0.027
      the same, Bonferroni                                 0.027
      honest split: choose on half, test on the other      0.053

**29.3%.** Nearly six times the rate you think you are running at. Not
because anyone did anything dishonest — that is simply what “I looked at
several and reported the interesting one” costs.

[`ilm_moderation()`](https://huttoncp.github.io/illume/reference/ilm_moderation.md)
exists to make that arithmetic visible and to charge for it.

## A worked example

Four hundred and eighty people, two arms, three sites, with age and a
baseline measure as covariates.

``` r

fit <- ilm_model(outcome ~ arm + site + age + baseline, data = trial)
ilm_moderation(fit, x = "arm")
#> Moderation of `arm`: 3 candidates, holm-adjusted
#>  moderator df statistic        p    p_adj
#>       site  2    15.544 2.90e-07 8.70e-07
#>   baseline  1     1.434 2.32e-01 4.64e-01
#>        age  1     0.009 9.26e-01 9.26e-01
#>
#>   `site`: the effect of `arm` is not the same for everyone.
```

`x = "arm"` is required and has no default. `arm:site` is the same term
whether you call `site` a moderator of `arm` or the reverse — the
algebra is symmetric and only the interpretation is not, so only you
know which is which. Passing an \[ilm_dag_model()\] is the one
exception: the graph has already named the exposure.

Three things in that output are doing work.

**`df` is what the candidate cost.** `site` has three levels, so its
interaction block is 2 df and the test is the **joint** test of that
block, not one coefficient. This matters more than it looks: on a
four-level moderator the three separate coefficients once gave p = 0.62,
0.0018 and 0.00002 for the same variable. There is one answer to “does
`site` moderate this”, and
[`ilm_anova()`](https://huttoncp.github.io/illume/reference/ilm_anova.md)
gives it.

**`p_adj` is what you report.** The raw `p` is there so you can see the
adjustment happening, not so you can quote it.

**The candidates were chosen for you**, defaulting to the model’s other
fixed effects — variables you had already judged worth adjusting for.
That keeps the multiplicity burden honest. Name `moderators` explicitly
to reach anything else in the data.

## Seeing it

``` r

m <- ilm_moderation(fit, x = "arm")
ilm_plot_moderation(m)
#>       xv        y       lo       hi        by
#>  central 50.60605 49.30866 51.90343   control
#>  central 52.66587 51.53020 53.80153 treatment
#>    north 49.24607 48.02000 50.47214   control
#>    north 52.15895 50.83507 53.48282 treatment
#>    south 50.22610 48.91922 51.53297   control
#>    south 58.92658 57.64629 60.20686 treatment
```

Two points of separation at central and north; nearly nine at south.
That is the moderation, and it is the kind of thing that a single
averaged treatment effect of 4.4 hides completely.

This is deliberately **not** the effect plot from
[`ilm_plot_model()`](https://huttoncp.github.io/illume/reference/ilm_plot_model.md).
That one varies a predictor with the others held at typical values,
which is right for “what does the model say about `arm`, all else equal”
and exactly wrong here: pinning `site` at its modal level is pinning the
thing you are trying to look at.

For a **continuous** exposure the picture changes shape, because the
quantity said to differ is a slope rather than a mean:

``` r

ilm_plot_moderation(m_age)        # slope of age within each site, with a
                                  # dashed line at zero
```

That comes from
[`ilm_trends()`](https://huttoncp.github.io/illume/reference/ilm_trends.md),
so the picture and the p-value are the same estimator on the same fit
rather than two code paths that happen to agree.

## What you do next

A moderation that survives is a **hypothesis**, not a finding. The
honest follow-up is to put it in the model and describe it properly:

``` r

fit2 <- ilm_model(outcome ~ arm * site + age + baseline, data = trial)
ilm_anova(fit2)
#>          Df F value    Pr(>F)
#> arm       1  68.576 1.268e-15 ***
#> site      2  18.898 1.273e-08 ***
#> age       1  11.577 0.0007245 ***
#> baseline  1  97.731 < 2.2e-16 ***
#> arm:site  2  15.544 2.900e-07 ***

ilm_emmeans(fit2, c("arm", "site"))
#>        arm    site estimate        se    lower    upper
#>    control central 50.60605 0.6602471 49.30866 51.90343
#>  treatment central 52.66587 0.5779471 51.53020 53.80153
#>    control   north 49.24607 0.6239534 48.02000 50.47214
#>  treatment   north 52.15895 0.6737258 50.83507 53.48282
#>    control   south 50.22610 0.6650750 48.91922 51.53297
#>  treatment   south 58.92658 0.6515423 57.64629 60.20686
```

Note that `arm:site` gives **15.544 on 2 df** in both places. It is the
same test;
[`ilm_moderation()`](https://huttoncp.github.io/illume/reference/ilm_moderation.md)
ran it three times and charged for the looking.

## A moderator you predicted in advance

If the interaction is already in the model, you hypothesised it a priori
— and it owes no multiplicity penalty, because you did not go searching.
Charging it one would make a pre-registered hypothesis *weaker* for
having been tested alongside exploratory ones.

[`ilm_moderation()`](https://huttoncp.github.io/illume/reference/ilm_moderation.md)
therefore takes it out of the search and says so:

``` r

ilm_moderation(fit2, x = "arm")
#> ilm_moderation(): `arm:site` is already in the model, so it was specified a
#> priori and is not charged a multiplicity penalty here. Test it with
#> ilm_anova().
#>
#> Moderation of `arm`: 2 candidates, holm-adjusted
#>  moderator df statistic     p p_adj
#>   baseline  1     1.434 0.232 0.464
#>        age  1     0.009 0.926 0.926
```

Two candidates now, not three, and `baseline`’s adjusted p is 0.464
rather than 0.348 — it is no longer sharing the correction with a term
that was never part of the search.

## Choosing the adjustment

`adjust` is passed straight to
[`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html), so every
method it knows is available: `"holm"` (the default), `"BH"`,
`"bonferroni"`, `"BY"`, `"hochberg"`, `"hommel"`, `"none"`.

Holm controls the family-wise error rate and is uniformly more powerful
than Bonferroni, which is why it is the default: if you are going to act
on the single strongest moderator, that is the quantity you want
controlled. If you intend to *report a set* of candidates for later
study, Benjamini-Hochberg controls the more appropriate thing and is a
little more generous:

``` r

ilm_moderation(fit, x = "arm", adjust = "BH")
#>  moderator df    statistic            p        p_adj
#>       site  2 15.543996217 2.899606e-07 8.698819e-07
#>   baseline  1  1.433625978 2.317733e-01 3.476599e-01
#>        age  1  0.008546022 9.263839e-01 9.263839e-01
```

With only three candidates the methods barely separate. They diverge
when there are many, or when several are real.

## The honest split, and why it is not the default

There is another way to be honest about a search: pick the moderator on
half the data and test it on the other half. The test half was never
used to choose, so nothing needs adjusting.

``` r

ilm_moderation(fit, x = "arm", split = TRUE)
#> Moderation of `arm`, honest split (240 chose, 240 tested)
#>   strongest of 3 candidates on the first half, tested on the second
#>  moderator df statistic        p    p_adj
#>       site  2     8.672 0.000233 0.000233
```

It is intellectually cleaner, and it is **better calibrated** — 0.053
against Holm’s conservative 0.027 in the simulation above. It is not the
default because of what it costs:

                                  Holm    honest split
      moderation by a numeric     0.547      0.320
      moderation by a 3-level     0.380      0.200

Roughly half the power. And the gap *widens* with degrees of freedom:
moving from a 1 df numeric moderator to a 2 df three-level factor costs
Holm 31% of its power and splitting 38%. Splitting is weakest precisely
in the three-level-factor case that experimental work is full of.

Splitting is genuinely **required** only when the hypotheses cannot be
counted — an open-ended tree search over covariates and split points,
where there is no adjustment to apply.
[`ilm_moderation()`](https://huttoncp.github.io/illume/reference/ilm_moderation.md)
searches a named, countable set, so there is.

## What a null here does and does not mean

If nothing survives adjustment, that is **not** evidence the effect is
constant. A moderation search is underpowered by construction: an
interaction needs far more data than a main effect, and a multi-df
candidate needs more still. In the simulation above, with a real
moderator present, Holm found it 54.7% of the time for a numeric
moderator and 38.0% for a three-level factor — at n = 500.

So the honest reading of an empty table is “this study could not detect
moderation of the sizes it was powered for”, which is a much weaker
claim than “the effect is the same for everyone”.
[`ilm_moderation()`](https://huttoncp.github.io/illume/reference/ilm_moderation.md)
says so in its own output rather than leaving you to remember it.

If heterogeneity matters to your conclusion, power for it deliberately —
[`ilm_power_design()`](https://huttoncp.github.io/illume/reference/ilm_power_design.md)
will take an interaction as its `term`.

## See also

[`vignette("effect-size-and-power")`](https://huttoncp.github.io/illume/articles/effect-size-and-power.md)
for powering a study that can detect this,
[`?ilm_anova`](https://huttoncp.github.io/illume/reference/ilm_anova.md)
for a moderator you specified in advance,
[`?ilm_trends`](https://huttoncp.github.io/illume/reference/ilm_trends.md)
for slopes within levels, and
[`?ilm_emmeans`](https://huttoncp.github.io/illume/reference/ilm_emmeans.md)
and
[`?ilm_contrast`](https://huttoncp.github.io/illume/reference/ilm_contrast.md)
for describing the moderation once you have decided it is real.
