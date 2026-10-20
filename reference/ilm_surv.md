# Right-censored follow-up, in the survival convention

A convenience wrapper for the commonest censoring there is: a study
where some subjects had the event and the rest were still event-free
when follow-up ended.

## Usage

``` r
ilm_surv(time, event, seed = 1L)
```

## Arguments

- time:

  Follow-up time, strictly positive.

- event:

  `1` or `TRUE` if the event was observed, `0` or `FALSE` if the subject
  was censored.

- seed:

  Random seed for the censoring times drawn for subjects who had the
  event; see the section below. Drawn once, so the result is
  reproducible and a replicate is censored the same way every time.

## Value

An `"ilm_censor"` object, to pass as `ilm_model(censor = )`.

## Details

The argument is `event`, in the convention
[`survival::Surv()`](https://rdrr.io/pkg/survival/man/Surv.html) uses:
`1` when the event was observed, `0` when the subject was censored. That
is the opposite of the code
[`ilm_censor()`](https://huttoncp.github.io/illume/reference/ilm_censor.md)
stores, which is why this wrapper exists – getting it backwards silently
fits the model to the wrong subjects, and the fit will not look wrong.

## How the diagnostics censor their replicates

Every envelope in this package simulates from the fit, and a replicate
has to be censored the way the study was or it is not the same process.
A subject censored in the data has a known censoring time, which is used
directly. A subject who had the event does not: theirs is only known to
exceed the time they failed at, so it is drawn from the censoring
distribution – the Kaplan-Meier with the roles of event and censoring
swapped – conditioned on exceeding it. Those draws are made once, here,
so they are reproducible and every replicate is censored consistently.

## See also

[`ilm_censor()`](https://huttoncp.github.io/illume/reference/ilm_censor.md)
for floors and ceilings,
[`ilm_model()`](https://huttoncp.github.io/illume/reference/ilm_model.md)
with `family = "weibull"`, `"lognormal"` or `"loglogistic"`.

## Examples

``` r
t <- c(5, 9, 12, 3, 20)
e <- c(1, 0, 1, 1, 0)
ilm_surv(t, e)
#> censoring: 3 observed, 0 left, 2 right (of 5)
#>   limits not recorded; diagnostics hold the pattern fixed
```
