# Mark censored observations

Describes which observations are known only as an interval: values at or
below a floor, at or above a ceiling, or both.

## Usage

``` r
ilm_censor(y, lower = NA, upper = NA, censored = NULL, event = NULL)
```

## Arguments

- y:

  The response.

- lower:

  Floor. Values at or below it are left-censored.

- upper:

  Ceiling. Values at or above it are right-censored.

- censored:

  Alternatively, the direction of censoring given directly: `-1`
  left-censored, `0` observed, `1` right-censored.

- event:

  Or an event indicator, in
  [`survival::Surv()`](https://rdrr.io/pkg/survival/man/Surv.html)'s
  convention: `1` (or `TRUE`) when the value was observed, `0` (or
  `FALSE`) when it is right-censored. Give `censored` or `event`, not
  both.

## Value

An integer vector of codes, classed `"ilm_censor"`, carrying the limits
it was built from.

## Details

Pass the result to `ilm_model(censor = )`. With `lower` or `upper`
given, the codes are derived from the response, and the limits are
carried along so that the simulation-based diagnostics can censor their
replicates the same way.

`censored` and `event` are for the case where the limits vary by
observation and only the outcome is known – right-censored follow-up in
a survival study, most commonly. Diagnostics then hold the censoring
pattern fixed across replicates rather than re-drawing it, which is
stated where it matters.

## Two codings, one meaning

The same data can be coded two opposite ways, and a fit given the wrong
one does not look wrong: it is fitted to the wrong rows.

|  |  |  |
|----|----|----|
|  | `censored =` | `event =` |
| says | the direction of censoring | whether the event was observed |
| an observed value | `0` | `1` |
| a value known only to be beyond it (right-censored) | `1` | `0` |
| a value known only to be below it (left-censored) | `-1` | – |
| the convention of | this function | [`survival::Surv()`](https://rdrr.io/pkg/survival/man/Surv.html) |

Give one or the other, never both. A `censored` vector of only 0s and 1s
is said once a session, since it is most often an event indicator given
to the wrong argument.

For time-to-event data,
[`ilm_surv()`](https://huttoncp.github.io/illume/reference/ilm_surv.md)
takes `event` in `Surv()`'s convention and also draws the censoring
times the diagnostics need, and `ilm_model(Surv(time, event) ~ ...)` is
read as exactly that.

## The Surv() trap

`survival::Surv(time, status)` reads its second argument as an EVENT
indicator: `1` means the event was observed. Passed here as `censored`,
the same column means the opposite – every event becomes a censored row
and every censored row an event. Use `event =` for that column, or
[`ilm_surv()`](https://huttoncp.github.io/illume/reference/ilm_surv.md).

## See also

[`ilm_surv()`](https://huttoncp.github.io/illume/reference/ilm_surv.md)
for time-to-event data;
[`ilm_model()`](https://huttoncp.github.io/illume/reference/ilm_model.md);
[`illumex::ilm_describe()`](https://huttoncp.github.io/illumex/reference/ilm_describe.html),
whose `p_zero` column is often the first sign of a floor.

## Examples

``` r
y <- c(0, 0, 1.4, 2.9, 5, 5)
ilm_censor(y, lower = 0, upper = 5)
#> censoring: 2 observed, 2 left, 2 right (of 6)
#>   limits: lower 0, upper 5
#>   67% censored; the fit rests mostly on interval probabilities
table(ilm_censor(y, lower = 0))
#> 
#> -1  0 
#>  2  4 

## the same three right-censored rows, coded both ways
t <- c(5, 9, 12, 3, 20)
ilm_censor(t, censored = c(0, 1, 0, 0, 1))   # 1 = censored
#> censored codes 1 as right-censored (not observed); for an event indicator (1 = event observed, as in Surv) use event = or ilm_surv()
#> censoring: 3 observed, 0 left, 2 right (of 5)
#>   limits not recorded; diagnostics hold the pattern fixed
ilm_censor(t, event = c(1, 0, 1, 1, 0))      # 1 = event observed, as in Surv()
#> censoring: 3 observed, 0 left, 2 right (of 5)
#>   limits not recorded; diagnostics hold the pattern fixed

## the trap: an event indicator given as `censored` flips every row
st <- c(1, 0, 1, 1, 0)                       # Surv()'s status, 1 = event
table(ilm_censor(t, censored = st))          # wrong: the events are censored
#> 
#> 0 1 
#> 2 3 
table(ilm_censor(t, event = st))             # right
#> 
#> 0 1 
#> 3 2 
```
