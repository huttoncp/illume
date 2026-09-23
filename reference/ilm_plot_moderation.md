# See a moderation

Draws the effect of `x` at each level, or several values, of a moderator
– the picture behind an
[`ilm_moderation()`](https://huttoncp.github.io/illume/reference/ilm_moderation.md)
row.

## Usage

``` r
ilm_plot_moderation(
  x,
  moderator = NULL,
  exposure = NULL,
  at = NULL,
  level = 0.95,
  data = NULL,
  main = NULL,
  colour = NULL,
  pch = NULL,
  ...
)
```

## Arguments

- x:

  An
  [`ilm_moderation()`](https://huttoncp.github.io/illume/reference/ilm_moderation.md)
  result, or a fitted
  [`ilm_model()`](https://huttoncp.github.io/illume/reference/ilm_model.md)
  that already contains the interaction.

- moderator:

  Which moderator to draw. Defaults to the strongest row.

- exposure:

  The treatment or exposure, when `x` is a plain model.

- at:

  For a continuous moderator, the values to evaluate at. Defaults to its
  quartiles, so the picture spans the data rather than a nominal range.

- level:

  Confidence level.

- data:

  Needed only if a refit is required and the moderator is not in the
  model frame.

- main, colour, pch:

  Passed through; `pch` takes a name.

- ...:

  Passed to
  [`tinyplot::tinyplot()`](https://grantmcdermott.com/tinyplot/man/tinyplot.html).

## Value

Invisibly, the data frame that was plotted.

## Details

For a categorical `x` this is the cell means, the interaction plot. For
a continuous `x` it is the **slope** of `x` within each level of the
moderator, from
[`ilm_trends()`](https://huttoncp.github.io/illume/reference/ilm_trends.md),
because a slope is the quantity being said to differ. Both carry
intervals.

This is deliberately not
[`ilm_plot_model()`](https://huttoncp.github.io/illume/reference/ilm_plot_model.md)'s
effect plot, which holds the other predictors at typical values –
pinning the moderator is exactly what would hide the moderation.

## See also

[`ilm_moderation()`](https://huttoncp.github.io/illume/reference/ilm_moderation.md),
[`ilm_trends()`](https://huttoncp.github.io/illume/reference/ilm_trends.md),
[`ilm_emmeans()`](https://huttoncp.github.io/illume/reference/ilm_emmeans.md).

## Examples

``` r
# \donttest{
set.seed(1)
d <- data.frame(tx = rbinom(300, 1, 0.5), age = rnorm(300),
                site = factor(sample(c("a", "b", "c"), 300, TRUE)))
d$y <- 0.4 * d$tx + 0.8 * d$tx * (d$site == "c") + rnorm(300)
fit <- ilm_model(y ~ tx + age + site, data = d, verbose = FALSE)
#> ilm_model(): family = "gaussian", inferred from `y`: continuous values from -2.72 to 4.16. Pass `family` to choose another.
m <- ilm_moderation(fit, x = "tx", progress = FALSE)
ilm_plot_moderation(m)

# }
```
