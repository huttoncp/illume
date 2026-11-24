# The rows a fit was made from

The fit's model frame, with the variables its transformed terms are
built from added back. The frame holds a transformed term as its own
column, such as `log(x)` or `poly(t, 2)`, and not the variable
underneath it. So code that places new rows on a fit, to predict at
changed values of `x` or `t` or to tell which rows lie beyond the data,
cannot rebuild the term from the frame alone. From these rows,
[`predict()`](https://rdrr.io/r/stats/predict.html) and
[`model.matrix()`](https://rdrr.io/r/stats/model.matrix.html) can.

## Usage

``` r
ilm_data(object)
```

## Arguments

- object:

  A fitted `"ilm_model"` object.

## Value

A data frame: the model frame, with the added columns after its own.
`NULL` for a fit that has no model frame.

## Details

The columns added are those of `data` that the model frame lacks and
that the model names: in its formula, in a dispersion or zero-inflation
formula, or, for a correlation over time given by its columns, as in
`ilm_ar1(~ t | g)`, its time and group. They are matched to the frame's
rows by row name, so after rows are dropped for missing values they are
the fit's own rows.

Nothing is added, and the frame is returned alone, when `data` was not a
data frame, when its row names do not all match the frame's, or when the
frame already holds every variable the model names, as when nothing in
it is transformed. A fit made before these columns were kept also
returns its frame alone.

## See also

[`model.frame.ilm_model()`](https://huttoncp.github.io/illume/reference/ilm_model-accessors.md)
for the frame alone.

## Examples

``` r
set.seed(1)
d <- data.frame(t = rep(1:12, 10), g = factor(rep(1:10, each = 12)))
d$y <- 0.2 * d$t + stats::rnorm(10)[d$g] + stats::rnorm(120)
f <- ilm_model(y ~ poly(t, 2) + (1 | g), data = d, family = "gaussian",
               verbose = FALSE)
names(model.frame(f))   # the frame holds poly(t, 2), not t
#> [1] "y"          "poly(t, 2)" "g"         
head(ilm_data(f))       # t is added back
#>            y poly(t, 2).1 poly(t, 2).2 g t
#> 1  1.0853274 -0.145443618  0.158691998 1 1
#> 2  0.1633894 -0.118999324  0.072132726 1 2
#> 3 -0.6476944 -0.092555030  0.002885309 1 3
#> 4 -2.0411537 -0.066110736 -0.049050254 1 4
#> 5  1.4984771 -0.039666441 -0.083673963 1 5
#> 6  0.5286126 -0.013222147 -0.100985817 1 6
```
