# Reduce a data frame's variables to a few dimensions

Takes the columns of a data frame and summarises them as a small number
of dimensions. Which method that means is decided by the column types:
PCA when they are all numeric, multiple correspondence analysis when
they are all categorical, and a mixed method when both are present. A
date or date-time column cannot be used by any of the three and is
dropped with a message; convert it to something numeric first if it
should count.

## Usage

``` r
ilm_reduce(data, cols = NULL, ndim = 5, method = c("pcamix", "glrm"), ...)
```

## Arguments

- data:

  A data frame.

- cols:

  Columns to use. A character vector of names, a regular expression, a
  predicate function such as `is.numeric`, or `NULL` for all of them –
  see
  [ilm_selection](https://huttoncp.github.io/illume/reference/ilm_selection.md).

- ndim:

  Number of dimensions to keep.

- method:

  `"pcamix"` (the default) for PCA, MCA or FAMD depending on the column
  types, in closed form. `"glrm"` fits a generalized low rank model
  instead, which uses a loss appropriate to each column's type rather
  than squared error on one-hot indicators, and reconstructs a category
  as a category. It costs an iterative fit, and on all-numeric data the
  two are the same model – see
  [`ilm_glrm()`](https://huttoncp.github.io/illume/reference/ilm_glrm.md)
  for when it is worth that.

- ...:

  Passed to
  [`ilm_glrm()`](https://huttoncp.github.io/illume/reference/ilm_glrm.md)
  when `method = "glrm"`.

## Value

An object of class `"ilm_reduce"`: `method` (`"pca"`, `"mca"` or
`"famd"`), `eig` (dimension, eigenvalue, percent of variance and its
cumulative total), `ind_coord` (`row_id` and one column per retained
dimension – this is what
[`ilm_cluster()`](https://huttoncp.github.io/illume/reference/ilm_cluster.md)
takes), `var_contrib` (`variable`, `dim`, `sqload`: how strongly each
original variable relates to each dimension, on a 0 to 1 scale, for
numeric and categorical variables alike), `n`, and `fit`, the underlying
[`PCAmixdata::PCAmix()`](https://rdrr.io/pkg/PCAmixdata/man/PCAmix.html)
object for anyone who wants to go past this wrapper.

## Details

The mixed method is Chavent et al.'s, which belongs to the same
generalised-PCA family as the FAMD of Pages without being a
reimplementation of it. The two were checked against each other directly
and agree for this purpose, matching on eigenvalues and on individual
coordinates.

## References

Chavent, M., Kuentz-Simonet, V., Labenne, A. and Saracco, J. (2014).
Multivariate analysis of mixed data: the PCAmixdata R package. arXiv
1411.4911.

## See also

[`ilm_cluster()`](https://huttoncp.github.io/illume/reference/ilm_cluster.md)
to group the rows,
[`ilm_profile()`](https://huttoncp.github.io/illume/reference/ilm_profile.md)
for the whole pipeline,
[`ilm_reduce_na()`](https://huttoncp.github.io/illume/reference/ilm_reduce_na.md)
for the same thing applied to missingness.

## Examples

``` r
r <- ilm_reduce(mtcars)
r
#> <ilm_reduce> method = pca, n = 32, 5 dimension(s) retained
#>   first 3 dimension(s) explain 89.9% of the variance
#> 
#>   strongest variable per dimension (squared loading)
#>     dim 1   cyl                      0.924
#>     dim 2   qsec                     0.569
#>     dim 3   carb                     0.175
#>     dim 4   drat                     0.197
#>     dim 5   vs                       0.080
head(r$var_contrib[order(-r$var_contrib$sqload), ])
#>   variable dim sqload
#> 1      cyl   1 0.9239
#> 2     disp   1 0.8958
#> 3      mpg   1 0.8685
#> 4       wt   1 0.7916
#> 5       hp   1 0.7199
#> 6       vs   1 0.6209
```
