# Reduce a data frame's missingness pattern to a few dimensions

The missingness counterpart to
[`ilm_reduce()`](https://craig-hutton.github.io/illume/reference/ilm_reduce.md).
Builds a present/missing marker for every column, drops any whose
missingness never varies, and reduces the markers – which, being
two-level categorical variables, always takes the MCA route. The
dimensions that come back describe which columns tend to go missing
*together*, which is what separates a block of variables lost to one
skipped section from values that went missing independently.

## Usage

``` r
ilm_reduce_na(data, cols = NULL, ndim = 5)
```

## Arguments

- data:

  A data frame.

- cols:

  Columns to use. A character vector of names, a regular expression, a
  predicate function such as `is.numeric`, or `NULL` for all of them –
  see
  [ilm_selection](https://craig-hutton.github.io/illume/reference/ilm_selection.md).

- ndim:

  Number of dimensions to keep.

## Value

An object of class `"ilm_reduce_na"`, which is also an `"ilm_reduce"` –
see
[`ilm_reduce()`](https://craig-hutton.github.io/illume/reference/ilm_reduce.md)
for the shared structure. In `var_contrib`, `sqload` means how strongly
a column's *missingness* relates to a dimension, not its values.

## See also

[`ilm_profile_na()`](https://craig-hutton.github.io/illume/reference/ilm_profile_na.md)
for the whole pipeline,
[`ilm_check_missing()`](https://craig-hutton.github.io/illume/reference/ilm_check_missing.md)
for whether any of it matters to your model.

## Examples

``` r
r <- ilm_reduce_na(airquality)
#> ilm_reduce_na(): dropping column(s) whose missingness never varies (always or never missing): Wind, Temp, Month, Day
r
#> <ilm_reduce_na> method = mca, n = 153, 2 dimension(s) retained
#>   first 2 dimension(s) explain 100.0% of the variance
#> 
#>   strongest variable per dimension (squared loading)
#>     dim 1   Ozone                    0.511
#>     dim 2   Ozone                    0.489
```
