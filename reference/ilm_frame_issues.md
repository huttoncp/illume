# Problems that belong to pairs of columns

Constant columns, identifier-like columns, duplicated columns and
near-perfect collinearity. These cannot live in a per-variable table
because they are properties of the data frame as a whole.

## Usage

``` r
ilm_frame_issues(data, cor_cut = 0.999)
```

## Arguments

- data:

  A data frame.

- cor_cut:

  Absolute correlation at or above which a numeric pair is reported as
  collinear.

## Value

A data frame with `issue`, `columns` and `detail`; zero rows when
nothing is found.

## Details

A continuous variable is unique per row by construction, so uniqueness
is reported as identifier-like only for discrete-valued columns.

## Examples

``` r
ilm_frame_issues(ilm_sim())
#>      issue columns                                 detail
#> 1 constant  cohort zero variance; breaks the model matrix
```
