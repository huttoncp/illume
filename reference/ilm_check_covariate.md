# Test whether residuals shift across a variable

Asks whether residuals are systematically higher or lower at some values
of a variable than others. If they are, the model is missing something
about that variable – it may need to be added, or added in a more
flexible form.

## Usage

``` r
ilm_check_covariate(
  object,
  x,
  name = NULL,
  nbin = 5L,
  B = 30L,
  ncores = 1L,
  seed = 1L,
  verbose = TRUE
)
```

## Arguments

- object:

  A fitted `"ilm_model"` object.

- x:

  A variable with one value per row of the data the model was fitted to.
  Rows the fit dropped for missing values are removed automatically, so
  the column can be passed straight from the original data frame.
  Numeric variables are binned; factors are used as-is.

- name:

  Character label for printing.

- nbin:

  Integer. Bins for a numeric variable.

- B:

  Integer. Simulated datasets.

- ncores:

  Integer. Worker processes.

- seed:

  Integer. Random seed.

- verbose:

  Logical. Print the result.

## Value

Invisibly, a list with `status`, the observed statistic and the test.

## Details

The variable does **not** have to be in the model; pointing the test at
variables you left out is exactly the useful case, and is the one
situation where residual testing has real power here. Broad summaries of
residuals detect almost nothing; a test aimed at specific structure
detects a great deal. See
[`ilm_rqr_test()`](https://huttoncp.github.io/illume/reference/ilm_rqr_test.md)
for why.

## See also

[`ilm_check_omitted()`](https://huttoncp.github.io/illume/reference/ilm_check_omitted.md).
