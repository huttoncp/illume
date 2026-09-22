# Screen variables the model does not use

Runs
[`ilm_check_covariate()`](https://huttoncp.github.io/illume/reference/ilm_check_covariate.md)
over every variable in a data frame that the model's formula does not
mention, and reports any whose residual pattern suggests it should have
been included.

## Usage

``` r
ilm_check_omitted(
  object,
  data,
  vars = NULL,
  B = 30L,
  ncores = 1L,
  seed = 1L,
  verbose = TRUE
)
```

## Arguments

- object:

  A fitted `"ilm_model"` object.

- data:

  A data frame containing candidate variables.

- vars:

  Optional character vector naming variables to test; defaults to
  everything not already in the model.

- B:

  Integer. Simulated datasets.

- ncores:

  Integer. Worker processes.

- seed:

  Integer. Random seed.

- verbose:

  Logical. Print results.

## Value

Invisibly, a named list of per-variable results.

## Details

Because one test is run per variable, some flags are expected by chance;
the output states how many. Treat a single borderline flag among many
variables with appropriate caution.
