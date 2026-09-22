# A simulated mixed-type dataset for testing and examples

Entirely synthetic, so it carries no third-party copyright. Every column
is present because some diagnostic should have something to say about
it: a fixture where nothing is wrong tests nothing.

## Usage

``` r
ilm_sim(n_id = 75L, n_period = 12L, seed = 2026L)
```

## Arguments

- n_id:

  Number of units.

- n_period:

  Observations per unit.

- seed:

  Random seed, so the fixture is reproducible.

## Value

A data frame with `n_id * n_period` rows and 13 columns.

## Details

- id:

  panel identifier, repeated across periods

- date:

  monthly, regular, repeated across units (long format)

- grp:

  factor with a deliberately rare level and an unused level

- site:

  character with case and whitespace variants, and empty strings

- flag:

  balanced logical

- consented:

  near-constant logical, a separation risk

- score:

  gaussian

- income:

  lognormal: right-skewed and bounded at zero

- visits:

  Poisson counts, dispersion near 1

- claims:

  negative binomial counts, dispersion above 1

- downtime:

  zero-inflated counts

- cohort:

  constant

- lab_value:

  gaussian with missing values

## Examples

``` r
d <- ilm_sim()
str(d)
#> 'data.frame':    900 obs. of  13 variables:
#>  $ id       : int  1 1 1 1 1 1 1 1 1 1 ...
#>  $ date     : Date, format: "2024-01-01" "2024-02-01" ...
#>  $ grp      : Factor w/ 5 levels "alpha","beta",..: 2 2 1 1 2 1 2 3 1 2 ...
#>  $ site     : chr  "East" "North" "South" "" ...
#>  $ flag     : logi  TRUE FALSE FALSE TRUE FALSE TRUE ...
#>  $ consented: logi  TRUE TRUE TRUE TRUE TRUE TRUE ...
#>  $ score    : num  50.3 49 43.3 49.1 50.6 ...
#>  $ income   : num  28146 49547 43454 62783 26390 ...
#>  $ visits   : int  6 2 6 3 0 4 3 1 1 3 ...
#>  $ claims   : num  0 3 1 0 0 2 7 3 2 5 ...
#>  $ downtime : int  0 1 0 5 0 4 0 0 4 0 ...
#>  $ cohort   : chr  "2024" "2024" "2024" "2024" ...
#>  $ lab_value: num  8.3 6.96 7.41 8.05 7.24 6.24 7.46 NA 7.6 7.69 ...
```
