# Recode a variable against a dictionary held as two vectors

The shape a lookup table takes when it comes out of a spreadsheet: one
vector of old codes, one of new, matched by position.

## Usage

``` r
ilm_translate(y, old, new, default = NA)
```

## Arguments

- y:

  A vector to recode.

- old:

  Values of the old coding scheme.

- new:

  Replacements, the same length as `old`.

- default:

  Value for entries not found in `old`. `NA` by default.

## Value

A vector of the same length as `y`.

## Details

Values with no entry in `old` become `NA` by default, which is
deliberate: a silent pass-through would hide codes the dictionary does
not cover. Use `default` to change that.

## See also

[`match()`](https://rdrr.io/r/base/match.html),
[`ilm_recode_errors()`](https://huttoncp.github.io/illume/reference/ilm_recode_errors.md).

## Examples

``` r
ilm_translate(c(1, 2, 3, 99), old = 1:3, new = c("low", "mid", "high"))
#> [1] "low"  "mid"  "high" NA    
ilm_translate(c(1, 2, 99), old = 1:2, new = c(10, 20), default = -1)
#> [1] 10 20 -1
```
