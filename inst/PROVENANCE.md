# Provenance

Where illume's methods come from: for each one, what it is, the code that
carries it out, the source it rests on and where that source came from, and
how it is checked. A source is a published paper or book, or, where no
publication gives what the code needs, a measurement of another package's
exported results, recorded here with its files.

This register is new and lists one entry so far. illume's other sources are
cited in its help pages (see `inst/REFERENCES.bib`) and are still to be
entered here.

## Entries

| What | Code | Source, and where it came from | Checked |
|---|---|---|---|
| The design of a factor-smooth interaction, `s(x, id, bs = "fs")`, at new covariate values: `(PredictMat(sm, newdata) %*% diag(trans.D))[, order(rind)]`, where `trans.D` and `rind` are elements of what `mgcv::smooth2random()` returns for the smooth. For this smooth it returns no `trans.U`, which the expression for other smooths uses. | `ilm_smooth_design()` (R/ilm_predict.R) | A measurement, not a publication. The expression was found by comparing what mgcv 1.9.1's exported functions `smooth2random()` and `PredictMat()` return for the same smooth on the same data; mgcv's source code was not read. The random blocks equal the expression to 2.2e-14 or better in eight set-ups (5, 30 and 60 levels; basis sizes 3 to 10; penalty orders 1 to 3; unbalanced groups; a cr basis; a two-dimensional smooth). It does not hold with a numeric `by` variable, which `ilm_model()` refuses for this smooth. The script and its output are `identity.R` (SHA-256 615b059d90814b4c8c79ff34e8ff095cbfe620f7bb8109e99a38c9d50c0924ac) and `identity-out.txt` (SHA-256 379b8e97b3eeff16eea75695912411d6861921b02b7781df6cff87535b17819e), run on 10 Nov 2026 with R 4.4.3 and kept with the project's working records, outside this repository. | tests/testthat/test-smooth-fs-predict.R: the expression against the fitted blocks, and `predict()` against `mgcv::predict.gam()` |
