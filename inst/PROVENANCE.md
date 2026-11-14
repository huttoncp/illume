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
| The design of a factor-smooth interaction, `s(x, id, bs = "fs")`, at new covariate values: `(PredictMat(sm, newdata) %*% diag(trans.D))[, order(rind)]`, where `trans.D` and `rind` are elements of what `mgcv::smooth2random()` returns for the smooth. For this smooth it returns no `trans.U`, which the expression for other smooths uses. | `ilm_smooth_design()` (R/ilm_predict.R) | A measurement, not a publication. The expression was found by comparing what mgcv 1.9.1's exported functions `smooth2random()` and `PredictMat()` return for the same smooth on the same data; neither mgcv's source code nor its help pages were read. The random blocks equal the expression to 2.2e-14 or better in eight set-ups (5, 30 and 60 levels; basis sizes 3 to 10; penalty orders 1 to 3; unbalanced groups; a cr basis; a two-dimensional smooth). It does not hold with a numeric `by` variable, which `ilm_model()` refuses for this smooth. The evidence is in `dev/provenance/fs-predict/`: the script `identity.R` (SHA-256 7986f52c6614718344aac3df99ea85f4284e4523d0078e92025fccabf07b1d6c), the simulation it sources, `sim.R` (SHA-256 931c19bd793cce908a4786db205de623ce9fa7a082f4050594b372c2819336e2), and its output, `identity-out.txt` (SHA-256 7cf4336d5e6d741f753f3535a85872af3727b41b87ec9d66a4f7a27397b8bcf9). The script runs from the repository root, as `Rscript dev/provenance/fs-predict/identity.R`; the output kept was made on 11 Nov 2026 with R 4.4.3, illume 0.0.8.9003 (commit 62e3e6a) and mgcv 1.9.1. | tests/testthat/test-smooth-fs-predict.R: the expression against the fitted blocks, and `predict()` against `mgcv::predict.gam()` |
