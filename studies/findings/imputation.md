# Imputation ablation

Which method, and when. Known cells are hidden, filled in, and scored against
what was really there -- on **reconstruction**, which is what everyone reports,
and on **coverage** of a downstream coefficient after Rubin pooling, which is
what an analysis actually needs.

Generated from `studies/runs/<version>/imputation/` by
`scripts/imputation_ablation.R`. Newest version first.

---

## 0.0.7.9000

40 replicates per cell, 10 imputations each, `beta_true = 0.6` on the first
predictor. `rmse` is per-column-standardised reconstruction error on hidden
numeric cells; `cat_acc` is the share of hidden categorical cells recovered;
`coverage` is of a nominal 95% interval. Monte Carlo error on coverage at 40
replicates is about 0.034, so differences under 0.07 are not differences.

```
-- n=200 p=8, true rank 3 --------------------------------------------------
     method  rmse cat_acc coverage width   bias
       mean 1.011     ---    0.850 0.275 -0.049
        fcs 0.706     ---    0.900 0.270 -0.016
    lowrank 1.477     ---    0.650 0.321 -0.141
 lowrank_pa 0.881     ---    0.825 0.272 -0.064
       glrm 1.113     ---    0.775 0.258 -0.070

-- n=200 p=8, FULL rank (no structure to find) -----------------------------
     method  rmse cat_acc coverage width   bias
       mean 0.996     ---    0.900 0.322  0.001
        fcs 1.424     ---    0.900 0.315 -0.014
    lowrank 3.355     ---    0.500 0.330 -0.204
 lowrank_pa 1.426     ---    0.875 0.326 -0.077
       glrm 1.355     ---    0.900 0.305 -0.058

-- n=400 p=12, rank 4, 30% missing -----------------------------------------
     method  rmse cat_acc coverage width   bias
       mean 0.995     ---    0.825 0.175 -0.031
        fcs 0.609     ---    0.975 0.154  0.006
    lowrank 0.725     ---    0.950 0.179 -0.020
 lowrank_pa 0.976     ---    0.475 0.177 -0.091
       glrm 1.152     ---    0.500 0.157 -0.081

-- n=60 p=80, rank 3 (wider than it is long) -------------------------------
     method  rmse cat_acc coverage width   bias
       mean 1.007     ---    0.925 0.395 -0.015
        fcs   ---     ---      ---   ---    ---
    lowrank 0.559     ---    0.975 0.384 -0.007
 lowrank_pa 0.876     ---    0.975 0.432 -0.043
       glrm 1.170     ---    0.975 0.378 -0.052

-- n=300, 6 numeric + 2 categorical, rank 2 --------------------------------
     method  rmse cat_acc coverage width   bias
       mean 1.002   0.508    0.900 0.265 -0.029
        fcs 0.815   0.634    0.900 0.257  0.001
    lowrank 1.449     ---    0.650 0.320 -0.114
 lowrank_pa 0.906     ---    0.800 0.257 -0.035
       glrm 1.051   0.548    0.675 0.248 -0.047

-- the same with 40% missing -----------------------------------------------
     method  rmse cat_acc coverage width   bias
       mean 1.007   0.502    0.775 0.303 -0.053
        fcs 0.889   0.615    0.950 0.299 -0.031
    lowrank 3.073     ---    0.500 0.326 -0.322
 lowrank_pa 1.065     ---    0.600 0.285 -0.091
       glrm 1.100   0.539    0.625 0.252 -0.091
```

### Chained equations wins wherever it can be fitted

On every design where FCS runs it has the best or joint-best coverage (0.900,
0.900, 0.975, 0.900, 0.950), the smallest bias, and the best reconstruction. It
is the default and the ablation gives no reason to change that.

It also recovers categories better than the generalized low rank model does --
0.634 against 0.548, and 0.615 against 0.539 -- which was not the expected
result and is reported because it was measured.

### The low-rank route exists for one case, and earns it there

`n = 60, p = 80` is the case: chained equations cannot be fitted at all, and the
low-rank reconstruction halves the error against filling in column means (0.559
against 1.007) while covering at 0.975. That is the whole justification for the
route and it holds.

Everywhere else the low-rank methods **under-cover badly** -- 0.50 to 0.65 in
four of the six cells. Reconstruction error does not reveal this: in the
full-rank cell `glrm` reconstructs no worse than `fcs` (1.355 against 1.424) and
covers at 0.900 while `lowrank` covers at 0.500. Coverage is the measure that
governs and the two do not agree.

### The rank selector: not the clean defect it looked like

An earlier single-design check suggested `ilm_lowrank_ncp()` -- which chooses
the rank by cross-validating held-out cells -- was simply picking badly. Across
six designs it is **design-dependent**, and neither criterion dominates:

| design | CV rank | parallel analysis |
|---|---|---|
| rank 3, n=200 p=8 | 0.650 | **0.825** |
| FULL rank | 0.500 | **0.875** |
| rank 4, n=400 p=12 | **0.950** | 0.475 |
| p > n | 0.975 | 0.975 |
| mixed, 20% missing | 0.650 | **0.800** |
| mixed, 40% missing | 0.500 | **0.600** |

(coverage; higher is better)

Cross-validation is worse in four cells and materially **better** in one -- the
design with the most data and the clearest structure, where it has enough
held-out cells to choose well. Parallel analysis is the safer default when
structure is weak or absent, and the wrong one when it is strong and
well-estimated.

So the selector is **not changed**. The earlier finding was real on the design
it was measured on and does not generalise, which is the reason for running six
designs rather than one. `ilm_anomaly()` continues to use parallel analysis,
for a different reason: it needs the number of directions that are real shared
structure, not the rank that best predicts a cell.

### Everything under-covers somewhat

Even mean-fill covers at 0.775 to 0.925 against a nominal 0.95, and the best
method reaches 0.975 in only two cells. With 40 replicates the Monte Carlo
error is 0.034, so 0.90 is a plausible 0.95; 0.50 is not. Read the table as
separating "acceptable" from "broken" rather than as ranking 0.90 against 0.95.

### What was learned about running it

Forcing chained equations past the point where it can be fitted does not fail
fast. With 80 predictors on 60 rows it attempts a rank-deficient regression per
column per cycle per imputation and runs for hours producing nothing. The
script now records that cell as inapplicable, which is what `method = "auto"`
does for a user anyway.
