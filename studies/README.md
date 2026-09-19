# Simulation studies

The evidence behind illume's validity claims. These are **not** part of the R
package — `.Rbuildignore` excludes this directory from `R CMD build`, so it does
not affect the package check. It is versioned here so the results are backed up
alongside the code that produced them.

Every script writes raw per-replicate results as `.rds` next to a summary
`.csv`, so any table can be recomputed or re-cut without refitting anything.
All were run on 2026-09-19 against illume 0.0.1.9000 (with the `lum_*` API;
the scripts have since been updated to `ilm_*`).

## Scripts

| script | what it answers | cost |
|---|---|---|
| `coverage_study.R` | Do Wald confidence intervals cover at 95%? | 22,000 fits, ~38 min |
| `power_study.R` | Is Type I error nominal, and what is the power curve? | 60,000 fits, ~35 min |
| `benchmark_freq.R` | Does illume agree with lme4 / glmmTMB / nnet, and how fast? | ~1,100 fits |
| `mclogit_compare.R` | illume (Laplace) vs mclogit (PQL) for multinomial mixed models | 2,000 fits |
| `brms_compare.R` | Does illume agree with an independent implementation of the same likelihood? | 3 fits + MCMC |

Each takes arguments, e.g.

```
Rscript coverage_study.R <nrep> <ncore> <outdir>
Rscript power_study.R    <nrep> <ncore> <outdir> [deltas] [cells]
Rscript mclogit_compare.R <nrep> <ncore> <outdir>
```

They require illume to be **installed** (not just loaded), because the parallel
workers call `library(illume)`.

## Results

- `coverage/` — one `cov_<cell>.rds` and `.csv` per design cell
- `power/` — one `pow_<cell>_d<delta>.rds` per cell and effect size, plus `power_summary.csv`
- `bench/`, `mclogit/`, `brms/` — summary CSVs and raw RDS

## Headline findings

**Coverage is nominal everywhere** — 0.947 to 0.954 across 11 designs, including
five families, AR(1), and deliberately under-powered cells.

**The latent budget governs convergence, not coverage.** At J=5 the convergence
rate ran 33% / 58% / 83% / 98% at 1.5 / 3 / 5 / 10 observations per latent
value, while coverage among converged fits stayed nominal. The approximation
fails loudly rather than returning overconfident intervals.

**Type I error is nominal** for both Wald and LRT, and null p-values are uniform
by Kolmogorov-Smirnov test at every level checked, not merely at 5%.

**mclogit's PQL silently under-covers.** With 100 clusters of 4 observations its
nominal 95% intervals cover 89.1%, against illume's 94.9%, because PQL shrinks
fixed effects to 73% of true magnitude and the bias does not show in its
standard errors. Caveat worth keeping attached: mclogit has *lower* RMSE
(0.336 vs 0.409), since shrinkage buys variance reduction — a reasonable trade
for prediction, but not for inference.

**Agreement with independent implementations** is at numerical tolerance
(lme4, glmmTMB, nnet: 1e-15 to 7e-5) and within 0.10 standard errors of brms,
at roughly 82x the speed.

## Caveats to carry into any write-up

- Coverage in cells with convergence failures is **conditional on convergence**.
  Failed fits are excluded, and if failure correlates with extreme estimates the
  surviving coverage is optimistic.
- The number of free covariance parameters matters independently of the latent
  budget: J=3 at 2.5 obs/latent converged 95% of the time, while J=5 at 3.0
  converged only 58%.
