# Simulation studies

The evidence behind illume's validity claims. These are **not** part of the R
package — `.Rbuildignore` excludes this directory from `R CMD build`, so it has
no effect on the package check or the distributed tarball. It is versioned here
so results stay with the code that produced them.

```
studies/
  scripts/     the study scripts, plus summarise_run.R and prune_runs.R
  runs/        <version>/<study>/  raw per-replicate .rds and summary .csv
               <version>/RUNINFO.dcf  version, date, R version
  findings/    <study>.md  permanent, generated record - newest version first
```

## Retention policy

**Raw per-replicate data is kept for the three most recent package versions.
For older versions only the findings are kept.**

The unit of retention is a version's whole **set** of studies, not each study
separately. Results are evidence about a particular build, so coverage from one
version and benchmarks from another are not a coherent picture; sets are kept or
dropped together. On version 7, the raw data for versions 5, 6 and 7 is retained
and everything earlier survives as findings only.

This ties retention to the `Version:` field in `DESCRIPTION`, so **bump the
development version when the engine changes materially** -- which is exactly
when the studies need re-running anyway. Re-running studies without bumping
overwrites that version's set rather than creating a new one, which is usually
what you want: one set of evidence per version.

Findings are *generated* from the raw data by `scripts/summarise_run.R`, never
hand-written, so the permanent record can never be lost to someone forgetting to
take notes. `scripts/prune_runs.R` enforces the policy and refuses to delete a
set unless every study in it has findings on disk.

```
Rscript scripts/prune_runs.R . 3 --dry-run   # show what would go
Rscript scripts/prune_runs.R . 3             # do it
```

Run it after each new study run, and commit the result.

### What this policy does not do

It prunes the **working tree, not git history**. Data that has been committed
stays in `.git` permanently, so the repository keeps growing even as the
checkout stays small — binary `.rds` files do not delta-compress, so each
version's set adds roughly its full size forever.

This policy keeps the directory comprehensible and fresh clones manageable. It
is not a way to cap repository size. If that becomes the goal, the real options
are to stop committing raw `.rds` at all, adopt git-lfs, or rewrite history —
each a larger decision than a retention rule.

## Scripts

| script | what it answers | cost |
|---|---|---|
| `coverage_study.R` | Do Wald confidence intervals cover at 95%? | 22,000 fits, ~38 min |
| `power_study.R` | Is Type I error nominal, and what is the power curve? | 60,000 fits, ~35 min |
| `benchmark_freq.R` | Does illume agree with lme4 / glmmTMB / nnet, and how fast? | ~1,100 fits |
| `mclogit_compare.R` | illume (Laplace) vs mclogit (PQL) for multinomial mixed models | 2,000 fits |
| `brms_compare.R` | Does illume agree with an independent implementation of the same likelihood? | 3 fits + MCMC |

Write output into `runs/<version>/<study>/`, which is the layout the summariser
and pruner expect, and record the run in `runs/<version>/RUNINFO.dcf`:

```
Rscript scripts/coverage_study.R  2000 10 runs/0.0.1.9000/coverage
Rscript scripts/power_study.R     2000 10 runs/0.0.1.9000/power
Rscript scripts/mclogit_compare.R  500 10 runs/0.0.1.9000/mclogit
```

They need illume **installed**, not merely loaded, because the parallel workers
call `library(illume)`.

## Keeping results honest

Committed results are a *snapshot*. Git preserves them faithfully but will never
re-run them, so after a material change to the fitting engine the repository
holds authoritative-looking numbers describing an older package — more dangerous
than having none, because nothing marks them stale.

When a change could move estimates or standard errors, re-run the affected study
and commit the new results in the same commit as the change. Re-run everything
before a release or paper submission.

Also check the scripts still use current function names. They broke silently at
the `lum_*` to `ilm_*` rename because nothing in the test suite covers this
directory.

## Findings

See `findings/*.md`. As of version 0.0.1.9000 (run 2026-09-19):

- **Coverage is nominal everywhere** — 0.947 to 0.954 across 11 designs.
- **The latent budget governs convergence, not coverage** — 33% / 58% / 83% /
  98% convergence at 1.5 / 3 / 5 / 10 observations per latent value, with
  coverage among converged fits staying nominal throughout.
- **Type I error is nominal** for Wald and LRT, with uniform null p-values.
- **mclogit's PQL silently under-covers** — 0.891 against illume's 0.949 at four
  observations per cluster.
- **Agreement with independent implementations** at numerical tolerance, and
  within 0.10 SE of brms at ~82x the speed.

Each findings file carries the caveats that belong with its result, including
that coverage in cells with convergence failures is conditional on convergence.
