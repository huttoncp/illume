# Handoff — illume 0.0.7.9000

Written 2026-09-22, at the end of a session that ran out of headroom. Craig is
relaunching on Opus 5.5. Everything below is uncommitted intent or in-flight
work; committed work is in the log.

Delete this file once the queue is empty.

---

## Standing constraints (do not re-derive, do not violate)

- **No GitHub repo creation, no pushing.** Craig reviews everything himself
  before anything is published. The two private remotes already exist:
  `https://github.com/huttoncp/illume`, `https://github.com/huttoncp/illuminator`.
- **No git history rewrites** without per-instance approval. One was approved
  on 2026-09-19; that approval covered only that one.
- **Do not list Claude as a package author or contributor.** The disclosures in
  README, the intro vignette and CONTRIBUTING are the transparency mechanism.
- **The "review in progress" paragraph in the AI-disclosure text is Craig's to
  update**, not ours.
- **Craig does the applied case study himself.**
- **illuminator is a separate repo.** Never develop it inside the illume tree.
- **Do not reference illuminator from illume** until Craig decides to publish
  it. illume must stand on its own.
- No new hard dependencies. Suggests behind `requireNamespace()` guards is fine.
  Prefer collapse over data.table.
- Never commit raw `.rds` from studies (`.gitignore` covers
  `studies/runs/**/*.rds`). Full findings go to the paper; a summary note goes
  to that version's `NEWS.md`.
- Every diagnostic names a remedy that exists in-package.
- `ilm_*` prefix, sum-to-zero contrasts, bespoke in-package diagnostics.
- Title, everywhere: **"A Unified Engine for Exploration and Frequentist
  Inference"** (DESCRIPTION and README both already updated).

## The lesson that keeps being confirmed

**Every real defect in this package was found from outside it** — by comparing
against lme4/nlme/survreg/mclogit, by Craig calling a function, or by
`R CMD check`. Never by the test suite. Twice in one session the full suite
passed while `R CMD check` caught real defects, because `devtools::test()`
never builds the package and `load_all()` exports everything regardless of
NAMESPACE. **Run `R CMD check`, not just the tests.**

Two traps that bit repeatedly this session, both now in memory:

- **Roxygen insertion trap.** A helper inserted between a documented
  function's `#'` block and its definition silently un-exports that function.
  Put helpers *above* the whole doc block.
- **Undocumented arguments.** Adding an argument without an `@param` is caught
  only by `R CMD check`. Verify package-wide with `tools::checkDocFiles(".")`.
- **Name collisions on new arguments.** Adding `keep =` to
  `ilm_anomaly_iforest()` silently collided with a local `keep <- rep(TRUE, n)`
  and the attribute came out a logical vector. Grep the function body for the
  name before adding a parameter.
- Bash heredocs eat backslashes and will silently corrupt R regex strings. Use
  the Write/Edit tools for anything containing `\`.

---

## Queue, in Craig's stated order

### 1. Messy-data study and benchmarking vignette  -- DONE

The run completed (400 reps x 6 regimes). `studies/findings/messy.md` is
written, `vignettes/benchmarking.Rmd` is complete and renders at 35 KB, and
NEWS has the summary note. Nothing left here except the two cross-references
below, which were not done:

- **The significance paraphrase still needs adding to the regression-models
  vignette** (not README), with a link to the benchmarking vignette. Craig
  approved the draft paraphrase in an earlier session.
- **README and the intro vignette do not yet reference the benchmarking
  vignette.** Check also that `benchmarking.Rmd` appears wherever the vignette
  index or pkgdown config lists them.
- `scripts/summarise_run.R` has no branch for this study, so `messy.md` was
  written from the csv by hand -- the file says so. Adding a branch would keep
  it consistent with the other findings files.

**The headline is not a clean win and the vignette says so.** illume covers
better in all six regimes (0.931 vs 0.853 with a sparse category; 0.945 vs
0.894 combined) and is 6-13x faster. But it *inflates* coefficients when a
category is sparse (attenuation 1.46, and mean absolute bias 0.394 against
mclogit's 0.271), it converged on only **63.3%** of replications when the true
random-effect sd was 0.05 against mclogit's 100%, and mclogit's RMSE is lower
in four of six regimes. Do not let a later edit quietly turn this into a
clean win.

A **brms agreement arm** against the messy regimes was never run. It would
need to be small -- 20-30 reps, one chain per worker, ~45s a fit. The existing
`studies/findings/brms.md` covers only the tidy case.

### 2. `ilm_plot_anomaly()` — approved, not started

Craig approved all four types this session ("yes, those all sound good to
me"). Design agreed: **one exported function**, `ilm_plot_anomaly(x, type =
c("scores", "drivers", "map", "row"), ...)`, default `"scores"`.

- **`"scores"`** — sorted score against rank, with the null envelope drawn as
  a ribbon and a vertical line at the number flagged. This is the point of the
  whole function: it is the only thing that distinguishes "five genuine
  outliers" from "the top 5% of a smooth continuum". When the curve is smooth,
  the remedy to name is: say plainly you are taking a fixed share, or accept
  the structure is not outlier structure and go to `ilm_profile()`.
- **`"drivers"`** — bar chart of `driver` among flagged rows. Nearly free,
  already computed. One column driving 90% of flags is a coding problem in
  that column, not a multivariate anomaly.
- **`"map"`** — first two `ilm_reduce()` dimensions, all rows grey, flagged
  highlighted and sized by score. Needs `attr(x, "data")`.
- **`"row"`** — for one flagged row, signed z-score *and* signed
  reconstruction residual per column, side by side. Small |z| with large
  |residual| means odd only in combination, which is the reconstruction
  method's actual selling point and is currently invisible.

**Groundwork is already committed.** `ilm_anomaly()` now stores `null_curve`,
an n × 3 matrix (`lo`, `mid`, `hi` = 2.5/50/97.5 percentiles of the per-rank
null score distribution across the B simulated datasets). Read it with
`attr(x, "null_curve")`.

**The honest constraint:** `method = "iforest"` sets `calibrated = FALSE` and
`p = NA`; it has no null, and its `alpha` is an assumed contamination rate.
`null_curve` is absent there. Do not draw the same reference band for both —
under iforest, label the line as a quantile rather than a test, and for
`type = "row"` say the combination view needs the reconstruction method.

Conventions to match (see `R/ilm_plot_profile.R`): tinyplot; validate the
class with an error naming the actual class; return `invisible(NULL)`; pass
`legend = list(title = ...)` **explicitly** — tinyplot deparses arguments for
legend titles and `do.call` has broken this twice. Use `ilm_pch()` for point
symbols. Test with realistic factor labels (`"north"/"central"/"south"`), not
`"a"/"b"/"c"` — the deparse bug hides behind short labels.

### 3. Missingness comparison study

Separate run, results into the **missing-data vignette**. Craig explicitly
wanted this kept apart from the messy-data study.

### 4. Final pass

- Documentation and `NEWS.md` current.
- `devtools::check()` (CRAN check).
- Reinstall on this machine.

### 5. Open, not scheduled

- macOS CI: `___kmpc_for_static_fini` is an LLVM/Intel OpenMP symbol;
  flat-namespace linking expects libomp pre-loaded. A macOS-only source-install
  step for TMB/RTMB was added to `.github/workflows/R-CMD-check.yaml`; unverified.
- Craig asked for "a succinct technical explanation plus accessible
  non-technical language of how the multinomial mixed effects model fitting
  works in illume". **This was only partially answered before an
  interruption and was never completed.** Worth finishing.

---

## State at handoff

Working tree clean apart from `studies/runs/messy/` (gitignored `.rds`; the
log and csv are not yet committed pending the run finishing).

Recent commits:

- `48c7725` Fix nested random effects, and let anomalies feed the explorers
- `668e57a` Keep a null reference curve on an anomaly scan; draft benchmarking vignette

Both features from Craig's lunch-break report are **done and tested**:

- `(1 | higher/lower)` now fits. Validated against lme4 — coefficients agree
  to 6.5e-07, RE sds to four decimals. `gapminder_unfiltered` fits 3,313 rows.
  `tests/testthat/test-nested-re.R`, 18 assertions.
- `ilm_anomaly()` keeps its data; `ilm_anomalous()` extracts the flagged rows;
  `ilm_profile()`, `ilm_cluster()`, `ilm_reduce()`, `ilm_describe()` and
  `ilm_describe_all()` all take the object directly.
  `tests/testthat/test-anomalous.R`, 28 assertions.
