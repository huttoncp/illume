# Handoff — illume 0.0.7.9000

First written 2026-09-22 at the end of a session that ran out of headroom;
rewritten at the end of the Opus 5.5 session that followed, the same day.
Items 1 and 2 of the original queue are done and committed (`8986280`); what
that session added, and what it left open, is below.

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

### 1. Messy-data study and benchmarking vignette -- DONE (`a7f5cdb`, `8986280`)

Including the three leftovers (the significance paraphrase in the
regression-models vignette, README/intro cross-references to the benchmarking
vignette, a `summarise_run.R` branch for the study).

**The headline is not a clean win and the vignette says so.** Do not let a
later edit turn it into one. The boundary work below changed how many fits in
the `flat_re` regime are *usable* (63.3% -> 97.5%), but the study's own tables
still report what was measured at the time; the next run of
`messy_compare.R` will count fits the new way.

### 2. `ilm_plot_anomaly()` -- DONE (`8986280`)

### 3. Missingness comparison study -- NOT STARTED

Separate run, results into the **missing-data vignette**. Craig explicitly
wanted this kept apart from the messy-data study.

### 4. brms vs illume on the messy suite -- REQUESTED, DESIGN NOT AGREED

Craig asked for it after the backtick/review/iml work, with a plan from Gemini
and "we can discuss if you have other ideas or reason to believe gemini is
wrong". Points to settle with him before building it:

- **"0 warmup" is not workable.** HMC adapts its step size and metric during
  warmup; without it the sampler is not calibrated and the comparison would be
  of a mistuned sampler, not of brms. Compiling once and reusing with
  `update(newdata = )` is right.
- **Flat priors on the variance components** are improper and are what makes
  brms struggle near a boundary -- which is exactly the regime of interest.
  Better: brms defaults (half-t on the SDs), stated as a difference between
  the methods, not hidden.
- **Parameterisation:** brms codes against a baseline category, illume sum to
  zero. Compare predicted probabilities, or transform coefficients with the
  contrast matrix -- never raw coefficients.
- **TOST** needs an equivalence margin agreed *before* the run.
- **Cost:** about 45 s a fit, so 20-30 replications per regime, one chain per
  worker. `studies/findings/brms.md` covers only the tidy case.

### 5. Plain-language explanation of multinomial mixed fitting -- OPEN

Asked for long ago ("a succinct technical explanation plus accessible
non-technical language"), partially answered, never finished.

### 6. Open, not scheduled

- macOS CI: `___kmpc_for_static_fini` (LLVM/Intel OpenMP). A macOS-only
  source-install step for TMB/RTMB is in `.github/workflows/R-CMD-check.yaml`;
  unverified.
- **Version bump** -- 0.0.7.9000 has grown a lot; Craig's call whether this is
  0.0.8.
- Found this session, not done: `ilm_scenario()` refuses a multinomial fit
  (it could report per-category probabilities); cluster-robust SEs are not
  available for a multinomial fit.
- **Recorded, worth a study:** with 60 clusters, the Wald test of a
  between-cluster effect in a multinomial mixed model rejected a true null
  5.8% (4 visits) and 7.8% (8 visits) of the time. The power functions count
  that test, and `ilm_interpret()` names `ilm_pb_lrt()` for fewer than 100
  groups. How size depends on the number of clusters is not yet measured.

---

## What the 2026-09-22 (Opus 5.5) session did

Committed in `8986280`: `ilm_plot_anomaly()`; column names needing backticks
work everywhere (mgcv stand-ins); one refit list for every refit (the ZI LRT
chi-square 68.8 vs 0.02 defect); boundary fits keep usable fixed effects
(63.3% -> 97.5% in `flat_re`); per-category `ilm_ame()`/`ilm_interpret()`;
`iml_*()` aliases (generated by `dev/make_iml_aliases.R`, enforced by a test);
`family = "auto"`.

After that commit (uncommitted until Craig says otherwise -- check `git
status`):

- **Power**: `ilm_power()` refits through the fit's own design (contrasts,
  zero part, dispersion, censoring, AR, rp all kept); counts the test the
  analysis reports (t/z, joint F/chi-square); multinomial and joint tests;
  matrix-normal RE draws; frequency weights; pre-fit design failures reported
  rather than counted. Scaffolds redraw the planned design per replicate,
  balanced by CELL; multinomial and ordinal scaffolds; latent-scale ICC.
  Validated: t-test power (0.387 vs 0.395 at n = 20), noncentral chi-square
  for multinomial, Superpower/simr (see the effect-size vignette).
- **rp refits** rebuilt the spline baseline from the new times (it kept the
  observed times' columns in every bootstrap/envelope/consistency refit).
- **Marginal means**: multinomial `ilm_emmeans()`/`ilm_contrast()` and
  `ilm_trends()` (emmeans/nnet agreement 1e-6), ordinal response-scale
  probabilities (emmeans/polr 1e-6).
- Ordinal `ilm_scores()` (with RPS) and `ilm_calibration()`; `ilm_pb_lrt()`
  verdict needs its Monte Carlo interval; clearer refusal messages;
  `set_coef()` takes fixed effects alone; simulation refits quieter.
- `ilm_interpret()`: inferred family, ordinal wording, BOUNDARY verdict,
  few-groups caveat; new methods for `ilm_power()` and `ilm_contrast()`.
- Docs: NEWS, README, regression-models (family auto; the stale "zero
  inflation is out of scope" line), effect-size-and-power (re-measured),
  package help title corrected to the standing title.
