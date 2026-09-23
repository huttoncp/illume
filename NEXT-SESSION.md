# Handoff — illume 0.0.7.9000

First written 2026-09-22 at the end of a session that ran out of
headroom; rewritten at the end of the Opus 5.5 session that followed,
the same day, and updated on 2026-09-23 after the CI failure on R
release (see “What the 2026-09-23 session did” at the end). Items 1 and
2 of the original queue are done and committed (`8986280`); the power
and marginal-means work is committed in `3f6e368`.

Delete this file once the queue is empty.

------------------------------------------------------------------------

## Standing constraints (do not re-derive, do not violate)

- **No GitHub repo creation, no pushing.** Craig reviews everything
  himself before anything is published. The two remotes already exist:
  `https://github.com/huttoncp/illume` (public since 2026-09-23, so its
  Actions pages can be read without signing in),
  `https://github.com/huttoncp/illuminator` (private).
- **No git history rewrites** without per-instance approval. One was
  approved on 2026-09-19; that approval covered only that one.
- **Do not list Claude as a package author or contributor.** The
  disclosures in README, the intro vignette and CONTRIBUTING are the
  transparency mechanism.
- **The “review in progress” paragraph in the AI-disclosure text is
  Craig’s to update**, not ours.
- **Craig does the applied case study himself.**
- **illuminator is a separate repo.** Never develop it inside the illume
  tree.
- **Do not reference illuminator from illume** until Craig decides to
  publish it. illume must stand on its own.
- No new hard dependencies. Suggests behind
  [`requireNamespace()`](https://rdrr.io/r/base/ns-load.html) guards is
  fine. Prefer collapse over data.table.
- Never commit raw `.rds` from studies (`.gitignore` covers
  `studies/runs/**/*.rds`). Full findings go to the paper; a summary
  note goes to that version’s `NEWS.md`.
- Every diagnostic names a remedy that exists in-package.
- `ilm_*` prefix, sum-to-zero contrasts, bespoke in-package diagnostics.
- Title, everywhere: **“A Unified Engine for Exploration and Frequentist
  Inference”** (DESCRIPTION and README both already updated).

## The lesson that keeps being confirmed

**Every real defect in this package was found from outside it** — by
comparing against lme4/nlme/survreg/mclogit, by Craig calling a
function, or by `R CMD check`. Never by the test suite. Twice in one
session the full suite passed while `R CMD check` caught real defects,
because `devtools::test()` never builds the package and `load_all()`
exports everything regardless of NAMESPACE. **Run `R CMD check`, not
just the tests.**

Two traps that bit repeatedly this session, both now in memory:

- **Roxygen insertion trap.** A helper inserted between a documented
  function’s `#'` block and its definition silently un-exports that
  function. Put helpers *above* the whole doc block.
- **Undocumented arguments.** Adding an argument without an `@param` is
  caught only by `R CMD check`. Verify package-wide with
  `tools::checkDocFiles(".")`.
- **Name collisions on new arguments.** Adding `keep =` to
  `ilm_anomaly_iforest()` silently collided with a local
  `keep <- rep(TRUE, n)` and the attribute came out a logical vector.
  Grep the function body for the name before adding a parameter.
- Bash heredocs eat backslashes and will silently corrupt R regex
  strings. Use the Write/Edit tools for anything containing `\`.

------------------------------------------------------------------------

## Queue, in Craig’s stated order

### 1. Messy-data study and benchmarking vignette – DONE (`a7f5cdb`, `8986280`)

Including the three leftovers (the significance paraphrase in the
regression-models vignette, README/intro cross-references to the
benchmarking vignette, a `summarise_run.R` branch for the study).

**The headline is not a clean win and the vignette says so.** Do not let
a later edit turn it into one. The boundary work below changed how many
fits in the `flat_re` regime are *usable* (63.3% -\> 97.5%), but the
study’s own tables still report what was measured at the time; the next
run of `messy_compare.R` will count fits the new way.

### 2. `ilm_plot_anomaly()` – DONE (`8986280`)

### 3. Missingness comparison study – NOT STARTED

Separate run, results into the **missing-data vignette**. Craig
explicitly wanted this kept apart from the messy-data study.

### 4. brms vs illume on the messy suite – REQUESTED, DESIGN MOSTLY AGREED

Craig asked for it after the backtick/review/iml work, with a plan from
Gemini. Settled on 2026-09-22: **no R-INLA arm** (Craig’s decision), and
the zero-warmup suggestion is dropped – warmup stays as normal. Still to
settle with him before building it:

- **Priors on the variance components.** Flat ones are improper and are
  what makes brms struggle near a boundary, the regime of interest. brms
  defaults (half-t on the SDs), stated as a difference between the
  methods.
- **Parameterisation:** brms codes against a baseline category, illume
  sum to zero. Compare predicted probabilities, or transform
  coefficients with the contrast matrix – never raw coefficients.
- **TOST** needs an equivalence margin agreed *before* the run.
- **Cost:** about 45 s a fit, so 20-30 replications per regime: an
  AGREEMENT study on the same datasets, not a coverage study.
  `studies/findings/brms.md` covers only the tidy case.

### 5. Plain-language explanation of multinomial mixed fitting – FOR THE PAPER

Craig, 2026-09-22: save it for the eventual paper rather than write it
now.

### 5a. Proposed, Craig interested – design before building

- **Multilevel imputation.**
  [`ilm_impute()`](https://huttoncp.github.io/illume/reference/ilm_impute.md)’s
  chained equations are single-level (checked: no grouping argument).
  Each conditional model could carry the analysis model’s random
  effects, for every family – few tools do nominal variables well, and
  illume fits multinomial mixed models natively.
- **An imputation argument on
  [`ilm_model()`](https://huttoncp.github.io/illume/reference/ilm_model.md)**
  (not
  [`ilm_fit()`](https://huttoncp.github.io/illume/reference/ilm_fit.md),
  which has no formula), so the analysis formula doubles as the
  imputation model with no extra work from the user. The principled form
  is substantive-model- compatible FCS (Bartlett et al. 2015, `smcfcs`):
  covariates imputed compatibly with the analysis model, interactions
  and random effects included; the outcome drawn from the model’s
  predictive; m fits pooled by Rubin’s rules. Open design question: what
  it returns, and which methods need pooling (summary and Wald tests
  from pooled coef/vcov; multi-df terms by the D1 statistic; emmeans is
  linear, so pooled coef/vcov suffice). Deterministic regression
  imputation (fill with the prediction) is the version NOT to build – it
  understates variance and overstates precision.
- **A boundary-avoiding penalty, opt-in – BUILT 2026-09-23,
  uncommitted.** `boundary = "avoid"` (Chung et al. 2013, 2015).
  Measured, not as expected: on the messy regimes (400 each, paired) it
  made every flat_re fit usable with coverage no worse (0.949 vs 0.946),
  but it inflates the SD (0.05 true: median 0.21 vs 0.10) and, through
  it, the fixed effects: x1.009-1.014 in four regimes, x1.11 with a rare
  category (coverage 0.930 -\> 0.916), x1.20 combined. So `"hold"` stays
  the default and the docs/message say what `"avoid"` costs. The pilot
  script and results are not in the repo (scratch). The same idea for
  fixed effects is Firth’s penalty (Kosmidis & Firth 2011, multinomial),
  the principled fix for the sparse-category inflation (1.46) the messy
  study found – still open, and the rare-category inflation above is one
  more reason for it.

### 6. Open, not scheduled

- macOS CI: `___kmpc_for_static_fini` (LLVM/Intel OpenMP). The
  macOS-only source-install step for TMB/RTMB in
  `.github/workflows/R-CMD-check.yaml` works: the 2026-09-23 run loaded
  RTMB and got as far as the tests, where it failed the same boundary
  test as Linux and Windows (fixed, see below).
- **Version bump** – 0.0.7.9000 has grown a lot; Craig’s call whether
  this is 0.0.8.
- Found this session, not done:
  [`ilm_scenario()`](https://huttoncp.github.io/illume/reference/ilm_scenario.md)
  refuses a multinomial fit (it could report per-category
  probabilities); cluster-robust SEs are not available for a multinomial
  fit.
- **Recorded, worth a study:** with 60 clusters, the Wald test of a
  between-cluster effect in a multinomial mixed model rejected a true
  null 5.8% (4 visits) and 7.8% (8 visits) of the time. The power
  functions count that test, and
  [`ilm_interpret()`](https://huttoncp.github.io/illume/reference/ilm_interpret.md)
  names
  [`ilm_pb_lrt()`](https://huttoncp.github.io/illume/reference/ilm_pb_lrt.md)
  for fewer than 100 groups. How size depends on the number of clusters
  is not yet measured.

------------------------------------------------------------------------

## What the 2026-09-22 (Opus 5.5) session did

Committed in `8986280`:
[`ilm_plot_anomaly()`](https://huttoncp.github.io/illume/reference/ilm_plot_anomaly.md);
column names needing backticks work everywhere (mgcv stand-ins); one
refit list for every refit (the ZI LRT chi-square 68.8 vs 0.02 defect);
boundary fits keep usable fixed effects (63.3% -\> 97.5% in `flat_re`);
per-category
[`ilm_ame()`](https://huttoncp.github.io/illume/reference/ilm_ame.md)/[`ilm_interpret()`](https://huttoncp.github.io/illume/reference/ilm_interpret.md);
`iml_*()` aliases (generated by `dev/make_iml_aliases.R`, enforced by a
test); `family = "auto"`.

After that commit, and committed since in `3f6e368`:

- **Power**:
  [`ilm_power()`](https://huttoncp.github.io/illume/reference/ilm_power.md)
  refits through the fit’s own design (contrasts, zero part, dispersion,
  censoring, AR, rp all kept); counts the test the analysis reports
  (t/z, joint F/chi-square); multinomial and joint tests; matrix-normal
  RE draws; frequency weights; pre-fit design failures reported rather
  than counted. Scaffolds redraw the planned design per replicate,
  balanced by CELL; multinomial and ordinal scaffolds; latent-scale ICC.
  Validated: t-test power (0.387 vs 0.395 at n = 20), noncentral
  chi-square for multinomial, Superpower/simr (see the effect-size
  vignette).
- **rp refits** rebuilt the spline baseline from the new times (it kept
  the observed times’ columns in every bootstrap/envelope/consistency
  refit).
- **Marginal means**: multinomial
  [`ilm_emmeans()`](https://huttoncp.github.io/illume/reference/ilm_emmeans.md)/[`ilm_contrast()`](https://huttoncp.github.io/illume/reference/ilm_contrast.md)
  and
  [`ilm_trends()`](https://huttoncp.github.io/illume/reference/ilm_trends.md)
  (emmeans/nnet agreement 1e-6), ordinal response-scale probabilities
  (emmeans/polr 1e-6).
- Ordinal
  [`ilm_scores()`](https://huttoncp.github.io/illume/reference/ilm_scores.md)
  (with RPS) and
  [`ilm_calibration()`](https://huttoncp.github.io/illume/reference/ilm_calibration.md);
  [`ilm_pb_lrt()`](https://huttoncp.github.io/illume/reference/ilm_pb_lrt.md)
  verdict needs its Monte Carlo interval; clearer refusal messages;
  `set_coef()` takes fixed effects alone; simulation refits quieter.
- [`ilm_interpret()`](https://huttoncp.github.io/illume/reference/ilm_interpret.md):
  inferred family, ordinal wording, BOUNDARY verdict, few-groups caveat;
  new methods for
  [`ilm_power()`](https://huttoncp.github.io/illume/reference/ilm_power.md)
  and
  [`ilm_contrast()`](https://huttoncp.github.io/illume/reference/ilm_contrast.md).
- Docs: NEWS, README, regression-models (family auto; the stale “zero
  inflation is out of scope” line), effect-size-and-power (re-measured),
  package help title corrected to the standing title.

------------------------------------------------------------------------

## What the 2026-09-23 session did (uncommitted until Craig says otherwise)

**CI on R release failed after `3f6e368`** (macOS, Windows, Ubuntu), all
in `test-boundary.R`: a multinomial `(1 | site)` fit at a correlation of
-1 ended in false convergence, gradient 1.08, standard errors 0. Release
is R 4.6.1 with **RTMB 2.0**, TMB 1.9.25, lme4 2.0-6; this machine has
RTMB 1.x, where the same fit stops at -9.9 and is fine.

- **Cause:** a correlation of +/-1 is a log-Cholesky diagonal at -Inf,
  and the RE density goes through `solve(t(L))` (centred
  parameterisation), so the inner Hessian’s Cholesky cancels e^(-2c)
  terms: the objective is noise below about -12 (C = 3) to -15 (C = 2),
  and RTMB 2.0’s optimiser walked in and found a spurious maximum.
  Scalar variances are exact to -20; rr is non-centred; diag has no
  correlation. So only `us` with C \>= 2 (multinomial only) and the AR
  innovation covariance.
- **Fix:** `ilm_logsd_floor = -10`, `ilm_floor_pos()`,
  `ilm_floor_refit()` in `R/ilm_fit.R`: a bounded refit only when the
  first pass went below it; the hold step respects it;
  `ilm_hess_recover()` refuses a non-stationary point (max \|grad\| \>
  1e-2, the gradient check’s FAIL line). Checked on the 2400 pilot fits:
  engaged in 144 (137 in flat_re!), coefficients moved \<= 2e-5, and 2
  fits that had SEs of exactly 0 *here too* are now fine. **RTMB 2.0 is
  not the trigger by itself.** With Craig’s go-ahead, RTMB 2.0 and TMB
  1.9.25 were compiled into a scratch library (gone with the session):
  the COMMITTED code passes there too, with the RTMB 1.9 answer to the
  digit. What still differs from CI is R 4.6.1 (local: 4.4.3) and Matrix
  1.7-5 (local: 1.7-2; TMB’s inner sparse Cholesky is Matrix’s CHOLMOD).
  So the fix is verified on mechanism – the hand-placed -18 test, seed
  8, the 144 pilot fits – but not yet under the conditions that failed:
  the next CI run is the test.
- **R CMD check** (as CI: `--no-manual --as-cran`), final sources: 0
  errors, 0 warnings, 2 NOTEs, both environmental (new submission;
  unable to verify the time). An earlier run caught `pen_re`/`pen_k`
  missing from
  [`globalVariables()`](https://rdrr.io/r/utils/globalVariables.html);
  fixed.
- **lme4 2.0** warns that `findbars`/`nobars` moved to reformulas: now
  `ilm_findbars()`/`ilm_nobars()` (reformulas first, lme4 fallback),
  reformulas in Suggests.
- **pkgdown’s warning** was the Node 20 deprecation:
  `actions/checkout@v6` and the Pages deploy action at v4.8.0 (pinned by
  SHA, as r-lib’s examples do). Its notice (ubuntu-latest -\> Ubuntu 26
  on 2026-10-19) needs nothing.
- **Hex sticker** (Craig’s artwork): `man/figures/logo.png`, 480 px
  wide, outside of the hex transparent and the ring’s edge un-matted so
  it has no dark fringe on white; made by `dev/logo/make_logo.R` from
  `dev/logo/illume-hex-source.jpg`. README header links it to the site;
  pkgdown puts it in the home-page header and uses it as og:image.
  `pkgdown/favicon/` came from pkgdown’s automatic favicon step (which
  sent the logo to realfavicongenerator.net – unintended, noted to
  Craig); keep it committed, and see make_logo.R for the two edits made
  to it. Local pkgdown is 2.0.9, whose template links the OLD favicon
  names; CI’s 2.2.1 links exactly the files present. Worth updating
  pkgdown locally.
- `boundary = "avoid"` finished: message at fit time naming its cost,
  docs with the pilot numbers, NEWS, and `test-boundary-avoid.R`
  (penalty against a brute-force exact penalised likelihood, unpenalised
  [`logLik()`](https://rdrr.io/r/stats/logLik.html), message, refits).
  `blme` was installed locally for the comparison only; it is not a
  dependency.

**Open, found this session:**

- **Non-centred parameterisation** (u = L z, as lme4 does) would remove
  the boundary numerics at the root instead of flooring them: the
  density of z never involves L’s inverse. A large change (objective,
  REML, the penalty, Sigma_d, AR, `set_coef()`), so only worth it if the
  floor proves not enough.
- **Sigma_d has the same structure** (`A <- solve(Ld)`): an
  intercept-slope correlation of +/-1 in ANY family reaches the same
  cancellation. No floor there, because Ld is in the slope’s units
  relative to the intercept’s, so any fixed floor is scale-dependent.
  Not seen failing yet.
- **tmb vs boundary is a coin toss for boundary fits.** Whether TMB’s
  own Hessian is positive definite at a boundary decides whether the
  fixed effects’ SEs carry the covariance’s uncertainty (“tmb”) or are
  conditional on it (“boundary”, as lme4). The floor flipped 42 of 144
  fits each way, and the SEs moved by up to 5% (10% in one). Worth
  deciding one rule.
- The pilot’s “usable” (convergence code 0) is stricter than illume’s
  `ok`, which grades code 8 with a small gradient WARN; the messy study
  should count fits the illume way when it is next run.
- blme’s `wishart(common.scale = FALSE)` fails inside blme
  (`repackageMerMod`, lme4 1.1.35.3); `common.scale = TRUE` plus a
  gamma(3) residual prior does not reproduce the same objective. Not
  worth more time.
