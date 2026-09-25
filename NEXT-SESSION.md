# Handoff — illume 0.0.8.9000

First written 2026-09-22. Rewritten on 2026-09-24, during the session
that added the random walk. The current state is in this first section.
Everything below “Standing constraints” is history.

Delete this file once the queue is empty.

## Where things stand (2026-09-24, evening)

- **Merged:** PRs \#6 to \#13 in illume and \#4 in illumex. They cover
  the 0.0.8.9000 bump, the interpretation work, the remedies and REML
  documentation, the joint-draw order (D1), averaging over every random
  term (D2 to D4), the quadrature and scenario follow-up, the random
  walk with
  [`ilm_cells()`](https://huttoncp.github.io/illume/reference/ilm_cells.md),
  and the padding in the prose numbers. Both packages were rebuilt from
  main (illume d793812, illumex c98ea7a) and installed into Craig’s R
  4.4 and R 4.6.1 libraries.
- **In review:** `kr-fix`, which fixes Kenward-Roger. The old one shrank
  the covariance, took Satterthwaite’s df, and its gate did not refuse
  what it said it refused. A review found it, and Craig put it ahead of
  the queue below.
- **Next, in the order agreed with Craig:**
  1.  The exports another agent asked for, which Craig approved:
      `ilm_matrices()` (new rows placed among the cells), `ilm_draws()`
      (joint draws of every parameter, matched by name), `ilm_dist()`
      (each family’s d/p/q/r in illume’s parameterisation),
      `ilm_normal_expect()` (the quadrature behind
      `predict(marginal = TRUE)`), and `ilm_ranef()` and
      `ilm_varcorr()`, registered on nlme’s generics. Also in that
      batch:
      - `fixef()` never reaches `fixef.ilm_model`, because the method is
        registered to illume’s own internal generic (`R/ilm_methods.R`),
        not nlme’s;
      - predict’s `marginal` help should say “set to zero (a typical
        group)”;
      - a sentence in `reml`’s help on few groups.
  2.  An outlier and influence check.
  3.  `ilm_plot_model_pdp()`, with its causal-language guard.
  4.  Reporting, version 1 (`ilm_session_log()` and `ilm_report()`,
      owned by illumex and extended by illume), as agreed with Craig on
      2026-09-24.
  5.  The multinomial slowdown, warm-started refits, and whether
      reformulas moves to Imports.
  6.  After that, for the same agent:
      [`ilm_scores()`](https://huttoncp.github.io/illume/reference/ilm_scores.md),
      [`ilm_calibration()`](https://huttoncp.github.io/illume/reference/ilm_calibration.md)
      and
      [`ilm_contrast()`](https://huttoncp.github.io/illume/reference/ilm_contrast.md)
      as S3 generics; `ilm_refit(fit, data)`; and optionally
      `ilm_support(fit, newdata)`.
- **Craig’s calls, still open:**
  - whether the latent-budget checks should FAIL at one observation per
    latent value for a non-gaussian family. A Poisson walk with one
    count per time point recovers its variance and is still marked
    unreliable.
  - the vocabulary for which random effects a prediction uses.
  - a formula interface for
    [`ilm_fit()`](https://huttoncp.github.io/illume/reference/ilm_fit.md),
    which was advised against.
- **At sprint end:** remind Craig about the open data case study, then
  the end-to-end illumex/illume experiment.

------------------------------------------------------------------------

## Standing constraints (do not re-derive, do not violate)

- **No GitHub repo creation, no pushing.** Craig reviews everything
  himself before anything is published. The illume remote already
  exists: `https://github.com/huttoncp/illume` (public since 2026-09-23,
  so its Actions pages can be read without signing in).
- **No git history rewrites** without per-instance approval. One was
  approved on 2026-09-19; that approval covered only that one.
- **Do not list Claude as a package author or contributor.** The
  disclosures in README, the intro vignette and CONTRIBUTING are the
  transparency mechanism.
- **The “review in progress” paragraph in the AI-disclosure text is
  Craig’s to update**, not ours.
- **Craig does the applied case study himself.**
- **Other projects on this computer stay out of illume until Craig
  publishes them:** never develop one inside the illume tree, never
  reference one from illume, and never name or describe one in anything
  pushed – commit messages, PR bodies, NEWS, these notes. illume must
  stand on its own.
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

### 0. Agreed on 2026-09-23, after the split (see the section at the end)

1.  **Publishing the split – in progress.** Craig created
    `https://github.com/huttoncp/illumex` (public, MIT licence commit)
    and asked for pull requests to review: one on illumex (the package),
    one on illume (the split). **Merge illumex’s first**: illume’s CI
    installs `Depends: illumex` through `Remotes: huttoncp/illumex`,
    i.e. from illumex’s `main`. Then enable GitHub Pages for illumex
    (gh-pages branch) as for illume.
2.  **Remedies the code can act on – BUILT on branch
    `actionable-remedies`,** not yet merged.
    `ilm_remedies(fit, dispersion =, zeros =, variance =)` lists a
    remedy for every check that is not OK, with a tier (numerical /
    structural / estimand) and the
    [`ilm_model()`](https://huttoncp.github.io/illume/reference/ilm_model.md)
    change that makes it; `ilm_apply_remedy(fit, remedies, id)` refits
    through the fit’s own call and keeps `fit$remedy_log`. The rules
    live in R/ilm_remedies.R and read the fit, not the suggestion text;
    `test-remedies.R` fails if a check has no rule. Writing the rules
    out also fixed: `re_struct` naming only some terms (it errored about
    `re_struct$NA`), `rr(1)` advice for one-category-dimension models (a
    no-op), “drop term ‘s(x)’” for a smooth shrunk to a line (which
    would drop the line too), and two remedies the package does not
    have.
3.  **A design document for the conductor (`ilm_analysis()`) –
    WRITTEN:** `dev/design-ilm-analysis.md` (build-ignored). Not to
    build now. It keeps the agreed shape (`ilm_plan()`, the remedy loop
    over the functions above with a revert rule and a refit budget, a
    decision log and sensitivity table, `ilm_report()`, validation on
    the messy regimes) and the rules that must survive: REML to report
    only where it is exact (gaussian), follow-ups declared and
    multiplicity-adjusted, only the exposure causal on the DAG route.
4.  **Gaps, all worth filling** (Craig agreed): Firth-type bias
    reduction for sparse multinomial categories; small-sample inference
    for GLMM Wald tests with few clusters; multinomial parity
    (cluster-robust SEs,
    [`ilm_scenario()`](https://huttoncp.github.io/illume/reference/ilm_scenario.md));
    the imputation argument on
    [`ilm_model()`](https://huttoncp.github.io/illume/reference/ilm_model.md)
    (5a below); unmeasured-confounding sensitivity (E-values) for
    [`ilm_dag_model()`](https://huttoncp.github.io/illume/reference/ilm_dag_model.md);
    a reporting function.

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

- **Python versions of both packages – agreed 2026-09-23, AFTER both are
  on CRAN.** The R packages are the reference implementation: every
  Python function gets an agreement test against them, run through rpy2
  in CI (as pyfixest does with fixest). illumex first (polars for
  collapse; it builds the cross-language harness), then a JAX prototype
  of illume’s engine – the multinomial GLMM first, benchmarked against R
  – before committing to a full port. The hard part is TMB’s sparse
  Laplace approximation (dense Hessians are fine up to a few thousand
  random effects; crossed designs need sparse Cholesky with hand-written
  derivatives); JAX’s vmap/jit would most speed up power, the parametric
  bootstrap and the envelopes. API: same formulas, names and arguments,
  `import pyillume as ilm` -\> `ilm.model()`. Names: `illume` is taken
  on PyPI (an unrelated placeholder); pyillume and pyillumex proposed,
  not yet chosen.
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

Committed in `8986280`: `ilm_plot_anomaly()`; column names needing
backticks work everywhere (mgcv stand-ins); one refit list for every
refit (the ZI LRT chi-square 68.8 vs 0.02 defect); boundary fits keep
usable fixed effects (63.3% -\> 97.5% in `flat_re`); per-category
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
  the next CI run is the test. **What CI then said** (read from the raw
  job log Craig supplied, R 4.6.1): the seed-1 test the floor was for
  PASSES; the failure moved to the seed-4 test, which asserted the hold
  path specifically. On R 4.6.1 that fit reaches the same optimum
  (logLik -311.3, correlation 1) but TMB’s own Hessian is positive
  definite there, so it takes the “tmb” route – the coin toss noted
  below, now seen in the wild. Fixed in the tests, not the package: the
  seed-4 test accepts either route and asserts what holds on both (SEs
  1.004-1.035 x the no-RE model’s on CI), and a new test forces the hold
  directly so it is exercised on every platform. oldrel (R 4.5) started
  failing at the same commit, presumably the same way; not seen in a
  log. **Lesson:** a WebFetch summary of the Actions page reported that
  failed run as a success, and it was relayed to Craig as a pass. Read
  CI status from the API (`/repos/{owner}/{repo}/actions/runs` and
  `.../jobs`), never from a page summary. Job logs need admin rights;
  Craig can supply the raw-log link.
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
- **tmb vs boundary was a coin toss for boundary fits – DECIDED AND
  BUILT.** Whether TMB’s own Hessian was positive definite at a boundary
  decided whether the fixed effects’ SEs carried the covariance’s
  uncertainty (“tmb”) or were conditional on it (“boundary”, as lme4).
  Three rules were measured on 4,000 fits
  (`studies/scripts/boundary_se.R`): A hold the whole term, B the status
  quo, and C hold only the flat directions. Craig chose C on 2026-09-23.
  `ilm_hess_recover()` now holds only the eigen-directions of the
  boundary terms’ block with curvature below `ilm_flat_rel` (1e-3) of
  its largest, whatever TMB said, which gives the reduced model’s
  standard errors. Next: the 0.0.8.9000 bump, with the full study set
  re-run (boundary_se included), then reinstall both packages locally.
- The pilot’s “usable” (convergence code 0) is stricter than illume’s
  `ok`, which grades code 8 with a small gradient WARN; the messy study
  should count fits the illume way when it is next run.
- blme’s `wishart(common.scale = FALSE)` fails inside blme
  (`repackageMerMod`, lme4 1.1.35.3); `common.scale = TRUE` plus a
  gamma(3) residual prior does not reproduce the same objective. Not
  worth more time.

------------------------------------------------------------------------

## The split into illume + illumex (2026-09-23)

Craig’s decision, after asking what a split would move. The exploratory
half is now **`illumex`**, a sibling folder:
`...\AI experiments\illumex` (the name went illumeda -\> illumex
mid-split; nothing of the first name survives). A full copy of illume
from just before the split, `.git` included, is in
`...\AI experiments\Illume-pre-split`.

- **The line**: data in (illumex) against model in (illume). Measured
  with a call graph before it was drawn: nothing that moved fits or
  reads a model, and the one call across is
  [`ilm_impute()`](https://huttoncp.github.io/illume/reference/ilm_impute.md)
  -\> `ilm_glrm()` (imported by name). 65 of 147 exports moved (the 64
  planned plus `ilm_sim()`), ~20% of the code.
- **Wiring**: illume `Depends: illumex`, so
  [`library(illume)`](https://github.com/huttoncp/illume) gives both and
  nothing downstream changed. `Remotes: huttoncp/illumex` until illumex
  is on CRAN; that field must go before illume’s CRAN submission, and
  illumex must be accepted first.
- **illumex must never depend on illume.** Its vignettes mention illume
  code in `eval = FALSE` chunks and its docs name illume functions as
  `illume::f()` text, never as links. A circular Suggests was considered
  and rejected.
- **Shared helpers are copied, not imported**: `%||%`, `ilm_bq()`,
  `ilm_wrap()`, `ilm_progress()` and the four `ilm_pch*()` helpers.
  illume’s `tests/testthat/test-shared-helpers.R` fails if a copy
  drifts. The `ilm_progress_arg` help page lives in illumex only.
- **Tests split by block**, not by file: the imputation blocks of
  test-missing.R became illume’s test-impute.R, the model-plot blocks of
  test-boot-and-plots.R became test-model-plots.R, the GLRM-imputation
  block became test-impute-glrm.R, and the anomaly-vs-imputation rank
  comparison became test-anomaly-rank.R (reaching
  `illumex:::ilm_anom_rank`).
- **Each package has its own** `dev/make_iml_aliases.R`,
  `pkgdown/make_reference_index.R` (navbars link each site to the
  other), workflows (illumex’s has no TMB step), README, NEWS,
  CONTRIBUTING and AI disclosure (the “review in progress” paragraph
  copied unchanged – Craig’s to update).
- **To work on illume locally, illumex has to be installed**
  (`R CMD INSTALL ../illumex`, or `devtools::install("../illumex")`),
  since `load_all()` resolves `Depends`. A change to illumex needs
  reinstalling before illume sees it.
- illumex’s git history starts fresh, on top of the licence commit
  GitHub made when Craig created the repo; the moved files’ history
  stays in illume’s. Its `LICENSE` is the two-line form CRAN needs, and
  GitHub’s full MIT text is `LICENSE.md` (build-ignored).
- Both packages were installed into Craig’s main library after the split
  (with vignettes), so
  [`library(illume)`](https://github.com/huttoncp/illume) there is the
  post-split version.
