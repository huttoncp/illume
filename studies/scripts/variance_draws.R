## Draws of a poorly determined variance (item 78, Arm B).
##
## PRE-REGISTERED DESIGN, committed before any code and before any run.
##
## ## The question
## A variance parameter that the data barely determine has a Gaussian
##   (Laplace) approximation on its log scale with a long upper tail. Draws of
##   it then explode: latent SDs to e^3.7 and beyond on 12 counts, and CRPS 1e7
##   against a log score of 1.7. No hold reaches such a parameter: it is not at
##   a boundary, only poorly determined. The study measures two things:
## 1. whether the standard error of its log SD, fixed as a flag by a rule set
##   in advance, picks out the fits where the approximation fails;
## 2. for each of the two remedies Craig will choose between, what it does to
##   interval coverage and forecast scores, and what it costs in time.
##
## ## Build and scope
## - **Build.** The runs use the hold pass as ruled, in a pinned library:
##   - σ at 0.2 × sd(y);
##   - the gaussian AR(1) SD and random-effect SD at 0.1 × sd(y), with the
##     flatness test;
##   - the both-flagged policy;
##   - the draws fixes already on main (0.0.8.9002).
##   So the arm studies only what the holds cannot fix. A parameter the build
##     holds is recorded as held and left out of every Arm B denominator, and
##     the count held is reported per cell.
## - **Parameters.** Every variance parameter the fit estimates:
##   - the gaussian residual SD (log σ);
##   - the AR(1) latent SD (lchol_ar), which is the MARGINAL SD. illume's
##     `Sigma$ar` is the stationary variance for AR(1) and CAR(1), and for a
##     random walk the variance per unit of time. So a correlation near 1 does
##     not inflate it, and ρ does not enter the label;
##   - each random-intercept SD (theta);
##   - the negative binomial's log k (logdisp), since its draws explode the
##     same way.
##   The AR correlation is recorded, not labelled; it is not a variance.
## - **Families.**
##   - In: gaussian, Poisson, negative binomial and binomial. Binomial goes
##     in B2 and B3: few-cluster logistic models are the classic case of a
##     poorly determined random-intercept SD.
##   - Binomial is out of B1: a single binary series of 12 to 48 points
##     almost never determines an AR latent, so its answer is known in advance
##     (nearly all unbounded) and would crowd out the informative cells.
##   - Out of the study: beta, ordinal, multinomial and survival. Their
##     variance parameters are the same latent SDs on a link scale, and the
##     rule chosen here extends to them by construction. Whether it holds there
##     is a later check on stored fits, not part of this pre-registration.
##
## ## Design
## Families: gaussian (noise SD 0.5), Poisson (mean about exp(1.5)), negative
##   binomial (size 2). 100 replicates per cell.
##
## | Arm | Model | Cells |
## |---|---|---|
## | B1 single series | y ~ 1 + AR(1), T in {12, 24, 48}, ρ in {0.3, 0.7}, latent SD in {0.3, 0.8} | 3 families × 12 = 36 |
## | B2 panel | 20 series, y ~ 1 + (1 \| s) + AR(1), T in {12, 24}, intercept SD in {0.2, 0.5}, ρ 0.7, latent SD 0.5 | 4 families × 4 = 16 |
## | B3 random intercept | y ~ x + (1 \| g), (L, m) in {(5, 10), (10, 5), (30, 2)}, intercept SD in {0.1, 0.5} | 4 × 6 = 24 |
##
## Binomial is binary, with a logit-scale intercept of 0. That is 76 cells
##   and 7,600 fits: 5,200 in B1 and B2, 2,400 in B3. B1 and B2 simulate T + 6
##   times and fit the first T, so 6 future times per series are held out for
##   forecasts at horizons 1 to 6. B3 has no horizon; it is scored on the
##   latent SD's coverage alone.
##
## ## The flag's signal, and the other package's signal
## - **s:** the SE of each variance parameter's log SD, from `vcov(fit, full
##   = TRUE)`.
## - **d:** the other package's signal, the SD over 400 joint draws
##   (`ilm_draws(nsim = 400, seed = 1)`) of the same log SD, computed exactly
##   as its cross-tab script does.
##
## Both are reported. How well they agree is a result: their correlation, and
##   agreement at the chosen cut-off.
##
## ## The label, from outside the flag: the profile likelihood
## For each unheld variance parameter:
## - A profile on a 25-point grid of its log SD, from the estimate − 4 to the
##   estimate + 6 (wider above, since that is where the tail is). Every other
##   fixed parameter is re-optimised at each point, and the random effects by
##   the Laplace inner step.
## - Its 95% interval, {ψ : 2(ℓ̂ − ℓp(ψ)) ≤ 3.84}, linearly interpolated.
## - The label has two classes, reported apart. **approx_fails** is either of
##   them:
##   - **unbounded:** the profile has no upper end within the grid. The data
##     cannot bound the SD, which is a different thing from a Wald tail that is
##     too long.
##   - **wald_too_wide:** the profile's upper end exists, and the Laplace
##     (Wald) upper end, exp(est + 1.96 s), exceeds it by more than a factor of
##     2.
##
## **Ridges.** In gaussian AR fits, σ and the AR SD can trade off: each alone
##   is flat while their sum is determined. In B2, a random-intercept SD and
##   the AR SD can trade off the same way. A coarse 9 × 9 joint profile of the
##   pair is taken on the fits where both one-parameter profiles fail. The
##   candidates are 1,600 gaussian fits for σ and the AR SD (B1 1,200 and B2
##   400), and B2's 1,600 fits over four families for the RI SD and the AR SD.
##   A fit is labelled **ridge** when both one-parameter profiles fail but the
##   profile of the total variance is bounded. Ridge fits are reported as a
##   class of their own. For the second remedy, the pair is drawn jointly from
##   the 9 × 9 profile, not separately. A random-intercept SD with an AR SD in
##   B2 is treated the same way.
##
## ## The flag's cut-off, by a rule fixed now
## Candidate cut-offs on s: 0.5, 0.75, 1, 1.5, 2, 3. **Rule:** the smallest
##   cut-off at which at most 5% of the parameters not labelled approx_fails
##   are flagged. That is specificity of at least 0.95, pooled over all cells,
##   and it is the cut-off that catches the most failures without flagging many
##   good fits. It is confirmed on fresh seeds (50 replicates per cell, a seed
##   offset), with sensitivity and specificity reported per family, per arm and
##   per label class.
##
## Its edges, stated now:
## - If no candidate reaches specificity 0.95, the flag is reported as not
##   viable, and remedy 1 fails with it.
## - If 0.5, the lowest candidate, already reaches it, that is reported as
##   meaning a lower cut-off might do better. The study does not search below
##   0.5.
##
## **"Calibrated"** means an interval's coverage lies within 2 Monte Carlo
##   standard errors of its nominal level.
##
## ## The two remedies, measured on the same fits
## - **Remedy 1, a check that flags.** The draws are unchanged, and a flag
##   appears in `summary()` and in forecasts. What it does to intervals is
##   measured two ways:
##   - (a) standard joint draws, for flagged against unflagged fits;
##   - (b) the fallback the flag would name: draws given θ (`given =
##     "theta"`, the variance held at its estimate). Only its forecast metrics
##     are reported. Holding the SD at its estimate makes its coverage of the
##     true SD degenerate by construction.
## - **Remedy 2, profile draws.** Each flagged parameter's log SD is drawn
##   from its profile, the normalised exp(ℓp) on the grid, interpolated. The
##   other parameters are drawn from the Gaussian conditional on it. Ridge
##   pairs are drawn from the joint 9 × 9 profile.
##   - **Cap, fixed now.** The grid's own edge, at +6, is arbitrary, and on
##     unbounded fits it would be all that stops the draws. So the upper end is
##     capped on the link scale's own terms: a latent SD may not exceed 5 ×
##     max(sd(g(y*)), 0.5), where g is the link and y* the response on its
##     scale (y for gaussian, log(y + 0.5) for counts, logit((y + 0.5) / 2) for
##     binary). A latent SD five times the response's own spread on the link
##     scale is beyond anything the data could have produced.
##   - Remedy 2 is reported separately for unbounded and wald_too_wide fits.
##   - **Two assumptions, stated in the findings:** a flat prior on the log
##     SD, and a conditional Gaussian for the other parameters, taken from the
##     covariance at the estimate.
##
## **Measured for each remedy, and for today's joint draws as the baseline:**
## - Coverage of the true latent SD (and σ, and the RI SD) by the 95%
##   interval from the draws.
## - For B1 and B2 forecasts at horizons 1, 3 and 6:
##   - coverage of 80% and 95% prediction intervals;
##   - CRPS, mean and median (the mean is what explodes);
##   - the share of forecast rows with CRPS above 100 times the CRPS of the
##     true-parameter predictive distribution. That reference is the one this
##     study can compute exactly; the other package's ETS reference is reported
##     where comparable.
##   - All with Monte Carlo SEs, per cell and pooled over flagged fits.
## - **Time:** seconds per fit for the profile or profiles and the draws,
##   against today's draws, on one core, as median and 90th percentile.
##
## ## Named fixtures
## Each is refitted on the build and reported by name. For each: s and d per
##   parameter, the label, the flag at the chosen cut-off, whether a hold took
##   it, and the three remedies' interval for the latent SD and forecast CRPS
##   at horizon 1.
## - The single series s36: gaussian, cell 2, rep 11, T 24.
## - The covariate panel's reps 4, 11 and 59: gaussian AR-only panel, T 24.
## - The 12-count Poisson series s14: cell 7, rep 1.
##
## Rows of the other package's per-fit CSV whose fits change class between
##   illume 9c29c0e (its run) and this build are listed. Those are the draws
##   fixes and the new holds at work.
##
## ## Sizes
## - Fits: 7,600 in the main run and 3,800 fresh.
## - Each unheld variance parameter needs a 25-point profile, each point a
##   re-optimisation of the other parameters. The ridge fits add 81 points for
##   their pair.
## - Time, MEASURED. A smoke run of 45 seconds, on illume 0.0.8.9002, timed a
##   25-point profile per parameter type in one fit of a representative cell
##   per arm, in seconds per fit:
##
##   | Arm | Family | Fit | AR SD | σ or log k | RI SD |
##   |---|---|---|---|---|---|
##   | B1 (T 24) | gaussian | 1.2 | 2.3 | 0.2 | |
##   | B1 | Poisson | 0.1 | 1.8 | | |
##   | B1 | NB | 0.1 | 4.4 | 3.0 | |
##   | B2 (T 12) | Poisson | 0.2 | 9.8 | | 5.4 |
##   | B2 | binomial | 0.3 | 12.4 | | 0.6 |
##   | B3 (30 × 2) | binomial | 0.1 | | | 3.1 |
##
##   With the ridge pairs' 81-point profiles (about 3 × a single profile,
##     assumed on a third of the candidate fits), and the draws and forecasts
##     at about 15% on top:
##   - B1: about 5 core-hours;
##   - B2: about 10;
##   - B3: about 2.
##   That is about 17 core-hours, **about 9 to 10 hours for the main run at 2
##     cores, and 5 for the fresh one.** The earlier 3.5 hours was a guess, and
##     it is withdrawn.
##
##   If that is too long, one option for the conductor, decided before any
##     run: a 13-point grid over the same range, which about halves it (5 hours
##     and 2.5). It costs resolution at the interval's ends: 0.83 on the log
##     scale between points, against 0.42.
## - Build: about a day for the script, the profile draws and the summariser.
##
## ## Decisions it gives Craig
## 1. The flag's cut-off, from the rule above.
## 2. Remedy 1, remedy 2 or both, on the measured coverage, scores and time.
## 3. Whether ridge fits need their own treatment, if they behave differently
##   from the rest.
##
## Usage: Rscript variance_draws.R <nrep> <ncore> <outdir> [offset]
## ---------------------------------------------------------------------------
## (the study's code follows in a later commit; this commit is the design)
