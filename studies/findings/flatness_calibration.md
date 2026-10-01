# flatness calibration -- findings log

The design is pre-registered in `../scripts/flatness_calibration.R`
(addendum A1 there too); the verdicts below are from
`../scripts/flatness_summarise.R` on the run's own files.

## 0.0.8.9003

Run on 2026-10-01 (item 269), illume main a387c93 with illumex ff5e579, on
two float paths: P1 R 4.4.3 / RTMB 1.9 and P2 R 4.6.1 / RTMB 2.0, one core,
51 and 57 minutes. 12,408 fits (6,204 per path): the saved phase-1 fits of
dispersion_limit.R (600), dispersion_limit_gaussian.R (500) and
re_sd_limit.R (5,100), rebuilt from their seeds, and the four fits behind
main's red CI. 2,919 fits per path are candidates (past their line on main).

**The rules.** New, as ruled: held if pushes of 1.5, 3 and 6 towards the
limit all move the objective by at most 0.05; unconverged if any lowers it by
more than 0.05. Old: one push of 3, held at a rise of at most 1e-3.

### The registered verdicts

| verdict | result |
|---|---|
| V1 robustness: the new decision the same on both paths for at least 99.5% of candidates, and no more flips than the old rule | **does not hold**: 99.66% the same (2,909 of 2,919), but 10 flips against the old rule's 3 |
| V2 no false holds against the labels | **does not hold as registered**: 4 (random effects, both paths), against the old rule's 2; see below -- the labels are wrong there |
| V2b missed holds (reported) | new 10, old 157 (random effects, each path); 0 in the negative binomial |
| V3 the band [0.025, 0.1] | dispersion 5 of 249, sigma 11 to 12 of 136 to 137 (8%), random-effect SD 70 of 2,526 (2.8%) |
| V4 unconverged | 29 flags: 9 beta fits (both paths) and 5 or 4 gaussian sigmas, and the curved CI case |
| V5 the four CI cases | k at its limit held (pushes ~4e-7); the flat beta held (largest 0.011); the curved beta unconverged (push at 6 of -261); the stalled sigma not held on P1 (+0.070 at 6) and held on P2 (+0.0048) |

### Where the failures sit

**Every flip is a gaussian sigma far below its line.** The random effects
(2,526 candidates) and the dispersions (252) flip 0 times under either rule.
The 9 sigma flips, and the stalled CI sigma, have sigma at 2e-5 to 3e-4 of
sd(y), against the hold-pass line of 0.2 of sd(y): already e^-6.5 to e^-9 below
the line. Pushing log sigma a further 6 down, to about 1e-8 of sd(y), moves
the objective by +-0.08 to 0.4 differently on the two paths: that is the
objective's arithmetic at a vanishing sigma, not its likelihood. The
unconverged flags are the same regime: gaussian sigmas at 2e-5 to 3e-4 of
sd(y), and beta fits with 1/sqrt(phi) of e^-10 to e^-15 (phi up to e^31),
where the beta density's arithmetic runs out (pushes of -53 to -7,184).

**The four V2 "false holds" are label failures.** Each is a binomial with an
AR(1) beside the random intercept. The dropped model's fit (the label) sits in
a worse local basin of the AR correlation, of the opposite sign (full fit
rho -1, -0.98, -0.59, -0.65; dropped 0.89, -0.56, 0.57, 0.50), and stays there
with 10 restarts from its own solution. The pushes, which hold rho at the full
fit's value, show the term flat: pushing the SD to about zero costs at most
0.0097, so the dropped model at the full fit's rho fits within that, and the
label's cost (0.05 to 1.7) is its local optimum's, not the model's. Two of the
four were the old rule's false holds as well. (Three of the four full fits put
rho at or near -1 for a truth of 0.8: the binary-AR observation in
illume-work notes.)

### The question for Craig: how far a push may go

Hold-pass sets the gaussian lines relative to sd(y) -- 0.2 for sigma, 0.1 for
an AR or random-effect SD -- and the dispersion line at 1/sqrt(phi) = 1e-2. A
term already far below its line is at its limit for every practical purpose:
its estimate is arbitrary and its SE meaningless. Pushing it further only
probes the objective's arithmetic, which differs between builds -- the whole
of V1's failure and nearly all of V4's flags.

**Post hoc, not a verdict** (chosen after the run): a floor at e^6 below the
line. A candidate whose estimate is already beyond it is held by value, with no
push; the others are judged by the pushes as ruled. On this run that holds 230
candidates by value (sigma 89, dispersion 25 per path, and the curved CI case),
and leaves 1 flip (the stalled CI sigma, which the floor would also hold) and 2
unconverged flags (one beta fit, both paths: push at 6 of -0.078, code 0).
Its cost: a fit that ran away past the floor -- the curved CI beta, phi e^28,
code 1, gradient 3.96 -- is held by value rather than reported unconverged;
its optimizer and gradient checks still FAIL.

**Caveats.** Two float paths on one machine stand in for CI's five
platforms. The labels exist only for the negative binomial and the random
effects. The floor variant was chosen after seeing the data.
