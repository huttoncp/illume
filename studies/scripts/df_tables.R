## Finite-sample degrees of freedom in the tables of a gaussian mixed model.
##
## PRE-REGISTERED DESIGN, committed before any code and before any run.
##
## ## The question
## A gaussian mixed model's tables -- summary()'s coefficients, ilm_anova()'s
##   Wald tests, ilm_emmeans()' intervals and ilm_contrast()'s rows -- use a
##   normal or chi-square reference today, and only ilm_trends() uses finite
##   degrees of freedom. With few clusters a between-cluster effect's z test
##   is anti-conservative. Craig ruled (D-DF1 to D-DF5): Satterthwaite's df by
##   default in all four tables, Kenward-Roger on request, ilm_anova() as F
##   with a denominator df (and `statistic = "Chisq"` for the old table), and
##   this study as the evidence for the default. The formulas themselves are
##   established by agreement with lmerTest and pbkrtest in the package tests;
##   the study measures what the reference does to the answers.
## It also answers one open question the tests alone cannot: what the df and
##   the coverage are when a variance component is held at its boundary under
##   rule C (the flat direction held, the rest of the covariance drawn from),
##   where Satterthwaite's gradient runs through a direction the fit holds.
##
## ## Arms and cells (gaussian, REML throughout, as KR requires)
## | Arm | Model | Cells |
## |---|---|---|
## | D1 between-cluster effect | y ~ trt + (1 \| g), trt constant within cluster, half the clusters treated | clusters G in {6, 10, 20} x cluster size m in {5, 20} x ICC in {0.05, 0.3} = 12 |
## | D2 within-cluster slope | y ~ x + (1 + x \| g), x varying within cluster, slope SD 0.3, correlation 0.3 | G in {6, 10, 20} x m in {5, 20} = 6 |
## | D3 held boundary | y ~ trt + (1 \| g), the intercept SD 0 in truth, so about half the fits hold it | G in {6, 10, 20}, m = 5 = 3 |
## | D4 multi-df F | y ~ f + (1 \| g), f a 3-level between-cluster factor, 2 clusters per level at least | G in {6, 12, 21}, m = 5, ICC 0.3 = 3 |
##
## 24 cells. Residual SD 1; ICC sets the intercept SD. The effect of interest
##   has a true value of 0.5 in D1 to D3 (coverage, and the size of a test of
##   the true value), and every level of f is equal in D4 (the size of the F
##   test). 1,000 replicates per cell, 24,000 fits.
##
## ## What is recorded, per fit
## - the estimate and its SE (and KR's adjusted SE);
## - df by each method: z (Inf), Satterthwaite, Kenward-Roger;
## - the 95% interval and the test of the true value by each, and in D4 the F
##   test's DenDF and p by Satterthwaite and KR and the chi-square's p;
## - whether a variance component was held, and the fit's checks;
## - time for each df method.
##
## ## What counts, fixed now
## - "Calibrated": coverage within 2 Monte Carlo SEs of 0.95 (the SE at 1,000
##   replicates is 0.0069, so 0.936 to 0.964); a test's size within 2 MC SEs
##   of 0.05.
## - The default is supported if Satterthwaite is calibrated in every D1 and
##   D2 cell with G >= 10, and closer to 0.95 than z in every cell with G = 6.
##   Where it is not calibrated at G = 6 it is reported, and whether KR is.
## - D3 answers the open question by the numbers: coverage and size in the
##   held fits and in the unheld fits separately, and whether any df is NaN,
##   infinite or below 1. If Satterthwaite is not calibrated in the held fits
##   while KR is, or the df are undefined there, that goes to Craig with a
##   proposal before the default ships for held fits.
## - D4: the F test's size by Satterthwaite and KR, and the chi-square's.
##
## ## Fresh seeds
## The whole design again at 500 replicates per cell with a seed offset, on
##   the same build: the verdicts above must hold there too, and any that does
##   not is reported as such.
##
## ## Build
## A pinned library built from the branch that carries the df change, shared
##   with the spatial validation study (item 115), whose gaussian cells depend
##   on the denominator df; versions and install times recorded before and
##   after the runs.
##
## ## Size
## About 24,000 fits and 12,000 fresh, each a few tenths of a second with KR
##   the largest part at G = 20, m = 20: about 3 to 5 core-hours in all.
##
## Usage: Rscript df_tables.R <nrep> <ncore> <outdir> [offset]
## ---------------------------------------------------------------------------
## (the study's code follows in a later commit; this commit is the design)
