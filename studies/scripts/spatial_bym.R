## Areal spatial terms: the BYM form, validated (item 115).
##
## PRE-REGISTERED DESIGN, committed before any code and before any run.
##
## ## The question
## The disease-surveillance case study fits a BYM-form areal term: an
##   intrinsic CAR field and an unstructured region effect, written in illume
##   as s(region, bs = "mrf", xt = list(nb = nb)) + (1 | region), the first
##   an mgcv Markov random field smooth (Wood 2017, 5.8.1), the second the
##   package's own random intercept. Nothing in illume's evidence covers it.
##   This study asks, for data generated from that model:
##   Q1 do the fixed effects' 95% intervals cover (the covariate's slope);
##   Q2 do the region-level relative risks' 95% intervals cover, the quantity
##      a disease map reports;
##   Q3 does illume agree with the reference fit of the same model;
##   Q4 how the fit's checks grade these models -- one count per region is
##      the rule in disease mapping, and obs_per_level[region] grades the
##      unstructured term FAIL there (seen in the design probe), so how many
##      fits are graded usable, and what coverage is among the ones that
##      are not;
##   Q5 whether ilm_variogram(coords = ) finds the spatial structure a model
##      without the term leaves behind, names the term that removes it, and
##      stops flagging once the term is in.
##
## ## The reference
## mgcv::gam() with the same Markov random field basis and the unstructured
##   effect as s(region, bs = "re"), method = "REML" (Laplace-approximate REML
##   for a count family). mgcv is a recommended package, already in the pinned
##   library lib-dfs (1.9-1); nothing is installed. It is the same model and
##   the same basis, so it tests illume's fitting, not the choice of prior.
## A disease-mapping reference proper -- an MCMC BYM fit, CARBayes::S.CARbym
##   on 100 replicates of two S1 cells, comparing the slope and the relative-
##   risk map -- is PENDING Craig's ruling (item 200). It needs a new package
##   and its dependencies in a pinned library, so nothing is installed until
##   he rules and the install is cleared. If it is added, it is specified in
##   an addendum committed before it runs; the arms below do not depend on
##   it.
##
## ## The maps
## Rook-adjacency lattices, 6 x 6 (R = 36) and 10 x 10 (R = 100), and the
##   centroids as coordinates. The structured field phi is drawn from the
##   intrinsic CAR with precision D - W on its non-null space and scaled to
##   its target SD; the unstructured theta ~ N(0, sd). A new field every
##   replicate. The covariate x ~ N(0, 1) per region, independent of phi, so
##   spatial confounding is not what is measured; slope 0.3.
##
## ## Arms and cells
## | Arm | Data | Cells |
## |---|---|---|
## | S1 counts | one Poisson count per region, y ~ Poisson(E exp(-0.2 + 0.3 x + phi + theta)), offset log(E), E ~ Gamma with mean Ebar | R in {36, 100} x Ebar in {5, 50} x total SD in {0.3, 0.6} x structured share rho in {0.2, 0.8} = 16 |
## | S2 gaussian | m observations per region, y = 1 + 0.3 x + phi + theta + N(0, 1), x per observation | R in {36, 100} x m in {3, 10} x rho in {0.2, 0.8}, total SD 0.6 = 8 |
## | S3 variogram | S1's counts at Ebar 50, rho 0.8: the model without the spatial term, y ~ x + offset(log(E)) + (1 | region), then the BYM model; ilm_variogram(coords = centroids, B = 39) on each | R in {36, 100} x total SD in {0.3, 0.6} = 4 |
##
## The total SD is sqrt(var(phi) + var(theta)) on the log (S1) or response
##   (S2) scale; rho = var(phi) / total. 500 replicates per S1 and S2 cell
##   (12,000 fits); the mgcv reference on the first 100 of each (2,400
##   fits); S3 at 100 replicates per cell, the BYM refits in its variogram
##   only at R = 36 (the R = 100 refits cost about 5 s each, times 39).
##
## ## What is recorded, per fit
## - the slope's estimate, SE, and 95% interval (z for S1; Satterthwaite t for
##   S2, the df tables' default) and whether it covers 0.3;
## - each region's relative risk (S1: exp of the linear predictor without the
##   offset; S2: the region mean) with its 95% interval from predict(interval
##   = "confidence", nsim = 500), and the share of regions covered;
## - the SDs of phi and theta as estimated, whether either was held at its
##   boundary, the checks' statuses (obs_per_level, latent_budget, hessian,
##   variance_boundary, optimizer, gradient) and ok;
## - each unheld variance parameter's SE on its log scale, from vcov(fit,
##   full = TRUE), and the fit's variance flag: any of them above 0.75, the
##   cut-off Arm B chose. The region intervals come from joint draws, and
##   this build predates the remedy for flagged fits (item 170); Arm B found
##   joint draws explode on them, so region coverage is reported by flag;
## - on the reference replicates, mgcv's slope, SE and the two SDs;
## - S3: the variogram's verdict, the bins flagged, and the advice line;
## - time for each fit.
##
## ## What counts, fixed now
## - "Calibrated": within 2 Monte Carlo SEs of 0.95. At 500 replicates the
##   slope's band is 0.931 to 0.969. The region-level coverage is a share of
##   R correlated regions per fit, so its MC SE is the SD of the per-fit
##   shares over replicates divided by sqrt(500), and its band is computed
##   from that.
## - Coverage is reported over ALL fits (the primary figure: a user sees a
##   fit whether or not a check fails) and separately over the fits graded
##   ok and not ok, held and unheld, and variance-flagged and not.
## - V1: the slope calibrated in every S1 cell with Ebar = 50 and every S2
##   cell. V2: the slope calibrated in every S1 cell with Ebar = 5 (reported
##   either way; low counts are where a Laplace fit is least sure).
## - V3: the region relative risks calibrated in every S1 and S2 cell.
## - V4: agreement with mgcv in every cell: the slope's median |illume -
##   mgcv| / mgcv's SE at most 0.1, and the median ratio of SEs within 0.9
##   to 1.1. Gaussian cells, where illume's REML and mgcv's are the same
##   criterion, fitted with reml = TRUE, are held to 1e-3 on the slope.
## - V5: S3 flags the model without the spatial term in at least 80% of
##   fits at total SD 0.6; V6: after the BYM term, flags in at most 10% (the
##   variogram's own verdict, at its default min_effect); V7: every flagged
##   fit's advice names the Markov random field term.
## - Q4: the share graded ok per cell, and the coverage in each grade. The
##   not-ok fits "cover as well as the ok ones" when, in every S1 cell where
##   both grades hold at least 50 fits, the slope's coverage and the region
##   coverage differ between the grades by less than 2 Monte Carlo SEs of
##   the difference (each grade's MC SE from its own count, combined in
##   quadrature). In a cell with fewer than 50 fits graded ok -- every fit
##   there, if the check fails them all, as the probe suggests -- the not-ok
##   fits are compared with nominal instead: they cover as well when they
##   are calibrated by V1's and V3's rules.
## - If V1 or V3 fails in a cell with Ebar = 50, or the not-ok fits cover as
##   well as the ok ones in that sense while the checks call them unusable,
##   that goes to Craig with a proposal before the case study relies on the
##   term or the check is changed.
##
## ## Fresh seeds
## S1 and S2 again at 250 replicates per cell (band 0.923 to 0.977) with a
##   seed offset, no reference, on the same build; the verdicts must hold
##   there too, and any that does not is reported as such.
##
## ## Build
## S1 and S2 need nothing beyond the df tables' build, and run on its pinned
##   library lib-dfs (illume acf93c2 with mgcv 1.9-1), counts and gaussian
##   alike -- the gaussian cells' t intervals are the df tables' default,
##   which that build carries. S3 tests the variogram's advice, which this
##   item rewrites to name the areal term (today it says illume fits no
##   spatial covariance), so it runs on a later pinned library built from
##   the branch with that change, with its own manifest.
##
## ## Size
## S1 and S2: 12,000 illume fits, about 1 s at R = 36 and 5 s at R = 100,
##   so about 10 core-hours with the predictions; 2,400 mgcv fits, about 3;
##   fresh seeds about 5; S3 about 3. Some 21 core-hours in all.
##
## Usage: Rscript spatial_bym.R <arm> <nrep> <ncore> <outdir> [offset]
## ---------------------------------------------------------------------------
## (the study's code follows in a later commit; this commit is the design)
