## Rescaling the fixed-effect columns inside the fit (item 6).
##
## PRE-REGISTERED DESIGN, committed before the change is measured.
##
## ## The problem, measured before anything was changed
## A zero-inflated negative binomial with a random intercept and income in
##   dollars (mean 33,000, SD 9,000) in both parts: 12 of 12 fits failed,
##   the Hessian not positive definite and three standard errors NaN, while
##   the checks advised simplifying the covariance. The same data with income
##   in thousands: 11 of 12 fits graded ok, every SE finite. A zero-inflated
##   Poisson with income only in the count part survived in dollars, but
##   only through the Hessian's recomputation, and its income SE was 0.7%
##   off the thousands fit's. Columns on very different scales make the
##   optimiser's steps and the Hessian's differencing badly conditioned; the
##   model is the same model either way.
##
## ## The change (Craig's ruling, item 6)
## Inside ilm_fit(), each non-constant column j of the fixed design (and of
##   the zero part's and the dispersion model's) is divided by s_j, its SD --
##   and, in the centring variant, has its mean m_j taken off first. The
##   optimisation and the Hessian run on those columns. At the end the
##   estimates and their covariance are converted back EXACTLY, by the linear
##   map beta = A beta_s (+ the intercept's shift when centring):
##     scale only    beta_j = beta_s,j / s_j
##     centring      beta_j = beta_s,j / s_j for j > 1, and the intercept
##                   beta_1 = beta_s,1 - sum_j m_j beta_s,j / s_j
##   V = A V_s A'. The joint precision's fixed block is mapped the same way.
##   Everything downstream sees the user's scale, as a REML fit's already
##   does after its own reshaping: opt$par and sdr$cov.fixed are rewritten
##   there, and a tape on the user's columns is built for code that
##   evaluates the objective at user-scale parameters (as obj_ml is for
##   REML). The remedies name rescaling first when a fit fails on columns
##   whose scales differ by more than a factor of 1e3.
## A model without an intercept is not centred (there is no intercept to
##   absorb the shift); nor is a column that is an indicator (0/1), whose
##   scale is already fine.
##
## ## Scope, stated (and said in the help)
## - Offsets are never rescaled: they enter the linear predictor with a
##   fixed coefficient of 1.
## - A smooth's penalised part is untouched. Its unpenalised (null-space)
##   columns are ordinary fixed columns, rescaled like any other; a test
##   shows the smooth's fitted values and edf unchanged against HAND.
## - Random-slope covariates are out of scope for this pass: their scale
##   sits in the random-effect covariance, not in a fixed column. A badly
##   scaled one keeps the remedy that names rescaling by hand.
## - The separation check's flat catch ("beyond 8 on the link scale") runs
##   on the internal standardised coefficients, where 8 means the same for
##   every column; a test shows the separation cases caught identically
##   before and after the change.
##
## ## The study: scale only against scale and centre
## | Arm | Model | Cells |
## |---|---|---|
## | R1 | the problem case: ZI negative binomial, (1 \| g), income in dollars in both parts | n in {300, 600} |
## | R2 | Poisson, (1 \| g), one covariate in dollars and one of SD 1e-4 | n in {300, 600} |
## | R3 | binomial, (1 \| g), income in dollars | n in {300, 600} |
## | R4 | well-scaled control: each family above with covariates of SD about 1 | 3 families |
## | R5 | a smooth beside a dollar-scale covariate: Poisson, y ~ income + s(x) + (1 \| g) | n = 600 |
## | R6 | out of scope, recorded: Poisson with a random slope on income in dollars, (1 + income \| g) | n = 600 |
##
## 11 cells, 100 replicates each. R5 checks exactness through the smooth
##   machinery (its fitted values and edf as well as the coefficients). R6
##   records the out-of-scope behaviour: the remedy naming rescaling by hand
##   fires, and nothing is silently wrong -- a fit graded ok there has finite
##   SEs that agree with HAND to E1's tolerance. Every data set is fitted four ways: the
##   current build on the raw columns (BEFORE), the current build on columns
##   the user rescaled by hand to SD about 1 (HAND), and the new build on the
##   raw columns in each variant (SCALE, CENTRE).
##
## What is recorded: whether the fit is graded ok; the count of non-finite
##   SEs; each coefficient and SE on the user's scale; the objective; the
##   checks that are not OK; time.
##
## ## What counts, fixed now
## - E1, exactness: wherever HAND and a new variant both converge (ok), their
##   user-scale coefficients and SEs agree to 1e-6 relative, the intercept
##   included under centring, and their objectives to 1e-8. This is the
##   stated test of the conversion back.
## - E2, the failures are fixed: in R1 to R3, each new variant is graded ok
##   in at least as many fits as HAND, less 2 (a variant that loses more
##   than 2 fits of 100 against the hand rescaling does not pass).
## - The noise floor, measured first: in R4, BEFORE refitted from a jittered
##   start (every parameter moved by N(0, 0.1)), against BEFORE from the
##   default start. The floor is the 99th percentile over fits of the
##   largest relative difference in a coefficient or SE.
## - E3, nothing well-scaled changes: in R4, each variant's coefficients and
##   SEs agree with BEFORE within twice that floor, and its grade is the same
##   in every fit. (E1 compares nearly identical problems and keeps 1e-6.)
## - R6, recorded: in every fit where the random slope's covariate is badly
##   scaled and the fit is not ok, the remedies name rescaling; no fit graded
##   ok differs from HAND beyond E1's tolerance.
## - The default is scale only if it passes E1 to E3, since it leaves the
##   intercept exactly where the user's model puts it; centring becomes the
##   default only if scale only fails E2 and centring passes all three. If
##   neither passes, that goes to Craig with the numbers before anything
##   ships.
##
## ## Build
## BEFORE and HAND on the pinned library of the branch the change is made on,
##   before the change (its parent commit); SCALE and CENTRE on a pinned
##   library of the change, the variant chosen by an internal option for the
##   study only. Manifests for both.
##
## ## Build and running
## One pinned library built from the branch's change commit, with its
##   manifest, and one of its parent for BEFORE and HAND; the runs as a
##   sequential chain on one core after the queued full checks, with each
##   fit's check lines recorded.
##
## ## Size
## 1,100 data sets, four fits each and the jittered refits of R4, about a
##   second a fit: about an hour and a half on one core.
##
## Usage: Rscript rescale.R <arm> <nrep> <outdir>
## ---------------------------------------------------------------------------
## (the study's code follows in a later commit; this commit is the design)
