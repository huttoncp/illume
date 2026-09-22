# messy study -- findings log

One section per package version, newest first. Raw per-replicate data
is kept only for the three most recent versions (see ../README.md);
this file is the permanent record. Unlike the other findings files it
is written from `messy_compare.csv` rather than by
`scripts/summarise_run.R`, which has no branch for this study yet.

## 0.0.7.9000

Run on 2026-09-22, R 4.4.3. 400 replications per regime, 10 cores.

illume (Laplace) against mclogit::mblogit (PQL) on data that misbehaves.
`mclogit_compare.md` asked the textbook question on tidy data; this varies
the four things most likely to break either approximation, one at a time
and then together. Nominal coverage 0.95; Monte Carlo error about 0.011.
Attenuation is the least-squares slope of estimate on truth, so 1 means
none, below 1 is shrinkage and above 1 is inflation.

Base design: 60 clusters x 8 observations, J = 3, random-effect sd 1.0.

| regime | cover_illume | cover_mclogit | atten_illume | atten_mclogit | bias_illume | bias_mclogit | rmse_illume | rmse_mclogit | conv_illume | conv_mclogit |
|---|---|---|---|---|---|---|---|---|---|---|
| clean | 0.950 | 0.937 | 1.011 | 0.864 | 0.008 | 0.102 | 0.321 | 0.289 | 1.000 | 1.000 |
| unbalanced | 0.943 | 0.932 | 1.006 | 0.881 | 0.024 | 0.104 | 0.349 | 0.308 | 1.000 | 1.000 |
| rare | 0.931 | 0.853 | 1.462 | 0.816 | 0.394 | 0.271 | 1.477 | 1.598 | 0.973 | 1.000 |
| flat_re | 0.948 | 0.947 | 1.018 | 1.008 | 0.012 | 0.008 | 0.249 | 0.243 | 0.633 | 1.000 |
| heavy_re | 0.936 | 0.902 | 0.982 | 0.846 | 0.059 | 0.130 | 0.326 | 0.301 | 1.000 | 1.000 |
| combined | 0.945 | 0.894 | 1.301 | 0.918 | 0.370 | 0.368 | 1.681 | 1.732 | 0.910 | 1.000 |

Variance-component recovery, and runtime:

| regime | true re_sd | re_sd illume | re_sd mclogit | secs illume | secs mclogit |
|---|---|---|---|---|---|
| clean | 1.00 | 0.931 | 1.166 | 0.55 | 4.27 |
| unbalanced | 1.00 | 0.946 | 1.241 | 0.57 | 5.59 |
| rare | 1.00 | 1.493 | 1.207 | 0.55 | 5.43 |
| flat_re | 0.05 | 0.099 | 0.004 | 0.60 | 3.43 |
| heavy_re | 1.00 | 0.873 | 1.227 | 0.52 | 4.16 |
| combined | 1.00 | 0.867 | 1.150 | 0.50 | 6.49 |

Realised rarest-category share: 0.259, 0.245, **0.036**, 0.180, 0.284,
**0.025**. The rare regime's intercepts were calibrated to land there rather
than guessed.

### What holds

illume's intervals cover better than mclogit's in all six regimes, and the
gap widens exactly where theory says PQL should struggle: 0.931 against
0.853 with a 3.6% category, 0.945 against 0.894 with everything combined.
mclogit shrinks fixed effects in five of six regimes (attenuation 0.82 to
0.92), and its reported standard errors never reflect it. illume is 6 to 13
times faster throughout.

Non-Gaussian random effects cost illume less than expected. Both methods
assume Gaussian REs, so `heavy_re` is misspecification for both, and illume's
attenuation stays at 0.982 against mclogit's 0.846.

### Three findings that go AGAINST illume, and must travel with the rest

**1. illume inflates coefficients when a category is sparse.** Attenuation
1.462 in `rare` and 1.301 in `combined` -- the opposite direction from PQL's
shrinkage, but still a bias, and in `rare` illume's mean absolute bias (0.394)
is *worse* than mclogit's (0.271). Its intervals are wide enough to cover
anyway, which is what honest uncertainty buys, but the point estimate in a
sparse-cell multinomial should not be trusted at face value. This is the
separation-adjacent regime: with ~17 of 480 observations in a category spread
over 60 clusters, the likelihood is close to flat and the MLE drifts outward.

**2. illume's convergence collapses when a variance component sits near the
boundary.** In `flat_re` (true sd 0.05) illume converged on 63.3% of
replications against mclogit's 100%. Its 0.948 coverage there describes only
those 63.3% and is conditional on convergence. Neither method recovers the
variance: illume reports a median 0.099 (twice the truth), mclogit 0.004
(driven almost to zero). Both are wrong; they are wrong in opposite
directions.

**3. mclogit's RMSE is lower in four of the six regimes** (clean, unbalanced,
flat_re, heavy_re), by 8 to 13%. Shrinkage buys variance reduction. illume's
RMSE is lower only in the two sparse-cell regimes, and there only narrowly.
The existing caveat from `mclogit.md` is confirmed, not overturned: if you
are predicting, mclogit may be the better choice, and the case for illume is
about inference rather than accuracy in every sense.

### Reading

The overall picture is that the two methods fail differently rather than one
dominating. PQL is reliable and biased: it always converges, always shrinks,
and never says so. The Laplace approximation is honest and fragile: its
intervals mean what they claim, and it declines to answer more often --
loudly, which is the design intent (see the latent budget in the regression
vignette), but declining is still a cost that a convergence-conditional table
hides. A user with a near-boundary variance component or a 3% category should
read illume's convergence status, not just its coefficients.
