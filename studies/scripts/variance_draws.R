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
## Usage: Rscript variance_draws.R <nrep> <ncore> <outdir> [offset] [cells]
## ---------------------------------------------------------------------------

args   <- commandArgs(trailingOnly = TRUE)
NREP   <- if (length(args) >= 1) as.integer(args[1]) else 100L
NCORE  <- if (length(args) >= 2) as.integer(args[2]) else 1L
sp     <- if (length(args) >= 3) args[3] else "."
## a seed offset, so the confirmation runs on data the cut-off was not chosen on
OFFSET <- if (length(args) >= 4) as.integer(args[4]) else 0L
## a subset of cells, for the smoke run only
ONLY   <- if (length(args) >= 5) as.integer(strsplit(args[5], ",")[[1]]) else NULL
dir.create(sp, showWarnings = FALSE, recursive = TRUE)

## ---- the cells -------------------------------------------------------------------
g <- function(...) expand.grid(..., stringsAsFactors = FALSE)
cells <- rbind(
  data.frame(arm = "B1", g(family = c("gaussian", "poisson", "nbinom"), Tn = c(12, 24, 48),
                           rho = c(0.3, 0.7), lat_sd = c(0.3, 0.8)), ri_sd = NA, L = NA, m = NA),
  data.frame(arm = "B2", g(family = c("gaussian", "poisson", "nbinom", "binomial"),
                           Tn = c(12, 24), ri_sd = c(0.2, 0.5)), rho = 0.7, lat_sd = 0.5,
             L = NA, m = NA),
  data.frame(arm = "B3", family = rep(c("gaussian", "poisson", "nbinom", "binomial"), 6),
             Tn = NA, rho = NA, lat_sd = NA,
             ri_sd = rep(rep(c(0.1, 0.5), each = 4), 3),
             L = rep(c(5, 10, 30), each = 8), m = rep(c(10, 5, 2), each = 8)))
cells$noise <- ifelse(cells$family == "gaussian", 0.5, NA)
cells$cell <- seq_len(nrow(cells))
stopifnot(nrow(cells) == 76L)
H <- 6L; HS <- c(1L, 3L, 6L)             # horizons simulated, and reported
NB_SIZE <- 2
MU0 <- c(gaussian = 5, poisson = 1.5, nbinom = 1.5, binomial = 0)
CUTS <- c(0.5, 0.75, 1, 1.5, 2, 3)       # the flag's candidate cut-offs on s
NDRAW <- 1000L

## ---- simulation: the data, the truth, and the future ------------------------------
ar_path <- function(n, rho, sd_marg)
  as.numeric(stats::arima.sim(list(ar = rho), n, sd = sd_marg * sqrt(1 - rho^2)))
## the response given its linear predictor; `disp` is sigma or the NB size
draw_y <- function(fam, eta, disp) switch(fam,
  gaussian = eta + stats::rnorm(length(eta), 0, disp),
  poisson  = stats::rpois(length(eta), exp(pmin(eta, 700))),
  nbinom   = stats::rnbinom(length(eta), mu = exp(pmin(eta, 700)), size = disp),
  binomial = stats::rbinom(length(eta), 1, stats::plogis(eta)))
true_disp <- function(ce) switch(ce$family, gaussian = ce$noise, nbinom = NB_SIZE, NULL)
## each row keeps its latent `a` and the rest of its true linear predictor `m`,
## so the true-parameter predictive can start from the truth at the last time
simulate <- function(ce) {
  mu <- MU0[[ce$family]]
  if (ce$arm == "B3") {
    d <- data.frame(g = factor(rep(seq_len(ce$L), each = ce$m)))
    d$x <- stats::rnorm(nrow(d))
    d$y <- draw_y(ce$family, mu + 0.3 * d$x + stats::rnorm(ce$L, 0, ce$ri_sd)[d$g],
                  true_disp(ce))
    return(list(d = d, fut = NULL, fml = y ~ x + (1 | g), ar = FALSE))
  }
  G <- if (ce$arm == "B1") 1L else 20L
  n <- ce$Tn + H
  a <- lapply(seq_len(G), function(i) ar_path(n, ce$rho, ce$lat_sd))
  b <- if (ce$arm == "B2") stats::rnorm(G, 0, ce$ri_sd) else rep(0, G)
  full <- do.call(rbind, lapply(seq_len(G), function(i)
    data.frame(g = i, t = seq_len(n), a = a[[i]], m = mu + b[i])))
  full$g <- factor(full$g)
  full$y <- draw_y(ce$family, full$m + full$a, true_disp(ce))
  list(d = full[full$t <= ce$Tn, ], fut = full[full$t > ce$Tn, ],
       fml = if (ce$arm == "B1") y ~ 1 else y ~ 1 + (1 | g), ar = TRUE)
}

## ---- the fit's variance parameters, on their log scales ---------------------------
## `id` indexes the fixed parameters (opt$par, vcov), `row` the draws' full
## layout: the same parameter, the first of its block
var_params <- function(f, ce, map) {
  pn <- names(f$opt$par); out <- list()
  add <- function(nm, block, truth, held) {
    out[[nm]] <<- list(id = which(pn == block)[1L],
                       row = if (is.null(map)) NA_integer_ else map$row[map$block == block][1L],
                       truth = truth, held = held)
  }
  if ("lchol_ar" %in% pn) add("ar_sd", "lchol_ar", log(ce$lat_sd), "ar")
  if ("theta" %in% pn) add("ri_sd", "theta", log(ce$ri_sd), "g")
  if ("logdisp" %in% pn && ce$family == "gaussian")
    add("sigma", "logdisp", log(ce$noise), "dispersion")
  if ("logdisp" %in% pn && ce$family == "nbinom")
    add("nb_logk", "logdisp", log(NB_SIZE), "dispersion")
  out
}
## the objective with some parameters pinned, every other fixed parameter
## re-optimised (the random effects by the Laplace inner step). A log SD the
## fit left near zero (below -5; a NB log k above 5) has no gradient to leave
## by, so the re-optimisation also starts once with those at -1 (+1), and the
## better of the two is kept
pinned_obj <- function(f, id, value) {
  p <- f$opt$par; pn <- names(p)
  fn <- function(q) { pp <- p; pp[-id] <- q; pp[id] <- value; f$obj$fn(pp) }
  run <- function(st) {
    o <- tryCatch(suppressWarnings(stats::nlminb(st, fn)), error = function(e) NULL)
    if (is.null(o) || !is.finite(o$objective)) Inf else o$objective
  }
  best <- run(p[-id])
  st <- p; sd_par <- pn %in% c("logdisp", "lchol_ar", "theta")
  st[sd_par & p < -5] <- -1; st[pn == "logdisp" & p > 5] <- 1
  if (any(st[-id] != p[-id])) best <- min(best, run(st[-id]))
  invisible(tryCatch(f$obj$fn(p), error = function(e) NULL))
  if (is.finite(best)) best else NA_real_
}
## a 25-point profile of one log SD and its 95% interval's upper end, and the
## label: unbounded, wald_too_wide (Wald's upper end over twice the profile's
## on the SD scale), ok, or profile_failed where a grid point's re-optimisation
## failed before the interval's end was reached
PGRID <- seq(-4, 6, length.out = 25L)
profile_one <- function(f, id, se) {
  est <- f$opt$par[[id]]; f0 <- f$obj$fn(f$opt$par)
  psi <- est + PGRID
  lp <- -vapply(psi, function(v) pinned_obj(f, id, v), 0)
  dev <- 2 * ((-f0) - lp)
  ## walk up from the estimate (deviance 0) to the first point past 3.84
  up <- c(est, psi[psi > est]); du <- c(0, dev[psi > est])
  k <- 1L
  while (k < length(up) && is.finite(du[k + 1L]) && du[k + 1L] <= 3.84) k <- k + 1L
  prof_up <- NA_real_
  label <- if (k == length(up)) "unbounded"
           else if (!is.finite(du[k + 1L])) "profile_failed"
           else {
             prof_up <- up[k] + (3.84 - du[k]) / (du[k + 1L] - du[k]) * (up[k + 1L] - up[k])
             if (is.finite(se) && est + 1.96 * se - prof_up > log(2)) "wald_too_wide" else "ok"
           }
  ## the grid with the estimate, the profile's peak, in its place
  o <- order(c(psi, est))
  list(psi = c(psi, est)[o], lp = c(lp, -f0)[o], est = est, prof_up = prof_up,
       wald_up = est + 1.96 * se, label = label)
}
## a 9 x 9 joint profile of two log SDs, others re-optimised, taken where both
## one-parameter profiles fail. RIDGE: the total variance is bounded, taken as
## the total variance staying within 4 times its value at the estimate over
## the pair's joint 95% region (deviance <= 5.99)
R9 <- seq(-4, 6, length.out = 9L)
ridge_pair <- function(f, id1, id2) {
  e <- f$opt$par[c(id1, id2)]; f0 <- f$obj$fn(f$opt$par)
  gr <- expand.grid(a = e[[1]] + R9, b = e[[2]] + R9)
  lp <- -vapply(seq_len(nrow(gr)), function(i)
    pinned_obj(f, c(id1, id2), c(gr$a[i], gr$b[i])), 0)
  dev <- 2 * ((-f0) - lp)
  V <- exp(2 * gr$a) + exp(2 * gr$b); V0 <- sum(exp(2 * e))
  inreg <- is.finite(dev) & dev <= 5.99
  list(grid = gr, lp = lp, ridge = any(inreg) && max(V[inreg]) <= 4 * V0)
}

## ---- remedy 2: draws from the profile, capped, others conditional ----------------
## the cap, fixed in the design: a latent SD at most 5 x max(sd(g(y*)), 0.5)
cap_of <- function(fam, y) {
  ys <- switch(fam, gaussian = y, poisson = , nbinom = log(y + 0.5),
               binomial = stats::qlogis((y + 0.5) / 2))
  log(5 * max(stats::sd(ys), 0.5))
}
## the normalised exp(profile), linearly interpolated on a fine grid below the
## cap (a spline overshoots wildly beside a steep point). `trim` leaves out
## that much above the grid's lower edge: the lower-plateau check only
sample_profile <- function(pr, n, cap, trim = 0) {
  ok <- is.finite(pr$lp) & pr$psi <= cap
  if (sum(ok) < 3L) return(rep(min(pr$est, cap), n))
  fine <- seq(min(pr$psi[ok]), max(pr$psi[ok]), length.out = 400L)
  l <- stats::approx(pr$psi[ok], pr$lp[ok], xout = fine)$y
  keep <- fine >= min(pr$psi[ok]) + trim
  if (!any(keep)) keep <- fine >= max(fine)
  sample(fine[keep], n, replace = TRUE, prob = exp(l[keep] - max(l[keep])))
}
## THE LOWER PLATEAU -- EXPLORATORY, POST-REGISTRATION (added after the smoke
## run, before any main-run data, at the conductor's request; it measures the
## registered grid, and changes no label, flag or registered remedy). Where the
## profile is flat at the grid's lower edge, remedy 2 spreads mass down to
## est - 4, so the edge acts as a floor. Recorded per parameter: whether the deviance
## changes by less than 0.1 over the lowest grid step, the deviance there, and
## the share of remedy 2's draws within one grid step of the edge. The forecasts
## add remedy 2 with that step left out ("profile_trim"), from the same seed as
## remedy 2, so the difference in CRPS is paired by fit and seed.
STEP <- diff(PGRID)[1]
low_edge <- function(pr, cap, draws) {
  ok <- which(is.finite(pr$lp) & pr$psi <= cap)
  if (length(ok) < 2L) return(c(low_plateau = NA, low_dev = NA, edge_mass = NA))
  top <- max(pr$lp[ok]); a <- ok[1L]; b <- ok[2L]
  c(low_plateau = abs(pr$lp[b] - pr$lp[a]) * 2 < 0.1, low_dev = 2 * (top - pr$lp[a]),
    edge_mass = mean(draws < pr$psi[a] + STEP))
}
## a pair from the 9 x 9 profile, uniform within the chosen grid cell
sample_ridge <- function(rp, n, cap) {
  ok <- is.finite(rp$lp) & rp$grid$a <= cap & rp$grid$b <= cap
  if (!any(ok)) return(NULL)
  k <- sample(which(ok), n, replace = TRUE, prob = exp(rp$lp[ok] - max(rp$lp[ok])))
  step <- diff(R9)[1]
  rbind(rp$grid$a[k] + stats::runif(n, -step / 2, step / 2),
        rp$grid$b[k] + stats::runif(n, -step / 2, step / 2))
}
## rows `rows` of a draw matrix replaced by `new` (rows x draws), every other
## row moved by its regression on them: the Gaussian conditional, from the
## joint draws' own covariance
condition_on <- function(D, rows, new) {
  S <- stats::cov(t(D)); others <- setdiff(seq_len(nrow(D)), rows)
  Srr <- S[rows, rows, drop = FALSE]
  if (any(!is.finite(Srr)) || min(eigen(Srr, TRUE, TRUE)$values) <= 1e-12) {
    D[rows, ] <- new; return(D)
  }
  B <- S[others, rows, drop = FALSE] %*% solve(Srr)
  D[others, ] <- D[others, ] + B %*% (new - D[rows, , drop = FALSE])
  D[rows, ] <- new
  D
}

## ---- forecasts ------------------------------------------------------------------------
## each draw carries its series' AR latent forward from its last cell (B_ar is
## added to the linear predictor as it stands; rho = tanh(rho_raw); lchol_ar is
## the log marginal SD), and the response is drawn at each horizon: scored
## against the future simulated with the data, and against the true-parameter
## predictive
crps_sample <- function(x, y) {
  x <- x[is.finite(x)]; if (length(x) < 2L) return(Inf)
  x <- sort(x); n <- length(x)
  mean(abs(x - y)) - sum((2 * seq_len(n) - n - 1) * x) / n^2
}
last_rows <- function(f, map) {
  ce <- ilm_cells(f)
  vapply(levels(droplevels(ce$group)), function(s) {
    k <- which(ce$group == s); k <- k[which.max(ce$time[k])]
    map$row[map$block == "B_ar" & map$cell == k][1L]
  }, 0L)
}
forecast_score <- function(D, map, ce, fut, ref, lastrow) {
  blk <- map$block
  b0 <- D[blk == "beta", , drop = FALSE][1L, ]
  sd_a <- exp(D[blk == "lchol_ar", ]); rho <- tanh(D[blk == "rho_raw", ])
  disp <- if (any(blk == "logdisp")) exp(D[blk == "logdisp", ]) else NULL
  rows <- list()
  for (s in names(lastrow)) {
    a <- D[lastrow[[s]], ]
    k <- which(blk == "bvec" & map$level == s)
    re <- if (length(k)) D[k[1L], ] else 0
    fs <- fut[fut$g == s, ]; fs <- fs[order(fs$t), ]
    for (h in seq_len(H)) {
      a <- rho * a + stats::rnorm(length(a), 0, sd_a * sqrt(pmax(1 - rho^2, 0)))
      x <- draw_y(ce$family, b0 + re + a, disp)
      if (!(h %in% HS)) next
      y <- fs$y[h]
      q <- stats::quantile(x, c(0.025, 0.1, 0.9, 0.975), na.rm = TRUE, names = FALSE)
      cr <- crps_sample(x, y)
      rows[[length(rows) + 1L]] <- data.frame(s = s, h = h, in80 = y >= q[2] & y <= q[3],
        in95 = y >= q[1] & y <= q[4], crps = cr, ref = ref[[paste(s, h)]],
        explosive = cr > 100 * ref[[paste(s, h)]])
    }
  }
  do.call(rbind, rows)
}
## the true-parameter predictive, from the true latent at the last fitted time
truth_ref <- function(ce, d, fut, n = NDRAW) {
  out <- list()
  for (s in levels(droplevels(fut$g))) {
    ds <- d[d$g == s, ]; ds <- ds[order(ds$t), ]
    fs <- fut[fut$g == s, ]; fs <- fs[order(fs$t), ]
    a <- rep(ds$a[nrow(ds)], n)
    for (h in seq_len(H)) {
      a <- ce$rho * a + stats::rnorm(n, 0, ce$lat_sd * sqrt(1 - ce$rho^2))
      out[[paste(s, h)]] <- crps_sample(draw_y(ce$family, fs$m[h] + a, true_disp(ce)), fs$y[h])
    }
  }
  out
}
summ_fc <- function(sc, remedy, cut = NA_real_) {
  if (is.null(sc) || !nrow(sc)) return(NULL)
  do.call(rbind, lapply(HS, function(h) { z <- sc[sc$h == h, ]
    data.frame(remedy = remedy, cut = cut, h = h, n_rows = nrow(z), cov80 = mean(z$in80),
               cov95 = mean(z$in95), crps_mean = mean(z$crps),
               crps_median = stats::median(z$crps), ref_mean = mean(z$ref),
               n_explosive = sum(z$explosive)) }))
}
qint <- function(x) stats::quantile(x, c(0.025, 0.975), names = FALSE, na.rm = TRUE)

## ---- one fit ----------------------------------------------------------------------------
one <- function(ce, fixture = NULL) {
  suppressMessages(library(illume))
  seed <- if (is.null(fixture)) 9973L * ce$cell + ce$rep + OFFSET else 1L
  set.seed(seed)
  s <- if (is.null(fixture)) simulate(ce) else fixture$sim
  key <- data.frame(cell = ce$cell, rep = ce$rep, seed = seed,
                    fixture = if (is.null(fixture)) NA_character_ else fixture$name)
  t0 <- proc.time()[["elapsed"]]
  f <- tryCatch(suppressMessages(suppressWarnings(ilm_model(s$fml, data = s$d,
        family = ce$family, ar = if (s$ar) ilm_ar1(~ t | g) else NULL, verbose = FALSE))),
        error = function(e) NULL)
  t_fit <- proc.time()[["elapsed"]] - t0
  if (is.null(f)) return(list(par = cbind(key, ok = FALSE), fc = NULL, time = NULL))
  held <- f$hessian_held
  ck <- f$checks
  nonconv <- any(ck$status[ck$check %in% c("gradient", "optimizer")] == "FAIL")
  V <- tryCatch(suppressWarnings(vcov(f, full = TRUE)), error = function(e) NULL)
  ## d: the other package's signal, as its cross-tab computes it
  dr <- tryCatch(suppressWarnings(ilm_draws(f, nsim = 400L, seed = 1L)),
                 error = function(e) NULL)
  ## today's joint draws, the baseline and remedy 1 (a)
  t1 <- proc.time()[["elapsed"]]
  base <- tryCatch(suppressWarnings(ilm_draws(f, nsim = NDRAW, seed = 2L, natural = FALSE)),
                   error = function(e) NULL)
  t_draws <- proc.time()[["elapsed"]] - t1
  map <- if (!is.null(base)) base$map else if (!is.null(dr)) dr$map else NULL
  vp <- var_params(f, ce, map)
  cap <- cap_of(ce$family, s$d$y)
  prs <- list(); prow <- list(); t_prof <- 0
  for (nm in names(vp)) {
    p <- vp[[nm]]; is_held <- p$held %in% held
    se <- if (!is.null(V)) sqrt(V[p$id, p$id]) else NA_real_
    dsd <- if (!is.null(dr)) stats::sd(dr$draws[p$row, ]) else NA_real_
    r <- data.frame(param = nm, held = is_held, est = f$opt$par[[p$id]], truth = p$truth,
                    s = se, d = dsd, label = NA_character_, prof_up = NA_real_,
                    wald_up = NA_real_, cap = cap, base_lo = NA_real_, base_hi = NA_real_,
                    r2_lo = NA_real_, r2_hi = NA_real_, low_plateau = NA,
                    low_dev = NA_real_, edge_mass = NA_real_)
    if (!is_held) {
      tp <- proc.time()[["elapsed"]]
      pr <- profile_one(f, p$id, se); prs[[nm]] <- pr
      t_prof <- t_prof + proc.time()[["elapsed"]] - tp
      r$label <- pr$label; r$prof_up <- pr$prof_up; r$wald_up <- pr$wald_up
      if (!is.null(base)) { qb <- qint(base$draws[p$row, ]); r$base_lo <- qb[1]; r$base_hi <- qb[2] }
      d2 <- sample_profile(pr, NDRAW, cap)
      q2 <- qint(d2); r$r2_lo <- q2[1]; r$r2_hi <- q2[2]
      le <- low_edge(pr, cap, d2)
      r$low_plateau <- as.logical(le[["low_plateau"]]); r$low_dev <- le[["low_dev"]]
      r$edge_mass <- le[["edge_mass"]]
    }
    prow[[nm]] <- r
  }
  par <- do.call(rbind, prow)
  if (is.null(par)) par <- data.frame(param = NA_character_)
  ## ridge pairs, where both one-parameter profiles failed
  par$ridge <- NA
  fails <- c("unbounded", "wald_too_wide")
  ridges <- list(); t_ridge <- 0
  for (pp in list(c("sigma", "ar_sd"), c("ri_sd", "ar_sd")))
    if (all(pp %in% names(prs)) && all(par$label[match(pp, par$param)] %in% fails)) {
      tr <- proc.time()[["elapsed"]]
      rp <- ridge_pair(f, vp[[pp[1]]]$id, vp[[pp[2]]]$id)
      t_ridge <- t_ridge + proc.time()[["elapsed"]] - tr
      par$ridge[match(pp, par$param)] <- rp$ridge
      ridges[[paste(pp, collapse = "+")]] <- list(pair = pp, rp = rp)
    }
  ## forecasts: today's draws, remedy 1 (b) given theta, and remedy 2 at each
  ## candidate cut-off
  fc <- NULL; t_r2 <- NA_real_
  if (s$ar && !is.null(base)) {
    lr <- last_rows(f, base$map)
    ref <- truth_ref(ce, s$d, s$fut)
    fc <- summ_fc(forecast_score(base$draws, base$map, ce, s$fut, ref, lr), "joint")
    gt <- tryCatch(suppressWarnings(ilm_draws(f, nsim = NDRAW, seed = 2L, given = "theta",
                                              natural = FALSE)), error = function(e) NULL)
    if (!is.null(gt))
      fc <- rbind(fc, summ_fc(forecast_score(gt$draws, gt$map, ce, s$fut, ref, lr),
                              "given_theta"))
    t2 <- proc.time()[["elapsed"]]
    for (k in seq_along(CUTS)) {
      cu <- CUTS[k]
      fl <- par$param[!par$held & is.finite(par$s) & par$s > cu & par$param %in% names(prs)]
      if (!length(fl)) next
      ## remedy 2, and the same with the lowest grid step left out (the
      ## lower-plateau check), each from the same seed
      for (tr in c(0, STEP)) {
        set.seed(seed + 7919L * k)
        D <- base$draws; done <- character(0)
        for (rn in names(ridges)) {
          pp <- ridges[[rn]]$pair
          if (isTRUE(ridges[[rn]]$rp$ridge) && all(pp %in% fl)) {
            nw <- sample_ridge(ridges[[rn]]$rp, NDRAW, cap)
            if (!is.null(nw)) {
              D <- condition_on(D, c(vp[[pp[1]]]$row, vp[[pp[2]]]$row), nw)
              done <- c(done, pp)
            }
          }
        }
        for (nm in setdiff(fl, done))
          D <- condition_on(D, vp[[nm]]$row,
                            rbind(sample_profile(prs[[nm]], NDRAW, cap, trim = tr)))
        fc <- rbind(fc, summ_fc(forecast_score(D, base$map, ce, s$fut, ref, lr),
                                if (tr == 0) "profile" else "profile_trim", cu))
      }
    }
    t_r2 <- proc.time()[["elapsed"]] - t2
  }
  list(par = cbind(key, ok = TRUE, nonconv = nonconv, par),
       fc = if (is.null(fc)) NULL else cbind(key, fc),
       time = cbind(key, t_fit = t_fit, t_draws = t_draws, t_prof = t_prof,
                    t_ridge = t_ridge, t_r2 = t_r2))
}

## ---- the named fixtures, rebuilt by their reporters' own generators ---------------------
## the forecasting package's simulate_panel() (40 series of T + 7), copied,
## keeping each row's latent and the rest of its linear predictor
fx_panel <- function(Tn, family, G = 40L) {
  n_t <- Tn + 3L + 4L
  rho <- rep(0.7, G)                                  # its "shared" truth
  sd_lat <- if (family == "gaussian") rep(0.8, G) else rep(0.5, G)
  sd_noise <- rep(0.5, G)
  mu <- if (family == "gaussian") stats::rnorm(G, 5, 1) else stats::rnorm(G, 1.5, 0.5)
  do.call(rbind, lapply(seq_len(G), function(k) {
    a <- as.numeric(stats::arima.sim(list(ar = rho[k]), n_t, sd = sd_lat[k] * sqrt(1 - rho[k]^2)))
    y <- if (family == "gaussian") mu[k] + a + stats::rnorm(n_t, 0, sd_noise[k])
         else stats::rpois(n_t, exp(mu[k] + a))
    data.frame(g = sprintf("s%02d", k), t = seq_len(n_t), a = a, m = mu[k], y = y)
  }))
}
fx_single <- function(name, cell, rep, Tn, family, series, lat_sd) {
  set.seed(5e6 + 1e4 * cell + rep)
  d <- fx_panel(Tn, family); d <- d[d$g == series, ]; d$g <- factor(d$g)
  list(name = name,
       ce = data.frame(arm = "B1", family = family, Tn = Tn, rho = 0.7, lat_sd = lat_sd,
                       ri_sd = NA, L = NA, m = NA,
                       noise = if (family == "gaussian") 0.5 else NA, cell = 0L, rep = 0L),
       sim = list(d = d[d$t <= Tn, ], fut = d[d$t > Tn & d$t <= Tn + H, ], fml = y ~ 1, ar = TRUE))
}
## the AR-only panel of its covariate reports, noise-free: set.seed(rep), 20
## series of 24 with AR 0.8 and marginal SD 1. Its generator stops at 24, so
## the six future times are drawn after it, from the true AR on from each
## series' last value, with the random stream carrying on
fx_covpanel <- function(rep) {
  set.seed(rep)
  d <- do.call(rbind, lapply(1:20, function(k)
    data.frame(g = k, t = 1:24,
               y = as.numeric(stats::arima.sim(list(ar = 0.8), 24, sd = sqrt(1 - 0.64))))))
  fut <- do.call(rbind, lapply(1:20, function(k) {
    a <- d$y[d$g == k & d$t == 24]; y <- numeric(H)
    for (h in seq_len(H)) { a <- 0.8 * a + stats::rnorm(1, 0, sqrt(1 - 0.64)); y[h] <- a }
    data.frame(g = k, t = 24 + seq_len(H), y = y)
  }))
  d$a <- d$y; fut$a <- fut$y; d$m <- 0; fut$m <- 0
  d$g <- factor(d$g); fut$g <- factor(fut$g, levels = levels(d$g))
  list(name = paste0("covpanel_rep", rep),
       ce = data.frame(arm = "B2", family = "gaussian", Tn = 24, rho = 0.8, lat_sd = 1,
                       ri_sd = NA, L = NA, m = NA, noise = 0, cell = 0L, rep = rep),
       sim = list(d = d, fut = fut, fml = y ~ 1, ar = TRUE))
}
fixtures <- function() c(
  list(fx_single("single_s36", 2, 11, 24L, "gaussian", "s36", 0.8)),
  lapply(c(4L, 11L, 59L), fx_covpanel),
  list(fx_single("poisson_s14", 7, 1, 12L, "poisson", "s14", 0.5)))

## ---- run ------------------------------------------------------------------------------------
t_all <- Sys.time()
use <- if (is.null(ONLY)) cells else cells[cells$cell %in% ONLY, ]
jobs <- merge(use, data.frame(rep = seq_len(NREP)))
jl <- split(jobs, seq_len(nrow(jobs)))
tag <- if (OFFSET) "fresh" else if (!is.null(ONLY)) "smoke" else "main"
## CHECKPOINTS (an I/O change, no effect on results): each fit's result is
## saved as it finishes, and a run started again reads the fits already done
## instead of refitting them. Every fit sets its own seed from its cell and
## replicate, so a resumed run gives the same numbers as one uninterrupted.
ck <- file.path(sp, paste0("checkpoints_", tag))
dir.create(ck, showWarnings = FALSE, recursive = TRUE)
one_ck <- function(job) {
  fp <- file.path(ck, sprintf("cell%02d_rep%03d.rds", job$cell, job$rep))
  if (file.exists(fp)) return(readRDS(fp))
  r <- one(job)
  saveRDS(r, paste0(fp, ".part")); file.rename(paste0(fp, ".part"), fp)
  r
}
res <- if (NCORE > 1L) {
  cl <- parallel::makeCluster(NCORE)
  parallel::clusterExport(cl, setdiff(ls(globalenv()), c("cl", "jobs", "jl")), envir = globalenv())
  ## one job at a time: by default the jobs go out in one contiguous chunk
  ## per worker, and on a resume the unfinished ones -- the last replicates --
  ## all fall in the last worker's chunk
  r <- parallel::parLapplyLB(cl, jl, one_ck, chunk.size = 1L)
  parallel::stopCluster(cl)
  r
} else lapply(jl, one_ck)
## rows of failed fits carry fewer columns; each is filled out with NA
bind <- function(x, k) {
  l <- Filter(Negate(is.null), lapply(x, `[[`, k))
  if (!length(l)) return(NULL)
  nm <- unique(unlist(lapply(l, names)))
  do.call(rbind, lapply(l, function(d) { for (n in setdiff(nm, names(d))) d[[n]] <- NA; d[nm] }))
}
put <- function(d, nm) if (!is.null(d))
  utils::write.csv(d, file.path(sp, paste0("variance_draws_", nm, ".csv")), row.names = FALSE)
for (k in c("par", "fc", "time")) put(bind(res, k), paste(tag, k, sep = "_"))
if (!OFFSET) {
  fx <- lapply(fixtures(), function(x) one(x$ce, fixture = x))
  for (k in c("par", "fc", "time")) put(bind(fx, k), paste("fixtures", k, sep = "_"))
}
cat("fits:", nrow(jobs), " minutes:",
    round(as.numeric(difftime(Sys.time(), t_all, units = "mins")), 1), "\n")
