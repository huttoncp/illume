## Calibration of the flatness rule (Craig's item 269).
##
## PRE-REGISTERED DESIGN, committed before any registered run. A pilot (2
## replicates of every cell, and the four CI cases, on illume main a387c93)
## ran first; it checked the rebuild and set the estimates below.
##
## ## The rule, as Craig ruled
## A variance or dispersion past its line is held when its log-scale
##   parameter, pushed 1.5, 3 and 6 further towards its limit (the other outer
##   parameters fixed, the random effects re-optimised), moves the objective
##   by at most 0.05 at every step. A push that LOWERS the objective by more
##   than 0.05 marks the fit unconverged, not held. summary() always shows a
##   held limit. Today's rule: one push of 3, held at a rise of at most 1e-3.
##
## ## What is measured
## The saved studies' phase-1 fits, rebuilt from their recorded seeds with the
##   studies' own generators (dispersion_limit.R: 600 fits, negative binomial
##   and beta; dispersion_limit_gaussian.R: 500, sigma; re_sd_limit.R: 5,100,
##   a random intercept's SD), plus the four fits behind main's red CI (item
##   243: a negative binomial k at its limit, a flat beta, a "curved" beta, a
##   gaussian sigma stalled near zero). Per fit: the three pushes, the
##   objective re-evaluated at the optimum after them (its noise), whether the
##   term is past its line on main (the candidate), the build's own hold, and
##   a label from OUTSIDE the rule where the studies had one -- the Poisson
##   refit for a negative binomial, and the model without the term for a
##   random effect -- as the log-likelihood difference.
## The pilot reproduced the saved studies: push at 3 against each study's
##   saved push, at most 6e-5 apart (dispersion), 0 (sigma), 6e-14 (SD).
##
## ## Two float paths
## Every fit on two builds of the same illume (main a387c93) with illumex
##   ff5e579: P1 R 4.4.3 with RTMB 1.9, and P2 R 4.6.1 with RTMB 2.0 (CI's
##   RTMB). Two paths on one machine are a proxy for CI's platforms, not the
##   same thing; the findings say so.
##
## ## Decisions, per candidate fit
## - new: held if all three pushes are in [-0.05, 0.05]; unconverged if any
##   is below -0.05; else not held.
## - old: held if the push at 3 is at most 1e-3.
## - label (labelled studies): the limit (or dropped) model fits as well, its
##   log-likelihood within 0.05 of the fit's (and, as the studies had it,
##   within 1e-3).
##
## ## What counts, fixed now
## V1 robustness: among candidates, the new decision is the same on P1 and P2
##    for at least 99.5% of fits, and flips no more often than the old one.
## V2 no false holds: in the labelled studies, the new rule holds no fit
##    whose label says the limit or dropped model is worse by more than 0.05
##    (the studies required no false hold of the old rule too); misses, fits
##    the label puts at the limit that the rule does not hold, are reported.
## V2b missed holds (addendum A1, committed before the run; reported, not a
##    gate): per study and path, the candidates whose label puts the limit or
##    dropped model within 0.05 that the new rule does not hold.
## V3 margin: the share of candidates whose largest |push| is within a factor
##    of 2 of 0.05 ([0.025, 0.1]), the band where a decision could flip, per
##    study and path.
## V4 unconverged: how many fits a push below -0.05 marks, each listed with
##    its pushes and the build's own convergence code.
## V5 the CI cases: the new decision for each of the four on both paths. The
##    pilot (P1) gives: k at its limit, pushes ~4e-7, held; the flat beta,
##    largest 0.011, held; the "curved" beta, push at 6 of -261, unconverged;
##    the stalled sigma (log sigma -9.3 of sd(y)), push at 6 of +0.070, NOT
##    held -- pushed to about e^-15 the gaussian objective may be arithmetic
##    rather than likelihood, so V5 reports it, and the findings say whether
##    the rule needs a floor on how far a push may go (a decision for Craig).
## If V1 or V2 fails, the numbers go to Craig before anything is built.
##
## ## Size and stops
## Piloted seconds per fit: dispersion 1.5, sigma 0.5, SD 0.4 (with the
##   label refits): about 51 minutes per path, 1.7 hours for both, one
##   process at a time on one core. A path stops at twice its estimate (1.7
##   hours) and reports the fits it finished.
##
## Usage: Rscript flatness_calibration.R <study> <rep_from> <rep_to> <out.csv>
##   study: disp | gauss | re (the saved studies' phase 1, rebuilt by seed)
##          | cases (the four CI fits)

a <- commandArgs(TRUE)
STUDY <- a[1]; R1 <- as.integer(a[2]); R2 <- as.integer(a[3]); OUT <- a[4]
suppressPackageStartupMessages(library(illume))
q <- function(e) tryCatch(suppressMessages(suppressWarnings(e)), error = function(e) NULL)
fn_at <- function(f, p) tryCatch(f$obj$fn(p), error = function(e) NA_real_)

## ---- the three studies' generators, exactly as their scripts ---------------
## dispersion_limit.R (seed 7919 cell + rep): negative binomial and beta
disp_cells <- rbind(
  data.frame(family = "nbinom", design = "ar1",  disp = c(Inf, 50, 5, 1)),
  data.frame(family = "nbinom", design = "ri",   disp = c(Inf, 50, 5)),
  data.frame(family = "beta",   design = "ar1",  disp = c(1e4, 200, 20)),
  data.frame(family = "beta",   design = "ri",   disp = c(200, 20)))
disp_cells$cell <- seq_len(nrow(disp_cells))
gen_disp <- function(job) {
  seed <- 7919L * job$cell + job$rep; set.seed(seed)
  if (job$design == "ar1") {
    G <- 20; Tn <- 16
    d <- expand.grid(t = seq_len(Tn), g = factor(seq_len(G)))
    lat <- unlist(lapply(seq_len(G), function(i) as.numeric(
      stats::arima.sim(list(ar = 0.7), Tn, sd = 0.8 * sqrt(1 - 0.49)))))
  } else {
    G <- 30; Tn <- 10
    d <- expand.grid(t = seq_len(Tn), g = factor(seq_len(G)))
    lat <- stats::rnorm(G, 0, 0.6)[d$g]
  }
  d$x <- stats::rnorm(nrow(d)); eta <- 0.5 + 0.3 * d$x + lat
  if (job$family == "nbinom") {
    mu <- exp(eta)
    d$y <- if (is.infinite(job$disp)) stats::rpois(nrow(d), mu)
           else stats::rnbinom(nrow(d), size = job$disp, mu = mu)
  } else {
    mu <- stats::plogis(eta - 0.5)
    d$y <- pmin(pmax(stats::rbeta(nrow(d), mu * job$disp, (1 - mu) * job$disp), 1e-6), 1 - 1e-6)
  }
  list(d = d, seed = seed, family = job$family, fml = y ~ x + (1 | g),
       ar = if (job$design == "ar1") ilm_ar1(~ t | g), par = "logdisp", side = 1,
       limit_family = if (job$family == "nbinom") "poisson")
}
## dispersion_limit_gaussian.R (seed 6007 cell + rep): sigma
gauss_cells <- rbind(data.frame(design = "ar1", noise = c(0.5, 0.2, 0.05)),
                     data.frame(design = "ri1", noise = 0.5),
                     data.frame(design = "ri",  noise = 0.5))
gauss_cells$cell <- seq_len(nrow(gauss_cells))
gen_gauss <- function(job) {
  seed <- 6007L * job$cell + job$rep; set.seed(seed)
  if (job$design == "ar1") {
    G <- 20; Tn <- 16
    d <- expand.grid(t = seq_len(Tn), g = factor(seq_len(G)))
    lat <- unlist(lapply(seq_len(G), function(i) as.numeric(
      stats::arima.sim(list(ar = 0.7), Tn, sd = 0.8 * sqrt(1 - 0.49)))))
  } else if (job$design == "ri1") {
    d <- data.frame(g = factor(seq_len(300))); lat <- stats::rnorm(300, 0, 0.6)
  } else {
    d <- data.frame(g = factor(rep(seq_len(30), each = 10))); lat <- stats::rnorm(30, 0, 0.6)[d$g]
  }
  d$x <- stats::rnorm(nrow(d))
  d$y <- 0.5 + 0.3 * d$x + lat + stats::rnorm(nrow(d), 0, job$noise)
  list(d = d, seed = seed, family = "gaussian", fml = y ~ x + (1 | g),
       ar = if (job$design == "ar1") ilm_ar1(~ t | g), par = "logdisp", side = -1,
       limit_family = NULL)
}
## re_sd_limit.R (seed 8191 cell + rep): a random intercept's SD
fams <- data.frame(family = c("gaussian", "poisson", "binomial"), b0 = c(0.5, log(3), 0),
                   stringsAsFactors = FALSE)
re_cells <- rbind(
  merge(fams, merge(data.frame(design = "ri_ar", size = c(12, 24)), data.frame(sd = c(0, 0.1, 0.3, 0.5)))),
  merge(fams, merge(data.frame(design = "ri", size = c(5, 10, 30)), data.frame(sd = c(0, 0.1, 0.3)))))
re_cells$cell <- seq_len(nrow(re_cells))
gen_re <- function(job) {
  seed <- 8191L * job$cell + job$rep; set.seed(seed)
  if (job$design == "ri_ar") {
    G <- 20L; Tn <- job$size
    d <- expand.grid(t = seq_len(Tn), g = factor(seq_len(G)))
    ar <- unlist(lapply(seq_len(G), function(i) as.numeric(
      stats::arima.sim(list(ar = 0.8), Tn, sd = 0.5 * sqrt(1 - 0.64)))))
    lat <- ar + stats::rnorm(G, 0, job$sd)[d$g]; arspec <- ilm_ar1(~ t | g)
  } else {
    d <- data.frame(g = factor(rep(seq_len(job$size), each = 10L)))
    lat <- stats::rnorm(job$size, 0, job$sd)[d$g]; arspec <- NULL
  }
  d$x <- stats::rnorm(nrow(d)); eta <- job$b0 + 0.3 * d$x + lat
  d$y <- switch(job$family, gaussian = eta + stats::rnorm(nrow(d), 0, 0.5),
                poisson = stats::rpois(nrow(d), exp(eta)),
                binomial = stats::rbinom(nrow(d), 1, stats::plogis(eta)))
  list(d = d, seed = seed, family = job$family, fml = y ~ x + (1 | g), ar = arspec,
       par = "theta", side = -1, drop_fml = y ~ x)
}

## the four fits behind main's red CI (item 243's probe, illume-work
## notes/ci_env_probe.R): a negative binomial k at its limit, a flat and a
## "curved" beta, and a gaussian residual SD stalled near zero
case_cells <- data.frame(case = c("nbinom_k_limit", "beta_flat_71272", "beta_curved_63362",
                                  "gaussian_stalled"), cell = 1:4)
gen_case <- function(job) {
  dl <- function(seed, phi = NULL) {
    set.seed(seed); G <- 20; Tn <- 16
    d <- expand.grid(t = seq_len(Tn), g = factor(seq_len(G)))
    if (is.null(phi)) {
      d$x <- stats::rnorm(nrow(d))
      lat <- unlist(lapply(seq_len(G), function(i) as.numeric(
        stats::arima.sim(list(ar = 0.7), Tn, sd = 0.8 * sqrt(1 - 0.49)))))
      d$y <- stats::rnbinom(nrow(d), size = 5, mu = exp(0.5 + 0.3 * d$x + lat))
    } else {
      lat <- unlist(lapply(seq_len(G), function(i) as.numeric(
        stats::arima.sim(list(ar = 0.7), Tn, sd = 0.8 * sqrt(1 - 0.49)))))
      d$x <- stats::rnorm(nrow(d)); mu <- stats::plogis(0.5 + 0.3 * d$x + lat - 0.5)
      d$y <- pmin(pmax(stats::rbeta(nrow(d), mu * phi, (1 - mu) * phi), 1e-6), 1 - 1e-6)
    }
    d
  }
  if (job$case == "gaussian_stalled") {
    set.seed(1293099); G <- 20L
    d <- expand.grid(t = 1:12, g = factor(seq_len(G)))
    lat <- unlist(lapply(seq_len(G), function(i)
      as.numeric(stats::arima.sim(list(ar = 0.7), 12, sd = 0.8 * sqrt(1 - 0.49)))))
    d$x <- stats::rnorm(nrow(d))
    d$y <- 1 + 0.3 * d$x + stats::rnorm(G, 0, 0.5)[d$g] + lat + stats::rnorm(nrow(d), 0, 0)
    return(list(d = d, seed = 1293099L, family = "gaussian", fml = y ~ x + (1 | g),
                ar = ilm_ar1(~ t | g), par = "logdisp", side = -1))
  }
  d <- switch(job$case, nbinom_k_limit = dl(5), beta_flat_71272 = dl(71272, 200),
              beta_curved_63362 = dl(63362, 1e4))
  list(d = d, seed = NA_integer_, family = if (job$case == "nbinom_k_limit") "nbinom" else "beta",
       fml = y ~ x + (1 | g), ar = ilm_ar1(~ t | g), par = "logdisp", side = 1,
       limit_family = if (job$case == "nbinom_k_limit") "poisson")
}

cells <- switch(STUDY, disp = disp_cells, gauss = gauss_cells, re = re_cells, cases = case_cells)
gen <- switch(STUDY, disp = gen_disp, gauss = gen_gauss, re = gen_re, cases = gen_case)
jobs <- merge(cells, data.frame(rep = R1:R2))

STEPS <- c(1.5, 3, 6)
for (i in seq_len(nrow(jobs))) {
  job <- jobs[i, ]; s <- gen(job)
  t0 <- proc.time()[["elapsed"]]
  f <- q(ilm_model(s$fml, data = s$d, family = s$family, ar = s$ar, verbose = FALSE))
  row <- data.frame(study = STUDY, cell = job$cell, rep = job$rep, seed = s$seed,
                    family = s$family, ok = !is.null(f), stringsAsFactors = FALSE)
  if (!is.null(f)) {
    pn <- names(f$opt$par); id <- which(pn == s$par)[1]
    p0 <- f$opt$par; f0 <- fn_at(f, p0)
    pushes <- vapply(STEPS, function(h) { p <- p0; p[id] <- p[id] + s$side * h; fn_at(f, p) - f0 }, 0)
    rep_diff <- fn_at(f, p0) - f0                       # the same point, after the pushes
    sd_y <- stats::sd(s$d$y)
    val <- unname(p0[id])
    ## the candidate line on main (1beb589 .. a387c93): dispersion 1/sqrt below
    ## 1e-2; sigma below 0.2 sd(y); a random effect's SD below 0.1 (of sd(y)
    ## for a gaussian response)
    cand <- switch(if (STUDY == "cases") (if (s$side > 0) "disp" else "gauss") else STUDY,
      disp = exp(-val / 2) < 1e-2,
      gauss = exp(val) / sd_y < 0.2,
      re = { sd <- sqrt(ilm_varcorr(f)$re$g[1, 1]); if (s$family == "gaussian") sd / sd_y < 0.1 else sd < 0.1 })
    ## labels from outside the rule: the limit model, or the model without the
    ## term, fits as well (its log-likelihood within 1e-3, the studies' own
    ## label, and within 0.05, the new rule's scale)
    ll <- as.numeric(logLik(f)); ll_alt <- NA_real_
    alt <- if (!is.null(s$limit_family)) q(ilm_model(s$fml, data = s$d, family = s$limit_family, ar = s$ar, verbose = FALSE))
           else if (!is.null(s$drop_fml)) q(ilm_model(s$drop_fml, data = s$d, family = s$family, ar = s$ar, verbose = FALSE))
    if (!is.null(alt)) ll_alt <- as.numeric(logLik(alt))
    row <- cbind(row, data.frame(value = val, candidate = cand,
      push_1.5 = pushes[1], push_3 = pushes[2], push_6 = pushes[3], repeat_diff = rep_diff,
      ll = ll, ll_alt = ll_alt, held_now = paste(f$hessian_held, collapse = "+"),
      how_now = f$hessian_how, conv = f$opt$convergence, secs = round(proc.time()[["elapsed"]] - t0, 2)))
  }
  utils::write.table(row, OUT, sep = ",", row.names = FALSE, col.names = !file.exists(OUT),
                     append = TRUE)
}
