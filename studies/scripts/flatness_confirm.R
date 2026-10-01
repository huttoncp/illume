## Confirmation of the flatness rule with a floor (Craig's item 288, the
## option "the floor, confirmed first").
##
## REGISTRATION. Drafted before Craig's ruling; he chose this option on
## 2026-10-01, and it is committed before any run, as drafted, with two
## additions: the paragraph on why the fair label is fair, and the optimizer
## and gradient checks' statuses in each row, which V6 needs.
##
## ## Why a confirmation
## The calibration (flatness_calibration.R, findings in
##   studies/findings/flatness_calibration.md) failed V1 and V2 as registered,
##   and a floor chosen AFTER that run -- a term already more than e^6 below its
##   line held by value, the rest judged by the pushes -- removed nearly all of
##   V1's flips and V4's flags. Because e^6 was chosen on that run's data, it
##   is confirmed here on fresh seeds, with the floor fixed in advance, and with
##   a label that cannot be fooled the way the calibration's was.
##
## ## The rule under test, fixed now
## A term is a candidate when it is past its line on main (dispersion
##   1/sqrt(phi) < 1e-2; gaussian sigma < 0.2 of sd(y); a random effect's SD
##   < 0.1, of sd(y) for a gaussian response). Its distance below the line is
##   measured on the log scale the pushes use (log of 1/sqrt(phi), sigma / sd(y)
##   or SD, minus log of the line).
## - Below the floor -- more than 6 below the line on that scale, e^-6 -- the
##   term is HELD BY VALUE, with no push.
## - Otherwise, the pushes of 1.5, 3 and 6 towards the limit (the other outer
##   parameters fixed, the random effects re-optimised): held if all three
##   move the objective by at most 0.05; unconverged if any lowers it by more
##   than 0.05; else not held.
## - Beside it, as in the calibration: the rule without the floor, and today's
##   rule (one push of 3, held at a rise of at most 1e-3).
##
## ## Data: fresh seeds
## The same three studies' designs and generators as the calibration, on
##   replicates the calibration never used: dispersion_limit.R reps 1001-1050
##   (600 fits), dispersion_limit_gaussian.R reps 1001-1100 (500),
##   re_sd_limit.R reps 1001-1100 (5,100), plus the four CI fits (fixed data,
##   reported again under the floor). Seeds as each study's generator makes
##   them from (cell, rep).
##
## ## The label, made fair
## The calibration's random-effect label (the model refitted without the term)
##   was fooled where the refit stopped in a worse local basin of an AR
##   correlation, rho of the opposite sign to the full fit's. Here the label
##   for a random effect, and for a negative binomial's k, is the BETTER of:
##   (a) the ordinary refit, without the term (or as a Poisson); and
##   (b) the profile at the limit: the full model's own objective with the
##       term's log parameter moved 20 beyond its estimate towards the limit,
##       minimised over every other outer parameter STARTING FROM THE FULL
##       FIT'S VALUES (so from its rho).
##   The label's cost is the better one's rise in the objective over the full
##   fit's; the limit or dropped model is "as good" when that cost is at most
##   0.05. Beta precisions and gaussian sigmas have no label: at their limits
##   the objective's arithmetic is the very thing under question.
## Why taking the lower cost is fair, not tilted towards holds: the label asks
##   how well the best model at the limit (or without the term) fits. Both (a)
##   and (b) are fits of that model, by minimisation from different starts, so
##   each is an upper bound on its best fit's cost, and the lower of the two
##   is the closer estimate of it. A false hold is a hold where the best fit at
##   the limit is worse by more than 0.05; taking the higher cost would count
##   a refit stuck in a worse local basin as evidence against the hold, which
##   is how the calibration's four "false holds" arose.
##
## ## Paths, and what counts
## Two float paths, as in the calibration: P1 R 4.4.3 / RTMB 1.9 and P2 R 4.6.1
##   / RTMB 2.0, the same illume build (main at the time of the run, recorded)
##   and illumex ff5e579. Judged on the floor rule:
## V1 robustness: among candidates, the same decision on P1 and P2 for at
##    least 99.5% of fits, and no more flips than today's rule.
## V2 no false holds: no fit held (by value or by the pushes) whose fair label
##    puts the dropped or limit model worse by more than 0.05.
## V2b missed holds (reported): candidates the fair label puts within 0.05 that
##    the floor rule does not hold, per study and path.
## V3 the band: the share of pushed candidates whose largest |push| is in
##    [0.025, 0.1].
## V4 unconverged: each flag listed, with the build's own convergence code.
## V5 the CI cases under the floor rule, both paths.
## V6 the floor's reach: how many candidates are held by value, per study and
##    path, and their distances below the line; and every fit held by value
##    whose optimizer or gradient check FAILs (the runaway case), listed.
## If V1 or V2 fails, the numbers go back to Craig before anything is built.
##
## ## Size and stops
## The calibration took 51 and 57 minutes per path; the profile label adds one
##   optimisation per labelled fit (about 5,450 per path, a fraction of a
##   second each): about 2.5 hours for both paths on one core. Each path stops
##   at twice its estimate.
##
## Usage: Rscript flatness_confirm.R <study> <rep_from> <rep_to> <out.csv>
##   study: disp | gauss | re | cases

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
    ## (b) the profile at the limit: the full objective with the term moved 20
    ## beyond its estimate, minimised over every other outer parameter from the
    ## full fit's own values (so from its rho)
    cost_prof <- NA_real_
    if (!is.null(s$limit_family) || !is.null(s$drop_fml)) {
      oth <- setdiff(seq_along(p0), id); pl <- p0; pl[id] <- p0[id] + s$side * 20
      gfn <- function(po) { p <- pl; p[oth] <- po; fn_at(f, p) }
      ggr <- function(po) { p <- pl; p[oth] <- po; tryCatch(f$obj$gr(p)[oth], error = function(e) rep(NA_real_, length(oth))) }
      op <- tryCatch(stats::nlminb(p0[oth], gfn, ggr), error = function(e) NULL)
      if (!is.null(op) && is.finite(op$objective)) cost_prof <- op$objective - f0
      invisible(fn_at(f, p0))
    }
    cost_refit <- if (is.na(ll_alt)) NA_real_ else ll - ll_alt
    cost_fair <- suppressWarnings(min(c(cost_refit, cost_prof), na.rm = TRUE)); if (!is.finite(cost_fair)) cost_fair <- NA_real_
    ## the distance below the line, on the log scale the pushes use
    below <- switch(if (STUDY == "cases") (if (s$side > 0) "disp" else "gauss") else STUDY,
      disp = (-val / 2) - log(1e-2), gauss = log(exp(val) / sd_y / 0.2),
      re = { sd <- sqrt(ilm_varcorr(f)$re$g[1, 1]); log(max(sd, 1e-300) / (if (s$family == "gaussian") 0.1 * sd_y else 0.1)) })
    row <- cbind(row, data.frame(value = val, candidate = cand,
      push_1.5 = pushes[1], push_3 = pushes[2], push_6 = pushes[3], repeat_diff = rep_diff,
      ll = ll, ll_alt = ll_alt, cost_refit = cost_refit, cost_prof = cost_prof, cost_fair = cost_fair,
      below = below, held_now = paste(f$hessian_held, collapse = "+"),
      how_now = f$hessian_how, conv = f$opt$convergence,
      check_optimizer = f$checks$status[f$checks$check == "optimizer"][1],
      check_gradient = f$checks$status[f$checks$check == "gradient"][1], secs = round(proc.time()[["elapsed"]] - t0, 2)))
  }
  utils::write.table(row, OUT, sep = ",", row.names = FALSE, col.names = !file.exists(OUT),
                     append = TRUE)
}
