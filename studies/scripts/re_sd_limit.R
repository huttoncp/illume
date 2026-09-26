## A random effect's SD at zero that the rule does not see.
##
## PRE-REGISTERED DESIGN, committed before any run.
##
## A random intercept beside an AR(1): the correlation over time takes up each
## series' level, the intercept's variance has its optimum at zero, and the
## optimiser stops short of it on a flat likelihood -- at an SD of 0.0204 in
## the reported case, above the 1e-3 line at which a variance is taken as
## zero. Nothing flags it, TMB calls its Hessian positive definite, the SE of
## log SD(g) is 535, and draws of it run to +/-1800. The same early stop the
## dispersion hold dealt with (dispersion_limit.R), on a covariance.
##
## The rule to be placed, as Craig ruled: a random effect's SD is at its
## boundary when it is below a LINE and the likelihood is FLAT below it --
## the objective, the random effects re-optimised, moves by at most 5e-3 when
## the log SD is pushed 3 lower. This study places the line.
##
## Labels, taken from OUTSIDE the rule, as the dispersion study did. A fit's
## term is AT ITS BOUNDARY when both hold:
##   drop_ok  the model without the term fits as well: its log-likelihood
##            is within 1e-3 of the full model's
##   push_ok  the objective rises by at most 1e-3 when the log SD is pushed
##            3 lower
##
## Candidate lines on the SD, on the linear predictor's scale (the scale of
## the existing 1e-3 variance line): 1e-3, 1e-2, 0.05, 0.1, 0.2, and the
## flatness test alone. For each: how many fits it holds, how many of those
## are NOT at the boundary by the labels (false holds, to be zero), and how
## many at the boundary it misses. The line adopted is the widest with no
## false hold, reported with the gap between the SDs at and not at the
## boundary, as for the dispersion.
##
## Design:
##   families  gaussian (noise SD 0.5); Poisson, mean about 3; binary, p
##             about 0.5
##   ri_ar     20 series x T times, T in {12, 24}; AR(1) rho 0.8, SD 0.5;
##             a random intercept per series with SD in {0, 0.1, 0.3, 0.5}
##   ri        L groups of 10 rows, L in {5, 10, 30}; random intercept SD in
##             {0, 0.1, 0.3}; no correlation over time
##   100 replicates per cell: 24 + 27 = 51 cells, 5,100 fits, each refitted
##   without the term for the label
## A slope x (coefficient 0.3) in every model, whose standard error the
## verify phase compares with the refit without the term.
##
## With `verify`, phase 2 on the built rule: the fits held have the refit's
## standard error for x, and draws of the held SD that stay finite.
##
## Usage: Rscript re_sd_limit.R <nrep> <ncore> <outdir> [verify]
## ---------------------------------------------------------------------------
args  <- commandArgs(trailingOnly = TRUE)
NREP  <- if (length(args) >= 1) as.integer(args[1]) else 100L
NCORE <- if (length(args) >= 2) as.integer(args[2]) else 2L
sp    <- if (length(args) >= 3) args[3] else "."
VERIFY <- length(args) >= 4 && identical(args[4], "verify")
## a seed offset, so the confirmation runs on data the line was not chosen on
OFFSET <- if (length(args) >= 5) as.integer(args[5]) else 0L
dir.create(sp, showWarnings = FALSE, recursive = TRUE)

fams <- data.frame(family = c("gaussian", "poisson", "binomial"),
                   b0 = c(0.5, log(3), 0), stringsAsFactors = FALSE)
cells <- rbind(
  merge(fams, merge(data.frame(design = "ri_ar", size = c(12, 24)),
                    data.frame(sd = c(0, 0.1, 0.3, 0.5)))),
  merge(fams, merge(data.frame(design = "ri", size = c(5, 10, 30)),
                    data.frame(sd = c(0, 0.1, 0.3)))))
cells$cell <- seq_len(nrow(cells))
jobs <- merge(cells, data.frame(rep = seq_len(NREP)))

one <- function(job) {
  suppressMessages(library(illume))
  seed <- 8191L * job$cell + job$rep + OFFSET
  set.seed(seed)
  if (job$design == "ri_ar") {
    G <- 20L; Tn <- job$size
    d <- expand.grid(t = seq_len(Tn), g = factor(seq_len(G)))
    ar <- unlist(lapply(seq_len(G), function(i) as.numeric(
      stats::arima.sim(list(ar = 0.8), Tn, sd = 0.5 * sqrt(1 - 0.64)))))
    lat <- ar + stats::rnorm(G, 0, job$sd)[d$g]
    arspec <- ilm_ar1(~ t | g)
  } else {
    d <- data.frame(g = factor(rep(seq_len(job$size), each = 10L)))
    lat <- stats::rnorm(job$size, 0, job$sd)[d$g]
    arspec <- NULL
  }
  d$x <- stats::rnorm(nrow(d))
  eta <- job$b0 + 0.3 * d$x + lat
  d$y <- switch(job$family,
    gaussian = eta + stats::rnorm(nrow(d), 0, 0.5),
    poisson  = stats::rpois(nrow(d), exp(eta)),
    binomial = stats::rbinom(nrow(d), 1, stats::plogis(eta)))
  fitit <- function(fml) tryCatch(suppressMessages(suppressWarnings(
    ilm_model(fml, data = d, family = job$family, ar = arspec,
              verbose = FALSE))), error = function(e) NULL)
  f <- fitit(y ~ x + (1 | g))
  out <- data.frame(job, seed = seed, ok = !is.null(f), sd_hat = NA_real_,
                    se_logsd = NA_real_, push = NA_real_, ll = NA_real_,
                    ll_drop = NA_real_, how = NA_character_,
                    held = NA_character_, boundary = NA_character_,
                    se_x = NA_real_, se_x_drop = NA_real_, draws_finite = NA,
                    stringsAsFactors = FALSE)
  if (is.null(f)) return(out)
  pn <- names(f$opt$par); it <- which(pn == "theta")
  out$sd_hat <- sqrt(ilm_varcorr(f)$re$g[1, 1])
  V <- tryCatch(suppressWarnings(vcov(f, full = TRUE)), error = function(e) NULL)
  if (!is.null(V)) out$se_logsd <- sqrt(V[it, it])
  out$ll <- as.numeric(logLik(f))
  out$how <- f$hessian_how
  out$held <- paste(f$hessian_held, collapse = ",")
  out$boundary <- paste(f$boundary_terms, collapse = ",")
  p2 <- f$opt$par; p2[it] <- p2[it] - 3
  out$push <- tryCatch(f$obj$fn(p2) - f$obj$fn(f$opt$par),
                       error = function(e) NA_real_)
  out$se_x <- tryCatch(sqrt(diag(vcov(f)))[["x"]], error = function(e) NA_real_)
  f0 <- fitit(y ~ x)
  if (!is.null(f0)) {
    out$ll_drop <- as.numeric(logLik(f0))
    out$se_x_drop <- tryCatch(sqrt(diag(vcov(f0)))[["x"]],
                              error = function(e) NA_real_)
  }
  if (VERIFY) {
    dr <- tryCatch(ilm_draws(f, nsim = 200, seed = 1, natural = FALSE),
                   error = function(e) NULL)
    if (!is.null(dr))
      out$draws_finite <- all(abs(dr$draws[rownames(dr$draws) == "theta", ]) < 50)
  }
  out
}

t0 <- Sys.time()
cl <- parallel::makeCluster(NCORE)
parallel::clusterExport(cl, c("VERIFY", "OFFSET"))
res <- do.call(rbind, parallel::parLapplyLB(cl, split(jobs, seq_len(nrow(jobs))), one))
parallel::stopCluster(cl)
res$drop_ok <- res$ll_drop >= res$ll - 1e-3
res$push_ok <- res$push <= 1e-3
res$at <- res$drop_ok & res$push_ok
tag <- if (VERIFY) "verify" else "phase1"
utils::write.csv(res, file.path(sp, paste0("re_sd_limit_", tag, ".csv")), row.names = FALSE)
cat("fits:", nrow(res), " failed:", sum(!res$ok), " minutes:",
    round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "\n")

r <- res[res$ok & !is.na(res$at), ]
cat("at the boundary by both labels:", sum(r$at), "of", nrow(r),
    "| labels disagree:", sum(r$drop_ok != r$push_ok), "\n")
q <- function(v) if (length(v)) signif(stats::quantile(v, c(0, .5, .99, 1), na.rm = TRUE), 3) else NA
cat("SD at the boundary (min, median, 99%, max):", q(r$sd_hat[r$at]), "\n")
cat("SD not at the boundary (min, median, 99%, max):", q(r$sd_hat[!r$at]), "\n")
## the pre-registered tolerance, 5e-3, and the tighter 1e-3 the first run
## pointed to
for (TOL in c(5e-3, 1e-3)) for (line in c(1e-3, 1e-2, 0.05, 0.1, 0.2, Inf)) {
  rl <- r$sd_hat < line & r$push <= TOL
  cat(sprintf("tol %.0e line %-6s + flat: holds %4d | at %4d | FALSE %3d | missed %3d\n",
              TOL, format(line), sum(rl), sum(rl & r$at), sum(rl & !r$at), sum(!rl & r$at)))
}
print(stats::aggregate(cbind(fits = 1, at = at, flagged_now = boundary != "") ~
                         family + design + size + sd, data = r, FUN = sum))
if (VERIFY) {
  h <- grepl("(^|,)g(,|$)", r$held) | grepl("(^|,)g(,|$)", r$boundary)
  cat("held:", sum(h), "| of them at the boundary:", sum(h & r$at),
      "| SE(x) ratio to the refit without the term (min, median, max):",
      q((r$se_x / r$se_x_drop)[h]), "| finite draws:",
      sum(r$draws_finite[h], na.rm = TRUE), "of", sum(h), "\n")
}
