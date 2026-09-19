## ---------------------------------------------------------------------------
## Monte Carlo coverage of illume's Wald confidence intervals for fixed effects.
##
## Design note.  The cells deliberately span the LATENT BUDGET -- observations
## per latent value, which governs how accurate the Laplace approximation is.
## Two multinomial cells sit below the ~3 obs/latent working rule, where
## coverage is expected to degrade.  That is deliberate: a coverage study that
## only samples the comfortable region cannot fail, and so proves nothing.
##
## The gaussian fixed-effects cell is the opposite check.  There the intervals
## are exact t intervals, so anything other than ~95% would indicate a bug
## rather than an approximation error.
##
## Raw per-replicate estimates are written to disk per cell, so a long run can
## be inspected while it is still going, and any summary can be recomputed
## later without refitting.
##
## Usage: Rscript coverage_study.R <nrep> <ncore> <outdir> [cell1,cell2,...]
## ---------------------------------------------------------------------------

suppressPackageStartupMessages(library(parallel))

args   <- commandArgs(trailingOnly = TRUE)
NREP   <- if (length(args) >= 1) as.integer(args[1]) else 10L
NCORE  <- if (length(args) >= 2) as.integer(args[2]) else 6L
OUTDIR <- if (length(args) >= 3) args[3] else "."
ONLY   <- if (length(args) >= 4) strsplit(args[4], ",")[[1]] else NULL

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

LEVEL <- 0.95

## ---- cell definitions ------------------------------------------------------
## nlat/obs_per_lat are recorded so the results table can be read against the
## latent budget without recomputing it.
cells <- list(
  list(name = "gauss_fixed", family = "gaussian",   J = NA, re = FALSE,
       ncl = NA, per = NA, N = 400, sd_re = NA),
  list(name = "gauss_mm",    family = "gaussian",   J = NA, re = TRUE,
       ncl = 30, per = 20, sd_re = 0.6),
  list(name = "binom_mm",    family = "binomial",   J = NA, re = TRUE,
       ncl = 30, per = 20, sd_re = 0.6),
  list(name = "pois_mm",     family = "poisson",    J = NA, re = TRUE,
       ncl = 30, per = 20, sd_re = 0.6),
  list(name = "mn_J3_rich",  family = "multinomial", J = 3, re = TRUE,
       ncl = 30, per = 20, sd_re = 0.6),
  list(name = "mn_J3_thin",  family = "multinomial", J = 3, re = TRUE,
       ncl = 60, per = 5,  sd_re = 0.6),
  ## A latent-budget dose-response series at J = 5 (C = 4).  Observations per
  ## latent value are 1.5, 3.0, 5.0 and 10.0, straddling the ~3 working rule, so
  ## the study can say where the approximation starts to fail rather than only
  ## that it does somewhere.
  list(name = "mn_J5_thin",  family = "multinomial", J = 5, re = TRUE,
       ncl = 40, per = 6,  sd_re = 0.6),
  list(name = "mn_J5_mid",   family = "multinomial", J = 5, re = TRUE,
       ncl = 30, per = 12, sd_re = 0.6),
  list(name = "mn_J5_rich",  family = "multinomial", J = 5, re = TRUE,
       ncl = 30, per = 20, sd_re = 0.6),
  list(name = "mn_J5_vrich", family = "multinomial", J = 5, re = TRUE,
       ncl = 30, per = 40, sd_re = 0.6),
  list(name = "mn_J3_ar",    family = "multinomial", J = 3, re = TRUE,
       ncl = 25, Tt = 10, rep_per = 8, sd_re = 0.5, rho = 0.5, ar = TRUE)
)
names(cells) <- vapply(cells, function(z) z$name, "")
if (!is.null(ONLY)) cells <- cells[ONLY]

## ---- true fixed effects ----------------------------------------------------
## Fixed once per cell, shared by every replicate, so "truth" is a constant and
## coverage is a property of the intervals rather than of a moving target.
truth_of <- function(cell) {
  C <- if (identical(cell$family, "multinomial")) cell$J - 1L else 1L
  p <- 4L                                   # (Intercept), x1, grpb, grpc
  set.seed(99L)
  matrix(round(runif(p * C, -0.7, 0.7), 3), p, C)
}

## ---- one simulated data set ------------------------------------------------
gen <- function(cell, seed) {
  ## truth_of() sets its own seed so that truth is constant across replicates.
  ## It must therefore be called BEFORE the replicate seed, or it resets the
  ## stream and every replicate simulates identical data.
  C  <- if (identical(cell$family, "multinomial")) cell$J - 1L else 1L
  Bt <- truth_of(cell)
  set.seed(seed)

  if (isTRUE(cell$ar)) {
    ## group x time x replicate: several rows share each (group, time) cell, so
    ## the AR latent values are actually identified
    ncl <- cell$ncl; Tt <- cell$Tt; rp <- cell$rep_per
    N <- ncl * Tt * rp
    dd <- data.frame(
      g    = factor(rep(seq_len(ncl), each = Tt * rp)),
      time = rep(rep(seq_len(Tt), each = rp), ncl))
    dd$x1  <- rnorm(N)
    dd$grp <- factor(sample(c("a", "b", "c"), N, TRUE))
    X <- model.matrix(~ x1 + grp, dd)
    ## AR(1) latent field, one row per (group, time), C columns
    nb <- ncl * Tt
    Bar <- matrix(0, nb, C)
    for (g in seq_len(ncl)) {
      i0 <- (g - 1L) * Tt
      Bar[i0 + 1L, ] <- rnorm(C, 0, cell$sd_re)
      for (t in 2:Tt)
        Bar[i0 + t, ] <- cell$rho * Bar[i0 + t - 1L, ] +
          rnorm(C, 0, cell$sd_re * sqrt(1 - cell$rho^2))
    }
    idx <- (as.integer(dd$g) - 1L) * Tt + dd$time
    eta <- X %*% Bt + Bar[idx, , drop = FALSE]
    attr(dd, "ar") <- list(idx = as.integer(idx), n_group = ncl, Tt = Tt)
  } else if (isTRUE(cell$re)) {
    ncl <- cell$ncl; per <- cell$per; N <- ncl * per
    dd <- data.frame(g = factor(rep(seq_len(ncl), each = per)))
    dd$x1  <- rnorm(N)
    dd$grp <- factor(sample(c("a", "b", "c"), N, TRUE))
    X <- model.matrix(~ x1 + grp, dd)
    b <- matrix(rnorm(ncl * C, 0, cell$sd_re), ncl, C)
    eta <- X %*% Bt + b[as.integer(dd$g), , drop = FALSE]
  } else {
    N <- cell$N
    dd <- data.frame(x1 = rnorm(N),
                     grp = factor(sample(c("a", "b", "c"), N, TRUE)))
    X <- model.matrix(~ x1 + grp, dd)
    eta <- X %*% Bt
  }

  dd$y <- switch(cell$family,
    gaussian = as.numeric(eta) + rnorm(nrow(dd), 0, 1),
    poisson  = rpois(nrow(dd), exp(pmin(as.numeric(eta), 5))),
    binomial = rbinom(nrow(dd), 1, 1 / (1 + exp(-as.numeric(eta)))),
    multinomial = {
      Tc <- contr.sum(cell$J)
      P  <- exp(eta %*% t(Tc)); P <- P / rowSums(P)
      labs <- paste0("c", seq_len(cell$J))
      factor(labs[apply(P, 1, function(pr) sample.int(cell$J, 1L, prob = pr))],
             levels = labs)
    })
  dd
}

## ---- one replicate ---------------------------------------------------------
run_rep <- function(i, cell) {
  dd <- gen(cell, seed = 10000L + i)
  fm <- if (isTRUE(cell$re)) y ~ x1 + grp + (1 | g) else y ~ x1 + grp
  f <- tryCatch(suppressWarnings(
         illume::ilm_model(fm, data = dd, family = cell$family,
                           ar = attr(dd, "ar"), verbose = FALSE)),
       error = function(e) structure(list(msg = conditionMessage(e)),
                                     class = "failed"))
  if (inherits(f, "failed"))
    return(list(ok = FALSE, err = f$msg, b = NULL, s = NULL))

  b <- coef(f)
  s <- suppressWarnings(sqrt(diag(vcov(f))))
  ok <- isTRUE(f$opt$convergence == 0L) && isTRUE(f$sdr$pdHess) &&
        length(s) == length(b) && all(is.finite(s)) && all(s > 0)
  crit <- if (isTRUE(f$exact_df)) qt(1 - (1 - LEVEL) / 2, f$resid_df)
          else qnorm(1 - (1 - LEVEL) / 2)
  list(ok = ok, err = NA_character_, b = unname(b), s = unname(s),
       crit = crit, nm = names(b))
}

## ---- run each cell ---------------------------------------------------------
cl <- makePSOCKcluster(NCORE)
on.exit(stopCluster(cl), add = TRUE)
invisible(clusterEvalQ(cl, suppressPackageStartupMessages(library(illume))))
clusterExport(cl, c("gen", "truth_of", "run_rep", "LEVEL"), envir = environment())

for (cell in cells) {
  t0 <- Sys.time()
  clusterExport(cl, "cell", envir = environment())
  res <- parLapply(cl, seq_len(NREP), function(i) run_rep(i, cell))
  el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  tv  <- as.vector(truth_of(cell))
  okv <- vapply(res, function(z) isTRUE(z$ok), TRUE)
  nb  <- length(tv)
  B <- t(vapply(res, function(z)
         if (is.null(z$b)) rep(NA_real_, nb) else z$b[seq_len(nb)], numeric(nb)))
  S <- t(vapply(res, function(z)
         if (is.null(z$s)) rep(NA_real_, nb) else z$s[seq_len(nb)], numeric(nb)))
  cr <- vapply(res, function(z) if (is.null(z$crit)) NA_real_ else z$crit, 1)
  nm <- res[[which(!vapply(res, function(z) is.null(z$nm), TRUE))[1]]]$nm
  nm <- if (is.null(nm)) paste0("b", seq_len(nb)) else nm[seq_len(nb)]

  keep <- okv & stats::complete.cases(B) & stats::complete.cases(S)
  Bk <- B[keep, , drop = FALSE]; Sk <- S[keep, , drop = FALSE]; ck <- cr[keep]
  lo <- Bk - ck * Sk; hi <- Bk + ck * Sk
  cov_i <- sweep(lo, 2, tv, "<=") & sweep(hi, 2, tv, ">=")

  out <- data.frame(
    cell      = cell$name,
    family    = cell$family,
    J         = if (is.null(cell$J) || is.na(cell$J)) NA_integer_ else cell$J,
    coef      = nm,
    truth     = tv,
    n_attempt = NREP,
    n_used    = sum(keep),
    coverage  = colMeans(cov_i),
    bias      = colMeans(Bk) - tv,
    emp_sd    = apply(Bk, 2, sd),
    mean_se   = colMeans(Sk),
    stringsAsFactors = FALSE)
  out$se_ratio <- out$mean_se / out$emp_sd
  out$mc_se    <- sqrt(out$coverage * (1 - out$coverage) / out$n_used)
  out$secs     <- el
  ## the latent budget this cell was built to sit at, carried into the results
  ## so the table can be read against it directly
  Cc <- if (identical(cell$family, "multinomial")) cell$J - 1L else 1L
  if (isTRUE(cell$ar)) {
    out$N <- cell$ncl * cell$Tt * cell$rep_per
    out$n_latent <- cell$ncl * cell$Tt * Cc
  } else if (isTRUE(cell$re)) {
    out$N <- cell$ncl * cell$per
    out$n_latent <- cell$ncl * Cc
  } else {
    out$N <- cell$N; out$n_latent <- 0L
  }
  out$obs_per_latent <- ifelse(out$n_latent > 0, out$N / out$n_latent, NA_real_)
  out$conv_rate <- sum(keep) / NREP

  saveRDS(list(cell = cell, raw = res, summary = out),
          file.path(OUTDIR, paste0("cov_", cell$name, ".rds")))
  write.csv(out, file.path(OUTDIR, paste0("cov_", cell$name, ".csv")),
            row.names = FALSE)
  ## Why replicates were dropped matters as much as how many.  An outright
  ## error and a non-positive-definite Hessian are different failures needing
  ## different explanations, and coverage below is conditional on survival.
  errs <- unlist(lapply(res, function(z) if (!is.na(z$err)) z$err else NULL))
  n_err <- length(errs); n_nonpd <- sum(!okv) - n_err
  cat(sprintf("%-12s  n=%3d/%3d  cover %.3f  se_ratio %.3f  %.0fs\n",
              cell$name, sum(keep), NREP, mean(out$coverage),
              mean(out$se_ratio), el))
  if (n_err || n_nonpd)
    cat(sprintf("               dropped: %d error(s), %d non-pd/non-finite%s\n",
                n_err, n_nonpd,
                if (n_err) paste0("  [", substr(names(sort(table(errs),
                  decreasing = TRUE))[1], 1, 70), "]") else ""))
  flush.console()
}
cat("done\n")
