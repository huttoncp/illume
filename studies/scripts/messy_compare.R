## ---------------------------------------------------------------------------
## illume (Laplace) against mclogit::mblogit (PQL) on data that misbehaves.
##
## mclogit_compare.R asked the textbook question: does PQL attenuate when
## clusters are small and the random-effect variance is large?  It does, and
## illume does not.  But both regimes there were tidy -- balanced clusters,
## three roughly equal outcome categories, Gaussian random effects.  Real data
## is none of those things, and a method that only wins on tidy data has not
## earned much.
##
## So this varies the four things most likely to break a Laplace approximation
## or a PQL one, ONE AT A TIME and then all together.  One-at-a-time says
## which feature is responsible; combined says whether they compound.
##
##   clean       60 clusters x 8, sd = 1.0, Gaussian REs, balanced categories
##   unbalanced  the same n spread over wildly uneven cluster sizes, many
##               singletons -- a singleton cluster carries a random effect no
##               row can identify, which is where PQL's working weights and
##               the Laplace approximation's mode both get thin
##   rare        one outcome category at roughly 3% -- the multinomial case of
##               a sparse cell, where the likelihood is flat in one direction
##   flat_re     sd = 0.05, a variance component effectively at the boundary.
##               The honest answer is "about zero"; the failure mode is a
##               variance driven to exactly zero with a standard error that
##               pretends otherwise
##   heavy_re    random effects from a two-component mixture, standardised to
##               the same sd -- most clusters alike, a few far out.  Both
##               methods ASSUME Gaussian REs, so this is misspecification
##               rather than difficulty, and neither should be let off
##   combined    all four at once
##
## Reported per regime: coverage of the nominal 95% interval, attenuation,
## bias, RMSE, the share of fits that converged, and runtime.  Convergence
## rate matters as much as accuracy here: an estimator that is excellent on
## the 60% of data sets it can handle is not better than one that is decent on
## all of them, and a comparison conditioned on convergence hides that.
##
## Comparison is on mclogit's baseline-category scale, since illume's
## sum-to-zero coefficients are not directly comparable.
##
## Usage: Rscript messy_compare.R <nrep> <ncore> <outdir>
## ---------------------------------------------------------------------------

suppressPackageStartupMessages(library(parallel))

args   <- commandArgs(trailingOnly = TRUE)
NREP   <- if (length(args) >= 1) as.integer(args[1]) else 400L
NCORE  <- if (length(args) >= 2) as.integer(args[2]) else 10L
OUTDIR <- if (length(args) >= 3) args[3] else "."
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

J <- 3L; C <- J - 1L; P4 <- 4L
NCL <- 60L; PER <- 8L; NTOT <- NCL * PER

regimes <- list(
  clean      = list(name = "clean",      unbal = FALSE, rare = FALSE,
                    sd_re = 1.00, re_dist = "gauss"),
  unbalanced = list(name = "unbalanced", unbal = TRUE,  rare = FALSE,
                    sd_re = 1.00, re_dist = "gauss"),
  rare       = list(name = "rare",       unbal = FALSE, rare = TRUE,
                    sd_re = 1.00, re_dist = "gauss"),
  flat_re    = list(name = "flat_re",    unbal = FALSE, rare = FALSE,
                    sd_re = 0.05, re_dist = "gauss"),
  heavy_re   = list(name = "heavy_re",   unbal = FALSE, rare = FALSE,
                    sd_re = 1.00, re_dist = "mixture"),
  combined   = list(name = "combined",   unbal = TRUE,  rare = TRUE,
                    sd_re = 1.00, re_dist = "mixture")
)

## The same fixed effects throughout, so regimes are comparable; only the
## intercepts move, and only to make a category rare.  With sum-to-zero
## contrasts the full effect for category j is B %*% t(contr.sum(J))[, j], so
## category 3's intercept is -(B[1, 1] + B[1, 2]) -- raising both free
## intercepts together pushes category 3 down without touching the slopes.
truth_B <- function(rare = FALSE) {
  set.seed(99L)
  B <- matrix(round(runif(P4 * C, -0.7, 0.7), 3), P4, C)
  ## calibrated, not guessed: 1.2 puts the rarest category at 3.6% on
  ## average and never below 1.5%, i.e. about 17 of 480 observations -- sparse
  ## enough to be the hard case, never empty enough to be a different model
  if (rare) B[1, ] <- c(1.2, 1.2)
  B
}

to_baseline <- function(B) {
  full <- B %*% t(contr.sum(J))
  full[, -1, drop = FALSE] - full[, 1]        # P4 x (J-1)
}

## Cluster sizes: balanced, or a long tail that keeps the total the same.  The
## tail is the realistic shape -- a few big clusters, a lot of singletons --
## and the total is held fixed so that any difference is the SHAPE of the
## design rather than simply less data.
cluster_sizes <- function(unbal, ncl = NCL, ntot = NTOT) {
  if (!unbal) return(rep(ntot %/% ncl, ncl))
  w <- c(rep(1, 24), rep(2, 12), rep(3, 8), rep(5, 8), rep(12, 5), rep(30, 3))
  w <- w[seq_len(ncl)]
  s <- pmax(1L, as.integer(round(w * ntot / sum(w))))
  ## round-off lands wherever; put it on the largest cluster, which can absorb
  ## it without changing the shape
  d <- ntot - sum(s)
  s[which.max(s)] <- s[which.max(s)] + d
  s
}

## Random effects.  The mixture is standardised to the requested sd so that
## "non-normal" is not silently confounded with "more variable": the second
## moment is held fixed and only the shape changes.
draw_re <- function(ncl, sd_re, dist) {
  if (dist == "gauss") return(matrix(rnorm(ncl * C, 0, sd_re), ncl, C))
  z <- matrix(NA_real_, ncl, C)
  for (cc in seq_len(C)) {
    g <- stats::rbinom(ncl, 1L, 0.15)
    v <- ifelse(g == 1L, rnorm(ncl, 2.2, 0.6), rnorm(ncl, -0.4, 0.5))
    z[, cc] <- (v - mean(v)) / stats::sd(v) * sd_re
  }
  z
}

gen <- function(rg, seed) {
  Bt <- truth_B(rg$rare)
  set.seed(seed)
  sz <- cluster_sizes(rg$unbal)
  n  <- sum(sz)
  dd <- data.frame(g = factor(rep(seq_along(sz), times = sz)))
  dd$x1  <- rnorm(n)
  dd$grp <- factor(sample(c("a", "b", "c"), n, TRUE))
  X <- model.matrix(~ x1 + grp, dd)
  b <- draw_re(length(sz), rg$sd_re, rg$re_dist)
  eta <- X %*% Bt + b[as.integer(dd$g), , drop = FALSE]
  Pm <- exp(eta %*% t(contr.sum(J))); Pm <- Pm / rowSums(Pm)
  dd$y <- factor(paste0("c", apply(Pm, 1, function(p) sample.int(J, 1L, prob = p))),
                 levels = paste0("c", seq_len(J)))
  dd
}

target_names <- function() {
  xn <- c("(Intercept)", "x1", "grpb", "grpc")
  as.vector(vapply(2:J, function(j) paste0("c", j, "~", xn), character(P4)))
}

run_rep <- function(i, rg) {
  dd <- gen(rg, seed = 70000L + i)
  nb <- P4 * C
  ## What the data set actually looked like, recorded rather than assumed:
  ## the rarest category's share and the smallest cluster are the two features
  ## the regimes are meant to vary, and a generator can always be wrong.
  tabs <- table(dd$y)
  rare_share <- min(tabs) / nrow(dd)
  ## A category absent altogether is not a hard case, it is a different model.
  ## It is counted, and those replications are reported separately.
  degenerate <- any(tabs == 0L)

  t0 <- Sys.time()
  f <- tryCatch(suppressWarnings(
         illume::ilm_model(y ~ x1 + grp + (1 | g), data = dd,
                           family = "multinomial", verbose = FALSE)),
       error = function(e) NULL)
  t_ill <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  b_ill <- rep(NA_real_, nb); s_ill <- rep(NA_real_, nb)
  fit_ill <- !is.null(f); ok_ill <- FALSE; sd_ill <- NA_real_
  if (fit_ill) {
    ok_ill <- isTRUE(f$opt$convergence == 0L) && isTRUE(f$sdr$pdHess)
    sd_ill <- tryCatch(sqrt(f$Sigma[[1]][1, 1]), error = function(e) NA_real_)
    if (ok_ill) {
      b_ill <- as.vector(to_baseline(matrix(coef(f), P4, C)))
      Tc <- contr.sum(J)
      M  <- Tc[-1, , drop = FALSE] - matrix(Tc[1, ], J - 1L, C, byrow = TRUE)
      L  <- kronecker(M, diag(P4))
      Vb <- L %*% suppressWarnings(vcov(f)) %*% t(L)
      s_ill <- sqrt(pmax(diag(Vb), 0))
      ok_ill <- all(is.finite(b_ill)) && all(is.finite(s_ill)) && all(s_ill > 0)
    }
  }

  t0 <- Sys.time()
  m <- tryCatch(suppressWarnings(suppressMessages(
         mclogit::mblogit(y ~ x1 + grp, random = ~ 1 | g, data = dd))),
       error = function(e) NULL)
  t_mc <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  b_mc <- rep(NA_real_, nb); s_mc <- rep(NA_real_, nb)
  fit_mc <- !is.null(m); ok_mc <- FALSE; sd_mc <- NA_real_
  if (fit_mc) {
    cf <- tryCatch(coef(m), error = function(e) NULL)
    se <- tryCatch(summary(m)$coefficients[, 2], error = function(e) NULL)
    sd_mc <- tryCatch(sqrt(m$VarCov[[1]][1, 1]), error = function(e) NA_real_)
    tn <- target_names()
    if (!is.null(cf) && !is.null(se) && all(tn %in% names(cf)) &&
        all(tn %in% names(se))) {
      b_mc <- unname(cf[tn]); s_mc <- unname(se[tn])
      ok_mc <- all(is.finite(b_mc)) && all(is.finite(s_mc)) && all(s_mc > 0)
    }
  }
  list(fit_ill = fit_ill, fit_mc = fit_mc, ok_ill = ok_ill, ok_mc = ok_mc,
       b_ill = b_ill, b_mc = b_mc, s_ill = s_ill, s_mc = s_mc,
       sd_ill = sd_ill, sd_mc = sd_mc, t_ill = t_ill, t_mc = t_mc,
       rare_share = rare_share, degenerate = degenerate, n = nrow(dd))
}

cl <- makePSOCKcluster(NCORE)
on.exit(stopCluster(cl), add = TRUE)
invisible(clusterEvalQ(cl, {
  suppressPackageStartupMessages(library(illume))
  suppressPackageStartupMessages(library(mclogit))
}))
clusterExport(cl, c("gen", "truth_B", "to_baseline", "run_rep", "target_names",
                    "cluster_sizes", "draw_re", "J", "C", "P4", "NCL", "PER",
                    "NTOT"), envir = environment())

allrows <- list()
for (rg in regimes) {
  clusterExport(cl, "rg", envir = environment())
  res <- parLapply(cl, seq_len(NREP), function(i) run_rep(i, rg))
  tv <- as.vector(to_baseline(truth_B(rg$rare)))
  nb <- length(tv)

  gm <- function(bf, sf, okf) {
    keep <- vapply(res, function(z) isTRUE(z[[okf]]), TRUE)
    M <- t(vapply(res[keep], function(z) z[[bf]], numeric(nb)))
    S <- t(vapply(res[keep], function(z) z[[sf]], numeric(nb)))
    cv <- if (sum(keep) > 1) {
      lo <- M - 1.959964 * S; hi <- M + 1.959964 * S
      colMeans(sweep(lo, 2, tv, "<=") & sweep(hi, 2, tv, ">="))
    } else rep(NA_real_, nb)
    list(n = sum(keep), M = M, S = S, cov = cv)
  }
  A <- gm("b_ill", "s_ill", "ok_ill"); Bb <- gm("b_mc", "s_mc", "ok_mc")
  att <- function(M) if (nrow(M) > 1) sum(colMeans(M) * tv) / sum(tv^2) else NA_real_
  pull <- function(f) vapply(res, function(z) z[[f]], if (f %in%
    c("fit_ill", "fit_mc", "degenerate")) TRUE else 1)

  r <- data.frame(
    regime = rg$name, sd_re = rg$sd_re, re_dist = rg$re_dist,
    unbal = rg$unbal, rare = rg$rare,
    coef = target_names(), truth = tv,
    nrep = NREP,
    ## three distinct things, kept apart: the fit threw an error; the fit
    ## returned but the optimiser or the Hessian said no; the fit is usable
    fit_illume  = mean(pull("fit_ill")),
    fit_mclogit = mean(pull("fit_mc")),
    conv_illume  = A$n  / NREP,
    conv_mclogit = Bb$n / NREP,
    n_illume = A$n, n_mclogit = Bb$n,
    degenerate = mean(pull("degenerate")),
    rare_share = mean(pull("rare_share")),
    n_obs = mean(pull("n")),
    bias_illume  = if (A$n  > 1) colMeans(A$M)  - tv else NA_real_,
    bias_mclogit = if (Bb$n > 1) colMeans(Bb$M) - tv else NA_real_,
    sd_illume    = if (A$n  > 1) apply(A$M, 2, sd) else NA_real_,
    sd_mclogit   = if (Bb$n > 1) apply(Bb$M, 2, sd) else NA_real_,
    cover_illume  = A$cov,
    cover_mclogit = Bb$cov,
    mean_se_illume  = if (A$n  > 1) colMeans(A$S)  else NA_real_,
    mean_se_mclogit = if (Bb$n > 1) colMeans(Bb$S) else NA_real_,
    att_illume  = att(A$M),
    att_mclogit = att(Bb$M),
    ## the variance component itself, which is what "flat_re" is about
    resd_illume  = median(pull("sd_ill"), na.rm = TRUE),
    resd_mclogit = median(pull("sd_mc"), na.rm = TRUE),
    t_illume  = median(pull("t_ill")),
    t_mclogit = median(pull("t_mc")),
    stringsAsFactors = FALSE)
  r$rmse_illume  <- sqrt(r$bias_illume^2  + r$sd_illume^2)
  r$rmse_mclogit <- sqrt(r$bias_mclogit^2 + r$sd_mclogit^2)
  allrows[[length(allrows) + 1L]] <- r

  ## Monte Carlo error on a coverage estimate, so that a difference can be
  ## read as a difference rather than as noise
  mcse <- function(p, n) if (is.finite(p) && n > 1) sqrt(p * (1 - p) / n) else NA_real_
  cat(sprintf("%-11s n=%.0f  rare share %.3f  usable: illume %d/%d  mclogit %d/%d\n",
              rg$name, r$n_obs[1], r$rare_share[1], A$n, NREP, Bb$n, NREP))
  cat(sprintf("            COVERAGE    illume %.3f (+-%.3f)   mclogit %.3f (+-%.3f)\n",
              mean(r$cover_illume), mcse(mean(r$cover_illume), A$n),
              mean(r$cover_mclogit), mcse(mean(r$cover_mclogit), Bb$n)))
  cat(sprintf("            attenuation illume %.3f            mclogit %.3f   (1 = none)\n",
              r$att_illume[1], r$att_mclogit[1]))
  cat(sprintf("            mean|bias|  illume %.4f            mclogit %.4f\n",
              mean(abs(r$bias_illume)), mean(abs(r$bias_mclogit))))
  cat(sprintf("            RMSE        illume %.4f            mclogit %.4f\n",
              mean(r$rmse_illume), mean(r$rmse_mclogit)))
  cat(sprintf("            RE sd       illume %.3f            mclogit %.3f   (true %.2f)\n",
              r$resd_illume[1], r$resd_mclogit[1], rg$sd_re))
  cat(sprintf("            median secs illume %.2f             mclogit %.2f\n",
              r$t_illume[1], r$t_mclogit[1]))
  flush.console()
}
out <- do.call(rbind, allrows)
write.csv(out, file.path(OUTDIR, "messy_compare.csv"), row.names = FALSE)
saveRDS(out, file.path(OUTDIR, "messy_compare.rds"))
cat("done\n")
