## ---------------------------------------------------------------------------
## illume (Laplace) against mclogit::mblogit (PQL) for multinomial mixed models.
##
## This is the only head-to-head frequentist comparison available: mclogit is
## the sole other R package fitting a random-effects baseline-category logit.
##
## The comparison is not arbitrary.  Penalised quasi-likelihood is known to
## attenuate fixed effects for discrete outcomes, worst when clusters are small
## and the random-effect variance is large (Breslow & Lin 1995; Rodriguez &
## Goldman 1995).  Two regimes are used:
##
##   "large"  40 clusters x 25 observations, sd = 0.7  -- PQL should be fine
##   "small" 100 clusters x  4 observations, sd = 1.2  -- PQL should attenuate
##
## If the literature is right, the two methods agree in the first regime and
## separate in the second, with illume closer to truth.  If illume is wrong,
## this is where it shows.
##
## Comparison is on mclogit's baseline-category scale, since illume's
## sum-to-zero coefficients are not directly comparable.
##
## Usage: Rscript mclogit_compare.R <nrep> <ncore> <outdir>
## ---------------------------------------------------------------------------

suppressPackageStartupMessages(library(parallel))

args   <- commandArgs(trailingOnly = TRUE)
NREP   <- if (length(args) >= 1) as.integer(args[1]) else 20L
NCORE  <- if (length(args) >= 2) as.integer(args[2]) else 6L
OUTDIR <- if (length(args) >= 3) args[3] else "."
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

J <- 3L; C <- J - 1L; P4 <- 4L

regimes <- list(
  large = list(name = "large", ncl = 40L, per = 25L, sd_re = 0.7),
  small = list(name = "small", ncl = 100L, per = 4L, sd_re = 1.2)
)

truth_B <- function() { set.seed(99L); matrix(round(runif(P4 * C, -0.7, 0.7), 3), P4, C) }

## sum-to-zero coefficients expressed as baseline-category contrasts vs c1
to_baseline <- function(B) {
  full <- B %*% t(contr.sum(J))
  full[, -1, drop = FALSE] - full[, 1]        # P4 x (J-1)
}

gen <- function(rg, seed) {
  Bt <- truth_B()
  set.seed(seed)
  n <- rg$ncl * rg$per
  dd <- data.frame(g = factor(rep(seq_len(rg$ncl), each = rg$per)))
  dd$x1  <- rnorm(n)
  dd$grp <- factor(sample(c("a", "b", "c"), n, TRUE))
  X <- model.matrix(~ x1 + grp, dd)
  b <- matrix(rnorm(rg$ncl * C, 0, rg$sd_re), rg$ncl, C)
  eta <- X %*% Bt + b[as.integer(dd$g), , drop = FALSE]
  Pm <- exp(eta %*% t(contr.sum(J))); Pm <- Pm / rowSums(Pm)
  dd$y <- factor(paste0("c", apply(Pm, 1, function(p) sample.int(J, 1L, prob = p))),
                 levels = paste0("c", seq_len(J)))
  dd
}

## mblogit names coefficients "c2~x1"; put them in illume's column-major order
target_names <- function() {
  xn <- c("(Intercept)", "x1", "grpb", "grpc")
  as.vector(vapply(2:J, function(j) paste0("c", j, "~", xn), character(P4)))
}

run_rep <- function(i, rg) {
  dd <- gen(rg, seed = 50000L + i)
  nb <- P4 * C
  t0 <- Sys.time()
  f <- tryCatch(suppressWarnings(
         illume::ilm_model(y ~ x1 + grp + (1 | g), data = dd,
                           family = "multinomial", verbose = FALSE)),
       error = function(e) NULL)
  t_ill <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  b_ill <- rep(NA_real_, nb); s_ill <- rep(NA_real_, nb); ok_ill <- FALSE
  if (!is.null(f)) {
    ## usable as the package says: a positive definite Hessian, or a
    ## covariance at its boundary with the unresolvable direction held,
    ## which has pdHess FALSE by design
    ok_ill <- isTRUE(f$opt$convergence == 0L) &&
              (isTRUE(f$sdr$pdHess) || length(f$hessian_held) > 0L)
    if (ok_ill) {
      b_ill <- as.vector(to_baseline(matrix(coef(f), P4, C)))
      ## The baseline contrasts are a linear map of the sum-to-zero
      ## coefficients acting on the category index only, so the covariance
      ## transforms as L V L'.  Comparing raw sum-to-zero SEs with mclogit's
      ## would be comparing different quantities.
      Tc <- contr.sum(J)
      M  <- Tc[-1, , drop = FALSE] - matrix(Tc[1, ], J - 1L, C, byrow = TRUE)
      L  <- kronecker(M, diag(P4))
      Vb <- L %*% suppressWarnings(vcov(f)) %*% t(L)
      s_ill <- sqrt(pmax(diag(Vb), 0))
      ok_ill <- all(is.finite(s_ill)) && all(s_ill > 0)
    }
  }

  t0 <- Sys.time()
  m <- tryCatch(suppressWarnings(suppressMessages(
         mclogit::mblogit(y ~ x1 + grp, random = ~ 1 | g, data = dd))),
       error = function(e) NULL)
  t_mc <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  b_mc <- rep(NA_real_, nb); s_mc <- rep(NA_real_, nb); ok_mc <- FALSE
  if (!is.null(m)) {
    cf <- tryCatch(coef(m), error = function(e) NULL)
    se <- tryCatch(summary(m)$coefficients[, 2], error = function(e) NULL)
    tn <- target_names()
    if (!is.null(cf) && !is.null(se) && all(tn %in% names(cf)) &&
        all(tn %in% names(se))) {
      b_mc <- unname(cf[tn]); s_mc <- unname(se[tn])
      ok_mc <- all(is.finite(b_mc)) && all(is.finite(s_mc)) && all(s_mc > 0)
    }
  }
  list(ok_ill = ok_ill, ok_mc = ok_mc, b_ill = b_ill, b_mc = b_mc,
       s_ill = s_ill, s_mc = s_mc, t_ill = t_ill, t_mc = t_mc)
}

cl <- makePSOCKcluster(NCORE)
on.exit(stopCluster(cl), add = TRUE)
invisible(clusterEvalQ(cl, {
  suppressPackageStartupMessages(library(illume))
  suppressPackageStartupMessages(library(mclogit))
}))
clusterExport(cl, c("gen", "truth_B", "to_baseline", "run_rep", "target_names",
                    "J", "C", "P4"), envir = environment())

allrows <- list()
for (rg in regimes) {
  clusterExport(cl, "rg", envir = environment())
  res <- parLapply(cl, seq_len(NREP), function(i) run_rep(i, rg))
  tv <- as.vector(to_baseline(truth_B()))
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
  ## Attenuation as the least-squares slope of estimate on truth.  A mean of
  ## per-coefficient ratios is useless here: some true values are near zero.
  att <- function(M) if (nrow(M) > 1) sum(colMeans(M) * tv) / sum(tv^2) else NA_real_

  r <- data.frame(
    regime = rg$name, ncl = rg$ncl, per = rg$per, sd_re = rg$sd_re,
    coef = target_names(), truth = tv,
    n_illume = A$n, n_mclogit = Bb$n,
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
    t_illume  = median(vapply(res, function(z) z$t_ill, 1)),
    t_mclogit = median(vapply(res, function(z) z$t_mc, 1)),
    stringsAsFactors = FALSE)
  ## RMSE captures bias and variance together, which is what actually matters
  r$rmse_illume  <- sqrt(r$bias_illume^2  + r$sd_illume^2)
  r$rmse_mclogit <- sqrt(r$bias_mclogit^2 + r$sd_mclogit^2)
  ## attenuation: ratio of fitted to true magnitude, averaged over slopes only
  allrows[[length(allrows) + 1L]] <- r

  cat(sprintf("%-6s ncl=%3d per=%2d  n=%d/%d\n", rg$name, rg$ncl, rg$per, A$n, Bb$n))
  cat(sprintf("        mean|bias|  illume %.4f   mclogit %.4f\n",
              mean(abs(r$bias_illume)), mean(abs(r$bias_mclogit))))
  cat(sprintf("        RMSE        illume %.4f   mclogit %.4f\n",
              mean(r$rmse_illume), mean(r$rmse_mclogit)))
  cat(sprintf("        attenuation illume %.3f    mclogit %.3f   (1 = none)\n",
              r$att_illume[1], r$att_mclogit[1]))
  cat(sprintf("        COVERAGE    illume %.3f    mclogit %.3f   (nominal 0.95)\n",
              mean(r$cover_illume), mean(r$cover_mclogit)))
  cat(sprintf("        median secs illume %.2f     mclogit %.2f\n",
              r$t_illume[1], r$t_mclogit[1]))
  flush.console()
}
out <- do.call(rbind, allrows)
write.csv(out, file.path(OUTDIR, "mclogit_compare.csv"), row.names = FALSE)
saveRDS(out, file.path(OUTDIR, "mclogit_compare.rds"))
cat("done\n")
