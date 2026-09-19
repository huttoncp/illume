## ---------------------------------------------------------------------------
## illume against the established frequentist packages.
##
## Two things are measured, and they matter for different reasons:
##
##   AGREEMENT is a validity check.  Where another package fits the same model,
##   illume should recover the same estimates.  lme4 and glmmTMB are independent
##   implementations, so agreement is evidence the engine is right, not merely
##   self-consistent.  nnet::multinom checks the multinomial likelihood itself,
##   on the fixed-effects case where a reference exists at all.
##
##   TIMING is a practicality check, and is the weaker of the two: these
##   packages are not solving identical problems (different optimisers,
##   convergence rules and parameterisations), so the numbers say whether
##   illume is in a usable range, not who "wins".
##
## nnet uses baseline-category coding and illume uses sum-to-zero, so the
## comparison converts illume's coefficients to the baseline scale first.
## Usage: Rscript benchmark_freq.R <nrep> <outdir>
## ---------------------------------------------------------------------------

args   <- commandArgs(trailingOnly = TRUE)
NREP   <- if (length(args) >= 1) as.integer(args[1]) else 3L
OUTDIR <- if (length(args) >= 2) args[2] else "."
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({
  library(illume); library(lme4); library(glmmTMB); library(nnet)
})

tm <- function(expr) {
  t0 <- Sys.time()
  v  <- tryCatch(suppressWarnings(suppressMessages(force(expr))),
                 error = function(e) structure(list(m = conditionMessage(e)),
                                               class = "bmfail"))
  list(secs = as.numeric(difftime(Sys.time(), t0, units = "secs")), val = v)
}
okv <- function(x) !inherits(x, "bmfail")

sim_uni <- function(family, ncl, per, seed) {
  set.seed(seed)
  N <- ncl * per
  dd <- data.frame(g = factor(rep(seq_len(ncl), each = per)))
  dd$x1 <- rnorm(N); dd$grp <- factor(sample(c("a", "b", "c"), N, TRUE))
  X <- model.matrix(~ x1 + grp, dd)
  bt <- c(0.3, 0.5, -0.4, 0.2)
  eta <- as.numeric(X %*% bt) + rnorm(ncl, 0, 0.6)[as.integer(dd$g)]
  dd$y <- switch(family,
    gaussian = eta + rnorm(N, 0, 1),
    poisson  = rpois(N, exp(pmin(eta, 5))),
    binomial = rbinom(N, 1, 1 / (1 + exp(-eta))))
  dd
}

sim_mn <- function(J, ncl, per, seed, re = TRUE) {
  set.seed(seed)
  C <- J - 1L; N <- ncl * per
  dd <- data.frame(g = factor(rep(seq_len(ncl), each = per)))
  dd$x1 <- rnorm(N); dd$grp <- factor(sample(c("a", "b", "c"), N, TRUE))
  X <- model.matrix(~ x1 + grp, dd)
  set.seed(99L); Bt <- matrix(round(runif(4 * C, -0.7, 0.7), 3), 4, C)
  set.seed(seed + 1L)
  eta <- X %*% Bt
  if (re) eta <- eta + matrix(rnorm(ncl * C, 0, 0.6), ncl, C)[as.integer(dd$g), , drop = FALSE]
  Tc <- contr.sum(J); P <- exp(eta %*% t(Tc)); P <- P / rowSums(P)
  labs <- paste0("c", seq_len(J))
  dd$y <- factor(labs[apply(P, 1, function(pr) sample.int(J, 1L, prob = pr))],
                 levels = labs)
  dd
}

## illume's sum-to-zero coefficients on nnet's baseline-category scale:
## eta_full[, j] = B %*% Tc[j, ], then differences from the first category
to_baseline <- function(B, J) {
  Tc <- contr.sum(J)
  full <- B %*% t(Tc)                      # p x J
  full[, -1, drop = FALSE] - full[, 1]     # p x (J-1), vs category 1
}

rows <- list()
add <- function(...) rows[[length(rows) + 1L]] <<- data.frame(..., stringsAsFactors = FALSE)

## ---- univariate families: lme4 and glmmTMB ---------------------------------
for (fam in c("gaussian", "poisson", "binomial")) {
  for (r in seq_len(NREP)) {
    dd <- sim_uni(fam, 40, 25, seed = 1000L + r)
    a <- tm(ilm_model(y ~ x1 + grp + (1 | g), data = dd, family = fam,
                      verbose = FALSE))
    b <- tm(if (fam == "gaussian") lmer(y ~ x1 + grp + (1 | g), data = dd,
                                        REML = FALSE)
            else glmer(y ~ x1 + grp + (1 | g), data = dd,
                       family = if (fam == "poisson") poisson() else binomial()))
    cc <- tm(glmmTMB(y ~ x1 + grp + (1 | g), data = dd,
                     family = if (fam == "gaussian") gaussian()
                              else if (fam == "poisson") poisson() else binomial()))
    ca <- if (okv(a$val)) coef(a$val) else rep(NA, 4)
    cb <- if (okv(b$val)) lme4::fixef(b$val) else rep(NA, 4)
    ccf <- if (okv(cc$val)) glmmTMB::fixef(cc$val)$cond else rep(NA, 4)
    add(task = paste0(fam, "_mm"), rep = r, N = nrow(dd),
        illume_s = a$secs, lme4_s = b$secs, glmmTMB_s = cc$secs,
        max_absdiff_lme4 = max(abs(ca - cb)),
        max_absdiff_glmmTMB = max(abs(ca - ccf)))
    cat(sprintf("%-12s rep %d  illume %.2fs  lme4 %.2fs  glmmTMB %.2fs  maxdiff %.2e / %.2e\n",
                fam, r, a$secs, b$secs, cc$secs,
                max(abs(ca - cb)), max(abs(ca - ccf)))); flush.console()
  }
}

## ---- multinomial, no random effects: nnet::multinom ------------------------
for (J in c(3L, 5L)) {
  for (r in seq_len(NREP)) {
    dd <- sim_mn(J, 40, 25, seed = 2000L + r, re = FALSE)
    a <- tm(ilm_model(y ~ x1 + grp, data = dd, family = "multinomial",
                      verbose = FALSE))
    b <- tm(multinom(y ~ x1 + grp, data = dd, trace = FALSE))
    d <- NA_real_
    if (okv(a$val) && okv(b$val)) {
      B <- matrix(coef(a$val), 4L, J - 1L)          # p x C, category-major
      d <- max(abs(to_baseline(B, J) - t(coef(b$val))))
    }
    add(task = paste0("multinom_fixed_J", J), rep = r, N = nrow(dd),
        illume_s = a$secs, nnet_s = b$secs, max_absdiff_nnet = d)
    cat(sprintf("mn fixed J=%d rep %d  illume %.2fs  nnet %.2fs  maxdiff %.2e\n",
                J, r, a$secs, b$secs, d)); flush.console()
  }
}

## ---- multinomial mixed: no frequentist competitor --------------------------
## Recorded for timing only; brms is handled separately because it compiles.
for (J in c(3L, 5L)) {
  for (r in seq_len(NREP)) {
    dd <- sim_mn(J, 40, 25, seed = 3000L + r, re = TRUE)
    a <- tm(ilm_model(y ~ x1 + grp + (1 | g), data = dd,
                      family = "multinomial", verbose = FALSE))
    add(task = paste0("multinom_mixed_J", J), rep = r, N = nrow(dd),
        illume_s = a$secs)
    cat(sprintf("mn mixed J=%d rep %d  illume %.2fs\n", J, r, a$secs))
    flush.console()
  }
}

## ---- how illume scales with n ---------------------------------------------
for (ncl in c(20L, 40L, 80L, 160L)) {
  dd <- sim_mn(3L, ncl, 25, seed = 4000L, re = TRUE)
  a <- tm(ilm_model(y ~ x1 + grp + (1 | g), data = dd,
                    family = "multinomial", verbose = FALSE))
  add(task = "scaling_J3", rep = ncl, N = nrow(dd), illume_s = a$secs)
  cat(sprintf("scaling  clusters %3d  N %5d  illume %.2fs\n",
              ncl, nrow(dd), a$secs)); flush.console()
}

## different tasks record different columns (there is no lme4 time for a
## multinomial fit), so pad to the union rather than merging on shared keys
allnm <- unique(unlist(lapply(rows, names)))
res <- do.call(rbind, lapply(rows, function(d) {
  miss <- setdiff(allnm, names(d))
  for (m in miss) d[[m]] <- NA
  d[, allnm, drop = FALSE]
}))
write.csv(res, file.path(OUTDIR, "benchmark_freq.csv"), row.names = FALSE)
saveRDS(rows, file.path(OUTDIR, "benchmark_freq.rds"))
cat("done\n")
