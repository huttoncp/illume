## Item 4: at a covariance boundary, which standard errors for the fixed
## effects? Every dataset is fitted ONCE (boundary = "hold", the default), and
## at that one optimum the fixed-effect covariance is computed three ways:
##   B  the fit's own standard errors -- when this was first run, the status
##      quo: TMB's Hessian inverted when it is positive definite ("tmb"), the
##      whole boundary term held when it is not
##   A  always hold the whole boundary term, as lme4 does
##   C  hold only the flat directions of the boundary term: the eigenvectors of
##      its block of the accurate Hessian with curvature below a threshold,
##      which is the reduced model's covariance (checked against a refit of
##      that model where the package can fit it)
## The paired design isolates the choice of standard error: the estimates are
## the same in all three.
##
## Written for item 4 of NEXT-SESSION.md (which standard errors at a
## boundary); first run 2026-09-23 on 0.0.7.9000, R 4.4.3, 400 replicates,
## 8.6 minutes on 10 cores. That run chose C, and C at the 1e-3 threshold
## (C3 below) is now what the package does (ilm_hess_recover()), so from then
## on B -- the fit's own standard errors -- IS C3, and a re-run checks that
## it is: B should equal C3 in every boundary fit, and A still shows what
## holding the whole term would cost.
## Not yet wired into summarise_run.R: the tables are printed at the end of
## this script.
##
## Usage: Rscript boundary_se.R <nrep> <ncore> <outdir> <package root>
## ---------------------------------------------------------------------------
args  <- commandArgs(trailingOnly = TRUE)
NREP  <- if (length(args) >= 1) as.integer(args[1]) else 400L
NCORE <- if (length(args) >= 2) as.integer(args[2]) else 10L
sp    <- if (length(args) >= 3) args[3] else "."
pkg   <- normalizePath(if (length(args) >= 4) args[4] else ".")
dir.create(sp, showWarnings = FALSE, recursive = TRUE)
J <- 3L; C <- J - 1L; P4 <- 4L; NCL <- 60L; PER <- 8L; NTOT <- NCL * PER

## ---- multinomial regimes: the messy-data study's generator and seeds -------
mn_regimes <- list(
  clean      = list(unbal = FALSE, rare = FALSE, sd_re = 1.00, re_dist = "gauss"),
  unbalanced = list(unbal = TRUE,  rare = FALSE, sd_re = 1.00, re_dist = "gauss"),
  rare       = list(unbal = FALSE, rare = TRUE,  sd_re = 1.00, re_dist = "gauss"),
  flat_re    = list(unbal = FALSE, rare = FALSE, sd_re = 0.05, re_dist = "gauss"),
  heavy_re   = list(unbal = FALSE, rare = FALSE, sd_re = 1.00, re_dist = "mixture"),
  combined   = list(unbal = TRUE,  rare = TRUE,  sd_re = 1.00, re_dist = "mixture"),
  ## new: the category covariance truly rank one (correlation +1, SDs 1)
  rank_one   = list(unbal = FALSE, rare = FALSE, sd_re = 1.00, re_dist = "rank1"))
truth_B <- function(rare = FALSE) {
  set.seed(99L); B <- matrix(round(runif(P4 * C, -0.7, 0.7), 3), P4, C)
  if (rare) B[1, ] <- c(1.2, 1.2); B }
to_baseline <- function(B) { full <- B %*% t(contr.sum(J)); full[, -1, drop = FALSE] - full[, 1] }
cluster_sizes <- function(unbal, ncl = NCL, ntot = NTOT) {
  if (!unbal) return(rep(ntot %/% ncl, ncl))
  w <- c(rep(1, 24), rep(2, 12), rep(3, 8), rep(5, 8), rep(12, 5), rep(30, 3))
  w <- w[seq_len(ncl)]; s <- pmax(1L, as.integer(round(w * ntot / sum(w))))
  d <- ntot - sum(s); s[which.max(s)] <- s[which.max(s)] + d; s }
draw_re <- function(ncl, sd_re, dist) {
  if (dist == "gauss") return(matrix(rnorm(ncl * C, 0, sd_re), ncl, C))
  if (dist == "rank1") { z <- rnorm(ncl, 0, sd_re); return(cbind(z, z)) }
  z <- matrix(NA_real_, ncl, C)
  for (cc in seq_len(C)) {
    g <- stats::rbinom(ncl, 1L, 0.15)
    v <- ifelse(g == 1L, rnorm(ncl, 2.2, 0.6), rnorm(ncl, -0.4, 0.5))
    z[, cc] <- (v - mean(v)) / stats::sd(v) * sd_re }
  z }
gen_mn <- function(rg, seed) {
  Bt <- truth_B(rg$rare); set.seed(seed)
  sz <- cluster_sizes(rg$unbal); n <- sum(sz)
  dd <- data.frame(g = factor(rep(seq_along(sz), times = sz)))
  dd$x1 <- rnorm(n); dd$grp <- factor(sample(c("a", "b", "c"), n, TRUE))
  X <- model.matrix(~ x1 + grp, dd)
  b <- draw_re(length(sz), rg$sd_re, rg$re_dist)
  eta <- X %*% Bt + b[as.integer(dd$g), , drop = FALSE]
  Pm <- exp(eta %*% t(contr.sum(J))); Pm <- Pm / rowSums(Pm)
  dd$y <- factor(paste0("c", apply(Pm, 1, function(p) sample.int(J, 1L, prob = p))),
                 levels = paste0("c", seq_len(J)))
  dd }

## ---- binomial random-slope regimes -----------------------------------------
sl_regimes <- list(
  slope_zero  = list(sd0 = 1, sd1 = 0.00, cor = 0),   # slope variance truly 0
  slope_small = list(sd0 = 1, sd1 = 0.25, cor = 0),   # small but real
  slope_cor1  = list(sd0 = 1, sd1 = 0.50, cor = 1))   # intercept-slope corr. 1
beta_sl <- c(-0.3, 0.5)
gen_sl <- function(rg, seed) {
  set.seed(seed); n <- NCL * PER
  dd <- data.frame(g = factor(rep(seq_len(NCL), each = PER)), x = rnorm(n))
  z0 <- rnorm(NCL); z1 <- rg$cor * z0 + sqrt(1 - rg$cor^2) * rnorm(NCL)
  b0 <- rg$sd0 * z0; b1 <- rg$sd1 * z1
  eta <- beta_sl[1] + beta_sl[2] * dd$x + b0[dd$g] + b1[dd$g] * dd$x
  dd$y <- rbinom(n, 1L, plogis(eta))
  dd }

regimes <- c(names(mn_regimes), names(sl_regimes))

## ---- the three covariances at one optimum ----------------------------------
## positive definite on the package's own scale-free test
pd <- function(M) {
  dg <- diag(M)
  if (!length(dg) || any(!is.finite(dg)) || any(dg <= 0)) return(FALSE)
  S <- M / sqrt(outer(dg, dg))
  ev <- tryCatch(eigen(S, symmetric = TRUE, only.values = TRUE)$values,
                 error = function(e) NA_real_)
  all(is.finite(ev)) && min(ev) > 1e-6
}
three_cov <- function(f) {
  pn <- names(f$opt$par); ib <- which(pn == "beta")
  cb <- ilm_cov_blocks(f$re, f$Sigma, f$Sigma_d, f$ty, f$rk, f$dk, f$toff,
                       f$ar, f$opt$par, pn)
  cf <- f$sdr$cov.fixed
  VB <- if (!is.null(cf) && all(dim(cf) >= max(ib))) cf[ib, ib, drop = FALSE] else NULL
  out <- list(flagged = cb$flagged, how = f$hessian_how, B = VB, A = VB, C = VB,
              C3 = VB, nflat = NA_integer_, nflat3 = NA_integer_, ev = NULL)
  if (!length(cb$flagged)) return(out)
  H <- tryCatch(ilm_hessian(function(p) as.numeric(f$obj$gr(p)), f$opt$par),
                error = function(e) NULL)
  g0 <- tryCatch(abs(as.numeric(f$obj$gr(f$opt$par))), error = function(e) NULL)
  invisible(tryCatch(f$obj$fn(f$opt$par), error = function(e) NULL))
  out$A <- NULL; out$C <- NULL; out$C3 <- NULL
  if (is.null(H) || !all(is.finite(H)) || is.null(g0)) return(out)
  d <- unlist(cb$blocks[cb$flagged], use.names = FALSE)
  kp <- setdiff(seq_along(pn), d)
  ## A: the whole term held -- as the package's "boundary" route does
  if (max(g0[kp]) <= 1e-2 && pd(H[kp, kp, drop = FALSE]))
    out$A <- solve(H[kp, kp, drop = FALSE])[match(ib, kp), match(ib, kp), drop = FALSE]
  ## C: only the flat directions of the term held. Flat is judged against
  ## the block's own largest curvature, so the test does not depend on how
  ## much data there is; two thresholds, to see whether the choice matters.
  e <- eigen(H[d, d, drop = FALSE], symmetric = TRUE)
  out$ev <- e$values
  gfull <- as.numeric(f$obj$gr(f$opt$par))
  invisible(tryCatch(f$obj$fn(f$opt$par), error = function(e) NULL))
  hold_flat <- function(rel) {
    flat <- e$values < rel * max(e$values[1], 1e-8)
    G <- e$vectors[, !flat, drop = FALSE]
    Tm <- matrix(0, length(pn), length(kp) + ncol(G))
    Tm[cbind(kp, seq_along(kp))] <- 1
    if (ncol(G)) Tm[d, length(kp) + seq_len(ncol(G))] <- G
    HC <- t(Tm) %*% H %*% Tm
    gC <- abs(as.numeric(t(Tm) %*% gfull))
    V <- if (max(gC) <= 1e-2 && pd(HC))
      solve(HC)[match(ib, kp), match(ib, kp), drop = FALSE] else NULL
    list(V = V, nflat = sum(flat))
  }
  c4 <- hold_flat(1e-4); c3 <- hold_flat(1e-3)
  out$C <- c4$V; out$C3 <- c3$V
  out$nflat <- c4$nflat; out$nflat3 <- c3$nflat
  out
}

one <- function(job) {
  r <- job$regime; i <- job$i
  mn <- r %in% names(mn_regimes)
  dd <- if (mn) gen_mn(mn_regimes[[r]], 70000L + i) else gen_sl(sl_regimes[[r]], 80000L + i)
  f <- tryCatch(suppressMessages(suppressWarnings(
    if (mn) ilm_model(y ~ x1 + grp + (1 | g), data = dd, family = "multinomial",
                      verbose = FALSE)
    else ilm_model(y ~ x + (1 + x | g), data = dd, family = "binomial",
                   verbose = FALSE))), error = function(e) NULL)
  if (is.null(f)) return(list(regime = r, i = i, ok = FALSE))
  tc <- tryCatch(three_cov(f), error = function(e) NULL)
  if (is.null(tc)) return(list(regime = r, i = i, ok = FALSE))
  if (mn) {
    Tc <- contr.sum(J)
    M <- Tc[-1, , drop = FALSE] - matrix(Tc[1, ], J - 1L, C, byrow = TRUE)
    L <- kronecker(M, diag(P4))
    b <- as.vector(to_baseline(matrix(coef(f), P4, C)))
  } else { L <- diag(2); b <- as.numeric(coef(f)) }
  se <- function(V) if (is.null(V)) rep(NA_real_, length(b)) else
    sqrt(pmax(diag(L %*% V %*% t(L)), 0))
  ## the reduced model, where the package can fit it: C's check
  se_red <- rep(NA_real_, length(b))
  if (length(tc$flagged) && mn) {
    ev <- eigen(f$Sigma$g, symmetric = TRUE, only.values = TRUE)$values
    rk <- sum(ev > 1e-6 * max(ev, .Machine$double.eps) & sqrt(pmax(ev, 0)) > 1e-3)
    fr <- tryCatch(suppressMessages(suppressWarnings(if (rk == 0L)
      ilm_model(y ~ x1 + grp, data = dd, family = "multinomial", verbose = FALSE)
      else ilm_model(y ~ x1 + grp + (1 | g), data = dd, family = "multinomial",
                     verbose = FALSE, re_struct = list(g = list(type = "rr", rank = rk))))),
      error = function(e) NULL)
    if (!is.null(fr) && isTRUE(fr$sdr$pdHess))
      se_red <- se(fr$sdr$cov.fixed[names(fr$opt$par) == "beta",
                                    names(fr$opt$par) == "beta", drop = FALSE])
  } else if (length(tc$flagged) && !mn) {
    ## slope at zero: the reduced model is the random intercept alone
    Sd <- f$Sigma_d$g; sv <- sqrt(pmax(diag(Sd), 0))
    if (sv[2] / sv[1] < 1e-3) {
      fr <- tryCatch(suppressMessages(suppressWarnings(
        ilm_model(y ~ x + (1 | g), data = dd, family = "binomial", verbose = FALSE))),
        error = function(e) NULL)
      if (!is.null(fr) && isTRUE(fr$sdr$pdHess))
        se_red <- se(fr$sdr$cov.fixed[names(fr$opt$par) == "beta",
                                      names(fr$opt$par) == "beta", drop = FALSE])
    }
  }
  conv <- f$opt$convergence == 0L
  list(regime = r, i = i, ok = TRUE, conv = conv, b = b,
       flagged = length(tc$flagged) > 0L, how = tc$how,
       se_B = se(tc$B), se_A = se(tc$A), se_C = se(tc$C), se_C3 = se(tc$C3),
       se_red = se_red, ev = tc$ev, nflat3 = tc$nflat3,
       usable_B = conv && f$hessian_how %in% c("tmb", "recomputed", "boundary") &&
         all(is.finite(se(tc$B))) && all(se(tc$B) > 0),
       nflat = tc$nflat,
       grad = max(abs(as.numeric(f$obj$gr(f$opt$par)))))
}

jobs <- do.call(rbind, lapply(regimes, function(r)
  data.frame(regime = r, i = seq_len(NREP), stringsAsFactors = FALSE)))
cl <- parallel::makeCluster(NCORE)
parallel::clusterExport(cl, c("pkg", "mn_regimes", "sl_regimes", "beta_sl", "truth_B",
                              "to_baseline", "cluster_sizes", "draw_re", "gen_mn", "gen_sl",
                              "J", "C", "P4", "NCL", "PER", "NTOT", "pd", "three_cov"))
invisible(parallel::clusterEvalQ(cl, suppressMessages(devtools::load_all(pkg, quiet = TRUE))))
t0 <- Sys.time()
res <- parallel::parLapply(cl, split(jobs, seq_len(nrow(jobs))), one)
parallel::stopCluster(cl)
saveRDS(res, file.path(sp, "boundary_se.rds"))
cat("elapsed", format(Sys.time() - t0), "\n")

## ---- what it found ----------------------------------------------------------
truth <- function(r) if (r %in% c("slope_zero", "slope_small", "slope_cor1")) c(-0.3, 0.5) else
  as.vector(to_baseline(truth_B(r %in% c("rare", "combined"))))
regs <- unique(vapply(res, `[[`, "", "regime"))
z <- qnorm(0.975)

cat("fits that errored:", sum(!vapply(res, `[[`, TRUE, "ok")), "of", length(res), "\n\n")
rows <- list()
for (r in regs) {
  rr <- Filter(function(x) x$regime == r && isTRUE(x$ok), res)
  tv <- truth(r)
  fl <- vapply(rr, `[[`, TRUE, "flagged")
  how <- vapply(rr, `[[`, "", "how")
  ## usable under each option: B as the package decides; A and C where their
  ## covariance exists (unflagged fits use B's for all three)
  use <- function(o) vapply(rr, function(x) {
    s <- x[[paste0("se_", o)]]
    isTRUE(x$conv) && all(is.finite(s)) && all(s > 0) &&
      (o != "B" || isTRUE(x$usable_B)) && (x$flagged || isTRUE(x$usable_B))
  }, TRUE)
  cov_of <- function(o, keep) {
    ok <- which(keep)
    if (!length(ok)) return(c(NA, 0))
    hit <- unlist(lapply(rr[ok], function(x) abs(x$b - tv) <= z * x[[paste0("se_", o)]]))
    c(mean(hit), length(ok))
  }
  for (o in c("B", "A", "C", "C3")) {
    u <- use(o)
    ## all usable fits, and the flagged (boundary) ones alone
    a <- cov_of(o, u); f <- cov_of(o, u & fl)
    ## SE calibration over all usable fits: mean SE / SD of the estimates
    ok <- which(u)
    Bm <- do.call(rbind, lapply(rr[ok], `[[`, "b"))
    Sm <- do.call(rbind, lapply(rr[ok], function(x) x[[paste0("se_", o)]]))
    rows[[length(rows) + 1L]] <- data.frame(
      regime = r, opt = o, n = length(rr), flagged = sum(fl),
      tmb = sum(fl & how == "tmb"), held = sum(fl & how == "boundary"),
      usable = sum(u), cover_all = round(a[1], 4), usable_flag = sum(u & fl),
      cover_flag = round(f[1], 4), se_sd = round(mean(colMeans(Sm) / apply(Bm, 2, sd)), 3))
  }
}
tab <- do.call(rbind, rows)
print(tab, row.names = FALSE)

cat("\n---- among flagged fits usable under A, B and C alike: SE ratios to C\n")
for (r in regs) {
  rr <- Filter(function(x) x$regime == r && isTRUE(x$ok) && x$flagged, res)
  if (!length(rr)) next
  all3 <- Filter(function(x) all(is.finite(c(x$se_A, x$se_B, x$se_C))) && isTRUE(x$usable_B), rr)
  if (!length(all3)) next
  rA <- unlist(lapply(all3, function(x) x$se_A / x$se_C))
  rB <- unlist(lapply(all3, function(x) x$se_B / x$se_C))
  rC3 <- unlist(lapply(Filter(function(x) all(is.finite(x$se_C3)), all3),
                       function(x) x$se_C3 / x$se_C))
  red <- Filter(function(x) all(is.finite(x$se_red)), all3)
  rR <- unlist(lapply(red, function(x) x$se_C / x$se_red))
  howv <- vapply(all3, `[[`, "", "how")
  rBt <- unlist(lapply(all3[howv == "tmb"], function(x) x$se_B / x$se_C))
  rBh <- unlist(lapply(all3[howv == "boundary"], function(x) x$se_B / x$se_C))
  q <- function(v) if (length(v)) sprintf("%.3f [%.3f, %.3f]", median(v), min(v), max(v)) else "-"
  cat(sprintf("%-11s n=%3d  A/C %s  B/C %s (tmb %s | held %s)  C3/C %s  C/reduced-refit %s (n=%d)\n",
              r, length(all3), q(rA), q(rB), q(rBt), q(rBh), q(rC3), q(rR), length(red)))
}

cat("\n---- how flat is flat: block eigenvalues of flagged fits\n")
for (r in regs) {
  rr <- Filter(function(x) x$regime == r && isTRUE(x$ok) && x$flagged && length(x$ev), res)
  if (!length(rr)) next
  flat_max <- unlist(lapply(rr, function(x) { e <- x$ev; m <- max(e[1], 1e-8)
    f <- e[e < 1e-4 * m]; if (length(f)) max(f / m) else NA }))
  id_min <- unlist(lapply(rr, function(x) { e <- x$ev; m <- max(e[1], 1e-8)
    k <- e[e >= 1e-4 * m]; if (length(k)) min(k / m) else NA }))
  gap <- unlist(lapply(rr, function(x) { e <- sort(x$ev / max(x$ev[1], 1e-8)); d <- diff(log10(pmax(e, 1e-12))); max(d) }))
  n4 <- vapply(rr, function(x) x$nflat, 1L); n3 <- vapply(rr, function(x) x$nflat3, 1L)
  cat(sprintf("%-11s n=%3d  largest flat/max %.1e  smallest kept/max %.1e  flat count disagrees (1e-4 vs 1e-3): %d\n",
              r, length(rr), max(flat_max, na.rm = TRUE), min(id_min, na.rm = TRUE), sum(n4 != n3)))
}

cat("\n---- paired coverage, flagged fits usable under both: C against A\n")
for (r in regs) {
  rr <- Filter(function(x) x$regime == r && isTRUE(x$ok) && x$flagged &&
                 all(is.finite(c(x$se_A, x$se_C))), res)
  if (!length(rr)) next
  tv <- truth(r)
  cA <- unlist(lapply(rr, function(x) abs(x$b - tv) <= z * x$se_A))
  cC <- unlist(lapply(rr, function(x) abs(x$b - tv) <= z * x$se_C))
  cat(sprintf("%-11s fits %3d  intervals %4d  A %.4f  C %.4f   C covers where A misses: %d; A where C misses: %d\n",
              r, length(rr), length(cA), mean(cA), mean(cC), sum(cC & !cA), sum(cA & !cC)))
}
