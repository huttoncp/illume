## ilm_survival()'s "typical" interval on a model with a smooth beside a
## random intercept, against reference intervals. For S(3) at x = 0.25, 0.5
## and 0.75, with the random intercept at zero, it prints:
##   ilm_survival()   its point and interval, with nsim = 200 and nsim = 4000
##   ref, joint       percentiles of S(3) over 4000 joint draws of the fixed
##                    effects and the smooth's penalised coefficients, from
##                    ilm_draws(blocks = c("beta", "bvec"))
##   ref, beta only   the same with the penalised coefficients held at their
##                    estimates
##   delta, joint     the delta method on the linear predictor, from the
##                    joint covariance of the fixed effects and the smooth's
##                    penalised coefficients: the inverse of the fit's joint
##                    precision, restricted to them (flexible-parametric fits)
## for a Weibull and a flexible-parametric fit, each by REML and by ML; then
## the flexible-parametric interval without a smooth. Uses illume's exported
## functions, and its own internals ilm_newX(), ilm_basis_new(),
## ilm_Bhat_term(), ilm_rcs(), ilm_joint_prec() and ilm_full_par() for the
## by-hand references. Simulated data.
##
## Run from the repository root, with the build to be measured installed:
##   Rscript dev/measurements/survival-smooth/interval.R <label for the build>
suppressPackageStartupMessages(library(illume))
lab <- commandArgs(trailingOnly = TRUE)
cat("build:", if (length(lab)) lab[1] else "(unlabelled)",
    "| illume", as.character(packageVersion("illume")),
    "| mgcv", as.character(packageVersion("mgcv")),
    "| R", paste0(R.version$major, ".", R.version$minor),
    "| date", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"), "\n")
q <- function(e) suppressMessages(suppressWarnings(e))
set.seed(31)
n <- 400
d <- data.frame(x = runif(n), h = factor(sample(letters[1:8], n, TRUE)))
hl <- rnorm(8, 0, 0.4)
eta <- 1 + 1.2 * sin(2 * pi * d$x) + hl[as.integer(d$h)]
tt <- exp(eta + 0.5 * log(rexp(n)))
ct <- rexp(n, 1 / (2 * median(tt)))
d$time <- pmin(tt, ct); d$event <- as.integer(tt <= ct)
nd <- data.frame(x = c(0.25, 0.5, 0.75), h = factor("a", levels = levels(d$h)))
tm <- 3
fit <- function(form, family, reml) q(ilm_model(form, data = d, family = family, reml = reml,
                                               censor = ilm_surv(d$time, d$event), verbose = FALSE))
## the linear predictor at nd and t = tm for given coefficients, the random
## intercept at zero; and its gradient in (beta, bvec) as a design row
parts <- function(f) {
  ndx <- illume:::ilm_newX(f, nd)
  keep <- if (is.null(f$rp)) seq_len(ncol(f$X)) else setdiff(seq_len(ncol(f$X)), f$rp$cols)
  ks <- which(vapply(f$re, function(r) identical(r$kind, "basis"), TRUE))
  Bs <- if (!is.null(f$rp)) illume:::ilm_rcs(log(tm), f$rp$knots)
  list(ndx = ndx, keep = keep, ks = ks, Bs = Bs)
}
S_of <- function(f, P, beta, bvec) {
  e <- as.vector(P$ndx$X %*% beta[P$keep])
  for (k in P$ks) e <- e + as.vector(illume:::ilm_basis_new(f, P$ndx, k) %*% illume:::ilm_Bhat_term(f, k, bvec))
  if (is.null(f$rp)) exp(-exp((log(tm) - e) / unname(f$dispersion[[1]])))
  else f$family$surv(as.vector(P$Bs %*% beta[f$rp$cols]) + e)
}
## the delta method from the joint covariance, for a flexible-parametric fit:
## eta = B(t) gamma + x beta + r u, linear in (beta, u)
delta_joint <- function(f, P) {
  Q <- as.matrix(illume:::ilm_joint_prec(f)); rn <- names(illume:::ilm_full_par(f))
  ib <- which(rn == "beta"); iv <- which(rn == "bvec")
  Sig <- solve(Q)[c(ib, iv), c(ib, iv)]
  nb <- length(ib); bv0 <- illume:::ilm_full_par(f)[iv]
  out <- matrix(NA_real_, nrow(nd), 2)
  for (i in seq_len(nrow(nd))) {
    zb <- numeric(nb); zb[f$rp$cols] <- P$Bs; zb[P$keep] <- P$ndx$X[i, ]
    zv <- numeric(length(iv))
    for (k in P$ks) {
      Xk <- illume:::ilm_basis_new(f, P$ndx, k)[i, , drop = FALSE]
      ## the term's coefficients sit at f$b_idx[[k]] within bvec
      zv[f$b_idx[[k]]] <- zv[f$b_idx[[k]]] + as.vector(Xk)
    }
    z <- c(zb, zv)
    e <- sum(zb * as.vector(f$beta)) + sum(zv * bv0)
    se <- sqrt(sum((z %*% Sig) * z))
    out[i, ] <- f$family$surv(c(e + 1.96 * se, e - 1.96 * se))
  }
  out
}
row_fmt <- function(x, lo, hi) sprintf("[%-9s, %-9s] width %s", signif(lo, 4), signif(hi, 4), signif(hi - lo, 3))
for (cfg in list(list("weibull", TRUE), list("weibull", FALSE), list("rp", TRUE), list("rp", FALSE))) {
  f <- fit(time ~ s(x, k = 8) + (1 | h), cfg[[1]], cfg[[2]])
  cat("\n==", cfg[[1]], "| reml =", cfg[[2]], "| time ~ s(x, k = 8) + (1 | h), S(3), h at zero ==\n")
  P <- parts(f)
  s200 <- q(ilm_survival(f, newdata = nd, times = c(2, tm), groups = "typical", nsim = 200L))
  s4k <- q(ilm_survival(f, newdata = nd, times = c(2, tm), groups = "typical", nsim = 4000L))
  s200 <- s200[s200$time == tm, ]; s4k <- s4k[s4k$time == tm, ]
  dr <- ilm_draws(f, nsim = 4000L, seed = 7L, blocks = c("beta", "bvec"), natural = FALSE)
  w <- rownames(dr$draws)
  Sj <- sapply(seq_len(ncol(dr$draws)), function(j) S_of(f, P, dr$draws[w == "beta", j], dr$draws[w == "bvec", j]))
  Sb <- sapply(seq_len(ncol(dr$draws)), function(j) S_of(f, P, dr$draws[w == "beta", j], NULL))
  qj <- apply(Sj, 1, quantile, c(0.025, 0.975)); qb <- apply(Sb, 1, quantile, c(0.025, 0.975))
  dj <- if (!is.null(f$rp)) delta_joint(f, P)
  for (i in seq_len(nrow(nd))) {
    cat(sprintf("x = %.2f  point %-9s inside (4000): %s\n", nd$x[i], signif(s4k$surv[i], 4),
                s4k$lower[i] <= s4k$surv[i] && s4k$surv[i] <= s4k$upper[i]))
    cat("          ilm_survival(), nsim = 200:  ", row_fmt(0, s200$lower[i], s200$upper[i]), "\n")
    cat("          ilm_survival(), nsim = 4000: ", row_fmt(0, s4k$lower[i], s4k$upper[i]), "\n")
    cat("          ref, joint (4000 draws):     ", row_fmt(0, qj[1, i], qj[2, i]), "\n")
    cat("          ref, beta only (4000 draws): ", row_fmt(0, qb[1, i], qb[2, i]), "\n")
    if (!is.null(dj))
      cat("          delta, joint covariance:     ", row_fmt(0, dj[i, 1], dj[i, 2]), "\n")
  }
}
cat("\n== control: rp, time ~ x + (1 | h), S(3) widths at the three rows ==\n")
for (reml in c(TRUE, FALSE)) {
  f0 <- fit(time ~ x + (1 | h), "rp", reml)
  s0 <- q(ilm_survival(f0, newdata = nd, times = c(2, tm), groups = "typical"))
  s0 <- s0[s0$time == tm, ]
  cat(sprintf("  reml = %s: %s\n", reml, paste(signif(s0$upper - s0$lower, 3), collapse = " ")))
}
