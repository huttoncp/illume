## ---------------------------------------------------------------------------
## Type I error and power of illume's fixed-effect tests (ilm_anova).
##
## The null cell (delta = 0) is the important one.  Power is only meaningful if
## the test holds its nominal size, so the study reports the rejection rate at
## delta = 0 first and treats the rest of the curve as conditional on it.
##
## The x1 row of the true coefficient matrix is zeroed and then set to delta in
## the FIRST category dimension only, so under the null x1 has no effect at all,
## and under the alternative the joint test has to find an effect that lives in
## one dimension out of C.
##
## Usage: Rscript power_study.R <nrep> <ncore> <outdir> [deltas] [cells]
## ---------------------------------------------------------------------------

suppressPackageStartupMessages(library(parallel))

args   <- commandArgs(trailingOnly = TRUE)
NREP   <- if (length(args) >= 1) as.integer(args[1]) else 10L
NCORE  <- if (length(args) >= 2) as.integer(args[2]) else 6L
OUTDIR <- if (length(args) >= 3) args[3] else "."
## NB: `else` on a new line is a parse error at top level in R, so this stays
## braced rather than split across lines.
DELTAS <- if (length(args) >= 4) {
  as.numeric(strsplit(args[4], ",")[[1]])
} else c(0, 0.1, 0.2, 0.3, 0.5)
ONLY   <- if (length(args) >= 5) strsplit(args[5], ",")[[1]] else NULL

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
ALPHA <- 0.05

cells <- list(
  list(name = "gauss_fixed", family = "gaussian",    J = NA, re = FALSE,
       N = 400),
  list(name = "gauss_mm",    family = "gaussian",    J = NA, re = TRUE,
       ncl = 30, per = 20, sd_re = 0.6),
  list(name = "binom_mm",    family = "binomial",    J = NA, re = TRUE,
       ncl = 30, per = 20, sd_re = 0.6),
  list(name = "pois_mm",     family = "poisson",     J = NA, re = TRUE,
       ncl = 30, per = 20, sd_re = 0.6),
  list(name = "mn_J3_rich",  family = "multinomial", J = 3,  re = TRUE,
       ncl = 30, per = 20, sd_re = 0.6),
  list(name = "mn_J3_thin",  family = "multinomial", J = 3,  re = TRUE,
       ncl = 60, per = 5,  sd_re = 0.6),

  ## ---- the model types added after the first run ---------------------------
  ## Size is the question for every one of these. A new family whose test
  ## rejects at 8% rather than 5% would be invisible to a coverage study of
  ## point estimates, and invisible to a test suite that only asks whether the
  ## fit ran.
  list(name = "nbinom_mm",  kind = "glm", family = "nbinom", re = TRUE,
       ncl = 30, per = 20, sd_re = 0.6, k = 2),
  list(name = "gauss_ar1",  kind = "ar",  family = "gaussian", re = TRUE,
       ncl = 40, Tt = 10, rep_per = 3, sd_re = 0.6, rho = 0.6),
  list(name = "gauss_car1", kind = "car", family = "gaussian", re = TRUE,
       ncl = 60, per = 6, tmax = 30, sd_re = 0.6, rho = 0.85, sd_e = 0.6),
  list(name = "tobit_ceil", kind = "censor", family = "gaussian", re = FALSE,
       N = 600, q = 0.8),
  list(name = "aft_weib",   kind = "aft", family = "weibull", re = FALSE,
       N = 600, scale = 0.7, cr = 0.35),
  list(name = "rp_flex",    kind = "rp",  family = "rp", re = FALSE,
       N = 800, rp_df = 5, cut = 3, h1 = 0.15, h2 = 0.9),
  list(name = "disp_group", kind = "disp", family = "gaussian", re = FALSE,
       N = 600)
)
names(cells) <- vapply(cells, function(z) z$name, "")
if (!is.null(ONLY)) cells <- cells[ONLY]

truth_of <- function(cell, delta) {
  C <- if (identical(cell$family, "multinomial")) cell$J - 1L else 1L
  p <- 4L
  set.seed(99L)
  Bt <- matrix(round(runif(p * C, -0.7, 0.7), 3), p, C)
  Bt[2, ] <- 0                 # x1 contributes nothing under the null
  Bt[2, 1] <- delta            # effect placed in one category dimension
  Bt
}

## The model types added later, each generating a response from the same linear
## predictor so that delta means the same thing throughout: the coefficient on
## x1, zero under the null.
gen_extra <- function(cell, delta, seed) {
  Bt <- truth_of(cell, delta)      # sets its own seed; must come first
  set.seed(seed)
  cov1 <- function(N) data.frame(x1 = rnorm(N),
                                 grp = factor(sample(c("a", "b", "c"), N, TRUE)))

  if (cell$kind == "ar") {
    ncl <- cell$ncl; Tt <- cell$Tt; rp <- cell$rep_per; N <- ncl * Tt * rp
    dd <- data.frame(g = factor(rep(seq_len(ncl), each = Tt * rp)),
                     time = rep(rep(seq_len(Tt), each = rp), ncl))
    dd <- cbind(dd, cov1(N))
    eta <- as.numeric(model.matrix(~ x1 + grp, dd) %*% Bt)
    Ba <- numeric(ncl * Tt)
    for (g in seq_len(ncl)) {
      i0 <- (g - 1L) * Tt
      Ba[i0 + 1L] <- rnorm(1, 0, cell$sd_re)
      for (t in 2:Tt)
        Ba[i0 + t] <- cell$rho * Ba[i0 + t - 1L] +
          rnorm(1, 0, cell$sd_re * sqrt(1 - cell$rho^2))
    }
    dd$y <- eta + Ba[(as.integer(dd$g) - 1L) * Tt + dd$time] + rnorm(N, 0, 1)
    attr(dd, "ar") <- illume::ilm_ar1(dd$time, dd$g, verbose = FALSE)
    return(dd)
  }
  if (cell$kind == "car") {
    ncl <- cell$ncl; per <- cell$per; N <- ncl * per
    tl <- lapply(seq_len(ncl), function(i) sort(sample.int(cell$tmax, per)))
    dd <- data.frame(g = factor(rep(seq_len(ncl), each = per)),
                     time = unlist(tl))
    dd <- cbind(dd, cov1(N))
    eta <- as.numeric(model.matrix(~ x1 + grp, dd) %*% Bt)
    u <- unlist(lapply(tl, function(tt) {
      z <- numeric(length(tt)); z[1] <- rnorm(1)
      for (k in seq_along(tt)[-1]) {
        ph <- cell$rho^(tt[k] - tt[k - 1])
        z[k] <- ph * z[k - 1] + rnorm(1, 0, sqrt(1 - ph^2))
      }
      z
    }))
    dd$y <- eta + cell$sd_re * u + rnorm(N, 0, cell$sd_e)
    attr(dd, "ar") <- illume::ilm_car1(dd$time, dd$g, verbose = FALSE)
    return(dd)
  }
  if (cell$kind == "censor") {
    dd <- cov1(cell$N)
    ys <- as.numeric(model.matrix(~ x1 + grp, dd) %*% Bt) + rnorm(cell$N, 0, 1)
    up <- unname(quantile(ys, cell$q)); dd$y <- pmin(ys, up)
    attr(dd, "censor") <- illume::ilm_censor(dd$y, upper = up)
    return(dd)
  }
  if (cell$kind == "disp") {
    dd <- cov1(cell$N)
    dd$s <- factor(sample(c("lo", "hi"), cell$N, TRUE))
    dd$y <- as.numeric(model.matrix(~ x1 + grp, dd) %*% Bt) +
      rnorm(cell$N, 0, ifelse(dd$s == "hi", 2, 0.5))
    return(dd)
  }
  if (cell$kind == "aft") {
    dd <- cov1(cell$N)
    eta <- as.numeric(model.matrix(~ x1 + grp, dd) %*% Bt)
    tt <- exp(eta + cell$scale * log(rexp(cell$N)))
    q <- unname(quantile(tt, 1 - cell$cr))
    ct <- pmin(q, rexp(cell$N, rate = 1 / (3 * q)))
    dd$y <- pmin(tt, ct)
    attr(dd, "censor") <- illume::ilm_surv(dd$y, as.integer(tt <= ct))
    return(dd)
  }
  if (cell$kind == "rp") {
    dd <- cov1(cell$N)
    lp <- as.numeric(model.matrix(~ x1 + grp, dd) %*% Bt) - Bt[1L]
    e <- -log(runif(cell$N)); H1 <- cell$h1 * cell$cut
    t1 <- e / (cell$h1 * exp(lp))
    tt <- pmax(ifelse(t1 <= cell$cut, t1,
                      cell$cut + (e - H1 * exp(lp)) / (cell$h2 * exp(lp))), 1e-4)
    ct <- rexp(cell$N, rate = 1 / (2 * median(tt)))
    dd$y <- pmin(tt, ct)
    attr(dd, "censor") <- illume::ilm_surv(dd$y, as.integer(tt <= ct))
    return(dd)
  }
  stop("unknown cell kind: ", cell$kind)
}

gen <- function(cell, delta, seed) {
  if (!is.null(cell$kind) && cell$kind != "glm")
    return(gen_extra(cell, delta, seed))
  C  <- if (identical(cell$family, "multinomial")) cell$J - 1L else 1L
  Bt <- truth_of(cell, delta)      # sets its own seed; must come first
  set.seed(seed)
  if (isTRUE(cell$re)) {
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
    nbinom   = rnbinom(nrow(dd), size = cell$k,
                       mu = exp(pmin(as.numeric(eta), 5))),
    multinomial = {
      Tc <- contr.sum(cell$J)
      P  <- exp(eta %*% t(Tc)); P <- P / rowSums(P)
      labs <- paste0("c", seq_len(cell$J))
      factor(labs[apply(P, 1, function(pr) sample.int(cell$J, 1L, prob = pr))],
             levels = labs)
    })
  dd
}

## p-value for the x1 term, from whichever column the table carries
pval_x1 <- function(a) {
  if (is.null(a) || !nrow(a) || !("x1" %in% rownames(a))) return(NA_real_)
  cn <- intersect(c("Pr(>Chisq)", "Pr(>F)"), names(a))
  if (!length(cn)) return(NA_real_)
  a[["x1", cn[1]]]
}

run_rep <- function(i, cell, delta, do_lrt) {
  dd <- gen(cell, delta, seed = 20000L + i)
  kind <- if (is.null(cell$kind)) "glm" else cell$kind
  fm <- if (isTRUE(cell$re)) y ~ x1 + grp + (1 | g) else y ~ x1 + grp
  args <- list(formula = fm, data = dd, family = cell$family,
               ar = attr(dd, "ar"), censor = attr(dd, "censor"),
               verbose = FALSE)
  if (kind == "disp") args$dispformula <- ~ s
  if (kind == "rp")   args$rp_df <- cell$rp_df
  f <- tryCatch(suppressWarnings(do.call(illume::ilm_model, args)),
       error = function(e) NULL)
  if (is.null(f)) return(list(ok = FALSE, wald = NA_real_, lrt = NA_real_))
  ## Judged on the first-order condition, not on nlminb's stopping code: its
  ## "false convergence (8)" means it could not verify a descent direction, and
  ## on flexible parametric fits a third of the replicates reported it while
  ## sitting at a small gradient and recovering the right coefficients.
  gst <- f$checks$status[f$checks$check == "gradient"]
  ok <- (!length(gst) || gst != "FAIL") && isTRUE(f$sdr$pdHess)
  w <- tryCatch(suppressWarnings(pval_x1(illume::ilm_anova(f, type = 3))),
                error = function(e) NA_real_)
  l <- if (do_lrt)
         tryCatch(suppressWarnings(
           pval_x1(illume::ilm_anova(f, type = 3, test = "LRT"))),
           error = function(e) NA_real_)
       else NA_real_
  list(ok = ok, wald = w, lrt = l)
}

cl <- makePSOCKcluster(NCORE)
on.exit(stopCluster(cl), add = TRUE)
invisible(clusterEvalQ(cl, suppressPackageStartupMessages(library(illume))))
clusterExport(cl, c("gen", "gen_extra", "truth_of", "run_rep", "pval_x1"),
              envir = environment())

rows <- list()
for (cell in cells) for (delta in DELTAS) {
  ## the LRT refits a reduced model for every replicate, so it is run only
  ## under the null, where size is the quantity of interest
  do_lrt <- isTRUE(all.equal(delta, 0))
  t0 <- Sys.time()
  clusterExport(cl, c("cell", "delta", "do_lrt"), envir = environment())
  res <- parLapply(cl, seq_len(NREP),
                   function(i) run_rep(i, cell, delta, do_lrt))
  el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  ok <- vapply(res, function(z) isTRUE(z$ok), TRUE)
  w  <- vapply(res, function(z) z$wald, 1)
  l  <- vapply(res, function(z) z$lrt,  1)
  usew <- ok & is.finite(w); usel <- ok & is.finite(l)

  r <- data.frame(
    cell = cell$name, family = cell$family,
    J = if (is.null(cell$J) || is.na(cell$J)) NA_integer_ else cell$J,
    delta = delta, n_attempt = NREP,
    n_wald = sum(usew), rej_wald = mean(w[usew] < ALPHA),
    n_lrt = sum(usel),
    rej_lrt = if (any(usel)) mean(l[usel] < ALPHA) else NA_real_,
    secs = el, stringsAsFactors = FALSE)
  r$mc_se_wald <- sqrt(r$rej_wald * (1 - r$rej_wald) / max(r$n_wald, 1))
  rows[[length(rows) + 1L]] <- r

  saveRDS(list(cell = cell, delta = delta, wald = w, lrt = l, ok = ok),
          file.path(OUTDIR, sprintf("pow_%s_d%s.rds", cell$name, delta)))
  cat(sprintf("%-12s d=%.2f  n=%4d/%4d  rej(Wald) %.4f%s  %.0fs\n",
              cell$name, delta, sum(usew), NREP, r$rej_wald,
              if (do_lrt) sprintf("  rej(LRT) %.4f", r$rej_lrt) else "", el))
  flush.console()
  write.csv(do.call(rbind, rows), file.path(OUTDIR, "power_summary.csv"),
            row.names = FALSE)
}
cat("done\n")
