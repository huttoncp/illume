## ---------------------------------------------------------------------------
## illume against brms::brm(family = categorical()).
##
## These do not answer the same question -- one reports a maximum likelihood
## estimate with a Wald interval, the other a posterior mean with a credible
## interval -- so this is not a contest.  With flat-ish priors and a reasonable
## sample the two should land in nearly the same place, and that agreement is
## the point: brms is a mature, independently written implementation of the
## same likelihood, so it is the strongest external check available for the
## multinomial mixed model.
##
## Timing is reported with compilation separated from sampling, because the
## compile cost is paid once per model form and would otherwise dominate.
##
## brms categorical() uses baseline-category coding against the first level,
## so illume's sum-to-zero coefficients are converted before comparison.
##
## Usage: Rscript brms_compare.R <nrep> <outdir>
## ---------------------------------------------------------------------------

args   <- commandArgs(trailingOnly = TRUE)
NREP   <- if (length(args) >= 1) as.integer(args[1]) else 3L
OUTDIR <- if (length(args) >= 2) args[2] else "."
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({ library(illume); library(brms) })

J <- 3L; C <- J - 1L; P4 <- 4L
NCL <- 30L; PER <- 25L                      # N = 750

truth_B <- function() { set.seed(99L); matrix(round(runif(P4 * C, -0.7, 0.7), 3), P4, C) }

to_baseline <- function(B) {
  full <- B %*% t(contr.sum(J))
  full[, -1, drop = FALSE] - full[, 1]
}

gen <- function(seed) {
  Bt <- truth_B(); set.seed(seed)
  n <- NCL * PER
  dd <- data.frame(g = factor(rep(seq_len(NCL), each = PER)))
  dd$x1  <- rnorm(n)
  dd$grp <- factor(sample(c("a", "b", "c"), n, TRUE))
  X <- model.matrix(~ x1 + grp, dd)
  b <- matrix(rnorm(NCL * C, 0, 0.7), NCL, C)
  eta <- X %*% Bt + b[as.integer(dd$g), , drop = FALSE]
  Pm <- exp(eta %*% t(contr.sum(J))); Pm <- Pm / rowSums(Pm)
  dd$y <- factor(paste0("c", apply(Pm, 1, function(p) sample.int(J, 1L, prob = p))),
                 levels = paste0("c", seq_len(J)))
  dd
}

## brms names fixed effects b_muc2_x1 etc.
brms_names <- function() {
  xn <- c("Intercept", "x1", "grpb", "grpc")
  as.vector(vapply(2:J, function(j) paste0("b_mu", "c", j, "_", xn), character(P4)))
}

tv <- as.vector(to_baseline(truth_B()))
rows <- list(); fit0 <- NULL

for (r in seq_len(NREP)) {
  dd <- gen(6000L + r)

  t0 <- Sys.time()
  f <- ilm_model(y ~ x1 + grp + (1 | g), data = dd, family = "multinomial",
                 verbose = FALSE)
  t_ill <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  b_ill <- as.vector(to_baseline(matrix(coef(f), P4, C)))
  Tc <- contr.sum(J)
  M  <- Tc[-1, , drop = FALSE] - matrix(Tc[1, ], J - 1L, C, byrow = TRUE)
  L  <- kronecker(M, diag(P4))
  s_ill <- sqrt(pmax(diag(L %*% suppressWarnings(vcov(f)) %*% t(L)), 0))

  t0 <- Sys.time()
  if (is.null(fit0)) {
    bf <- suppressWarnings(suppressMessages(
      brm(y ~ x1 + grp + (1 | g), data = dd, family = categorical(),
          chains = 4, iter = 2000, warmup = 1000, cores = 4,
          refresh = 0, silent = 2, seed = 1L)))
    fit0 <- bf
    compiled <- TRUE
  } else {
    bf <- suppressWarnings(suppressMessages(
      update(fit0, newdata = dd, recompile = FALSE,
             chains = 4, iter = 2000, warmup = 1000, cores = 4,
             refresh = 0, silent = 2, seed = 1L)))
    compiled <- FALSE
  }
  t_brm <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  ps <- posterior::as_draws_matrix(bf)
  bn <- brms_names()
  have <- bn %in% colnames(ps)
  b_brm <- rep(NA_real_, length(bn)); s_brm <- rep(NA_real_, length(bn))
  if (all(have)) {
    b_brm <- apply(ps[, bn, drop = FALSE], 2, mean)
    s_brm <- apply(ps[, bn, drop = FALSE], 2, sd)
  }
  rh <- suppressWarnings(max(brms::rhat(bf), na.rm = TRUE))

  rows[[r]] <- data.frame(
    rep = r, N = nrow(dd), coef = bn, truth = tv,
    illume = b_ill, brms = unname(b_brm),
    se_illume = s_ill, sd_brms = unname(s_brm),
    t_illume = t_ill, t_brms = t_brm, brms_compiled = compiled,
    max_rhat = rh, stringsAsFactors = FALSE)

  cat(sprintf("rep %d  max|illume-brms| %.4f  mean|est-truth| ill %.4f brms %.4f  |  SE ratio %.3f  |  %.2fs vs %.1fs%s  max Rhat %.3f\n",
      r, max(abs(b_ill - b_brm)), mean(abs(b_ill - tv)), mean(abs(b_brm - tv)),
      mean(s_ill / s_brm), t_ill, t_brm,
      if (compiled) " (incl. compile)" else "", rh))
  flush.console()
  write.csv(do.call(rbind, rows), file.path(OUTDIR, "brms_compare.csv"),
            row.names = FALSE)
}
cat("done\n")
