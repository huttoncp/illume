## The typical group's interval in ilm_survival(), for each family and kind of
## fit: how far it is from the delta method on the linear predictor worked out
## by hand from the whole dense inverse of the joint precision (or vcov() when no
## smooth is in the curve), and whether nsim and seed change it.
##
## The reference builds the interval as S evaluated at eta -/+ crit * se, with
## se from a' V a, `a` being a row's fixed design (and the baseline spline's,
## at each time, for a flexible model) and its smooths' design. A build that
## forms the interval that way should agree to rounding; a build that took it
## from draws, or from the complementary log-log value shifted by crit * se /
## scale, will not.
##
## Run from the repository root, with the build to be measured installed:
##   Rscript dev/measurements/typical-interval-delta/measure.R <label for the build> [<file.rds to save the intervals in>]
suppressPackageStartupMessages(library(illume))
args <- commandArgs(trailingOnly = TRUE)
cat("build:", if (length(args)) args[1] else "(unlabelled)",
    "| illume", as.character(packageVersion("illume")),
    "| R", paste0(R.version$major, ".", R.version$minor),
    "| date", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"), "\n")
q <- function(e) suppressMessages(suppressWarnings(e))

set.seed(41)
n <- 250
d <- data.frame(x = runif(n), z = rnorm(n), h = factor(sample(letters[1:8], n, TRUE)))
hl <- rnorm(8, 0, 0.4)
eta <- 1 + 0.5 * d$z + 1.2 * sin(2 * pi * d$x) + hl[as.integer(d$h)]
tt <- exp(eta + 0.5 * log(rexp(n)))
ct <- rexp(n, 1 / (2 * median(tt)))
d$time <- pmin(tt, ct); d$event <- as.integer(tt <= ct)
nd <- data.frame(x = c(0.25, 0.5, 0.75), z = c(-1, 0, 1), h = factor("a", levels = levels(d$h)))
tg <- c(0.5, 2, 6)
fit <- function(form, family, reml = FALSE)
  q(ilm_model(form, data = d, family = family, censor = ilm_surv(d$time, d$event),
              reml = reml, verbose = FALSE))

aft_S <- function(nm, z)
  switch(nm, weibull = exp(-exp(z)), lognormal = pnorm(-z), loglogistic = 1 / (1 + exp(z)))
ref <- function(f) {
  crit <- qnorm(0.975)
  ndx <- illume:::ilm_newX(f, nd, fs_fill = TRUE)
  basis <- vapply(f$re, function(e) identical(e$kind, "basis"), TRUE)
  used <- which(basis & !illume:::ilm_fs_terms(f))
  rp <- f$rp; p <- ncol(f$X)
  keep <- if (is.null(rp)) seq_len(p) else setdiff(seq_len(p), rp$cols)
  rn <- names(f$obj$env$par); pb <- which(rn == "beta"); pv <- which(rn == "bvec")
  V <- if (length(used)) solve(as.matrix(illume:::ilm_joint_prec(f))) else as.matrix(vcov(f))
  se_of <- function(ab, sm) {
    if (length(used)) {
      a <- numeric(nrow(V)); a[pb] <- ab
      for (k in used) a[pv[f$b_idx[[k]]]] <- sm[[as.character(k)]]
      sqrt(as.numeric(t(a) %*% V %*% a))
    } else sqrt(as.numeric(t(ab) %*% V %*% ab))
  }
  do.call(rbind, lapply(seq_len(nrow(nd)), function(i) {
    sm <- lapply(used, function(k) as.matrix(illume:::ilm_basis_new(f, ndx, k))[i, ])
    names(sm) <- as.character(used)
    off <- sum(ndx$X[i, ] * f$beta[keep, 1L]) +
      sum(vapply(used, function(k) sum(sm[[as.character(k)]] * as.vector(illume:::ilm_Bhat_term(f, k))), 0))
    if (is.null(rp)) {
      se <- se_of(ndx$X[i, ], sm); sc <- unname(f$dispersion[[1]]); nm <- f$family$name
      data.frame(lower = aft_S(nm, (log(tg) - (off - crit * se)) / sc),
                 upper = aft_S(nm, (log(tg) - (off + crit * se)) / sc))
    } else {
      Bs <- illume:::ilm_rcs(log(tg), rp$knots)
      e <- as.vector(Bs %*% f$beta[rp$cols, 1L]) + off
      se <- vapply(seq_along(tg), function(j) {
        ab <- numeric(p); ab[rp$cols] <- Bs[j, ]; ab[keep] <- ndx$X[i, ]; se_of(ab, sm)
      }, 0)
      data.frame(lower = f$family$surv(e + crit * se), upper = f$family$surv(e - crit * se))
    }
  }))
}

forms <- list(plain = time ~ z, re = time ~ z + (1 | h), smooth = time ~ s(x, k = 8) + z + (1 | h))
saved <- list()
cat(sprintf("\n%-12s %-7s %-9s %-34s %s\n", "family", "fit", "reml", "max |limit - by hand| (lower, upper)",
            "max |change| with nsim 10/seed 1 -> nsim 2000/seed 7"))
for (fam in c("weibull", "lognormal", "loglogistic", "rp")) for (k in names(forms)) {
  reml <- k == "smooth"
  f <- fit(forms[[k]], fam, reml)
  a <- q(ilm_survival(f, newdata = nd, times = tg, groups = "typical", nsim = 10L, seed = 1L))
  b <- q(ilm_survival(f, newdata = nd, times = tg, groups = "typical", nsim = 2000L, seed = 7L))
  r <- ref(f)
  ## what a user sees with the defaults (nsim = 200, seed = 1), kept for compare.R
  saved[[paste(fam, k)]] <- q(ilm_survival(f, newdata = nd, times = tg, groups = "typical"))
  cat(sprintf("%-12s %-7s %-9s %-34s %s\n", fam, k, reml,
              paste(format(c(max(abs(a$lower - r$lower)), max(abs(a$upper - r$upper))), digits = 3), collapse = ", "),
              format(max(abs(a$lower - b$lower), abs(a$upper - b$upper)), digits = 3)))
}
if (length(args) > 1L) { saveRDS(saved, args[2]); cat("\nintervals saved to", args[2], "\n") }
