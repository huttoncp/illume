## A dispersion at its unbounded limit: where is the line?
##
## A negative binomial's size k, or a beta's precision phi, can run off to
## infinity: the data show no overdispersion beyond what the model's other
## terms -- a correlation over time, a random intercept -- already carry, and
## the maximum-likelihood fit is the limit model (for the negative binomial,
## the Poisson). The optimiser stops wherever the gradient is small, so the
## estimate lands somewhere large and arbitrary, with a standard error to
## match, and draws of it span hundreds on the log scale.
##
## The rule proposed mirrors the covariance one. A variance is at its boundary
## when its SD is below 1e-3 on the linear predictor's scale. A negative
## binomial is a Poisson with gamma noise whose SD on the log scale is about
## 1 / sqrt(k), and a beta's noise on the logit scale shrinks as
## 1 / sqrt(phi); so the dispersion is at its limit when 1 / sqrt(k) (or
## phi) < 1e-3. This study measures whether that line separates the fits
## whose optimum IS the limit from those with a finite one, and by how much.
##
## Each fit is labelled from outside the rule, twice:
##   limit_ll   the limit model fits as well: for the negative binomial, the
##              Poisson refit's log-likelihood is within 1e-3 of it; for the
##              beta, the objective with log(phi) pushed up by 3 (the other
##              parameters held, the random effects re-optimised) is no worse
##   push_down  the objective still falls, or stays level within 1e-3, when
##              the log dispersion is pushed up by 3 -- a one-sided profile
## and the fitted 1 / sqrt(dispersion) is recorded for each.
##
## Phase 1 (this script, before the hold): the distribution of
## 1 / sqrt(dispersion) in each label, and the gap between them.
## With `verify`, phase 2: the fits the rule holds have the Poisson refit's
## fixed-effect standard errors, and their draws are finite.
##
## Usage: Rscript dispersion_limit.R <nrep> <ncore> <outdir> [verify]
## ---------------------------------------------------------------------------
args  <- commandArgs(trailingOnly = TRUE)
NREP  <- if (length(args) >= 1) as.integer(args[1]) else 50L
NCORE <- if (length(args) >= 2) as.integer(args[2]) else 2L
sp    <- if (length(args) >= 3) args[3] else "."
VERIFY <- length(args) >= 4 && identical(args[4], "verify")
dir.create(sp, showWarnings = FALSE, recursive = TRUE)

cells <- rbind(
  data.frame(family = "nbinom", design = "ar1",  disp = c(Inf, 50, 5, 1)),
  data.frame(family = "nbinom", design = "ri",   disp = c(Inf, 50, 5)),
  data.frame(family = "beta",   design = "ar1",  disp = c(1e4, 200, 20)),
  data.frame(family = "beta",   design = "ri",   disp = c(200, 20)))
cells$cell <- seq_len(nrow(cells))
jobs <- merge(cells, data.frame(rep = seq_len(NREP)))

one <- function(job) {
  suppressMessages(library(illume))
  seed <- 7919L * job$cell + job$rep
  set.seed(seed)
  if (job$design == "ar1") {
    G <- 20; Tn <- 16
    d <- expand.grid(t = seq_len(Tn), g = factor(seq_len(G)))
    lat <- unlist(lapply(seq_len(G), function(i) as.numeric(
      stats::arima.sim(list(ar = 0.7), Tn, sd = 0.8 * sqrt(1 - 0.49)))))
  } else {
    G <- 30; Tn <- 10
    d <- expand.grid(t = seq_len(Tn), g = factor(seq_len(G)))
    lat <- stats::rnorm(G, 0, 0.6)[d$g]
  }
  d$x <- stats::rnorm(nrow(d))
  eta <- 0.5 + 0.3 * d$x + lat
  if (job$family == "nbinom") {
    mu <- exp(eta)
    d$y <- if (is.infinite(job$disp)) stats::rpois(nrow(d), mu)
           else stats::rnbinom(nrow(d), size = job$disp, mu = mu)
  } else {
    mu <- stats::plogis(eta - 0.5)
    d$y <- stats::rbeta(nrow(d), mu * job$disp, (1 - mu) * job$disp)
    d$y <- pmin(pmax(d$y, 1e-6), 1 - 1e-6)
  }
  ar <- if (job$design == "ar1") ilm_ar1(~ t | g) else NULL
  fml <- if (job$design == "ar1") y ~ x + (1 | g) else y ~ x + (1 | g)
  fit <- function(fam) tryCatch(suppressMessages(suppressWarnings(
    ilm_model(fml, data = d, family = fam, ar = ar, verbose = FALSE))),
    error = function(e) NULL)
  f <- fit(job$family)
  out <- data.frame(job, seed = seed, ok = !is.null(f), logdisp = NA_real_,
                    inv_sqrt = NA_real_, se_logdisp = NA_real_,
                    ll = NA_real_, ll_limit = NA_real_, push = NA_real_,
                    held = NA_character_, how = NA_character_,
                    se_ratio_x = NA_real_, draws_finite = NA)
  if (is.null(f)) return(out)
  pn <- names(f$opt$par); il <- which(pn == "logdisp")
  ld <- unname(f$opt$par[il])
  out$logdisp <- ld; out$inv_sqrt <- exp(-ld / 2)
  V <- tryCatch(vcov(f, full = TRUE), error = function(e) NULL)
  if (!is.null(V)) out$se_logdisp <- sqrt(V[il, il])
  out$ll <- as.numeric(logLik(f))
  out$held <- paste(f$hessian_held, collapse = ",")
  out$how <- f$hessian_how
  ## one-sided profile: log dispersion up by 3, the rest held, the random
  ## effects re-optimised by the objective itself; REML-free, so opt$par is
  ## the vector the objective takes
  p2 <- f$opt$par; p2[il] <- p2[il] + 3
  nll0 <- tryCatch(f$obj$fn(f$opt$par), error = function(e) NA_real_)
  nll1 <- tryCatch(f$obj$fn(p2), error = function(e) NA_real_)
  out$push <- nll1 - nll0
  if (job$family == "nbinom") {
    fp <- fit("poisson")
    if (!is.null(fp)) {
      out$ll_limit <- as.numeric(logLik(fp))
      if (VERIFY) {
        s1 <- sqrt(diag(vcov(f)))[["x"]]; s0 <- sqrt(diag(vcov(fp)))[["x"]]
        out$se_ratio_x <- s1 / s0
      }
    }
  }
  if (VERIFY) {
    dr <- tryCatch(ilm_draws(f, nsim = 200, seed = 1, natural = FALSE),
                   error = function(e) NULL)
    if (!is.null(dr))
      out$draws_finite <- all(abs(dr$draws[rownames(dr$draws) == "logdisp", ]) < 50)
  }
  out
}

t0 <- Sys.time()
cl <- parallel::makeCluster(NCORE)
parallel::clusterExport(cl, "VERIFY")
res <- do.call(rbind, parallel::parLapplyLB(cl, split(jobs, seq_len(nrow(jobs))), one))
parallel::stopCluster(cl)
res$limit_ll <- ifelse(res$family == "nbinom", res$ll_limit >= res$ll - 1e-3,
                       res$push <= 1e-3)
res$push_down <- res$push <= 1e-3
tag <- if (VERIFY) "verify" else "phase1"
utils::write.csv(res, file.path(sp, paste0("dispersion_limit_", tag, ".csv")),
                 row.names = FALSE)
cat("fits:", nrow(res), " failed:", sum(!res$ok), " minutes:",
    round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "\n")

## the line, and the gap either side of it
r <- res[res$ok, ]
for (fam in unique(r$family)) {
  x <- r[r$family == fam, ]
  x <- x[!is.na(x$limit_ll) & !is.na(x$push_down), ]
  cat("\n==", fam, "==\n")
  at <- x$limit_ll & x$push_down
  cat(sprintf("at the limit (both labels): %d of %d\n", sum(at), nrow(x)))
  q <- function(v) if (length(v)) signif(stats::quantile(v, c(0, .5, 1), na.rm = TRUE), 3) else NA
  cat("1/sqrt(disp), at the limit:     ", q(x$inv_sqrt[at]), "\n")
  cat("1/sqrt(disp), not at the limit: ", q(x$inv_sqrt[!at]), "\n")
  cat("labels disagree:", sum(x$limit_ll != x$push_down), "\n")
  ## the rule as the package has it: past the line, and flat beyond it
  rl <- x$inv_sqrt < 1e-2 & x$push <= 5e-3
  cat(sprintf("rule (1/sqrt < 1e-2, push <= 5e-3) flags: %d; at the limit: %d; not at it: %d; missed at the limit: %d\n",
              sum(rl), sum(rl & at), sum(rl & !at), sum(!rl & at)))
  for (line in c(1e-3, 1e-2))
    cat(sprintf("line %.0e alone flags: %d; not at the limit among them: %d\n",
                line, sum(x$inv_sqrt < line), sum(x$inv_sqrt < line & !at)))
  if (length(x$inv_sqrt[at]) && length(x$inv_sqrt[!at]))
    cat(sprintf("gap: largest at the limit %.3g, smallest not at it %.3g (%.1f decades)\n",
                max(x$inv_sqrt[at]), min(x$inv_sqrt[!at]),
                log10(min(x$inv_sqrt[!at]) / max(x$inv_sqrt[at]))))
  print(stats::aggregate(cbind(at = at, flagged = rl) ~ design + disp,
                         data = transform(x, at = at, rl = rl), FUN = sum))
  if (VERIFY) {
    h <- grepl("dispersion", x$held)
    cat("held by the fit:", sum(h), "| of them at the limit:", sum(h & at),
        "| at the limit but not held:", sum(!h & at), "\n")
    cat("SE(x) ratio to the Poisson refit, held fits (min, median, max):",
        q(x$se_ratio_x[h]), "\n")
    cat("held fits whose log-dispersion draws stayed within 50:",
        sum(x$draws_finite[h], na.rm = TRUE), "of", sum(h), "\n")
  }
}
