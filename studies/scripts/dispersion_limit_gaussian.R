## A gaussian residual SD at zero: the same edge as a negative binomial's k
## at infinity, at the other end of the log scale.
##
## With one observation per cell of a correlation over time, the latent
## process can take up all the noise, and the residual SD runs to zero: one
## fit in 500 of another agent's study had it at 0.0004, and draws of log
## sigma spanned -529 to +466. The rule adapted from dispersion_limit.R holds
## sigma when it is below 1e-3 of the response's SD -- as a random effect's
## SD below 1e-3 is taken as zero -- and the objective moves by at most 1e-3
## when log sigma is pushed 3 lower. This measures whether that line
## separates the fits whose optimum is sigma = 0 from the rest.
##
## Each fit is labelled from outside the rule: at the limit when the
## objective does not rise by more than 1e-3 as log sigma is pushed 3 lower
## (the other parameters held, the random effects re-optimised). Controls:
## a random intercept with one observation per level, where sigma and the
## intercept's variance trade along a ridge but neither need reach zero, and
## an ordinary random intercept.
##
## Usage: Rscript dispersion_limit_gaussian.R <nrep> <ncore> <outdir> [verify]
## ---------------------------------------------------------------------------
args  <- commandArgs(trailingOnly = TRUE)
NREP  <- if (length(args) >= 1) as.integer(args[1]) else 100L
NCORE <- if (length(args) >= 2) as.integer(args[2]) else 2L
sp    <- if (length(args) >= 3) args[3] else "."
VERIFY <- length(args) >= 4 && identical(args[4], "verify")
dir.create(sp, showWarnings = FALSE, recursive = TRUE)

cells <- rbind(
  data.frame(design = "ar1", noise = c(0.5, 0.2, 0.05)),
  data.frame(design = "ri1", noise = 0.5),
  data.frame(design = "ri",  noise = 0.5))
cells$cell <- seq_len(nrow(cells))
jobs <- merge(cells, data.frame(rep = seq_len(NREP)))

one <- function(job) {
  suppressMessages(library(illume))
  seed <- 6007L * job$cell + job$rep
  set.seed(seed)
  if (job$design == "ar1") {
    G <- 20; Tn <- 16
    d <- expand.grid(t = seq_len(Tn), g = factor(seq_len(G)))
    lat <- unlist(lapply(seq_len(G), function(i) as.numeric(
      stats::arima.sim(list(ar = 0.7), Tn, sd = 0.8 * sqrt(1 - 0.49)))))
  } else if (job$design == "ri1") {
    d <- data.frame(g = factor(seq_len(300)))
    lat <- stats::rnorm(300, 0, 0.6)
  } else {
    d <- data.frame(g = factor(rep(seq_len(30), each = 10)))
    lat <- stats::rnorm(30, 0, 0.6)[d$g]
  }
  d$x <- stats::rnorm(nrow(d))
  d$y <- 0.5 + 0.3 * d$x + lat + stats::rnorm(nrow(d), 0, job$noise)
  ar <- if (job$design == "ar1") ilm_ar1(~ t | g) else NULL
  fml <- if (job$design == "ar1") y ~ x + (1 | g) else y ~ x + (1 | g)
  f <- tryCatch(suppressMessages(suppressWarnings(
    ilm_model(fml, data = d, family = "gaussian", ar = ar, verbose = FALSE))),
    error = function(e) NULL)
  out <- data.frame(job, seed = seed, ok = !is.null(f), sigma = NA_real_,
                    rel = NA_real_, se_logsig = NA_real_, push = NA_real_,
                    held = NA_character_, how = NA_character_, draws_finite = NA)
  if (is.null(f)) return(out)
  pn <- names(f$opt$par); il <- which(pn == "logdisp")
  ls <- unname(f$opt$par[il])
  out$sigma <- exp(ls); out$rel <- exp(ls) / stats::sd(d$y)
  V <- tryCatch(suppressWarnings(vcov(f, full = TRUE)), error = function(e) NULL)
  if (!is.null(V)) out$se_logsig <- sqrt(V[il, il])
  out$held <- paste(f$hessian_held, collapse = ",")
  out$how <- f$hessian_how
  p2 <- f$opt$par; p2[il] <- p2[il] - 3
  nll0 <- tryCatch(f$obj$fn(f$opt$par), error = function(e) NA_real_)
  nll1 <- tryCatch(f$obj$fn(p2), error = function(e) NA_real_)
  out$push <- nll1 - nll0
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
res$at <- res$push <= 1e-3
res$rule <- res$rel < 1e-3 & res$push <= 1e-3
tag <- if (VERIFY) "verify" else "phase1"
utils::write.csv(res, file.path(sp, paste0("dispersion_limit_gaussian_", tag, ".csv")),
                 row.names = FALSE)
cat("fits:", nrow(res), " failed:", sum(!res$ok), " minutes:",
    round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "\n")
r <- res[res$ok & !is.na(res$at), ]
print(stats::aggregate(cbind(fits = 1, at_limit = at, rule = rule,
                             rel_below_1e3 = rel < 1e-3) ~ design + noise,
                       data = r, FUN = sum))
q <- function(v) if (length(v)) signif(stats::quantile(v, c(0, .5, .99, 1), na.rm = TRUE), 3) else NA
cat("sigma / sd(y) at the limit (min, median, 99%, max):", q(r$rel[r$at]), "\n")
cat("sigma / sd(y) not at the limit (min, median, 99%, max):", q(r$rel[!r$at]), "\n")
cat(sprintf("rule flags %d; at the limit %d; not at it %d; missed at the limit %d\n",
            sum(r$rule), sum(r$rule & r$at), sum(r$rule & !r$at), sum(!r$rule & r$at)))
if (VERIFY) {
  h <- grepl("dispersion", r$held)
  cat("held by the fit:", sum(h), "| at the limit:", sum(h & r$at),
      "| draws of log sigma within 50:", sum(r$draws_finite[h], na.rm = TRUE), "of", sum(h), "\n")
}
