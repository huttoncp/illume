## What a typical curve's interval costs, on a model with many groups: a Weibull
## and a flexible-parametric fit of time ~ s(x, k = 8) + z + (1 | h) with G groups
## by REML, so the joint precision has m = G + a handful of parameters rows.
## Two steps, in separate processes so that the peak memory of the second is the
## curve's alone and not the fit's:
##   Rscript dev/measurements/typical-interval-delta/cost.R fit <G> <folder>
##   Rscript dev/measurements/typical-interval-delta/cost.R run <folder> <label for the build> <weibull or rp>
## `fit` saves the fits to the folder; `run` reads one back (a process for each
## family, since the peak is the process's whole life), times ilm_survival()
## for a typical group (3 rows, 100 times, the defaults) three times, and reports
## the process's peak working set before and after, from Windows.
suppressPackageStartupMessages(library(illume))
args <- commandArgs(trailingOnly = TRUE)
q <- function(e) suppressMessages(suppressWarnings(e))
peak_mb <- function()
  as.numeric(system(sprintf("powershell -NoProfile -Command \"(Get-Process -Id %d).PeakWorkingSet64\"",
                            Sys.getpid()), intern = TRUE)) / 2^20
if (args[1] == "fit") {
  G <- as.integer(args[2]); dir.create(args[3], showWarnings = FALSE)
  set.seed(7)
  n <- 3L * G
  d <- data.frame(x = runif(n), z = rnorm(n), h = factor(sample(seq_len(G), n, TRUE)))
  hl <- rnorm(G, 0, 0.4)
  eta <- 1 + 0.5 * d$z + 1.2 * sin(2 * pi * d$x) + hl[as.integer(d$h)]
  tt <- exp(eta + 0.5 * log(rexp(n)))
  ct <- rexp(n, 1 / (2 * median(tt)))
  d$time <- pmin(tt, ct); d$event <- as.integer(tt <= ct)
  for (fam in c("weibull", "rp")) {
    tm <- system.time(f <- q(ilm_model(time ~ s(x, k = 8) + z + (1 | h), data = d, family = fam,
                                      censor = ilm_surv(d$time, d$event), reml = TRUE, verbose = FALSE)))
    cat(fam, "fit:", format(tm[["elapsed"]], digits = 4), "s elapsed | rows of the joint precision:",
        nrow(illume:::ilm_joint_prec(f)), "\n")
    saveRDS(list(fit = f, newdata = data.frame(x = c(0.25, 0.5, 0.75), z = c(-1, 0, 1),
                                               h = factor(levels(d$h)[1], levels = levels(d$h)))),
            file.path(args[3], paste0(fam, ".rds")))
  }
} else {
  cat("build:", args[3], "| illume", as.character(packageVersion("illume")), "| R",
      paste0(R.version$major, ".", R.version$minor), "| date",
      format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"), "\n")
  for (fam in args[4]) {
    x <- readRDS(file.path(args[2], paste0(fam, ".rds")))
    f <- x$fit
    cat(sprintf("\n%s, joint precision %d x %d\n", fam, nrow(illume:::ilm_joint_prec(f)),
                ncol(illume:::ilm_joint_prec(f))))
    cat("  peak working set after reading the fit:", format(peak_mb(), digits = 5), "MB\n")
    for (r in 1:3) {
      tm <- system.time(s <- q(ilm_survival(f, newdata = x$newdata, groups = "typical")))
      cat(sprintf("  run %d: %.2f s elapsed; mean interval width over the %d points: %.5f\n", r,
                  tm[["elapsed"]], nrow(s), mean(s$upper - s$lower)))
    }
    cat("  peak working set after the curves:", format(peak_mb(), digits = 5), "MB\n")
  }
}
