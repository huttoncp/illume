## Can sdreport() hand back a finite covariance of the fixed parameters when the
## Hessian is not positive definite? A saddle, x^2 - y^2 with a quartic to hold it,
## evaluated at its saddle point. Plain RTMB; no illume code.
## Run: Rscript dev/measurements/typical-interval-delta/pdhess-probe.R
suppressPackageStartupMessages(library(RTMB))
cat("RTMB", as.character(packageVersion("RTMB")), "| TMB", as.character(packageVersion("TMB")),
    "| R", paste0(R.version$major, ".", R.version$minor), "| date",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"), "\n")
f <- function(p) p$x^2 - p$y^2 + 0.01 * (p$x^4 + p$y^4)
obj <- MakeADFun(f, list(x = 0, y = 0), silent = TRUE)
sdr <- sdreport(obj, par.fixed = c(0, 0))
cat("pdHess:", sdr$pdHess, "\n")
cat("cov.fixed (finite entries: ", all(is.finite(sdr$cov.fixed)), "):\n", sep = "")
print(sdr$cov.fixed)
