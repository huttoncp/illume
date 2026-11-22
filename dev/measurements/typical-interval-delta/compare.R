## Two saved sets of intervals from measure.R, one per build: which are
## identical(), and how far apart the others are, in either limit.
##   Rscript dev/measurements/typical-interval-delta/compare.R <before.rds> <after.rds>
args <- commandArgs(trailingOnly = TRUE)
a <- readRDS(args[1]); b <- readRDS(args[2])
stopifnot(identical(names(a), names(b)))
cat(sprintf("%-18s %-10s %s\n", "fit", "identical", "max |difference| in either limit"))
for (k in names(a))
  cat(sprintf("%-18s %-10s %s\n", k, identical(a[[k]], b[[k]]),
              format(max(abs(a[[k]]$lower - b[[k]]$lower), abs(a[[k]]$upper - b[[k]]$upper)), digits = 3)))
