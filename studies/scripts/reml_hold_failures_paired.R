## A reading made AFTER the hold-rule study ran (2026-09-29), not part of its
## pre-registration: REML's failures against ML's on the same data, paired,
## by exact McNemar, over all fits and per cell.
## Usage: Rscript reml_hold_failures_paired.R <run dir> <out csv>
a <- commandArgs(TRUE)
rows <- list()
for (t in c("main", "fresh")) {
  r <- read.csv(file.path(a[1], sprintf("reml_hold_%s.csv", t)))
  fr <- !(r$reml_ok %in% TRUE) | r$reml_failed %in% TRUE
  fm <- !(r$ml_ok %in% TRUE) | r$ml_failed %in% TRUE
  one <- function(z, cell) {
    x <- sum(fr[z] & !fm[z]); y <- sum(!fr[z] & fm[z])
    data.frame(run = t, cell = cell, reml_fails = sum(fr[z]), ml_fails = sum(fm[z]),
               reml_only = x, ml_only = y,
               p = if (x + y) binom.test(x, x + y)$p.value else 1)
  }
  rows[[length(rows) + 1L]] <- one(rep(TRUE, nrow(r)), "all")
  for (cl in sort(unique(r$cell))) rows[[length(rows) + 1L]] <- one(r$cell == cl, as.character(cl))
}
out <- do.call(rbind, rows)
utils::write.csv(out, a[2], row.names = FALSE)
print(out)
