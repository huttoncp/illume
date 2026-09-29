## ---------------------------------------------------------------------------
## df tables, addendum A1 (the ML arm): ML against REML on the same
## replicates, PAIRED. Written after the ML arm ran (2026-09-29), so a
## post-hoc reading beside the pre-registered verdicts, never a verdict.
##
## The ML arm used the main run's seeds, so each replicate was fitted to the
## same data by both methods. Per cell and per method (z, Satterthwaite), the
## exact two-sided McNemar test on the replicates both arms fitted: whether
## the interval covered (D1 to D3) or the test rejected at 0.05 (D4), counting
## the replicates where only REML's did and where only ML's did.
##
## Usage: Rscript df_tables_ml_paired.R <ML per-fit csv> <REML main per-fit csv> <out csv>
## ---------------------------------------------------------------------------
a <- commandArgs(trailingOnly = TRUE)
ml <- utils::read.csv(a[1]); rm <- utils::read.csv(a[2])
ml <- ml[ml$ok %in% TRUE, ]; rm <- rm[rm$ok %in% TRUE, ]
m <- merge(rm, ml, by = c("cell", "rep"), suffixes = c("_reml", "_ml"))
mcnemar <- function(a, b) {
  ok <- !is.na(a) & !is.na(b); a <- a[ok]; b <- b[ok]
  n10 <- sum(a & !b); n01 <- sum(!a & b)
  p <- if (n10 + n01 == 0) 1 else stats::binom.test(n10, n10 + n01)$p.value
  c(n = length(a), reml_only = n10, ml_only = n01, p = p)
}
out <- do.call(rbind, lapply(sort(unique(m$cell)), function(cl) {
  z <- m[m$cell == cl, ]; arm <- z$arm_reml[1]
  ev <- function(k, sfx) if (arm == "D4") z[[paste0("p_", k, "_", sfx)]] < 0.05
                         else as.logical(z[[paste0("cover_", k, "_", sfx)]])
  do.call(rbind, lapply(c("z", "s"), function(k) {
    r <- mcnemar(ev(k, "reml"), ev(k, "ml"))
    data.frame(cell = cl, arm = arm, G = z$G_reml[1], m = z$m_reml[1], icc = z$icc_reml[1],
               what = if (arm == "D4") "rejects" else "covers", method = k,
               reml = mean(ev(k, "reml"), na.rm = TRUE), ml = mean(ev(k, "ml"), na.rm = TRUE),
               n = r[["n"]], reml_only = r[["reml_only"]], ml_only = r[["ml_only"]],
               p = r[["p"]])
  }))
}))
utils::write.csv(out, a[3], row.names = FALSE)
print(out, row.names = FALSE, digits = 3)
