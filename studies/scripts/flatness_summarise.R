## The flatness calibration's registered analysis (flatness_calibration.R,
## "What counts" and addendum A1). Usage:
##   Rscript flatness_summarise.R <rundir> <saved gaussian phase-1 csv>
a <- commandArgs(TRUE); O <- a[1]; GS <- a[2]
rd <- function(path, study) {
  d <- utils::read.csv(file.path(O, sprintf("cal_%s_%s.csv", path, study)))
  d$path <- path; d
}
studies <- c("disp", "gauss", "re", "cases")
d <- do.call(rbind, lapply(c("P1", "P2"), function(p) do.call(rbind, lapply(studies, function(s) {
  x <- rd(p, s); x[, setdiff(names(x), character(0))] }))))
d <- d[d$ok, ]
pmin3 <- pmin(d$push_1.5, d$push_3, d$push_6); pmax3 <- pmax(abs(d$push_1.5), abs(d$push_3), abs(d$push_6))
## a push the objective could not be evaluated at is its own category
nopush <- is.na(d$push_1.5) | is.na(d$push_3) | is.na(d$push_6)
d$new <- ifelse(!d$candidate, "not candidate", ifelse(nopush, "no push",
         ifelse(pmin3 < -0.05, "unconverged", ifelse(pmax3 <= 0.05, "held", "not held"))))
d$old <- ifelse(!d$candidate, "not candidate", ifelse(is.na(d$push_3), "no push",
         ifelse(d$push_3 <= 1e-3, "held", "not held")))
cat("No push (objective not evaluable at a pushed point), candidates:
")
print(table(d$study[d$candidate & nopush], d$path[d$candidate & nopush]))
d$maxabs <- pmax3
d$label <- ifelse(is.na(d$ll_alt), NA, d$ll_alt >= d$ll - 0.05)       # the limit or dropped model as good
d$label_1e3 <- ifelse(is.na(d$ll_alt), NA, d$ll_alt >= d$ll - 1e-3)
key <- function(x) paste(x$study, x$cell, x$rep)
cat("Fits:", nrow(d), " (P1", sum(d$path == "P1"), ", P2", sum(d$path == "P2"), ")\n")
cat("Candidates:", sum(d$candidate & d$path == "P1"), "(P1),", sum(d$candidate & d$path == "P2"), "(P2)\n\n")

## V1: robustness across the paths, among candidates on either path
p1 <- d[d$path == "P1", ]; p2 <- d[d$path == "P2", ]
m <- merge(p1[, c("study", "cell", "rep", "candidate", "new", "old")],
           p2[, c("study", "cell", "rep", "candidate", "new", "old")],
           by = c("study", "cell", "rep"), suffixes = c(".1", ".2"))
m <- m[m$candidate.1 | m$candidate.2, ]
v1 <- do.call(rbind, lapply(split(m, m$study), function(x) data.frame(study = x$study[1], candidates = nrow(x),
  new_same = sum(x$new.1 == x$new.2), new_flips = sum(x$new.1 != x$new.2),
  old_same = sum(x$old.1 == x$old.2), old_flips = sum(x$old.1 != x$old.2))))
v1 <- rbind(v1, data.frame(study = "all", candidates = sum(v1$candidates), new_same = sum(v1$new_same),
  new_flips = sum(v1$new_flips), old_same = sum(v1$old_same), old_flips = sum(v1$old_flips)))
v1$new_same_pct <- round(100 * v1$new_same / v1$candidates, 2)
cat("V1 robustness (candidates on either path; decisions the same on P1 and P2):\n"); print(v1, row.names = FALSE)
cat("V1 verdict:", if (v1$new_same_pct[v1$study == "all"] >= 99.5 &&
      v1$new_flips[v1$study == "all"] <= v1$old_flips[v1$study == "all"]) "HOLDS" else "DOES NOT HOLD", "\n")
fl <- m[m$new.1 != m$new.2, ]
if (nrow(fl)) { cat("  the new rule's flips:\n"); print(fl, row.names = FALSE) }

## V2 / V2b: against the labels (labelled rows: negative binomial, random effect)
lab <- d[d$candidate & !is.na(d$label), ]
v2 <- do.call(rbind, lapply(split(lab, list(lab$study, lab$path), drop = TRUE), function(x) data.frame(
  study = x$study[1], path = x$path[1], labelled = nrow(x),
  new_false_holds = sum(x$new == "held" & !x$label), old_false_holds = sum(x$old == "held" & !x$label),
  new_missed = sum(x$label & x$new != "held"), old_missed = sum(x$label & x$old != "held"))))
cat("\nV2 false holds (held while the limit/dropped model is worse by > 0.05) and V2b missed holds:\n")
print(v2, row.names = FALSE)
cat("V2 verdict:", if (sum(v2$new_false_holds) == 0) "HOLDS" else "DOES NOT HOLD", "\n")
mi <- lab[lab$label & lab$new != "held", ]
if (nrow(mi)) { cat("  V2b missed, by cell (P1):\n"); print(table(mi$study[mi$path == "P1"], mi$cell[mi$path == "P1"])) }

## V3: the band where a decision could flip
cand <- d[d$candidate, ]
v3 <- do.call(rbind, lapply(split(cand, list(cand$study, cand$path), drop = TRUE), function(x) data.frame(
  study = x$study[1], path = x$path[1], candidates = nrow(x),
  in_band = sum(x$maxabs >= 0.025 & x$maxabs <= 0.1), pct = round(100 * mean(x$maxabs >= 0.025 & x$maxabs <= 0.1), 2))))
cat("\nV3 the fragile band (largest |push| in [0.025, 0.1]):\n"); print(v3, row.names = FALSE)

## V4: unconverged
un <- d[d$new == "unconverged", c("path", "study", "cell", "rep", "family", "value", "push_1.5", "push_3", "push_6", "conv", "how_now")]
cat("\nV4 unconverged (a push below -0.05):", nrow(un), "\n"); if (nrow(un)) print(un, row.names = FALSE, digits = 4)

## V5: the CI cases
cs <- d[d$study == "cases", c("path", "cell", "value", "push_1.5", "push_3", "push_6", "new", "old", "conv", "how_now")]
cs$case <- c("nbinom_k_limit", "beta_flat_71272", "beta_curved_63362", "gaussian_stalled")[cs$cell]
cat("\nV5 the four CI cases:\n"); print(cs[order(cs$cell, cs$path), ], row.names = FALSE, digits = 4)

## The gaussian sigma question: where the not-held candidates sit against the
## relative line (0.2 of sd(y)), from the saved study's own sigma / sd(y)
g <- d[d$study == "gauss" & d$candidate, ]
sv <- utils::read.csv(GS)[, c("seed", "rel")]
g <- merge(g, sv, by = "seed")
cat("\nGaussian sigma candidates by the new decision, sigma / sd(y):\n")
print(do.call(rbind, lapply(split(g, list(g$new, g$path), drop = TRUE), function(x) data.frame(
  decision = x$new[1], path = x$path[1], n = nrow(x),
  rel_min = signif(min(x$rel), 3), rel_median = signif(stats::median(x$rel), 3), rel_max = signif(max(x$rel), 3),
  push6_median = signif(stats::median(x$push_6), 3)))), row.names = FALSE)
utils::write.csv(d, file.path(O, "calibration_all.csv"), row.names = FALSE)
