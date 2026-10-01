## The flatness confirmation's registered analysis (flatness_confirm.R, "The
## rule under test" and "Paths, and what counts"), written before the run.
## Usage: Rscript flatness_confirm_summarise.R <rundir>
a <- commandArgs(TRUE); O <- a[1]
studies <- c("disp", "gauss", "re", "cases")
d <- do.call(rbind, lapply(c("P1", "P2"), function(p) do.call(rbind, lapply(studies, function(s) {
  fn <- file.path(O, sprintf("conf_%s_%s.csv", p, s))
  if (!file.exists(fn)) { cat("missing:", fn, "\n"); return(NULL) }
  x <- utils::read.csv(fn); x$path <- p; x }))))
cat("Fits:", nrow(d), " not fitted:", sum(!d$ok), "\n")
d <- d[d$ok, ]
pmin3 <- pmin(d$push_1.5, d$push_3, d$push_6)
d$maxabs <- pmax(abs(d$push_1.5), abs(d$push_3), abs(d$push_6))
nopush <- is.na(d$push_1.5) | is.na(d$push_3) | is.na(d$push_6)
## the pushes as ruled, without the floor
d$pushes <- ifelse(!d$candidate, "not candidate", ifelse(nopush, "no push",
            ifelse(pmin3 < -0.05, "unconverged", ifelse(d$maxabs <= 0.05, "held", "not held"))))
## the rule under test: more than 6 below the line is held by value
d$by_value <- d$candidate & d$below < -6
d$floor <- ifelse(d$by_value, "held", d$pushes)
## today's rule: one push of 3, held at a rise of at most 1e-3
d$today <- ifelse(!d$candidate, "not candidate", ifelse(is.na(d$push_3), "no push",
           ifelse(d$push_3 <= 1e-3, "held", "not held")))
d$as_good <- ifelse(is.na(d$cost_fair), NA, d$cost_fair <= 0.05)
cat("Candidates:", sum(d$candidate & d$path == "P1"), "(P1),", sum(d$candidate & d$path == "P2"), "(P2)\n")
cat("No push (objective not evaluable at a pushed point), candidates not held by value:\n")
print(table(d$study[d$candidate & nopush & !d$by_value], d$path[d$candidate & nopush & !d$by_value]))

## V1 robustness, among candidates on either path
k <- c("study", "cell", "rep")
m <- merge(d[d$path == "P1", c(k, "candidate", "floor", "pushes", "today")],
           d[d$path == "P2", c(k, "candidate", "floor", "pushes", "today")],
           by = k, suffixes = c(".1", ".2"))
m <- m[m$candidate.1 | m$candidate.2, ]
v1 <- do.call(rbind, lapply(split(m, m$study), function(x) data.frame(study = x$study[1],
  candidates = nrow(x), floor_flips = sum(x$floor.1 != x$floor.2),
  pushes_flips = sum(x$pushes.1 != x$pushes.2), today_flips = sum(x$today.1 != x$today.2))))
v1 <- rbind(v1, data.frame(study = "all", candidates = sum(v1$candidates),
  floor_flips = sum(v1$floor_flips), pushes_flips = sum(v1$pushes_flips), today_flips = sum(v1$today_flips)))
v1$floor_same_pct <- round(100 * (1 - v1$floor_flips / v1$candidates), 2)
cat("\nV1 robustness (candidates on either path; decisions that differ between P1 and P2):\n")
print(v1, row.names = FALSE)
al <- v1[v1$study == "all", ]
cat("V1 verdict:", if (al$floor_same_pct >= 99.5 && al$floor_flips <= al$today_flips) "HOLDS" else "DOES NOT HOLD", "\n")
fl <- m[m$floor.1 != m$floor.2, ]
if (nrow(fl)) { cat("  the floor rule's flips:\n"); print(fl, row.names = FALSE) }

## V2 / V2b against the fair label (negative binomial k, random-effect SD)
lab <- d[d$candidate & !is.na(d$as_good), ]
v2 <- do.call(rbind, lapply(split(lab, list(lab$study, lab$path), drop = TRUE), function(x) data.frame(
  study = x$study[1], path = x$path[1], labelled = nrow(x),
  floor_false_holds = sum(x$floor == "held" & !x$as_good),
  of_them_by_value = sum(x$by_value & !x$as_good),
  pushes_false_holds = sum(x$pushes == "held" & !x$as_good),
  today_false_holds = sum(x$today == "held" & !x$as_good),
  floor_missed = sum(x$as_good & x$floor != "held"),
  today_missed = sum(x$as_good & x$today != "held"))))
cat("\nV2 false holds (held while the fair label puts the limit or dropped model worse by > 0.05),",
    "and V2b missed holds:\n")
print(v2, row.names = FALSE)
cat("V2 verdict:", if (sum(v2$floor_false_holds) == 0) "HOLDS" else "DOES NOT HOLD", "\n")
fh <- lab[lab$floor == "held" & !lab$as_good, c("path", k, "family", "value", "below", "by_value",
          "maxabs", "cost_refit", "cost_prof", "cost_fair")]
if (nrow(fh)) { cat("  the false holds:\n"); print(fh, row.names = FALSE, digits = 4) }
cat("  which label was the better (labelled candidates): refit", sum(lab$cost_fair == lab$cost_refit, na.rm = TRUE),
    ", profile", sum(lab$cost_fair == lab$cost_prof & lab$cost_fair != lab$cost_refit, na.rm = TRUE), "\n")

## V3 the band, among pushed candidates
pc <- d[d$candidate & !d$by_value & !nopush, ]
v3 <- do.call(rbind, lapply(split(pc, list(pc$study, pc$path), drop = TRUE), function(x) data.frame(
  study = x$study[1], path = x$path[1], pushed = nrow(x),
  in_band = sum(x$maxabs >= 0.025 & x$maxabs <= 0.1),
  pct = round(100 * mean(x$maxabs >= 0.025 & x$maxabs <= 0.1), 2))))
cat("\nV3 the fragile band (largest |push| in [0.025, 0.1]), pushed candidates:\n"); print(v3, row.names = FALSE)

## V4 unconverged under the floor rule
un <- d[d$floor == "unconverged", c("path", k, "family", "value", "below", "push_1.5", "push_3", "push_6",
        "conv", "check_optimizer", "check_gradient")]
cat("\nV4 unconverged (a push below -0.05):", nrow(un), "\n"); if (nrow(un)) print(un, row.names = FALSE, digits = 4)

## V5 the four CI cases
cs <- d[d$study == "cases", c("path", "cell", "value", "below", "by_value", "push_1.5", "push_3", "push_6",
        "floor", "today", "conv", "check_optimizer", "check_gradient")]
cs$case <- c("nbinom_k_limit", "beta_flat_71272", "beta_curved_63362", "gaussian_stalled")[cs$cell]
cat("\nV5 the four CI cases:\n"); print(cs[order(cs$cell, cs$path), ], row.names = FALSE, digits = 4)

## V6 the floor's reach
bv <- d[d$by_value, ]
v6 <- do.call(rbind, lapply(split(bv, list(bv$study, bv$path), drop = TRUE), function(x) data.frame(
  study = x$study[1], path = x$path[1], by_value = nrow(x),
  below_max = signif(max(x$below), 3), below_median = signif(stats::median(x$below), 3),
  below_min = signif(min(x$below), 3))))
cat("\nV6 held by value (more than 6 below the line), per study and path:\n"); print(v6, row.names = FALSE)
rw <- bv[bv$check_optimizer %in% "FAIL" | bv$check_gradient %in% "FAIL",
         c("path", k, "family", "value", "below", "conv", "check_optimizer", "check_gradient")]
cat("V6 held by value with a FAILed optimizer or gradient check (the runaway case):", nrow(rw), "\n")
if (nrow(rw)) print(rw, row.names = FALSE, digits = 4)

utils::write.csv(d, file.path(O, "confirm_all.csv"), row.names = FALSE)
