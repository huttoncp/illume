## The hold rule under REML, for a gaussian mixed model (Craig's item 249).
##
## PRE-REGISTERED DESIGN, committed before any run.
##
## ## The question
## Craig ruled (item 249) that a gaussian mixed model is fitted by REML by
##   default. When a covariance sits at its boundary, the fixed effects'
##   standard errors come from rule C: hold only the flat directions of the
##   covariance and estimate everything else (item 4, chosen 2026-09-23). Rule
##   C was measured on 4,000 multinomial and binomial fits, all by ML: its SEs
##   equalled the reduced model's (median ratio 1.000) and covered best. Under
##   REML the fit's outer parameters are the variances alone, the coefficients'
##   covariance is a Schur complement (7f57fed), and a fit whose every variance
##   is at its boundary is held as the plain model (P5, 2d54943). None of that
##   was in the item 4 study. This study asks whether rule C gives the reduced
##   model's SEs, and calibrated intervals, on the path that becomes the default.
##
## ## Arms and cells (gaussian; the effect of interest has true value 0.5)
## | Arm | Model | Truth at the boundary | Cells |
## |---|---|---|---|
## | H1 intercept at zero | y ~ trt + (1 \| g), trt constant within a cluster, half the clusters treated | intercept SD 0 | clusters G in {8, 20}, cluster size m = 5 |
## | H2 slope at zero | y ~ x + (1 + x \| g), x varying within a cluster | intercept SD 0.5, slope SD 0 | G in {8, 20} x m in {5, 10} |
## | H3 correlation at one | y ~ x + (1 + x \| g) | intercept SD 0.5, slope SD 0.3, correlation 1 | G in {8, 20}, m = 10 |
## | H4 crossed, one at zero | y ~ x + (1 \| a) + (1 \| b), every (a, b) pair once | a's SD 0.5, b's SD 0 | G_a = G_b in {8, 20} |
##
## 10 cells, residual SD 1, 1,000 replicates per cell.
##
## ## Per replicate
## - The REML fit (reml = TRUE), and the ML fit (reml = FALSE), of the arm's
##   model to the same data. For each: whether a direction was held and which
##   terms; the effect's estimate and SE; its Satterthwaite df and 95% interval
##   (the default table's), and its z interval; whether it covers 0.5.
## - The reference: the reduced model the boundary implies, refitted by REML --
##   H1 y ~ trt (no random term: an exact-t fit), H2 y ~ x + (1 | g), H4
##   y ~ x + (1 | a). A rank-one covariance of an intercept and a slope has no
##   form in illume, so H3's reference is lme4::lmer(REML = TRUE) at its own
##   fit, which is not graded (its boundary SEs are not rule C's).
## - Failures: a fit that errors, or whose optimiser or gradient check FAILs.
##
## ## What counts, fixed now
## - V1, the hold rule under REML: in H1, H2 and H4, among the REML fits held
##   at the term that is zero in truth, the median of (rule C's SE / the
##   reduced refit's SE) is within [0.99, 1.01] in every cell, and at least 95%
##   of those fits are within [0.97, 1.03].
## - V2, the intervals: the Satterthwaite coverage of the REML held fits is
##   within 2 Monte Carlo SEs of 0.95, at the number of held fits, in every
##   cell with G = 20.
## - V3, nothing blows up: no REML fit's SE exceeds 10 times its reference's
##   (the reduced refit, or lme4's in H3), and in no cell do REML fits fail
##   more often than ML fits.
## - If V1 or V3 does not hold, that goes to Craig with a proposal before the
##   REML default ships for fits held at a boundary.
## - Beside the verdicts, pre-registered and not graded: per cell, the paired
##   exact McNemar test of coverage, REML against ML, on the same data; the
##   share of fits held by each; H3's SE against lme4's.
##
## ## Fresh seeds
## The whole design again at 500 replicates per cell with seeds offset by
##   1,000,000, on the same build: V1 to V3 must hold there too, and any that
##   holds in one run and not the other is reported as unresolved.
##
## ## Build and size
## A pinned library built from this branch at the commit that pre-registers
##   the study, with lme4 from the user library; versions recorded before and
##   after the runs. About 30,000 fits and 15,000 fresh, a few tenths of a
##   second each: about 4 to 6 core-hours in all.
##
## ## Addendum A1 (2026-09-29, before any main run): the build
## A smoke run (3 replicates a cell) found every table's df infinite: the
##   branch at the pre-registration (f5fb8e1) did not yet carry the
##   Satterthwaite tables, which are on the df-tables branch. The default path
##   Craig ruled on is REML with those tables, so df-tables is merged into this
##   branch (0d792a7) and the pinned library is built from the commit that
##   records this addendum. Nothing else changes; the smoke run's numbers are
##   not results.
##
## Usage: Rscript reml_hold.R <nrep> <ncore> <outdir> [offset] [cells|all]
##        Rscript reml_hold.R summarise <outdir>
## ---------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
SUMMARISE <- length(args) >= 1 && identical(args[1], "summarise")
if (!SUMMARISE) {
  NREP   <- if (length(args) >= 1) as.integer(args[1]) else 100L
  NCORE  <- if (length(args) >= 2) as.integer(args[2]) else 1L
  sp     <- if (length(args) >= 3) args[3] else "."
  OFFSET <- if (length(args) >= 4) as.integer(args[4]) else 0L
  ## a subset of cells, for the smoke run only ("all" for every cell)
  ONLY   <- if (length(args) >= 5 && !identical(args[5], "all"))
    as.integer(strsplit(args[5], ",")[[1]]) else NULL
} else sp <- if (length(args) >= 2) args[2] else "."
dir.create(sp, showWarnings = FALSE, recursive = TRUE)
suppressPackageStartupMessages(library(illume))
`%||%` <- function(a, b) if (is.null(a)) b else a

## ---- the cells ----------------------------------------------------------------
cells <- rbind(
  data.frame(arm = "H1", G = c(8, 20), m = 5),
  data.frame(arm = "H2", expand.grid(G = c(8, 20), m = c(5, 10))),
  data.frame(arm = "H3", G = c(8, 20), m = 10),
  data.frame(arm = "H4", G = c(8, 20), m = NA))
cells$cell <- seq_len(nrow(cells))
stopifnot(nrow(cells) == 10L)
TRUTH <- 0.5
LEVEL <- 0.95
fml <- list(H1 = y ~ trt + (1 | g), H2 = y ~ x + (1 + x | g), H3 = y ~ x + (1 + x | g),
            H4 = y ~ x + (1 | a) + (1 | b))
reduced <- list(H1 = y ~ trt, H2 = y ~ x + (1 | g), H4 = y ~ x + (1 | a))
coef_of <- c(H1 = "trt", H2 = "x", H3 = "x", H4 = "x")
zero_term <- c(H1 = "g", H2 = "g", H3 = "g", H4 = "b")

## ---- the data -----------------------------------------------------------------
gen <- function(ce, seed) {
  set.seed(seed)
  if (ce$arm == "H4") {
    G <- ce$G
    d <- expand.grid(a = factor(seq_len(G)), b = factor(seq_len(G)))
    d$x <- rnorm(nrow(d))
    d$y <- 1 + TRUTH * d$x + rnorm(G, 0, 0.5)[d$a] + rnorm(nrow(d))
    return(d)
  }
  G <- ce$G; m <- ce$m
  d <- data.frame(g = factor(rep(seq_len(G), each = m)))
  if (ce$arm == "H1") {
    d$trt <- rep(rep(0:1, length.out = G)[sample.int(G)], each = m)
    d$y <- 1 + TRUTH * d$trt + rnorm(nrow(d))
  } else {
    d$x <- rnorm(nrow(d))
    b0 <- rnorm(G, 0, 0.5)
    b1 <- if (ce$arm == "H2") rep(0, G) else 0.3 * b0 / 0.5   # correlation 1
    d$y <- 1 + b0[d$g] + (TRUTH + b1[d$g]) * d$x + rnorm(nrow(d))
  }
  d
}

## ---- one replicate ----------------------------------------------------------------
fit_row <- function(f, nm, pre) {
  out <- list()
  if (is.null(f)) { out[[paste0(pre, "ok")]] <- FALSE; return(out) }
  ck <- f$checks
  ct <- tryCatch(suppressMessages(suppressWarnings(ilm_coef_table(f))), error = function(e) NULL)
  row <- if (is.data.frame(ct) && nm %in% rownames(ct)) ct[nm, ] else NULL
  est <- if (is.null(row)) NA_real_ else row$Estimate
  se <- if (is.null(row)) NA_real_ else row[["Std. Error"]]
  df <- if (!is.null(row) && "df" %in% names(row)) row$df else Inf
  q <- if (is.finite(df) && df > 0) qt(1 - (1 - LEVEL) / 2, df) else qnorm(1 - (1 - LEVEL) / 2)
  out[[paste0(pre, "ok")]] <- TRUE
  out[[paste0(pre, "failed")]] <- any(ck$status[ck$check %in% c("optimizer", "gradient")] == "FAIL")
  out[[paste0(pre, "held")]] <- length(f$hessian_held) > 0L
  out[[paste0(pre, "held_terms")]] <- paste(f$hessian_held, collapse = ",")
  out[[paste0(pre, "how")]] <- f$hessian_how %||% NA_character_
  out[[paste0(pre, "est")]] <- est
  out[[paste0(pre, "se")]] <- se
  out[[paste0(pre, "df")]] <- df
  out[[paste0(pre, "cover_s")]] <- abs(est - TRUTH) <= q * se
  out[[paste0(pre, "cover_z")]] <- abs(est - TRUTH) <= qnorm(1 - (1 - LEVEL) / 2) * se
  out
}
quiet <- function(expr) tryCatch(suppressMessages(suppressWarnings(expr)), error = function(e) NULL)

one <- function(job) {
  ce <- cells[cells$cell == job$cell, ]
  seed <- 104729L * ce$cell + job$rep + OFFSET
  d <- gen(ce, seed)
  nm <- coef_of[[ce$arm]]
  key <- list(cell = ce$cell, arm = ce$arm, G = ce$G, m = ce$m, rep = job$rep, seed = seed)
  fr <- quiet(ilm_model(fml[[ce$arm]], data = d, family = "gaussian", reml = TRUE, verbose = FALSE))
  fm <- quiet(ilm_model(fml[[ce$arm]], data = d, family = "gaussian", reml = FALSE, verbose = FALSE))
  ref_se <- NA_real_
  if (ce$arm %in% names(reduced)) {
    fx <- quiet(ilm_model(reduced[[ce$arm]], data = d, family = "gaussian", reml = TRUE,
                          verbose = FALSE))
    ct <- if (!is.null(fx)) quiet(ilm_coef_table(fx))
    if (is.data.frame(ct) && nm %in% rownames(ct)) ref_se <- ct[nm, "Std. Error"]
  } else if (requireNamespace("lme4", quietly = TRUE)) {
    lm4 <- quiet(lme4::lmer(fml[[ce$arm]], data = d, REML = TRUE))
    if (!is.null(lm4)) ref_se <- sqrt(diag(as.matrix(stats::vcov(lm4))))[[nm]]
  }
  as.data.frame(c(key, fit_row(fr, nm, "reml_"), fit_row(fm, nm, "ml_"),
                  list(ref_se = ref_se, zero_term = zero_term[[ce$arm]])),
                stringsAsFactors = FALSE)
}

## ---- summaries ----------------------------------------------------------------
mcnemar <- function(a, b) {
  ok <- !is.na(a) & !is.na(b); a <- a[ok]; b <- b[ok]
  n10 <- sum(a & !b); n01 <- sum(!a & b)
  c(reml_only = n10, ml_only = n01,
    p = if (n10 + n01 == 0) 1 else stats::binom.test(n10, n10 + n01)$p.value)
}
summarise_fits <- function(r) {
  do.call(rbind, lapply(sort(unique(r$cell)), function(cl) {
    z <- r[r$cell == cl, ]; arm <- z$arm[1]
    ok <- z$reml_ok %in% TRUE
    ## held at the term that is zero in truth
    hz <- ok & z$reml_held %in% TRUE &
      vapply(strsplit(z$reml_held_terms, ","), function(h) z$zero_term[1] %in% h, TRUE)
    ratio <- z$reml_se / z$ref_se
    mc <- mcnemar(z$reml_cover_s, z$ml_cover_s)
    nh <- sum(hz)
    covh <- mean(z$reml_cover_s[hz], na.rm = TRUE)
    data.frame(cell = cl, arm = arm, G = z$G[1], m = z$m[1], n = nrow(z),
               reml_ok = sum(ok), ml_ok = sum(z$ml_ok %in% TRUE),
               reml_failed = sum(!ok | z$reml_failed %in% TRUE),
               ml_failed = sum(!(z$ml_ok %in% TRUE) | z$ml_failed %in% TRUE),
               reml_held = mean(z$reml_held[ok] %in% TRUE),
               ml_held = mean(z$ml_held[z$ml_ok %in% TRUE] %in% TRUE),
               held_zero = nh,
               ratio_median = if (nh) stats::median(ratio[hz], na.rm = TRUE) else NA_real_,
               ratio_within = if (nh) mean(ratio[hz] >= 0.97 & ratio[hz] <= 1.03, na.rm = TRUE) else NA_real_,
               ratio_max = suppressWarnings(max(ratio[ok], na.rm = TRUE)),
               cover_held_s = covh,
               mcse_held = if (nh) sqrt(0.95 * 0.05 / nh) else NA_real_,
               cover_reml_s = mean(z$reml_cover_s[ok], na.rm = TRUE),
               cover_ml_s = mean(z$ml_cover_s[z$ml_ok %in% TRUE], na.rm = TRUE),
               cover_reml_z = mean(z$reml_cover_z[ok], na.rm = TRUE),
               cover_ml_z = mean(z$ml_cover_z[z$ml_ok %in% TRUE], na.rm = TRUE),
               reml_only = mc[["reml_only"]], ml_only = mc[["ml_only"]], p_paired = mc[["p"]],
               stringsAsFactors = FALSE)
  }))
}
verdicts <- function(s) {
  g <- s[s$arm %in% c("H1", "H2", "H4"), ]
  v1 <- all(!is.na(g$ratio_median) & abs(g$ratio_median - 1) <= 0.01 & g$ratio_within >= 0.95)
  g20 <- s[s$G == 20, ]
  v2 <- all(!is.na(g20$cover_held_s) & abs(g20$cover_held_s - 0.95) <= 2 * g20$mcse_held)
  v3 <- all(is.finite(s$ratio_max) & s$ratio_max <= 10) && all(s$reml_failed <= s$ml_failed)
  data.frame(verdict = c(
    "V1: rule C's SE / the reduced refit's, held fits: median within [0.99, 1.01], 95% within [0.97, 1.03], in every H1, H2, H4 cell",
    "V2: Satterthwaite coverage of REML held fits calibrated in every cell with G = 20",
    "V3: no REML SE above 10 times its reference, and REML fails no more often than ML in any cell"),
    holds = c(v1, v2, v3), stringsAsFactors = FALSE)
}

if (SUMMARISE) {
  for (tag in c("main", "fresh")) {
    fp <- file.path(sp, paste0("reml_hold_", tag, ".csv"))
    if (!file.exists(fp)) next
    s <- summarise_fits(utils::read.csv(fp, stringsAsFactors = FALSE))
    utils::write.csv(s, file.path(sp, paste0("reml_hold_", tag, "_summary.csv")), row.names = FALSE)
    v <- verdicts(s)
    utils::write.csv(v, file.path(sp, paste0("reml_hold_", tag, "_verdicts.csv")), row.names = FALSE)
    cat("==", tag, "==\n"); print(s, row.names = FALSE, digits = 3); print(v, row.names = FALSE)
  }
  quit(save = "no")
}

## ---- run ----------------------------------------------------------------------
t_all <- Sys.time()
tag <- if (OFFSET == 0L) "main" else "fresh"
ck <- file.path(sp, paste0("checkpoints_", tag))
dir.create(ck, showWarnings = FALSE)
jobs <- expand.grid(rep = seq_len(NREP), cell = if (is.null(ONLY)) cells$cell else ONLY)
jl <- split(jobs, seq_len(nrow(jobs)))
one_ck <- function(job) {
  f <- file.path(ck, sprintf("c%02d_r%04d.rds", job$cell, job$rep))
  if (file.exists(f)) return(readRDS(f))
  r <- one(job); saveRDS(r, f); r
}
if (NCORE > 1L) {
  cl <- parallel::makeCluster(NCORE)
  invisible(parallel::clusterEvalQ(cl, suppressPackageStartupMessages(library(illume))))
  parallel::clusterExport(cl, setdiff(ls(globalenv()), c("cl", "jobs", "jl")), envir = globalenv())
  r <- parallel::parLapplyLB(cl, jl, one_ck, chunk.size = 1L)
  parallel::stopCluster(cl)
} else r <- lapply(jl, one_ck)
res <- do.call(rbind, r)
utils::write.csv(res, file.path(sp, paste0("reml_hold_", tag, ".csv")), row.names = FALSE)
cat("fits:", nrow(res), " minutes:", round(as.numeric(difftime(Sys.time(), t_all, units = "mins"))), "\n")
