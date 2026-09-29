## Finite-sample degrees of freedom in the tables of a gaussian mixed model.
##
## PRE-REGISTERED DESIGN, committed before any code and before any run.
##
## ## The question
## A gaussian mixed model's tables -- summary()'s coefficients, ilm_anova()'s
##   Wald tests, ilm_emmeans()' intervals and ilm_contrast()'s rows -- use a
##   normal or chi-square reference today, and only ilm_trends() uses finite
##   degrees of freedom. With few clusters a between-cluster effect's z test
##   is anti-conservative. Craig ruled (D-DF1 to D-DF5): Satterthwaite's df by
##   default in all four tables, Kenward-Roger on request, ilm_anova() as F
##   with a denominator df (and `statistic = "Chisq"` for the old table), and
##   this study as the evidence for the default. The formulas themselves are
##   established by agreement with lmerTest and pbkrtest in the package tests;
##   the study measures what the reference does to the answers.
## It also answers one open question the tests alone cannot: what the df and
##   the coverage are when a variance component is held at its boundary under
##   rule C (the flat direction held, the rest of the covariance drawn from),
##   where Satterthwaite's gradient runs through a direction the fit holds.
##
## ## Arms and cells (gaussian, REML throughout, as KR requires)
## | Arm | Model | Cells |
## |---|---|---|
## | D1 between-cluster effect | y ~ trt + (1 \| g), trt constant within cluster, half the clusters treated | clusters G in {6, 10, 20} x cluster size m in {5, 20} x ICC in {0.05, 0.3} = 12 |
## | D2 within-cluster slope | y ~ x + (1 + x \| g), x varying within cluster, slope SD 0.3, correlation 0.3 | G in {6, 10, 20} x m in {5, 20} = 6 |
## | D3 held boundary | y ~ trt + (1 \| g), the intercept SD 0 in truth, so about half the fits hold it | G in {6, 10, 20}, m = 5 = 3 |
## | D4 multi-df F | y ~ f + (1 \| g), f a 3-level between-cluster factor, 2 clusters per level at least | G in {6, 12, 21}, m = 5, ICC 0.3 = 3 |
##
## 24 cells. Residual SD 1; ICC sets the intercept SD. The effect of interest
##   has a true value of 0.5 in D1 to D3 (coverage, and the size of a test of
##   the true value), and every level of f is equal in D4 (the size of the F
##   test). 1,000 replicates per cell, 24,000 fits.
##
## ## What is recorded, per fit
## - the estimate and its SE (and KR's adjusted SE);
## - df by each method: z (Inf), Satterthwaite, Kenward-Roger;
## - the 95% interval and the test of the true value by each, and in D4 the F
##   test's DenDF and p by Satterthwaite and KR and the chi-square's p;
## - whether a variance component was held, and the fit's checks;
## - time for each df method.
##
## ## What counts, fixed now
## - "Calibrated": coverage within 2 Monte Carlo SEs of 0.95 (the SE at 1,000
##   replicates is 0.0069, so 0.936 to 0.964); a test's size within 2 MC SEs
##   of 0.05.
## - The default is supported if Satterthwaite is calibrated in every D1 and
##   D2 cell with G >= 10, and closer to 0.95 than z in every cell with G = 6.
##   Where it is not calibrated at G = 6 it is reported, and whether KR is.
## - D3 answers the open question by the numbers: coverage and size in the
##   held fits and in the unheld fits separately, and whether any df is NaN,
##   infinite or below 1. If Satterthwaite is not calibrated in the held fits
##   while KR is, or the df are undefined there, that goes to Craig with a
##   proposal before the default ships for held fits.
## - D4: the F test's size by Satterthwaite and KR, and the chi-square's.
##
## ## Fresh seeds
## The whole design again at 500 replicates per cell with a seed offset, on
##   the same build: the verdicts above must hold there too, and any that does
##   not is reported as such.
##
## ## Build
## A pinned library built from the branch that carries the df change, shared
##   with the spatial validation study (item 115), whose gaussian cells depend
##   on the denominator df; versions and install times recorded before and
##   after the runs.
##
## ## Size
## About 24,000 fits and 12,000 fresh, each a few tenths of a second with KR
##   the largest part at G = 20, m = 20: about 3 to 5 core-hours in all.
##
## Two things the design leaves to the code, fixed here before any run: the
##   intercept SD in D2 is set by an ICC of 0.3, as in D4, and the intercept
##   of every arm is 1. A fit's "held" is `length(f$hessian_held) > 0`.
##
## ## Addendum A1 (2026-09-29, committed before any of its fits): the ML arm
## Every fit above was REML, as KR requires, and the study validated
##   Satterthwaite on REML fits. But ilm_model()'s default is maximum
##   likelihood, and a default gaussian mixed fit gets Satterthwaite df too,
##   computed at the ML variance estimates, which are smaller with few
##   clusters. Measured before this addendum: 6 clusters of 8, a
##   between-cluster slope, df 6.0 and SE 0.359 by ML against 4.0 and 0.440 by
##   REML. No study covered that path. So:
## - Arm ML: the same 24 cells, the MAIN seeds (offset 0), 1,000 replicates
##   per cell, reml = FALSE; everything else as above. Kenward-Roger is
##   derived for REML and refuses an ML fit, so its columns are NA, and the
##   refusal is recorded.
## - The same bands and the same verdict lines, applied to the ML fits. The
##   ML path's default is supported if the first two hold there: Satterthwaite
##   calibrated in every D1 and D2 cell with G >= 10, and closer to 0.95 than z
##   in every cell with G = 6. If they do not hold, that goes to Craig with a
##   proposal (REML by default for gaussian mixed models, item 249 option c;
##   or df from REML variance estimates for an ML fit). Until this arm
##   reports, the NEWS says the validated path is REML.
## - Beside the verdicts, not graded: per cell, ML against REML on the same
##   replicates (the same data, since the seeds are the main run's): the
##   coverage (and D4's size) by z and Satterthwaite, the median
##   Satterthwaite df, the mean SE, and the share of fits held.
## - The build: the main run's pinned library (lib-dfs), unchanged.
## - Size: 24,000 fits without KR, about 2 to 3 core-hours on one core.
##
## Usage: Rscript df_tables.R <nrep> <ncore> <outdir> [offset] [cells|all] [reml|ml]
##        Rscript df_tables.R summarise <outdir> [REML main csv, for the ML comparison]
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
  ## Addendum A1: the estimation method, REML unless "ml"
  REML   <- !(length(args) >= 6 && identical(args[6], "ml"))
} else sp <- if (length(args) >= 2) args[2] else "."
dir.create(sp, showWarnings = FALSE, recursive = TRUE)
suppressPackageStartupMessages(library(illume))

## ---- the cells ----------------------------------------------------------------
g <- function(...) expand.grid(..., stringsAsFactors = FALSE)
cells <- rbind(
  data.frame(arm = "D1", g(G = c(6, 10, 20), m = c(5, 20), icc = c(0.05, 0.3))),
  data.frame(arm = "D2", g(G = c(6, 10, 20), m = c(5, 20)), icc = 0.3),
  data.frame(arm = "D3", G = c(6, 10, 20), m = 5, icc = 0),
  data.frame(arm = "D4", G = c(6, 12, 21), m = 5, icc = 0.3))
cells$cell <- seq_len(nrow(cells))
stopifnot(nrow(cells) == 24L)
TRUTH <- 0.5
LEVEL <- 0.95

## ---- the data -----------------------------------------------------------------
gen <- function(ce, seed) {
  set.seed(seed)
  G <- ce$G; m <- ce$m
  d <- data.frame(g = factor(rep(seq_len(G), each = m)))
  sd0 <- sqrt(ce$icc / (1 - ce$icc))
  b0 <- rnorm(G, 0, sd0)[d$g]
  if (ce$arm %in% c("D1", "D3")) {
    ## half the clusters treated, the treatment constant within a cluster
    d$trt <- rep(rep(0:1, length.out = G)[sample.int(G)], each = m)
    d$y <- 1 + TRUTH * d$trt + b0 + rnorm(nrow(d))
  } else if (ce$arm == "D2") {
    ## a correlated random intercept and slope: slope SD 0.3, correlation 0.3
    S <- matrix(c(sd0^2, 0.3 * sd0 * 0.3, 0.3 * sd0 * 0.3, 0.3^2), 2)
    u <- matrix(rnorm(2 * G), G) %*% chol(S)
    d$x <- rnorm(nrow(d))
    d$y <- 1 + u[d$g, 1] + (TRUTH + u[d$g, 2]) * d$x + rnorm(nrow(d))
  } else {
    ## a 3-level between-cluster factor, G / 3 clusters a level, no effect
    d$f <- factor(rep(rep(c("a", "b", "c"), length.out = G)[sample.int(G)], each = m))
    d$y <- 1 + b0 + rnorm(nrow(d))
  }
  d
}
fml <- list(D1 = y ~ trt + (1 | g), D2 = y ~ x + (1 + x | g),
            D3 = y ~ trt + (1 | g), D4 = y ~ f + (1 | g))
coef_of <- c(D1 = "trt", D2 = "x", D3 = "trt")

## ---- one fit ------------------------------------------------------------------
## Every call is timed and its warnings counted, so a fallback to the normal is
## seen rather than read as a df.
timed <- function(expr) {
  nw <- 0L; msg <- character(0)
  t0 <- proc.time()[["elapsed"]]
  v <- tryCatch(withCallingHandlers(suppressMessages(expr), warning = function(w) {
    nw <<- nw + 1L; msg <<- c(msg, conditionMessage(w)); invokeRestart("muffleWarning")
  }), error = function(e) structure(list(), err = conditionMessage(e)))
  list(v = v, t = proc.time()[["elapsed"]] - t0, nw = nw,
       msg = paste(unique(msg), collapse = " | "), err = attr(v, "err") %||% NA_character_)
}
`%||%` <- function(a, b) if (is.null(a)) b else a

one <- function(job) {
  ce <- cells[cells$cell == job$cell, ]
  seed <- 7919L * ce$cell + job$rep + OFFSET
  d <- gen(ce, seed)
  key <- data.frame(cell = ce$cell, arm = ce$arm, G = ce$G, m = ce$m, icc = ce$icc,
                    rep = job$rep, seed = seed)
  t0 <- proc.time()[["elapsed"]]
  f <- tryCatch(suppressMessages(suppressWarnings(ilm_model(fml[[ce$arm]], data = d,
         family = "gaussian", reml = REML, verbose = FALSE))), error = function(e) NULL)
  t_fit <- proc.time()[["elapsed"]] - t0
  if (is.null(f)) return(cbind(key, ok = FALSE))
  ck <- f$checks
  key$ok <- TRUE
  key$held <- length(f$hessian_held) > 0L
  key$nonconv <- any(ck$status[ck$check %in% c("gradient", "optimizer")] == "FAIL")
  key$t_fit <- t_fit
  meths <- c(z = "asymptotic", s = "satterthwaite", kr = "kenward-roger")
  if (ce$arm != "D4") {
    nm <- coef_of[[ce$arm]]
    for (k in names(meths)) {
      r <- timed(ilm_coef_table(f, df = meths[[k]]))
      ct <- r$v
      row <- if (is.data.frame(ct) && nm %in% rownames(ct)) ct[nm, ] else NULL
      est <- if (is.null(row)) NA_real_ else row$Estimate
      se  <- if (is.null(row)) NA_real_ else row[["Std. Error"]]
      df  <- if (is.null(row)) NA_real_ else if (k == "z") Inf else row$df
      fb  <- is.data.frame(ct) && !is.null(attr(ct, "fallback")) &&
             match(nm, rownames(ct)) %in% attr(ct, "fallback")$rows
      q <- if (is.finite(df) && df > 0) qt(1 - (1 - LEVEL) / 2, df) else qnorm(1 - (1 - LEVEL) / 2)
      tv <- (est - TRUTH) / se
      p <- if (is.finite(df) && df > 0) 2 * pt(-abs(tv), df) else 2 * pnorm(-abs(tv))
      key[[paste0("est_", k)]] <- est
      key[[paste0("se_", k)]] <- se
      key[[paste0("df_", k)]] <- df
      key[[paste0("cover_", k)]] <- abs(est - TRUTH) <= q * se
      key[[paste0("p_", k)]] <- p
      key[[paste0("fallback_", k)]] <- fb
      key[[paste0("nwarn_", k)]] <- r$nw
      key[[paste0("warn_", k)]] <- r$msg
      key[[paste0("err_", k)]] <- r$err
      key[[paste0("t_", k)]] <- r$t
    }
  } else {
    ## the F test of f, by each df, and the Wald chi-square
    for (k in c("s", "kr")) {
      r <- timed(ilm_anova(f, df = meths[[k]]))
      a <- r$v
      ok <- is.data.frame(a) && "f" %in% rownames(a)
      key[[paste0("numdf_", k)]] <- if (ok) a["f", "NumDF"] else NA_real_
      key[[paste0("dendf_", k)]] <- if (ok) a["f", "DenDF"] else NA_real_
      key[[paste0("F_", k)]] <- if (ok) a["f", "F value"] else NA_real_
      key[[paste0("p_", k)]] <- if (ok) a["f", "Pr(>F)"] else NA_real_
      key[[paste0("nwarn_", k)]] <- r$nw
      key[[paste0("warn_", k)]] <- r$msg
      key[[paste0("err_", k)]] <- r$err
      key[[paste0("t_", k)]] <- r$t
    }
    r <- timed(ilm_anova(f, statistic = "Chisq"))
    a <- r$v
    ok <- is.data.frame(a) && "f" %in% rownames(a)
    key$chisq_z <- if (ok) a["f", "Chisq"] else NA_real_
    key$p_z <- if (ok) a["f", "Pr(>Chisq)"] else NA_real_
    key$nwarn_z <- r$nw; key$err_z <- r$err; key$t_z <- r$t
  }
  key
}

## ---- summaries ----------------------------------------------------------------
## Coverage (D1 to D3) and the size of the test (D4) by method, with the Monte
## Carlo SE and the pre-registered verdict: calibrated when within 2 MC SEs
## of the target. D3 is also split by whether the fit held a component.
mcse <- function(p, n) sqrt(p * (1 - p) / n)
summarise_fits <- function(x) {
  x <- x[x$ok %in% TRUE, ]
  rows <- list()
  add <- function(sub, cell, split) {
    if (!nrow(sub)) return(NULL)
    ce <- sub[1, c("cell", "arm", "G", "m", "icc")]
    for (k in c("z", "s", "kr")) {
      if (ce$arm != "D4") {
        v <- sub[[paste0("cover_", k)]]; target <- LEVEL; what <- "coverage"
      } else {
        v <- sub[[paste0("p_", k)]] < 0.05; target <- 0.05; what <- "size"
      }
      n <- sum(!is.na(v)); est <- mean(v, na.rm = TRUE)
      dfv <- if (k == "z") rep(Inf, nrow(sub)) else
        sub[[if (ce$arm == "D4") paste0("dendf_", k) else paste0("df_", k)]]
      rows[[length(rows) + 1L]] <<- data.frame(ce, split = split, method = k,
        what = what, n = n, estimate = est, mcse = mcse(target, n),
        calibrated = abs(est - target) <= 2 * mcse(target, n),
        gap = abs(est - target),
        df_median = if (k == "z") Inf else stats::median(dfv, na.rm = TRUE),
        df_bad = if (k == "z") 0L else sum(!is.finite(dfv) | dfv < 1, na.rm = TRUE) + sum(is.na(dfv)),
        errors = sum(!is.na(sub[[paste0("err_", k)]])),
        warned = sum(sub[[paste0("nwarn_", k)]] > 0, na.rm = TRUE),
        n_held = sum(sub$held), n_nonconv = sum(sub$nonconv),
        time_mean = mean(sub[[paste0("t_", k)]], na.rm = TRUE))
    }
  }
  for (cl in sort(unique(x$cell))) {
    sub <- x[x$cell == cl, ]
    add(sub, cl, "all")
    if (sub$arm[1] == "D3") {
      add(sub[sub$held, ], cl, "held")
      add(sub[!sub$held, ], cl, "unheld")
    }
  }
  do.call(rbind, rows)
}
## The pre-registered verdicts, one line each; NA where no fit bears on one
## (no held fits in D3, say), never a vacuous TRUE.
verdicts <- function(s) {
  all <- function(x) if (length(x)) base::all(x) else NA
  any <- function(x) if (length(x)) base::any(x) else NA
  a <- s[s$split == "all", ]
  at <- function(arm, G, meth) a[a$arm %in% arm & a$G %in% G & a$method == meth, ]
  s_big <- at(c("D1", "D2"), c(10, 20), "s")
  z6 <- at(c("D1", "D2"), 6, "z"); s6 <- at(c("D1", "D2"), 6, "s")
  closer <- merge(z6[, c("cell", "gap")], s6[, c("cell", "gap")], by = "cell",
                  suffixes = c("_z", "_s"))
  d3 <- s[s$arm == "D3" & s$split != "all", ]
  data.frame(
    verdict = c("Satterthwaite calibrated in every D1/D2 cell with G >= 10",
                "Satterthwaite closer to 0.95 than z in every D1/D2 cell with G = 6",
                "D3: Satterthwaite calibrated in the held fits",
                "D3: KR calibrated in the held fits",
                "D3: every Satterthwaite df defined (finite and at least 1)",
                "D4: Satterthwaite's F size calibrated in every cell",
                "D4: KR's F size calibrated in every cell",
                "D4: chi-square size calibrated in every cell"),
    holds = c(all(s_big$calibrated), all(closer$gap_s < closer$gap_z),
              all(d3$calibrated[d3$split == "held" & d3$method == "s"]),
              all(d3$calibrated[d3$split == "held" & d3$method == "kr"]),
              !any(a$df_bad[a$arm == "D3" & a$method == "s"] > 0),
              all(at("D4", c(6, 12, 21), "s")$calibrated),
              all(at("D4", c(6, 12, 21), "kr")$calibrated),
              all(at("D4", c(6, 12, 21), "z")$calibrated)))
}

## Addendum A1: ML against REML on the same replicates, per cell
ml_vs_reml <- function(ml, reml) {
  ml <- ml[ml$ok %in% TRUE, ]; reml <- reml[reml$ok %in% TRUE, ]
  m <- merge(reml, ml, by = c("cell", "rep"), suffixes = c("_reml", "_ml"))
  do.call(rbind, lapply(sort(unique(m$cell)), function(cl) {
    z <- m[m$cell == cl, ]; arm <- z$arm_reml[1]
    cov <- function(k, sfx) if (arm == "D4") mean(z[[paste0("p_", k, "_", sfx)]] < 0.05, na.rm = TRUE)
                            else mean(z[[paste0("cover_", k, "_", sfx)]], na.rm = TRUE)
    dfs <- function(sfx) stats::median(z[[paste0(if (arm == "D4") "dendf_s_" else "df_s_", sfx)]], na.rm = TRUE)
    se <- function(sfx) if (arm == "D4") NA_real_ else mean(z[[paste0("se_s_", sfx)]], na.rm = TRUE)
    data.frame(cell = cl, arm = arm, G = z$G_reml[1], m = z$m_reml[1], icc = z$icc_reml[1],
               pairs = nrow(z), what = if (arm == "D4") "size" else "coverage",
               z_reml = cov("z", "reml"), z_ml = cov("z", "ml"),
               s_reml = cov("s", "reml"), s_ml = cov("s", "ml"),
               df_s_reml = dfs("reml"), df_s_ml = dfs("ml"),
               se_reml = se("reml"), se_ml = se("ml"),
               held_reml = mean(z$held_reml), held_ml = mean(z$held_ml))
  }))
}

if (SUMMARISE) {
  for (tag in c("main", "fresh", "ml")) {
    fp <- file.path(sp, paste0("df_tables_", tag, ".csv"))
    if (!file.exists(fp)) next
    s <- summarise_fits(utils::read.csv(fp))
    utils::write.csv(s, file.path(sp, paste0("df_tables_", tag, "_summary.csv")), row.names = FALSE)
    v <- verdicts(s)
    utils::write.csv(v, file.path(sp, paste0("df_tables_", tag, "_verdicts.csv")), row.names = FALSE)
    cat("==", tag, "==\n"); print(v, row.names = FALSE)
  }
  fm <- file.path(sp, "df_tables_ml.csv")
  fr <- if (length(args) >= 3) args[3] else file.path(sp, "df_tables_main.csv")
  if (file.exists(fm) && file.exists(fr)) {
    cmp <- ml_vs_reml(utils::read.csv(fm), utils::read.csv(fr))
    utils::write.csv(cmp, file.path(sp, "df_tables_ml_vs_reml.csv"), row.names = FALSE)
    cat("== ML against REML, same replicates ==\n"); print(cmp, row.names = FALSE, digits = 3)
  }
  quit(save = "no")
}

## ---- run ----------------------------------------------------------------------
t_all <- Sys.time()
use <- if (is.null(ONLY)) cells else cells[cells$cell %in% ONLY, ]
jobs <- merge(use[, "cell", drop = FALSE], data.frame(rep = seq_len(NREP)))
jl <- split(jobs, seq_len(nrow(jobs)))
tag <- if (OFFSET) "fresh" else if (!is.null(ONLY)) "smoke" else "main"
if (!REML) tag <- paste0("ml", if (identical(tag, "main")) "" else paste0("_", tag))
## CHECKPOINTS (an I/O change, no effect on results): each fit's result is
## saved as it finishes, and a run started again reads the fits already done
## instead of refitting them. Every fit sets its own seed from its cell and
## replicate, so a resumed run gives the same numbers as one uninterrupted.
ck <- file.path(sp, paste0("checkpoints_", tag))
dir.create(ck, showWarnings = FALSE, recursive = TRUE)
one_ck <- function(job) {
  fp <- file.path(ck, sprintf("cell%02d_rep%04d.rds", job$cell, job$rep))
  if (file.exists(fp)) return(readRDS(fp))
  r <- one(job)
  saveRDS(r, paste0(fp, ".part")); file.rename(paste0(fp, ".part"), fp)
  r
}
res <- if (NCORE > 1L) {
  cl <- parallel::makeCluster(NCORE)
  invisible(parallel::clusterEvalQ(cl, suppressPackageStartupMessages(library(illume))))
  parallel::clusterExport(cl, setdiff(ls(globalenv()), c("cl", "jobs", "jl")), envir = globalenv())
  ## one job at a time, so a resume's unfinished fits do not all fall to one worker
  r <- parallel::parLapplyLB(cl, jl, one_ck, chunk.size = 1L)
  parallel::stopCluster(cl)
  r
} else lapply(jl, one_ck)
## failed fits carry fewer columns; each is filled out with NA
nm <- unique(unlist(lapply(res, names)))
out <- do.call(rbind, lapply(res, function(d) { for (n in setdiff(nm, names(d))) d[[n]] <- NA; d[nm] }))
utils::write.csv(out, file.path(sp, paste0("df_tables_", tag, ".csv")), row.names = FALSE)
cat("illume", format(utils::packageVersion("illume")), "from", find.package("illume"), "\n")
cat("fits:", nrow(jobs), " minutes:",
    round(as.numeric(difftime(Sys.time(), t_all, units = "mins")), 1), "\n")
