## ---------------------------------------------------------------------------
## Which imputation method, and when.
##
## Known cells are hidden, filled in, and scored against what was really there.
## Four methods, plus a baseline:
##
##   mean      column means -- the floor any method has to beat
##   fcs       chained equations, ilm_impute(method = "fcs")
##   lowrank   regularised iterative PCA, MIPCA (Josse & Husson)
##   lowrank_pa   the same, with the rank from PARALLEL ANALYSIS rather than
##                cross-validation, which is the open question this run settles
##   glrm      a generalized low rank model, a loss per column type
##
## TWO MEASURES, and the second is the one that governs.
##
## RECONSTRUCTION says how close the filled values are. It is what everyone
## reports and it is not what an analysis needs.
##
## COVERAGE says whether the interval around a coefficient estimated AFTER
## imputation still contains the truth 95% of the time. A method can reconstruct
## well and under-cover badly, because filling in a value and then treating it
## as data understates the uncertainty -- which is the whole reason for
## multiple imputation rather than single. Both are reported, and where they
## disagree the coverage decides.
##
## Usage: Rscript imputation_ablation.R <nrep> <ncore> <outdir> [cells]
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(parallel)
})

args   <- commandArgs(trailingOnly = TRUE)
NREP   <- if (length(args) >= 1) as.integer(args[1]) else 10L
NCORE  <- if (length(args) >= 2) as.integer(args[2]) else 4L
OUTDIR <- if (length(args) >= 3) args[3] else "."
ONLY   <- if (length(args) >= 4) strsplit(args[4], ",")[[1]] else NULL
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

M_IMP <- 10L        # imputations per replicate
ALPHA <- 0.05

## The parallel-analysis rank is internal to illumex, where it moved from
## illume in the split. Resolved here, before any replicate: inside one, a
## missing function is an error that run() records as NA, and at 0.0.8.9000
## that emptied lowrank_pa in the four cells with enough complete rows to
## call it -- 0 of 40 fits where 0.0.7.9000 had 40.
invisible(getFromNamespace("ilm_anom_rank", "illumex"))

## ---- the designs -----------------------------------------------------------
## Chosen so the methods disagree. A low-rank structure is where the low-rank
## methods should win; full rank is where imposing one should hurt; p > n is
## where chained equations cannot be fitted at all; and a mixed-type frame is
## where the loss function starts to matter.
cells <- list(
  list(name = "lowrank_n200_p8",   n = 200L, p = 8L,  rank = 3L, miss = 0.20,
       cat = 0L, desc = "n=200 p=8, true rank 3"),
  list(name = "fullrank_n200_p8",  n = 200L, p = 8L,  rank = NA, miss = 0.20,
       cat = 0L, desc = "n=200 p=8, FULL rank -- no structure to find"),
  list(name = "lowrank_n400_p12",  n = 400L, p = 12L, rank = 4L, miss = 0.30,
       cat = 0L, desc = "n=400 p=12, rank 4, 30% missing"),
  list(name = "wide_n60_p80",      n = 60L,  p = 80L, rank = 3L, miss = 0.20,
       cat = 0L, desc = "n=60 p=80, rank 3 -- wider than it is long"),
  list(name = "mixed_n300_p6c2",   n = 300L, p = 6L,  rank = 2L, miss = 0.20,
       cat = 2L, desc = "n=300, 6 numeric + 2 categorical, rank 2"),
  list(name = "mixed_n300_heavy",  n = 300L, p = 6L,  rank = 2L, miss = 0.40,
       cat = 2L, desc = "the same with 40% missing")
)
if (!is.null(ONLY)) cells <- Filter(function(z) z$name %in% ONLY, cells)

## ---- one replicate ---------------------------------------------------------
one_rep <- function(rep_id, cell) {
  suppressPackageStartupMessages(library(illume))
  set.seed(1000L * rep_id + nchar(cell$name))
  n <- cell$n; p <- cell$p; k <- cell$rank

  ## the complete data
  if (is.na(k)) {
    X <- matrix(rnorm(n * p), n, p)                    # full rank: no structure
  } else {
    F <- matrix(rnorm(n * k), n, k)
    L <- matrix(rnorm(k * p), k, p)
    X <- F %*% L + matrix(rnorm(n * p, 0, 0.5), n, p)
  }
  d <- as.data.frame(X)
  names(d) <- paste0("v", seq_len(p))
  ## categorical columns driven by the same factors, so they carry information
  if (cell$cat > 0L) {
    for (j in seq_len(cell$cat)) {
      eta <- if (is.na(k)) rnorm(n) else F[, min(j, ncol(F))]
      d[[paste0("c", j)]] <- factor(ifelse(plogis(1.5 * eta) > runif(n),
                                           "hi", "lo"))
    }
  }
  ## a downstream regression whose coefficient is the thing inference is about
  beta_true <- 0.6
  d$y <- 1 + beta_true * d$v1 + 0.3 * d$v2 + rnorm(n)

  ## hide cells at random among the PREDICTORS; the response stays complete,
  ## which is the usual situation and the one where imputation helps
  pred <- setdiff(names(d), "y")
  full <- d
  mi <- matrix(runif(n * length(pred)) < cell$miss, n, length(pred))
  for (j in seq_along(pred)) d[mi[, j], pred[j]] <- NA

  num_pred <- pred[vapply(full[pred], is.numeric, TRUE)]
  cat_pred <- setdiff(pred, num_pred)
  numi <- mi[, match(num_pred, pred), drop = FALSE]
  cati <- if (length(cat_pred))
    mi[, match(cat_pred, pred), drop = FALSE] else NULL

  ## reconstruction error on the hidden NUMERIC cells, standardised per column
  ## so no one column's units decide the answer
  scl <- vapply(full[num_pred], stats::sd, 0)
  scl[!is.finite(scl) | scl <= 0] <- 1
  rmse <- function(fill) {
    a <- as.matrix(fill[num_pred]); b <- as.matrix(full[num_pred])
    sqrt(mean(((a[numi] - b[numi]) / rep(scl, each = nrow(a))[numi])^2,
              na.rm = TRUE))
  }
  cacc <- function(fill) {
    if (is.null(cati)) return(NA_real_)
    ok <- unlist(lapply(seq_along(cat_pred), function(j) {
      i <- which(cati[, j])
      if (!length(i)) return(logical(0))
      as.character(fill[[cat_pred[j]]][i]) ==
        as.character(full[[cat_pred[j]]][i])
    }))
    if (!length(ok)) NA_real_ else mean(ok)
  }

  fo <- stats::as.formula(paste("y ~", paste(num_pred[1:2], collapse = " + ")))
  ## coverage of beta_true from the pooled fit
  pooled <- function(imps) {
    fits <- lapply(imps, function(z)
      try(suppressWarnings(ilm_model(fo, data = z, family = "gaussian",
                                     verbose = FALSE)), silent = TRUE))
    fits <- Filter(function(z) !inherits(z, "try-error"), fits)
    if (length(fits) < 2L) return(c(NA_real_, NA_real_, NA_real_))
    pl <- try(ilm_mi_pool(fits), silent = TRUE)
    if (inherits(pl, "try-error")) return(c(NA_real_, NA_real_, NA_real_))
    r <- pl[pl$term == num_pred[1L], , drop = FALSE]
    if (!nrow(r)) return(c(NA_real_, NA_real_, NA_real_))
    c(as.numeric(r$lower[1] < beta_true && beta_true < r$upper[1]),
      r$upper[1] - r$lower[1], r$estimate[1])
  }

  out <- list()
  ## ---- mean fill (single, so no pooling: this is the floor) ----------------
  mfill <- d
  for (v in num_pred) mfill[[v]][is.na(mfill[[v]])] <-
    mean(mfill[[v]], na.rm = TRUE)
  for (v in cat_pred) {
    tb <- sort(table(mfill[[v]]), decreasing = TRUE)
    mfill[[v]][is.na(mfill[[v]])] <- names(tb)[1L]
  }
  fm <- try(suppressWarnings(ilm_model(fo, data = mfill, family = "gaussian",
                                       verbose = FALSE)), silent = TRUE)
  cv <- if (inherits(fm, "try-error")) c(NA, NA, NA) else {
    ci <- stats::confint(fm)
    row <- match(num_pred[1L], rownames(ci))
    c(as.numeric(ci[row, 1] < beta_true && beta_true < ci[row, 2]),
      ci[row, 2] - ci[row, 1], unname(stats::coef(fm)[num_pred[1L]]))
  }
  out[["mean"]] <- c(rmse(mfill), cacc(mfill), cv)

  ## ---- the three model-based routes ----------------------------------------
  run <- function(label, expr) {
    im <- try(suppressWarnings(suppressMessages(expr)), silent = TRUE)
    if (inherits(im, "try-error")) return(c(NA, NA, NA, NA, NA))
    imps <- im$imputations
    ## a route that leaves a column untouched cannot be scored on it
    r <- tryCatch(rmse(imps[[1L]]), error = function(e) NA_real_)
    a <- tryCatch(cacc(imps[[1L]]), error = function(e) NA_real_)
    c(r, a, pooled(imps))
  }
  ## Chained equations needs more rows than predictors. Forcing it past that
  ## point does not fail fast -- it attempts a rank-deficient regression per
  ## column per cycle per imputation and takes hours to produce nothing -- so
  ## the cell where it cannot be fitted is recorded as such, which is what
  ## method = "auto" does for a user anyway.
  out[["fcs"]] <- if (length(num_pred) >= n - 2L) rep(NA_real_, 5L) else
    run("fcs", ilm_impute(d, m = M_IMP, method = "fcs", verbose = FALSE,
                          progress = FALSE))
  out[["lowrank"]] <- run("lowrank",
    ilm_impute(d, m = M_IMP, method = "lowrank", verbose = FALSE,
               progress = FALSE))
  ## the same machinery with the rank from parallel analysis instead
  out[["lowrank_pa"]] <- run("lowrank_pa", {
    Xc <- as.matrix(d[num_pred])
    cc <- stats::complete.cases(Xc)
    kk <- if (sum(cc) > 10L)
      illumex:::ilm_anom_rank(scale(Xc[cc, , drop = FALSE])) else 2L
    ilm_impute(d, m = M_IMP, method = "lowrank", ncp = kk, verbose = FALSE,
               progress = FALSE)
  })
  out[["glrm"]] <- run("glrm",
    ilm_impute(d, m = M_IMP, method = "glrm", verbose = FALSE,
               progress = FALSE))

  do.call(rbind, lapply(names(out), function(nm)
    data.frame(cell = cell$name, rep = rep_id, method = nm,
               rmse = out[[nm]][1], cat_acc = out[[nm]][2],
               covered = out[[nm]][3], width = out[[nm]][4],
               estimate = out[[nm]][5], row.names = NULL)))
}

## ---- run -------------------------------------------------------------------
for (cell in cells) {
  cat("== ", cell$name, ": ", cell$desc, "\n", sep = "")
  t0 <- proc.time()[3]
  res <- if (NCORE > 1L) {
    cl <- makeCluster(NCORE)
    on.exit(try(stopCluster(cl), silent = TRUE), add = TRUE)
    clusterExport(cl, c("one_rep", "cell", "M_IMP"), envir = environment())
    r <- parLapply(cl, seq_len(NREP), function(i) one_rep(i, cell))
    stopCluster(cl); on.exit()
    r
  } else lapply(seq_len(NREP), function(i) one_rep(i, cell))
  df <- do.call(rbind, res)
  saveRDS(df, file.path(OUTDIR, paste0(cell$name, ".rds")))
  cat("   ", NREP, " replicates in ", round(proc.time()[3] - t0, 1), "s\n",
      sep = "")
}

## ---- summary ---------------------------------------------------------------
files <- list.files(OUTDIR, pattern = "\\.rds$", full.names = TRUE)
all <- do.call(rbind, lapply(files, readRDS))
agg <- do.call(rbind, lapply(split(all, list(all$cell, all$method),
                                   drop = TRUE), function(z)
  data.frame(cell = z$cell[1L], method = z$method[1L],
             n_ok = sum(!is.na(z$rmse)),
             rmse = mean(z$rmse, na.rm = TRUE),
             cat_acc = mean(z$cat_acc, na.rm = TRUE),
             coverage = mean(z$covered, na.rm = TRUE),
             width = mean(z$width, na.rm = TRUE),
             bias = mean(z$estimate, na.rm = TRUE) - 0.6,
             row.names = NULL)))
agg <- agg[order(agg$cell, agg$method), ]
utils::write.csv(agg, file.path(OUTDIR, "summary.csv"), row.names = FALSE)
print(agg, row.names = FALSE)
