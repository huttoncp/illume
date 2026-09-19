## ---------------------------------------------------------------------------
## Turn a study run's raw output into a permanent findings record.
##
## The retention policy keeps raw per-replicate data for the three most recent
## runs of each study and discards the rest (see ../README.md).  What survives
## is whatever this script writes, so it must capture everything a reader would
## otherwise have gone back to the raw data for.  Findings are GENERATED rather
## than hand-written, so pruning never depends on someone having remembered to
## take notes.
##
## Usage: Rscript summarise_run.R <studies_dir> <study> <run>
##        Rscript summarise_run.R <studies_dir>          # every run found
## ---------------------------------------------------------------------------

num <- function(x, d = 3) formatC(x, format = "f", digits = d)

md_table <- function(df, digits = 3) {
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  for (j in seq_along(df)) {
    if (!is.numeric(df[[j]])) next
    v <- df[[j]]
    ## a count column rendered as 240.000 is just noise; decide per column, so
    ## a proportion that happens to reach exactly 1 keeps its decimals
    whole <- all(is.na(v) | v == round(v)) && any(abs(v) >= 2, na.rm = TRUE)
    df[[j]] <- if (whole) formatC(v, format = "d") else num(v, digits)
  }
  out <- c(paste0("| ", paste(names(df), collapse = " | "), " |"),
           paste0("|", paste(rep("---", ncol(df)), collapse = "|"), "|"))
  for (i in seq_len(nrow(df)))
    out <- c(out, paste0("| ", paste(unlist(df[i, ]), collapse = " | "), " |"))
  out
}

read_all <- function(dir, pat) {
  fs <- list.files(dir, pattern = pat, full.names = TRUE)
  if (!length(fs)) return(NULL)
  do.call(rbind, lapply(fs, utils::read.csv, stringsAsFactors = FALSE))
}

## ---- one block of findings per study type ---------------------------------
summarise_study <- function(dir, study) {
  if (study == "coverage") {
    d <- read_all(dir, "^cov_.*csv$"); if (is.null(d)) return(NULL)
    cells <- unique(d$cell)
    s <- do.call(rbind, lapply(cells, function(cl) {
      x <- d[d$cell == cl, ]
      data.frame(cell = cl, N = x$N[1],
                 obs_per_latent = x$obs_per_latent[1],
                 conv_rate = x$conv_rate[1],
                 coverage = mean(x$coverage),
                 worst_coef = min(x$coverage),
                 se_ratio = mean(x$se_ratio),
                 max_abs_bias = max(abs(x$bias)),
                 stringsAsFactors = FALSE)
    }))
    s <- s[order(s$conv_rate, s$cell), ]
    c(paste0("Replicates attempted per cell: ", max(d$n_attempt),
             ". Nominal coverage 0.95."),
      "", md_table(s), "",
      paste0("Coverage range across cells: ", num(min(s$coverage)), " to ",
             num(max(s$coverage)), ". Worst single coefficient: ",
             num(min(s$worst_coef)), " (", s$cell[which.min(s$worst_coef)], ")."),
      paste0("Lowest convergence rate: ", num(min(s$conv_rate)), " (",
             s$cell[which.min(s$conv_rate)], ")."),
      "",
      "Coverage in cells with convergence failures is CONDITIONAL ON CONVERGENCE:",
      "failed fits are excluded, so if failure correlates with extreme estimates",
      "the surviving coverage is optimistic.")

  } else if (study == "power") {
    d <- read_all(dir, "^power_summary.csv$"); if (is.null(d)) return(NULL)
    nul <- d[d$delta == 0, c("cell", "n_wald", "rej_wald", "rej_lrt")]
    names(nul) <- c("cell", "n", "TypeI_Wald", "TypeI_LRT")
    pw <- d[d$delta > 0, c("cell", "delta", "rej_wald")]
    w <- reshape(pw, idvar = "cell", timevar = "delta", direction = "wide")
    names(w) <- sub("^rej_wald\\.", "d=", names(w))
    c(paste0("Replicates per cell/effect size: ", max(d$n_attempt),
             ". Nominal alpha 0.05."),
      "", "**Type I error (delta = 0)** -- the primary result; power is only",
      "meaningful if size is correct.", "", md_table(nul),
      "", "**Power (rejection rate at 5%)**", "", md_table(w),
      "",
      paste0("Type I error range: ", num(min(nul$TypeI_Wald)), " to ",
             num(max(nul$TypeI_Wald)), " (Wald)."))

  } else if (study == "bench") {
    d <- read_all(dir, "^benchmark_freq.csv$"); if (is.null(d)) return(NULL)
    tk <- unique(d$task)
    s <- do.call(rbind, lapply(tk, function(t) {
      x <- d[d$task == t, ]
      g <- function(n) if (n %in% names(x) && any(is.finite(x[[n]])))
        stats::median(x[[n]], na.rm = TRUE) else NA_real_
      data.frame(task = t, N = x$N[1], illume_s = g("illume_s"),
                 lme4_s = g("lme4_s"), glmmTMB_s = g("glmmTMB_s"),
                 nnet_s = g("nnet_s"),
                 max_diff = suppressWarnings(max(c(
                   x$max_absdiff_lme4, x$max_absdiff_glmmTMB,
                   x$max_absdiff_nnet), na.rm = TRUE)),
                 stringsAsFactors = FALSE)
    }))
    s$max_diff[!is.finite(s$max_diff)] <- NA_real_
    c("Agreement is the validity check; timing is context, not a contest.",
      "", md_table(s, 4), "",
      paste0("Largest disagreement with any independent implementation: ",
             format(max(s$max_diff, na.rm = TRUE), digits = 3, scientific = TRUE)))

  } else if (study == "mclogit") {
    d <- read_all(dir, "^mclogit_compare.csv$"); if (is.null(d)) return(NULL)
    rg <- unique(d$regime)
    s <- do.call(rbind, lapply(rg, function(r) {
      x <- d[d$regime == r, ]
      data.frame(regime = r, clusters = x$ncl[1], per_cluster = x$per[1],
                 cover_illume = mean(x$cover_illume),
                 cover_mclogit = mean(x$cover_mclogit),
                 attenuation_illume = x$att_illume[1],
                 attenuation_mclogit = x$att_mclogit[1],
                 bias_illume = mean(abs(x$bias_illume)),
                 bias_mclogit = mean(abs(x$bias_mclogit)),
                 rmse_illume = mean(x$rmse_illume),
                 rmse_mclogit = mean(x$rmse_mclogit),
                 stringsAsFactors = FALSE)
    }))
    c("illume (Laplace) against mclogit::mblogit (PQL). Nominal coverage 0.95.",
      "", md_table(s), "",
      "CAVEAT that must travel with this result: mclogit's RMSE is LOWER,",
      "because shrinkage buys variance reduction. That is a defensible trade",
      "for prediction but not for inference, where the attenuation bias is",
      "invisible in the reported standard errors.")

  } else if (study == "brms") {
    d <- read_all(dir, "^brms_compare.csv$"); if (is.null(d)) return(NULL)
    d$diff_se <- abs(d$illume - d$brms) / d$se_illume
    samp <- d$t_brms[!d$brms_compiled]
    c(paste0("Datasets: ", length(unique(d$rep)), "; coefficients: ",
             length(unique(d$coef)), "; max Rhat ", num(max(d$max_rhat), 4), "."),
      "",
      paste0("- agreement with brms: max ", num(max(d$diff_se)),
             " SE, mean ", num(mean(d$diff_se)), " SE"),
      paste0("- SE ratio illume/brms: ", num(mean(d$se_illume / d$sd_brms))),
      paste0("- illume median ", num(stats::median(d$t_illume), 2), "s; brms ",
             "sampling-only median ",
             num(if (length(samp)) stats::median(samp) else NA_real_, 1), "s"),
      "",
      "These answer different questions (MLE with a Wald interval vs posterior",
      "mean with a credible interval), so exact agreement is not expected; the",
      "residual gap is prior shrinkage.")
  } else NULL
}

## ---- write one run's block into findings/<study>.md ------------------------
write_findings <- function(studies, study, run) {
  dir <- file.path(studies, "runs", study, run)
  body <- summarise_study(dir, study)
  if (is.null(body)) {
    message("  no summarisable output in ", dir); return(invisible(FALSE))
  }
  fp <- file.path(studies, "findings", paste0(study, ".md"))
  dir.create(dirname(fp), showWarnings = FALSE, recursive = TRUE)
  head <- c(paste0("# ", study, " study -- findings log"), "",
            "Newest first. Raw per-replicate data is kept only for the three most",
            "recent runs (see ../README.md); this file is the permanent record and",
            "is generated by scripts/summarise_run.R, never hand-written.", "")
  old <- if (file.exists(fp)) readLines(fp, warn = FALSE) else character(0)
  ## Parse the file into run blocks, replace or add this one, and rewrite in
  ## sorted order.  Inserting at the top instead would make the ordering depend
  ## on the order runs happen to be processed in, which is not something the
  ## caller should have to think about.
  starts <- which(startsWith(old, "## "))
  blocks <- list()
  if (length(starts)) {
    ends <- c(starts[-1] - 1L, length(old))
    for (i in seq_along(starts))
      blocks[[trimws(sub("^## ", "", old[starts[i]]))]] <- old[starts[i]:ends[i]]
  }
  blocks[[run]] <- c(paste0("## ", run), "", body, "")
  labs <- sort(names(blocks), decreasing = TRUE)     # newest first
  writeLines(c(head, unlist(blocks[labs], use.names = FALSE)), fp)
  message("  wrote findings for ", study, " / ", run)
  invisible(TRUE)
}

if (sys.nframe() == 0L) {
  a <- commandArgs(trailingOnly = TRUE)
  studies <- if (length(a) >= 1) a[1] else "."
  if (length(a) >= 3) {
    write_findings(studies, a[2], a[3])
  } else {
    for (s in list.dirs(file.path(studies, "runs"), recursive = FALSE))
      for (r in list.dirs(s, recursive = FALSE))
        write_findings(studies, basename(s), basename(r))
  }
}
