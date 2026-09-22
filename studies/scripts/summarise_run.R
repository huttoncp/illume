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
## A "run" is identified by the PACKAGE VERSION the studies were run against,
## because that is what the results are evidence about.
##
## Usage: Rscript summarise_run.R <studies_dir> <study> <version>
##        Rscript summarise_run.R <studies_dir>          # every run found
## ---------------------------------------------------------------------------

num <- function(x, d = 3) formatC(x, format = "f", digits = d)

## Newest package version first.  A plain string sort is wrong here: it puts
## 0.0.1.9010 before 0.0.1.999, so compare as versions and fall back to string
## order only for labels R cannot parse.
ver_sort <- function(v) {
  ok <- vapply(v, function(x)
    !inherits(try(numeric_version(x), silent = TRUE), "try-error"), TRUE)
  c(v[ok][order(numeric_version(v[ok]), decreasing = TRUE)],
    sort(v[!ok], decreasing = TRUE))
}

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

  } else if (study == "imputation") {
    d <- read_all(dir, "^summary.csv$"); if (is.null(d)) return(NULL)
    cells <- unique(d$cell)
    meths <- c("mean", "fcs", "lowrank", "lowrank_pa", "glrm")
    ## coverage as cells x methods, so "who wins where" is a row operation
    pull <- function(m) vapply(cells, function(cl) {
      z <- d$coverage[d$cell == cl & d$method == m]
      if (length(z)) z[1] else NA_real_
    }, numeric(1))
    M <- vapply(meths, pull, numeric(length(cells)))
    if (is.null(dim(M))) M <- matrix(M, 1, dimnames = list(cells, meths))
    best <- apply(M, 1, function(r) if (all(is.na(r))) NA_character_ else
      paste(meths[which(r == max(r, na.rm = TRUE))], collapse = "/"))

    nrep <- max(d$n_ok, na.rm = TRUE)
    mcse <- sqrt(0.95 * 0.05 / nrep)          # MC error at the nominal rate
    same <- 2 * mcse                          # below this is not a difference

    fcs_ok <- is.finite(M[, "fcs"])
    fcs_best <- sum(fcs_ok & grepl("fcs", best, fixed = TRUE))

    ## the question the study was built to settle: how the rank is chosen
    cv <- M[, "lowrank"]; pa <- M[, "lowrank_pa"]
    ok <- is.finite(cv) & is.finite(pa)
    pa_win <- sum(pa[ok] - cv[ok] > same)
    cv_win <- sum(cv[ok] - pa[ok] > same)

    c(paste0("Five methods over ", length(cells), " designs, ", nrep,
             " replicates each, scored on reconstruction of hidden cells AND"),
      paste0("on coverage of a downstream coefficient after Rubin pooling.",
             " Monte Carlo error on"),
      paste0("coverage is ", num(mcse), ", so a gap under ", num(same, 2),
             " is not a gap."),
      "", md_table(d), "",
      paste0("Best coverage per design: ",
             paste0(cells, " -> ", best, collapse = "; "), "."),
      "",
      paste0("Chained equations could not be fitted in ", sum(!fcs_ok), " of ",
             length(cells), " designs. Where it could, it takes or"),
      paste0("shares the best coverage in ", fcs_best, " of ", sum(fcs_ok),
             ", so it remains the default."),
      "",
      paste0("Rank selection -- cross-validation (lowrank) against parallel",
             " analysis (lowrank_pa):"),
      paste0("parallel analysis better in ", pa_win, " designs,",
             " cross-validation better in ", cv_win, ", the rest"),
      paste0("indistinguishable. ",
             if (pa_win > 0 && cv_win > 0)
               "Neither dominates, so ilm_lowrank_ncp() is UNCHANGED."
             else "One dominates; revisit ilm_lowrank_ncp()."),
      "",
      "CAVEAT that must travel with this result: reconstruction error and",
      "coverage disagree, and only coverage is what an analysis needs. A",
      "method can reconstruct hidden cells as well as the best one and still",
      "cover at half the nominal rate.")

  } else if (study == "messy") {
    ## Every sentence below is computed from the csv, the ones that go against
    ## illume included, so a re-run at a later version says what THAT run
    ## found. The first version of this file was written by hand, and got the
    ## range of mclogit's RMSE advantage wrong (8 to 13% where it was 3 to
    ## 14%) -- which is the reason findings are generated.
    d <- read_all(dir, "^messy_compare.csv$"); if (is.null(d)) return(NULL)
    rg <- unique(d$regime)
    first <- function(x, col) x[[col]][1]
    s <- do.call(rbind, lapply(rg, function(r) {
      x <- d[d$regime == r, ]
      data.frame(regime = r,
                 cover_illume = mean(x$cover_illume),
                 cover_mclogit = mean(x$cover_mclogit),
                 atten_illume = first(x, "att_illume"),
                 atten_mclogit = first(x, "att_mclogit"),
                 bias_illume = mean(abs(x$bias_illume)),
                 bias_mclogit = mean(abs(x$bias_mclogit)),
                 rmse_illume = mean(x$rmse_illume),
                 rmse_mclogit = mean(x$rmse_mclogit),
                 conv_illume = first(x, "conv_illume"),
                 conv_mclogit = first(x, "conv_mclogit"),
                 stringsAsFactors = FALSE)
    }))
    v <- do.call(rbind, lapply(rg, function(r) {
      x <- d[d$regime == r, ]
      data.frame(regime = r, true_re_sd = first(x, "sd_re"),
                 re_sd_illume = first(x, "resd_illume"),
                 re_sd_mclogit = first(x, "resd_mclogit"),
                 secs_illume = first(x, "t_illume"),
                 secs_mclogit = first(x, "t_mclogit"),
                 stringsAsFactors = FALSE)
    }))
    share <- vapply(rg, function(r) first(d[d$regime == r, ], "rare_share"), 0)
    nrep <- max(d$nrep)
    mcse <- sqrt(0.95 * 0.05 / nrep)
    lst <- function(i, a, b, dg = 3)
      paste0(s$regime[i], " (", num(s[[a]][i], dg), " against ",
             num(s[[b]][i], dg), ")", collapse = ", ")

    gap <- s$cover_illume - s$cover_mclogit
    better <- which(gap > 0)
    top <- utils::head(order(-gap), 2L)
    shrink <- which(s$atten_mclogit < 0.95)
    speed <- v$secs_mclogit / v$secs_illume
    inflate <- which(s$atten_illume > 1.1)
    worse_bias <- which(s$bias_illume > s$bias_mclogit)
    less_conv <- which(s$conv_illume < s$conv_mclogit)
    rm_lower <- which(s$rmse_mclogit < s$rmse_illume)
    rel <- (s$rmse_illume - s$rmse_mclogit) / s$rmse_mclogit
    edge <- which(v$true_re_sd < 0.2)

    c(paste0("illume (Laplace) against mclogit::mblogit (PQL) on data that ",
             "misbehave, ", nrep, " replications per regime."),
      paste0("Nominal coverage 0.95; Monte Carlo error about ", num(mcse),
             ". Attenuation is the least-squares slope of estimate on truth:"),
      "1 means none, below 1 is shrinkage and above 1 is inflation.",
      "", md_table(s), "",
      "Variance-component recovery, and median seconds per fit:",
      "", md_table(v), "",
      paste0("Realised rarest-category share: ",
             paste(paste0(rg, " ", num(share)), collapse = ", "), "."),
      "",
      "### What holds",
      "",
      paste0("illume's intervals cover better than mclogit's in ",
             length(better), " of ", length(rg), " regimes; the widest gaps ",
             "are ", lst(top, "cover_illume", "cover_mclogit"), "."),
      paste0("mclogit shrinks fixed effects (attenuation below 0.95) in ",
             length(shrink), " of ", length(rg), " regimes",
             if (length(shrink)) paste0(", from ", num(min(s$atten_mclogit[shrink])),
                                        " to ", num(max(s$atten_mclogit[shrink])))
             else "", "."),
      paste0("illume is ", round(min(speed)), " to ", round(max(speed)),
             " times faster, by median seconds per fit."),
      "",
      "### Findings that go AGAINST illume, and must travel with the rest",
      "",
      if (length(inflate))
        paste0("- illume INFLATES coefficients (attenuation above 1.1) in ",
               lst(inflate, "atten_illume", "atten_mclogit"), ".")
      else "- illume inflates coefficients in no regime.",
      if (length(worse_bias))
        paste0("- illume's mean absolute bias is higher than mclogit's in ",
               lst(worse_bias, "bias_illume", "bias_mclogit"), ".")
      else "- illume's mean absolute bias is no higher than mclogit's anywhere.",
      if (length(less_conv))
        paste0("- illume's fits were usable less often than mclogit's in ",
               lst(less_conv, "conv_illume", "conv_mclogit"),
               "; its coverage there is CONDITIONAL on the usable fits.")
      else "- illume's fits were usable as often as mclogit's everywhere.",
      if (length(edge))
        paste0("- With a true random-effect SD of ",
               paste(num(v$true_re_sd[edge], 2), collapse = ", "),
               " neither method recovers it: illume reports ",
               paste(num(v$re_sd_illume[edge]), collapse = ", "),
               " and mclogit ", paste(num(v$re_sd_mclogit[edge]), collapse = ", "), "."),
      if (length(rm_lower))
        paste0("- mclogit's RMSE is LOWER in ", length(rm_lower), " of ",
               length(rg), " regimes (", paste(s$regime[rm_lower], collapse = ", "),
               "), by ", round(100 * min(rel[rm_lower])), " to ",
               round(100 * max(rel[rm_lower])), "%.")
      else "- mclogit's RMSE is lower in no regime.",
      "",
      "CAVEAT that must travel with this result: the two methods fail",
      "differently rather than one dominating. PQL converges, shrinks, and does",
      "not say so. The Laplace approximation's intervals mean what they claim,",
      "and it declines to answer more often -- loudly, which is the intent, but",
      "a table conditional on the fits that were usable hides that cost.")
  } else NULL
}

## metadata written beside each version's results
run_info <- function(studies, version) {
  fp <- file.path(studies, "runs", version, "RUNINFO.dcf")
  if (!file.exists(fp)) return(NULL)
  tryCatch(as.list(read.dcf(fp)[1, ]), error = function(e) NULL)
}

## ---- write one run's block into findings/<study>.md ------------------------
write_findings <- function(studies, study, version) {
  dir <- file.path(studies, "runs", version, study)
  body <- summarise_study(dir, study)
  if (is.null(body)) {
    message("  no summarisable output in ", dir); return(invisible(FALSE))
  }
  ri <- run_info(studies, version)
  ## a study run on a different day from the rest of its set carries its own
  ## date as `Date.<study>` in the same RUNINFO
  when <- if (is.null(ri)) NULL else
    if (!is.null(ri[[paste0("Date.", study)]])) ri[[paste0("Date.", study)]]
    else ri$Date
  stamp <- if (!is.null(when))
    paste0("Run on ", when, if (!is.null(ri$R)) paste0(", R ", ri$R) else "", ".")
  else NULL
  body <- c(stamp, if (!is.null(stamp)) "", body)
  fp <- file.path(studies, "findings", paste0(study, ".md"))
  dir.create(dirname(fp), showWarnings = FALSE, recursive = TRUE)
  head <- c(paste0("# ", study, " study -- findings log"), "",
            "One section per package version, newest first. Raw per-replicate data",
            "is kept only for the three most recent versions (see ../README.md);",
            "this file is the permanent record and is generated by",
            "scripts/summarise_run.R, never hand-written.", "")
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
  blocks[[version]] <- c(paste0("## ", version), "", body, "")
  labs <- ver_sort(names(blocks))                    # newest version first
  writeLines(c(head, unlist(blocks[labs], use.names = FALSE)), fp)
  message("  wrote findings for ", study, " @ ", version)
  invisible(TRUE)
}

if (sys.nframe() == 0L) {
  a <- commandArgs(trailingOnly = TRUE)
  studies <- if (length(a) >= 1) a[1] else "."
  if (length(a) >= 3) {
    write_findings(studies, a[2], a[3])
  } else {
    for (v in list.dirs(file.path(studies, "runs"), recursive = FALSE))
      for (s in list.dirs(v, recursive = FALSE))
        write_findings(studies, basename(s), basename(v))
  }
}
