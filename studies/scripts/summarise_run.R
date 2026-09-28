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
  ds <- lapply(fs, utils::read.csv, stringsAsFactors = FALSE)
  ## a cell re-run by a later version of its script can carry columns the
  ## others do not; those are NA where they were not recorded
  nm <- unique(unlist(lapply(ds, names)))
  do.call(rbind, lapply(ds, function(x) { x[setdiff(nm, names(x))] <- NA; x[nm] }))
}

## The fits used with a covariance direction held at its boundary, which the
## package calls usable and the coverage study counts from 0.0.8.9000 on,
## with their coverage on their own. NULL where the run did not record them.
held_line <- function(d, order) {
  if (!"n_held" %in% names(d)) return(NULL)
  h <- do.call(rbind, lapply(order, function(cl) {
    x <- d[d$cell == cl, ]
    if (is.na(x$n_held[1]) || x$n_held[1] == 0) return(NULL)
    data.frame(cell = cl, n = x$n_held[1], used = x$n_used[1],
               cover = mean(x$coverage_held))
  }))
  if (is.null(h)) return(NULL)
  paste0("Fits with a covariance direction held at its boundary have usable ",
         "fixed-effect standard errors, as the package says, and count as ",
         "converged: ", paste(sprintf("%s %d of %d used (coverage %s)", h$cell,
                                      h$n, h$used, num(h$cover)), collapse = "; "),
         ".")
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
      held_line(d, s$cell),
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

  } else if (study == "re_sd_limit") {
    ## Where a random effect's SD short of zero is held: below a line, and
    ## the objective flat when its log SD is pushed 3 lower. phase1 placed
    ## the line; verify is the same design on fresh seeds, on the built rule.
    one <- function(file) {
      d <- read_all(dir, paste0("^", file, "$")); if (is.null(d)) return(NULL)
      d <- d[d$ok & !is.na(d$ll_drop) & !is.na(d$push), ]
      d$at <- d$ll_drop >= d$ll - 1e-3 & d$push <= 1e-3
      tab <- do.call(rbind, lapply(c(5e-3, 1e-3), function(tol)
        do.call(rbind, lapply(c(1e-3, 1e-2, 0.05, 0.1, 0.2, Inf), function(line) {
          rl <- d$sd_hat < line & d$push <= tol
          data.frame(tolerance = format(tol), line = format(line),
                     holds = sum(rl), false_holds = sum(rl & !d$at),
                     missed = sum(!rl & d$at), stringsAsFactors = FALSE)
        }))))
      list(d = d, tab = tab)
    }
    p1 <- one("re_sd_limit_phase1.csv"); v <- one("re_sd_limit_verify.csv")
    if (is.null(p1) && is.null(v)) return(NULL)
    lines <- c(paste0("A random effect's SD short of zero. A term is AT ITS BOUNDARY, by ",
                      "labels from outside the rule, when the model without it fits as ",
                      "well (log-likelihood within 1e-3) and the objective rises by at ",
                      "most 1e-3 with its log SD pushed 3 lower. The rule holds an SD ",
                      "below the line when the objective moves by at most the tolerance ",
                      "there; the pre-registered tolerance was the dispersion's 5e-3."))
    for (nm in c("phase1", "verify")) {
      x <- if (nm == "phase1") p1 else v
      if (is.null(x)) next
      lines <- c(lines, "", paste0("**", if (nm == "phase1") "Placing the line" else
                                     "Fresh seeds, on the built rule", "** (",
                                   nrow(x$d), " fits, ", sum(x$d$at),
                                   " at the boundary):"), "", md_table(x$tab))
    }
    if (!is.null(v) && "draws_finite" %in% names(v$d)) {
      h <- grepl("(^|,)g(,|$)", v$d$held) | grepl("(^|,)g(,|$)", v$d$boundary)
      r <- v$d$se_x[h] / v$d$se_x_drop[h]
      lines <- c(lines, "", paste0(
        "On the built rule (line 0.1, tolerance 1e-3): ", sum(h), " fits held, ",
        sum(h & v$d$at), " of them at the boundary by the labels; the standard ",
        "error of the slope runs ", num(min(r, na.rm = TRUE)), " to ",
        num(max(r, na.rm = TRUE)), " times the refit without the term (median ",
        num(stats::median(r, na.rm = TRUE)), ")."))
    }
    c(lines, "",
      "CAVEAT that must travel with this result: the tolerance of 1e-3 was chosen",
      "after the pre-registered run, which used 5e-3 and showed that no line above",
      "1e-2 avoided new false holds with it; the fresh-seed run is the check on",
      "that choice. The false holds that remain are those the old 1e-3 line",
      "already made, where the refit without the term landed at a worse optimum.")

  } else if (study == "dispersion_limit") {
    ## Where a negative binomial's k or a beta's phi is at its unbounded
    ## limit, and whether the package's rule -- 1 / sqrt(dispersion) below
    ## 1e-2 and the objective flat beyond it -- holds those fits and no
    ## others. phase1 is the run before the hold, verify the one after it.
    d1 <- read_all(dir, "^dispersion_limit_phase1.csv$")
    d2 <- read_all(dir, "^dispersion_limit_verify.csv$")
    d <- if (!is.null(d2)) d2 else d1
    if (is.null(d)) return(NULL)
    d <- d[d$ok & !is.na(d$limit_ll) & !is.na(d$push_down), ]
    d$at <- d$limit_ll & d$push_down
    d$rule <- d$inv_sqrt < 1e-2 & d$push <= 1e-3
    tab <- do.call(rbind, lapply(split(d, list(d$family, d$design, d$disp),
                                       drop = TRUE), function(x)
      data.frame(family = x$family[1], design = x$design[1],
                 true_disp = format(x$disp[1]), fits = nrow(x), at_limit = sum(x$at),
                 rule_flags = sum(x$rule),
                 held = if ("held" %in% names(x)) sum(grepl("dispersion", x$held))
                        else NA_integer_,
                 stringsAsFactors = FALSE)))
    tab <- tab[order(tab$family, tab$design, as.numeric(tab$true_disp)), ]
    fam_line <- function(f) {
      x <- d[d$family == f, ]
      paste0(f, ": ", sum(x$at), " of ", nrow(x), " fits at the limit; the rule ",
             "flags ", sum(x$rule), ", ", sum(x$rule & !x$at), " of them not at it, ",
             "and misses ", sum(!x$rule & x$at), ". At the limit 1 / sqrt(dispersion) ",
             "had a 99th percentile of ", signif(stats::quantile(x$inv_sqrt[x$at], .99), 3),
             "; not at it, a minimum of ", signif(min(x$inv_sqrt[!x$at]), 3), ".")
    }
    nb <- d[d$family == "nbinom", ]
    h <- grepl("dispersion", nb$held)
    c(paste0("A dispersion at its unbounded limit. Fits are labelled from outside ",
             "the rule: at the limit when a model at the limit fits as well (for ",
             "the negative binomial, the Poisson refit's log-likelihood within 1e-3) ",
             "AND the objective does not rise by more than 1e-3 when the log ",
             "dispersion is pushed 3 further. The rule is the package's: ",
             "1 / sqrt(dispersion) below 1e-2 and the objective rising by at most ",
             "1e-3 there (5e-3 before the random-effect study; see re_sd_limit)."),
      "", md_table(tab), "",
      fam_line("nbinom"), fam_line("beta"),
      paste0("The line alone, without the flatness test, flags ",
             sum(d$family == "beta" & d$inv_sqrt < 1e-2 & !d$at), " beta fits ",
             "not at the limit: curved optima at phi up to ",
             signif(max(exp(d$logdisp[d$family == "beta" & d$inv_sqrt < 1e-2 & !d$at])), 2),
             ", where a correlation over time interpolates nearly noiseless data."),
      if (!is.null(d2) && any(h))
        paste0("After the hold: ", sum(h), " negative binomial fits held; their ",
               "fixed-effect standard error of x runs ",
               num(min(nb$se_ratio_x[h], na.rm = TRUE)), " to ",
               num(max(nb$se_ratio_x[h], na.rm = TRUE)), " times the Poisson refit's ",
               "(median ", num(stats::median(nb$se_ratio_x[h], na.rm = TRUE)), "), and ",
               sum(d$draws_finite[grepl("dispersion", d$held)], na.rm = TRUE), " of ",
               sum(grepl("dispersion", d$held)), " held fits have draws of the log ",
               "dispersion within 50."),
      {
        ## the gaussian residual SD at zero, from dispersion_limit_gaussian.R
        g <- read_all(dir, "^dispersion_limit_gaussian_verify.csv$")
        if (is.null(g)) g <- read_all(dir, "^dispersion_limit_gaussian_phase1.csv$")
        if (is.null(g)) NULL else {
          g <- g[g$ok & !is.na(g$at), ]
          g$rule <- g$rel < 1e-3 & g$push <= 1e-3   # the rule as it stands
          hg <- if ("held" %in% names(g)) grepl("dispersion", g$held) else rep(FALSE, nrow(g))
          c("",
            paste0("Gaussian residual SD at zero (", nrow(g), " fits; AR(1) at one ",
                   "observation per cell with noise SD 0.5, 0.2 and 0.05, and two ",
                   "random-intercept controls): at the limit when the objective ",
                   "does not rise by more than 1e-3 as log sigma is pushed 3 lower. ",
                   sum(g$at), " fits were at the limit; the rule -- sigma below 1e-3 ",
                   "of the response's SD, and flat within 1e-3 -- flags ", sum(g$rule),
                   ", ", sum(g$rule & !g$at), " of them not at it, and misses ",
                   sum(!g$rule & g$at), ", whose sigma was ",
                   signif(min(g$rel[!g$rule & g$at]), 2), " to ",
                   signif(max(g$rel[!g$rule & g$at]), 2), " of the response's SD. ",
                   "Not at the limit, sigma was at least ", signif(min(g$rel[!g$at]), 2),
                   " of it, so no line separates those few: the flatness test does. ",
                   if (any(hg)) paste0("After the hold, ", sum(hg), " fits held, all ",
                                       "with draws of log sigma within 50: ",
                                       sum(g$draws_finite[hg], na.rm = TRUE), " of ",
                                       sum(hg), ".") else ""))
        }
      },
      "",
      "CAVEAT that must travel with this result: a fit at the limit that the",
      "optimiser left unconverged -- the objective still falling steeply as the",
      "dispersion grows -- is not held; it stays a failure, as it should.")

  } else if (study == "boundary_se") {
    ## Which standard errors the fixed effects get when a random-effect
    ## covariance sits at its boundary. Every dataset is fitted once and its
    ## covariance computed under each rule at the same optimum, so the rules
    ## differ in their intervals and nothing else. B is the fit's own: the
    ## rule the package used at that version. C3 is the one adopted at
    ## 0.0.7.9000 -- only the flat directions held, at the 1e-3 threshold.
    cv <- read_all(dir, "^boundary_se_cover.csv$"); if (is.null(cv)) return(NULL)
    rt <- read_all(dir, "^boundary_se_ratio.csv$")
    pr <- read_all(dir, "^boundary_se_paired.csv$")
    at <- unique(cv$regime[cv$boundary > 0])
    pick <- function(rule, col) vapply(at, function(r)
      cv[[col]][cv$regime == r & cv$rule == rule][1], 0)
    tab <- data.frame(regime = at, fits = pick("B", "n"),
                      at_boundary = pick("B", "boundary"),
                      cover_A = pick("A", "cover_boundary"),
                      cover_B = pick("B", "cover_boundary"),
                      cover_C3 = pick("C3", "cover_boundary"),
                      se_sd_A = pick("A", "se_sd"), se_sd_C3 = pick("C3", "se_sd"),
                      stringsAsFactors = FALSE)
    nb <- sum(tab$at_boundary)
    best <- sum(abs(tab$cover_C3 - 0.95) <= abs(tab$cover_A - 0.95))
    low <- tab$regime[which.min(tab$cover_C3)]
    same_B <- !is.null(rt) && all(abs(c(rt$B_C3_min, rt$B_C3_max) - 1) < 1e-3)
    ref <- if (is.null(rt)) NULL else rt[rt$refits > 0, ]
    c(paste0("Standard errors of the fixed effects when a random-effect covariance ",
             "sits at its boundary, under three rules computed at the same optimum:"),
      "A holds the whole boundary term; B is the fit's own, the rule the package",
      "used at this version; C3 holds only the directions with curvature below",
      "1e-3 of the term's largest. Coverage of 95% intervals among the fits at a",
      "boundary; se_sd is mean SE over the SD of the estimates across all usable fits.",
      "", md_table(tab), "",
      paste0(nb, " of ", sum(tab$fits), " fits ended at a boundary, in the ",
             nrow(tab), " of ", length(unique(cv$regime)), " regimes where any did. ",
             "C3's coverage is at least as close to 0.95 as A's in ", best, " of those ",
             nrow(tab), "; its lowest is ", num(min(tab$cover_C3)), " (", low, ")."),
      if (!is.null(rt))
        paste0("Against C3, A's standard errors run a median ",
               num(min(rt$A_C3_median)), " to ", num(max(rt$A_C3_median)),
               " times as large by regime, and as little as ", num(min(rt$A_C3_min)),
               " in the worst fit."),
      if (!is.null(rt))
        if (same_B) "B equals C3 in every boundary fit, to three decimals: the rule in use is the rule measured."
        else paste0("B, the rule in use at this version, runs ", num(min(rt$B_C3_min)),
                    " to ", num(max(rt$B_C3_max)), " times C3's."),
      if (!is.null(ref) && nrow(ref))
        paste0("Where the reduced model could be refitted (", sum(ref$refits),
               " fits), C3 matched it to a ratio of ", num(min(ref$C3_refit_min)),
               " to ", num(max(ref$C3_refit_max)), "."),
      if (!is.null(pr))
        paste0("Paired, C3 covered ", sum(pr$C3_not_A), " intervals that A missed, ",
               "and A covered ", sum(pr$A_not_C3), " that C3 missed."),
      if (!is.null(rt))
        paste0("Every boundary fit had a gap of at least ", num(min(rt$gap_decades_min), 1),
               " orders of magnitude between its flat and curved directions; ",
               "thresholds of 1e-3 and 1e-4 chose differently in ",
               sum(rt$threshold_disagree), " of ", nb, " fits."),
      "",
      "CAVEAT that must travel with this result: a regime can under-cover under",
      "every rule, and when it does the cause is not the choice of standard",
      "error -- in `combined` it is the rare category's bias. AR(1) and CAR(1)",
      "terms are held the same way but are not part of this study.")

  } else if (study == "variance_draws") {
    ## Draws of a poorly determined variance (item 78, Arm B): the flag on the
    ## SE of a log SD, its cut-off by the pre-registered rule, and the two
    ## remedies' coverage, forecast scores and time.
    par <- read_all(dir, "^variance_draws_main_par.csv$")
    if (is.null(par)) return(NULL)
    fc <- read_all(dir, "^variance_draws_main_fc.csv$")
    tm <- read_all(dir, "^variance_draws_main_time.csv$")
    fpar <- read_all(dir, "^variance_draws_fresh_par.csv$")
    fxp <- read_all(dir, "^variance_draws_fixtures_par.csv$")
    fxf <- read_all(dir, "^variance_draws_fixtures_fc.csv$")
    cutsv <- c(0.5, 0.75, 1, 1.5, 2, 3)
    fails <- c("unbounded", "wald_too_wide")
    ## the cells, as the script builds them, for the arm and family of each row
    cl <- rbind(
      data.frame(arm = "B1", expand.grid(family = c("gaussian", "poisson", "nbinom"),
                 Tn = c(12, 24, 48), rho = c(0.3, 0.7), lat_sd = c(0.3, 0.8),
                 stringsAsFactors = FALSE)[, "family", drop = FALSE]),
      data.frame(arm = "B2", expand.grid(family = c("gaussian", "poisson", "nbinom", "binomial"),
                 Tn = c(12, 24), ri_sd = c(0.2, 0.5),
                 stringsAsFactors = FALSE)[, "family", drop = FALSE]),
      data.frame(arm = "B3", family = rep(c("gaussian", "poisson", "nbinom", "binomial"), 6)))
    cl$cell <- seq_len(nrow(cl))
    prep <- function(p) {
      p <- merge(p[p$ok %in% TRUE & !is.na(p$param), ], cl, by = "cell")
      p$fails <- p$label %in% fails
      p
    }
    par <- prep(par)
    mcse <- function(x) sqrt(mean(x) * (1 - mean(x)) / length(x))
    pc <- function(x) paste0(num(mean(x), 3), " (", num(mcse(x), 3), ")")
    lines <- c(paste0(
      "Draws of a poorly determined variance. ", length(unique(paste(par$cell, par$rep))),
      " fits and ", nrow(par), " variance parameters; s is the SE of a parameter's log ",
      "scale from vcov(fit, full = TRUE), d the SD of its 400 joint draws (the other ",
      "package's signal). The label is the profile likelihood's, from outside the flag. ",
      "Parameters the build holds are counted here and left out of every rate below. ",
      "Monte Carlo SEs in brackets."))
    ## held, and the labels, per arm and family
    lt <- do.call(rbind, lapply(split(par, list(par$arm, par$family, par$param), drop = TRUE),
      function(x) data.frame(arm = x$arm[1], family = x$family[1], param = x$param[1],
        n = nrow(x), held = sum(x$held), ok = sum(x$label %in% "ok"),
        wald_too_wide = sum(x$label %in% "wald_too_wide"),
        unbounded = sum(x$label %in% "unbounded"),
        profile_failed = sum(x$label %in% "profile_failed"),
        ridge = sum(x$ridge %in% TRUE), stringsAsFactors = FALSE)))
    lt <- lt[order(lt$arm, lt$family, lt$param), ]
    lines <- c(lines, "", "**Holds and labels** (parameters):", "", md_table(lt))
    u <- par[!par$held & par$label %in% c("ok", fails) & is.finite(par$s), ]
    ## s against d
    okd <- is.finite(u$d) & u$d > 0 & u$s > 0
    lines <- c(lines, "", paste0(
      "**s against d.** Spearman correlation ", num(stats::cor(u$s[okd], u$d[okd],
      method = "spearman"), 3), " over ", sum(okd), " unheld parameters."))
    ## the cut-off, by the rule: the smallest candidate flagging at most 5% of
    ## the parameters not labelled approx_fails
    ct <- do.call(rbind, lapply(cutsv, function(k) {
      fl <- u$s > k
      data.frame(cut = k, specificity = mean(!fl[!u$fails]), sensitivity = mean(fl[u$fails]),
                 sens_unbounded = mean(fl[u$label == "unbounded"]),
                 sens_wald_too_wide = mean(fl[u$label == "wald_too_wide"]),
                 agree_d = mean(fl == (u$d > k), na.rm = TRUE))
    }))
    pick <- ct$cut[ct$specificity >= 0.95]
    chosen <- if (length(pick)) min(pick) else NA_real_
    lines <- c(lines, "", paste0("**The cut-off** (", nrow(u), " unheld parameters, ",
      sum(u$fails), " approx_fails; agree_d is the share where d gives the same flag):"),
      "", md_table(ct), "",
      if (is.na(chosen)) paste0("No candidate reaches specificity 0.95: by the rule the ",
                               "flag is NOT VIABLE, and remedy 1 fails with it.")
      else paste0("By the rule the cut-off is **", chosen, "**.",
        if (chosen == min(cutsv)) paste0(" It is the lowest candidate, so a lower cut-off ",
          "might do better; the study does not search below it.") else ""))
    if (is.na(chosen)) chosen <- max(cutsv)
    ## fresh seeds at the chosen cut-off
    if (!is.null(fpar)) {
      fp <- prep(fpar)
      fu <- fp[!fp$held & fp$label %in% c("ok", fails) & is.finite(fp$s), ]
      fu$fl <- fu$s > chosen
      by <- function(g) do.call(rbind, lapply(split(fu, fu[[g]]), function(x)
        data.frame(group = paste(g, x[[g]][1]), n = nrow(x), failing = sum(x$fails),
                   specificity = mean(!x$fl[!x$fails]), sensitivity = mean(x$fl[x$fails]))))
      ft <- rbind(by("family"), by("arm"), by("label"))
      lines <- c(lines, "", paste0("**Fresh seeds at the cut-off ", chosen, "** (",
        nrow(fu), " unheld parameters):"), "", md_table(ft))
    }
    ## coverage of the true SD
    u$fl <- u$s > chosen
    u$cb <- u$truth >= u$base_lo & u$truth <= u$base_hi
    ## remedy 2 redraws only the flagged; the others keep their joint draws
    u$c2 <- ifelse(u$fl, u$truth >= u$r2_lo & u$truth <= u$r2_hi, u$cb)
    cv <- function(x, nm) data.frame(group = nm, n = nrow(x), joint_draws = pc(x$cb),
                                     profile_draws = pc(x$c2), stringsAsFactors = FALSE)
    cvt <- rbind(cv(u[!u$fl, ], "not flagged"), cv(u[u$fl, ], "flagged"),
                 cv(u[u$fl & u$label == "unbounded", ], "flagged, unbounded"),
                 cv(u[u$fl & u$label == "wald_too_wide", ], "flagged, wald_too_wide"),
                 cv(u[u$fl & u$label == "ok", ], "flagged, ok"),
                 cv(u[u$fl & u$ridge %in% TRUE, ], "flagged, ridge"))
    lines <- c(lines, "", paste0("**Coverage of the true SD by the 95% interval** ",
      "(nominal 0.95; calibrated means within 2 Monte Carlo SEs of it; the fallback ",
      "given theta holds the SD at its estimate, so has none):"), "", md_table(cvt))
    ## the lower plateau: how often remedy 2's grid edge acts as a floor, and
    ## whether leaving its lowest step out moves the forecasts (paired by fit
    ## and seed, so the SE is of the per-fit difference)
    if ("low_plateau" %in% names(u)) {
      fu <- u[u$fl, ]
      lines <- c(lines, "", paste0(
        "**The lower plateau, a property of the registered grid** (estimate - 4 to + 6). ",
        "Of ", nrow(fu), " flagged parameters, ", sum(fu$low_plateau %in% TRUE),
        " have a profile flat (deviance change < 0.1) over the grid's lowest step. ",
        "Remedy 2's draws within one grid step of the lower edge: mean share ",
        num(mean(fu$edge_mass, na.rm = TRUE), 3), ", over 5% in ",
        sum(fu$edge_mass > 0.05, na.rm = TRUE), " and over 20% in ",
        sum(fu$edge_mass > 0.2, na.rm = TRUE), "."))
      if (!is.null(fc) && any(fc$remedy == "profile_trim")) {
        a <- fc[fc$remedy == "profile" & fc$cut %in% chosen, ]
        b <- fc[fc$remedy == "profile_trim" & fc$cut %in% chosen, ]
        m <- merge(a, b, by = c("cell", "rep", "h"), suffixes = c("", "_trim"))
        if (nrow(m)) {
          pt <- do.call(rbind, lapply(split(m, m$h), function(z) {
            dd <- z$crps_mean - z$crps_mean_trim
            data.frame(h = z$h[1], fits = nrow(z), crps_profile = signif(mean(z$crps_mean), 3),
                       crps_trimmed = signif(mean(z$crps_mean_trim), 3),
                       difference = signif(mean(dd), 3),
                       mc_se = signif(stats::sd(dd) / sqrt(nrow(z)), 3),
                       beyond_2se = abs(mean(dd)) > 2 * stats::sd(dd) / sqrt(nrow(z)))
          }))
          lines <- c(lines, "", paste0("Remedy 2 at ", chosen, " against the same with the ",
            "lowest grid step left out, mean CRPS on the fits it redraws:"), "", md_table(pt),
            "", if (any(pt$beyond_2se)) paste0("The floor MOVES the forecasts beyond Monte ",
              "Carlo error: an input to the remedy's next design (a wider or adaptive lower ",
              "grid).") else paste0("The difference is within Monte Carlo error at every ",
              "horizon: the floor is harmless in practice."))
        }
      }
    }
    ## forecasts, per fit, flagged where any unheld parameter is
    if (!is.null(fc)) {
      ff <- stats::aggregate(fl ~ cell + rep, u, any)
      fc <- merge(fc, ff, by = c("cell", "rep"), all.x = TRUE)
      fc$fl[is.na(fc$fl)] <- FALSE
      ## remedy 2 at the chosen cut-off; its fits not flagged there keep the
      ## joint draws, which is what remedy 2 would do to them
      pr <- fc[fc$remedy == "profile" & fc$cut %in% chosen, ]
      jn <- fc[fc$remedy == "joint", ]
      r2 <- rbind(pr, transform(jn[!paste(jn$cell, jn$rep) %in% paste(pr$cell, pr$rep), ],
                                remedy = "profile"))
      sets <- list(joint = jn, given_theta = fc[fc$remedy == "given_theta", ], profile = r2)
      ## per-fit rates, averaged over fits: the SE is over fits, since a
      ## panel's 20 series share one fit
      fsum <- function(x, nm) do.call(rbind, lapply(split(x, x$h), function(z) {
        w <- z$n_rows
        data.frame(remedy = nm, h = z$h[1], fits = nrow(z),
          cov80 = paste0(num(stats::weighted.mean(z$cov80, w), 3), " (",
                         num(stats::sd(z$cov80) / sqrt(nrow(z)), 3), ")"),
          cov95 = paste0(num(stats::weighted.mean(z$cov95, w), 3), " (",
                         num(stats::sd(z$cov95) / sqrt(nrow(z)), 3), ")"),
          crps_mean = signif(stats::weighted.mean(z$crps_mean, w), 3),
          crps_median = signif(stats::median(z$crps_median), 3),
          true_crps = signif(stats::weighted.mean(z$ref_mean, w), 3),
          explosive = num(sum(z$n_explosive) / sum(w), 4), stringsAsFactors = FALSE)
      }))
      for (grp in c("flagged", "not flagged")) {
        tab <- do.call(rbind, lapply(names(sets), function(nm) {
          x <- sets[[nm]]; fsum(x[x$fl == (grp == "flagged"), ], nm) }))
        lines <- c(lines, "", paste0("**Forecasts, fits ", grp, " at ", chosen,
          "** (B1 and B2; explosive is the share of rows with CRPS over 100 times the ",
          "true-parameter predictive's):"), "", md_table(tab))
      }
    }
    ## time
    if (!is.null(tm)) {
      tm$remedy2 <- tm$t_prof + tm$t_ridge + ifelse(is.na(tm$t_r2), 0, tm$t_r2)
      q <- function(x) paste0(num(stats::median(x, na.rm = TRUE), 2), " / ",
                              num(stats::quantile(x, 0.9, na.rm = TRUE, names = FALSE), 2))
      tm <- merge(tm, cl, by = "cell")
      tt <- do.call(rbind, lapply(split(tm, tm$arm), function(x)
        data.frame(arm = x$arm[1], fits = nrow(x), fit = q(x$t_fit), draws = q(x$t_draws),
                   profiles = q(x$t_prof), ridge = q(x$t_ridge), remedy2_all = q(x$remedy2))))
      lines <- c(lines, "", paste0("**Time**, seconds per fit on one core, median / 90th ",
        "percentile (remedy2_all is the profiles, the ridge profiles and the profile ",
        "draws at all six candidate cut-offs, so an upper bound for one):"), "", md_table(tt))
    }
    ## the named fixtures
    if (!is.null(fxp)) {
      fxp$flag <- ifelse(fxp$held, "held", ifelse(fxp$s > chosen, "flagged", "-"))
      sd_i <- function(lo, hi) ifelse(is.na(lo), "", paste0(signif(exp(lo), 3), " to ",
                                                            signif(exp(hi), 3)))
      ft <- data.frame(fixture = fxp$fixture, param = fxp$param, s = signif(fxp$s, 3),
        d = signif(fxp$d, 3), label = fxp$label, flag = fxp$flag,
        sd_true = signif(exp(fxp$truth), 3), sd_est = signif(exp(fxp$est), 3),
        joint_95 = sd_i(fxp$base_lo, fxp$base_hi), profile_95 = sd_i(fxp$r2_lo, fxp$r2_hi),
        stringsAsFactors = FALSE)
      lines <- c(lines, "", "**The named fixtures** (SDs on their own scale):", "",
                 md_table(ft))
      if (!is.null(fxf)) {
        h1 <- fxf[fxf$h == 1 & (fxf$remedy != "profile" | fxf$cut %in% chosen), ]
        lines <- c(lines, "", "Their forecast CRPS at horizon 1:", "",
          md_table(data.frame(fixture = h1$fixture, remedy = h1$remedy,
                              crps_mean = signif(h1$crps_mean, 3),
                              true_crps = signif(h1$ref_mean, 3))))
      }
    }
    c(lines, "",
      "CAVEATS that must travel with this result: remedy 2 assumes a flat prior on",
      "the log SD and a Gaussian for the other parameters given it, from the draws'",
      "covariance at the estimate. Its cap, 5 x max(sd(g(y*)), 0.5) on the link scale,",
      "was fixed before the run. The covariate panel fixtures' six future times were",
      "drawn after their generator's 24, from the true AR.")
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
