## A gaussian residual SD, or an AR latent's SD, left short of zero, that the
## rules do not see.
##
## PRE-REGISTERED DESIGN, committed before any run.
##
## With one observation per cell of a correlation over time, or one row per
## level of a random intercept, the latent terms can take up all the noise:
## the residual SD's optimum is at zero and the likelihood is flat towards it.
## The optimiser can stop short of it. Reported by another agent, reproduced
## on main (8700876):
##   - a single 24-point series, y ~ 1 + AR(1): sigma = 0.034 x sd(y), SE of
##     log sigma 207, nothing held, every check OK, half of the sigma draws
##     above 10 x sd(y) and the largest 5e279;
##   - a 20 x 24 panel of a noise-free AR(1), y ~ 1 + AR(1): sigma = 0.023 x
##     sd, SE 27, the 90th percentile of the sigma draws 1.8e14 (and two more
##     replicates of the same design at 0.041 and 0.051 x sd; replicates that
##     stopped at 0.16 to 0.19 x sd had SEs of 0.3 to 0.4 and sane draws).
## The same agent's per-fit table shows the mirror case: the AR latent's SD
## at about zero, with a log-SD spread of 11.8 over the draws, while sigma is
## fine (gaussian series with weak correlation, T = 48).
## Today's rules: sigma is held below 1e-3 x sd(y) (ilm_sigma_limit), with
## the likelihood flat below it (ilm_disp_flat, pushed 3 lower, tolerance
## 1e-3); the AR block's SD only below 1e-3, with no flatness test.
##
## The rules to be placed, as for the random-effect SD (re_sd_limit.R), with
## the tolerance kept at the ruled 1e-3. Each SD is at its boundary when it is
## below a LINE and the objective, the random effects re-optimised and the
## other fixed parameters held, moves by at most 1e-3 when its log is pushed
## 3 lower. This study places both lines: sigma's on sigma / sd(y); the AR
## SD's on AR SD / sd(y), with the absolute SD (the random-effect rule's
## scale) reported beside it.
##
## Labels, taken from OUTSIDE the rules. An SD is AT ITS BOUNDARY when both
## hold, with every other fixed parameter re-optimised:
##   floor_ok  pinned at log(1e-4 x sd(y)), the fit is as good: its
##             log-likelihood within 1e-3 of the full fit's
##   push_ok   pinned 3 lower on the log scale, the objective rises by at most
##             1e-3 -- stricter than the rule's own push
## Each fit is classed by the labels: sigma flat, AR SD flat, BOTH, or
## neither; and joint_flat records whether pushing both 3 lower together
## (the others re-optimised) also moves the objective by at most 1e-3.
##
## Candidate lines, on each scale: 1e-3 (today's), 1e-2, 0.05, 0.1, 0.2, and
## the flatness test alone. For each: holds, false holds (held but not at the
## boundary by the labels, to be zero), misses. The line adopted for each SD
## is the widest with no false hold, confirmed on fresh seeds (a seed
## offset).
##
## BOTH FLAGGED, pre-registered policy. Where the noise and the latent are
## nearly interchangeable (weak correlation), either can absorb the other's
## variance, so each can look flat on its own while their sum is determined.
## Holding both would set the total variance to zero, which the data refute;
## holding neither leaves the explosive draws. The rule will therefore HOLD
## ONLY THE ONE WHOSE OWN PUSH MOVES THE OBJECTIVE LEAST, and not the other.
## Expected: in fits flat on both, joint_flat is FALSE, and holding the one
## leaves the other curved (its SE ordinary) with finite draws. Evidence that
## would overturn it:
##   - joint_flat TRUE in most both-flat fits: then both are at zero together
##     (only fixed effects remain) and both should be held;
##   - in the verify phase, the held fit's first fixed-effect SE more than 5%
##     from the floor fit's, or the unheld SD's draws exploding: then the
##     choice by smaller push picks the wrong one, and neither is held,
##     with the fit flagged instead.
## The per-cell cross-tab of the four classes is reported in every run.
##
## Design, all gaussian (63 cells, 100 replicates each, 6,300 fits):
##   single     one series, y ~ 1 + AR(1); T in {12, 24, 48}; rho in {0.1,
##              0.3, 0.7}; latent SD 0.8; noise SD in {0, 0.1, 0.3, 0.5}
##                                                                     36 cells
##   panel_ri   20 series, y ~ x + (1 | g) + AR(1); T in {12, 24}; AR rho
##              0.7, SD 0.8; intercept SD 0.5; x coefficient 0.3; noise SD
##              in {0, 0.1, 0.5}                                        6 cells
##   ri         L groups of m rows, y ~ x + (1 | g); (L, m) in {(10, 1),
##              (10, 2), (30, 2)}; intercept SD 0.8; noise SD in {0, 0.1,
##              0.5}                                                    9 cells
##   ar_panel   20 series, y ~ 1 + AR(1), no intercept term; T in {12, 24};
##              rho in {0.5, 0.8}; marginal SD 1; noise SD in {0, 0.1,
##              0.5}                                                   12 cells
##
## With `verify`, phase 2 on the built rules:
##   - held fits: the held SDs' draws finite and below 10 x sd(y) in 99% of
##     draws;
##   - the first fixed effect's SE within 5% of the floor fit's (from the
##     Hessian of the floor objective);
##   - the both-flat fits: which was held, against the policy;
##   - the reported cases (the 24-point series, and the panel's replicates 4,
##     11 and 59) are held.
##
## Usage: Rscript sigma_limit.R <nrep> <ncore> <outdir> [verify] [offset]
## ---------------------------------------------------------------------------
args  <- commandArgs(trailingOnly = TRUE)
NREP  <- if (length(args) >= 1) as.integer(args[1]) else 100L
NCORE <- if (length(args) >= 2) as.integer(args[2]) else 2L
sp    <- if (length(args) >= 3) args[3] else "."
VERIFY <- length(args) >= 4 && identical(args[4], "verify")
## a seed offset, so the confirmation runs on data the line was not chosen on
OFFSET <- if (length(args) >= 5) as.integer(args[5]) else 0L
dir.create(sp, showWarnings = FALSE, recursive = TRUE)

g <- function(...) expand.grid(..., stringsAsFactors = FALSE)
cells <- rbind(
  data.frame(design = "single",   g(Tn = c(12, 24, 48), rho = c(0.1, 0.3, 0.7), noise = c(0, 0.1, 0.3, 0.5)), L = NA, m = NA),
  data.frame(design = "panel_ri", g(Tn = c(12, 24), rho = 0.7, noise = c(0, 0.1, 0.5)), L = NA, m = NA),
  data.frame(design = "ri",       g(Tn = NA, rho = NA, noise = c(0, 0.1, 0.5), Lm = c("10x1", "10x2", "30x2"))[, c("Tn", "rho", "noise")],
             L = rep(c(10, 10, 30), each = 3), m = rep(c(1, 2, 2), each = 3)),
  data.frame(design = "ar_panel", g(Tn = c(12, 24), rho = c(0.5, 0.8), noise = c(0, 0.1, 0.5)), L = NA, m = NA))
cells$cell <- seq_len(nrow(cells))
stopifnot(nrow(cells) == 63L)
jobs <- merge(cells, data.frame(rep = seq_len(NREP)))

ar_series <- function(Tn, rho, sd_marg)
  as.numeric(stats::arima.sim(list(ar = rho), Tn, sd = sd_marg * sqrt(1 - rho^2)))

## the data and the model for one job
simulate <- function(job) {
  if (job$design == "single") {
    d <- data.frame(g = factor("s1"), t = seq_len(job$Tn))
    d$y <- 5 + ar_series(job$Tn, job$rho, 0.8) + stats::rnorm(job$Tn, 0, job$noise)
    list(d = d, fml = y ~ 1, ar = TRUE)
  } else if (job$design == "panel_ri") {
    G <- 20L
    d <- expand.grid(t = seq_len(job$Tn), g = factor(seq_len(G)))
    lat <- unlist(lapply(seq_len(G), function(i) ar_series(job$Tn, job$rho, 0.8)))
    d$x <- stats::rnorm(nrow(d))
    d$y <- 1 + 0.3 * d$x + stats::rnorm(G, 0, 0.5)[d$g] + lat +
      stats::rnorm(nrow(d), 0, job$noise)
    list(d = d, fml = y ~ x + (1 | g), ar = TRUE)
  } else if (job$design == "ri") {
    d <- data.frame(g = factor(rep(seq_len(job$L), each = job$m)))
    d$x <- stats::rnorm(nrow(d))
    d$y <- 1 + 0.3 * d$x + stats::rnorm(job$L, 0, 0.8)[d$g] +
      stats::rnorm(nrow(d), 0, job$noise)
    list(d = d, fml = y ~ x + (1 | g), ar = FALSE)
  } else {
    G <- 20L
    d <- expand.grid(t = seq_len(job$Tn), g = factor(seq_len(G)))
    d$y <- unlist(lapply(seq_len(G), function(i) ar_series(job$Tn, job$rho, 1))) +
      stats::rnorm(nrow(d), 0, job$noise)
    list(d = d, fml = y ~ 1, ar = TRUE)
  }
}

## The reported cases, rebuilt as their reporter's scripts drew them.
reported <- function() {
  single <- local({
    set.seed(5e6 + 1e4 * 2 + 11)
    G <- 40L; n_t <- 24L + 7L
    rho <- rep(0.7, G); sd_lat <- rep(0.8, G); sd_noise <- rep(0.5, G)
    mu <- stats::rnorm(G, 5, 1)
    d <- do.call(rbind, lapply(seq_len(G), function(k) {
      a <- ar_series(n_t, rho[k], sd_lat[k])
      data.frame(g = sprintf("s%02d", k), t = seq_len(n_t),
                 y = mu[k] + a + stats::rnorm(n_t, 0, sd_noise[k]))
    }))
    d <- d[d$g == "s36" & d$t <= 24, ]; d$g <- factor(d$g); d
  })
  panel <- lapply(c(4L, 11L, 59L), function(rep) {
    set.seed(3e6 + 1e4 * 2 + rep)
    n_t <- 30L
    full <- expand.grid(t = seq_len(n_t), g = factor(sprintf("s%02d", 1:20)))
    full$y <- unlist(lapply(1:20, function(k) ar_series(n_t, 0.8, 1)))
    full[full$t <= 24, ]
  })
  c(list(single_s36 = single), stats::setNames(panel, c("panel_4", "panel_11", "panel_59")))
}

## the objective with some parameters pinned, every other fixed parameter
## free: `id` their positions, `value` their pinned values
pinned <- function(f, id, value) {
  p <- f$opt$par
  fn <- function(q) { pp <- p; pp[-id] <- q; pp[id] <- value; f$obj$fn(pp) }
  o <- tryCatch(stats::nlminb(p[-id], fn), error = function(e) NULL)
  invisible(tryCatch(f$obj$fn(p), error = function(e) NULL))   # tape back at the optimum
  o
}
## the objective moved by pinning, with the others held at the optimum: the
## rule's own push
nudge <- function(f, id, value) {
  p <- f$opt$par; p2 <- p; p2[id] <- value
  out <- tryCatch(f$obj$fn(p2) - f$obj$fn(p), error = function(e) NA_real_)
  invisible(tryCatch(f$obj$fn(p), error = function(e) NULL))
  out
}

measure <- function(f, d, VERIFY) {
  num <- c("sigma", "ysd", "se_log", "push_rule", "ll", "ll_floor", "push_full",
           "ar_sd", "se_log_ar", "push_rule_ar", "ll_floor_ar", "push_full_ar",
           "push_rule_both", "push_full_both", "se_b", "se_b_floor")
  out <- c(stats::setNames(as.list(rep(NA_real_, length(num))), num),
           list(held = NA_character_, boundary = NA_character_, draws_ok = NA))
  if (is.null(f)) return(out)
  p <- f$opt$par; pn <- names(p)
  id <- which(pn == "logdisp"); ia <- which(pn == "lchol_ar")
  if (length(id) != 1L) return(out)
  V <- tryCatch(suppressWarnings(vcov(f, full = TRUE)), error = function(e) NULL)
  se <- function(i) if (!is.null(V) && length(i) == 1L && nrow(V) >= i) sqrt(V[i, i]) else NA_real_
  out$ysd <- stats::sd(d$y)
  f0 <- f$obj$fn(p); out$ll <- -f0
  ## the residual SD
  out$sigma <- exp(p[[id]]); out$se_log <- se(id)
  out$push_rule <- nudge(f, id, p[[id]] - 3)
  fl <- pinned(f, id, log(1e-4 * out$ysd))
  if (!is.null(fl)) out$ll_floor <- -fl$objective
  pu <- pinned(f, id, p[[id]] - 3)
  if (!is.null(pu)) out$push_full <- pu$objective - f0
  ## the AR latent's SD (lchol_ar is its log, for one AR term), on the same
  ## two labels, and the two pushed together
  if (length(ia) == 1L) {
    out$ar_sd <- exp(p[[ia]]); out$se_log_ar <- se(ia)
    out$push_rule_ar <- nudge(f, ia, p[[ia]] - 3)
    fa <- pinned(f, ia, log(1e-4 * out$ysd))
    if (!is.null(fa)) out$ll_floor_ar <- -fa$objective
    pa <- pinned(f, ia, p[[ia]] - 3)
    if (!is.null(pa)) out$push_full_ar <- pa$objective - f0
    out$push_rule_both <- nudge(f, c(id, ia), p[c(id, ia)] - 3)
    pb <- pinned(f, c(id, ia), p[c(id, ia)] - 3)
    if (!is.null(pb)) out$push_full_both <- pb$objective - f0
  }
  out$held <- paste(f$hessian_held, collapse = ",")
  out$boundary <- paste(f$boundary_terms, collapse = ",")
  if (VERIFY) {
    out$se_b <- tryCatch(sqrt(diag(vcov(f)))[[1L]], error = function(e) NA_real_)
    ## the floor fit of whichever the built rule held; sigma's when both are
    hid <- if (grepl("dispersion", out$held) || length(ia) != 1L) id
           else if (grepl("(^|,)ar(,|$)", out$held)) ia else id
    fl2 <- pinned(f, hid, log(1e-4 * out$ysd))
    if (!is.null(fl2)) {
      fn <- function(q) { pp <- p; pp[-hid] <- q; pp[hid] <- log(1e-4 * out$ysd); f$obj$fn(pp) }
      H <- tryCatch(stats::optimHess(fl2$par, fn), error = function(e) NULL)
      invisible(tryCatch(f$obj$fn(p), error = function(e) NULL))
      if (!is.null(H)) out$se_b_floor <- tryCatch(sqrt(solve(H)[1L, 1L]), error = function(e) NA_real_)
    }
    dr <- tryCatch(ilm_draws(f, nsim = 200, seed = 1), error = function(e) NULL)
    if (!is.null(dr)) {
      ld <- dr$draws[dr$map$block %in% c("logdisp", "lchol_ar"), , drop = FALSE]
      out$draws_ok <- all(is.finite(ld)) &&
        all(apply(exp(ld), 1L, stats::quantile, 0.99) < 10 * out$ysd)
    }
  }
  out
}

one <- function(job) {
  suppressMessages(library(illume))
  seed <- 7919L * job$cell + job$rep + OFFSET
  set.seed(seed)
  s <- simulate(job)
  f <- tryCatch(suppressMessages(suppressWarnings(
    ilm_model(s$fml, data = s$d, family = "gaussian",
              ar = if (s$ar) ilm_ar1(~ t | g) else NULL, verbose = FALSE))),
    error = function(e) NULL)
  base <- data.frame(job, seed = seed, ok = !is.null(f), stringsAsFactors = FALSE)
  cbind(base, as.data.frame(measure(f, s$d, VERIFY), stringsAsFactors = FALSE))
}

t0 <- Sys.time()
cl <- parallel::makeCluster(NCORE)
parallel::clusterExport(cl, c("VERIFY", "OFFSET", "simulate", "ar_series", "pinned", "nudge", "measure"))
res <- do.call(rbind, parallel::parLapplyLB(cl, split(jobs, seq_len(nrow(jobs))), one))
parallel::stopCluster(cl)
TOL <- 1e-3
res$at_sigma <- res$ll_floor >= res$ll - TOL & res$push_full <= TOL
res$at_ar <- res$ll_floor_ar >= res$ll - TOL & res$push_full_ar <= TOL
res$joint_flat <- res$push_full_both <= TOL
res$ratio <- res$sigma / res$ysd
res$ratio_ar <- res$ar_sd / res$ysd
res$class <- ifelse(is.na(res$at_sigma), NA,
             ifelse(res$at_sigma & res$at_ar %in% TRUE, "both",
             ifelse(res$at_sigma, "sigma", ifelse(res$at_ar %in% TRUE, "ar", "neither"))))
tag <- if (VERIFY) "verify" else "phase1"
utils::write.csv(res, file.path(sp, paste0("sigma_limit_", tag, ".csv")), row.names = FALSE)
cat("fits:", nrow(res), " failed:", sum(!res$ok), " minutes:",
    round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "\n")

r <- res[res$ok & !is.na(res$at_sigma), ]
q <- function(v) if (length(v)) signif(stats::quantile(v, c(0, .5, .99, 1), na.rm = TRUE), 3) else NA
cat("\n## the residual SD\n")
cat("sigma/sd(y) at the boundary (min, median, 99%, max):", q(r$ratio[r$at_sigma]), "\n")
cat("sigma/sd(y) not at the boundary:", q(r$ratio[!r$at_sigma]), "\n")
for (line in c(1e-3, 1e-2, 0.05, 0.1, 0.2, Inf)) {
  rl <- r$ratio < line & r$push_rule <= TOL
  cat(sprintf("line %-6s + flat: holds %4d | at %4d | FALSE %3d | missed %3d\n",
              format(line), sum(rl), sum(rl & r$at_sigma), sum(rl & !r$at_sigma), sum(!rl & r$at_sigma)))
}
ra <- r[!is.na(r$at_ar), ]
cat("\n## the AR latent's SD, as a share of sd(y) and absolute (the RI rule's scale)\n")
cat("AR SD/sd(y) at the boundary:", q(ra$ratio_ar[ra$at_ar]), "| not:", q(ra$ratio_ar[!ra$at_ar]), "\n")
for (sc in c("ratio_ar", "ar_sd")) for (line in c(1e-3, 1e-2, 0.05, 0.1, 0.2, Inf)) {
  rl <- ra[[sc]] < line & ra$push_rule_ar <= TOL
  cat(sprintf("%-8s line %-6s + flat: holds %4d | at %4d | FALSE %3d | missed %3d\n", sc,
              format(line), sum(rl), sum(rl & ra$at_ar), sum(rl & !ra$at_ar), sum(!rl & ra$at_ar)))
}
cat("\n## sigma and the AR SD together, by the outside labels\n")
print(table(class = ra$class, joint_flat = ra$joint_flat, useNA = "ifany"))
cat("per cell:\n")
print(stats::xtabs(~ paste(design, Tn, rho, noise) + class, data = ra))
## the pre-registered policy for both flagged: hold the one whose push moves
## the objective least, and not the other
bb <- ra[ra$class %in% "both", ]
if (nrow(bb)) {
  pick <- ifelse(bb$push_rule <= bb$push_rule_ar, "sigma", "ar")
  cat("both flat:", nrow(bb), "| joint push flat too:", sum(bb$joint_flat, na.rm = TRUE),
      "| the policy picks sigma:", sum(pick == "sigma"), " ar:", sum(pick == "ar"), "\n")
}

if (VERIFY) {
  hs <- grepl("dispersion", r$held) | grepl("dispersion", r$boundary)
  ha <- grepl("(^|,)ar(,|$)", r$held) | grepl("(^|,)ar(,|$)", r$boundary)
  cat("\nheld sigma:", sum(hs), "(at the boundary:", sum(hs & r$at_sigma), ")",
      "| held AR SD:", sum(ha), "(at:", sum(ha & r$at_ar %in% TRUE), ")",
      "| held both:", sum(hs & ha), "\n")
  cat("of fits both flat by the labels:", sum(r$class %in% "both"), "| held both:",
      sum(r$class %in% "both" & hs & ha), "| held one:", sum(r$class %in% "both" & xor(hs, ha)),
      "| held neither:", sum(r$class %in% "both" & !hs & !ha), "\n")
  h <- hs | ha
  cat("draws within 10 x sd(y) in held fits:", sum(r$draws_ok[h], na.rm = TRUE), "of", sum(h), "\n")
  ratio_se <- r$se_b[h] / r$se_b_floor[h]
  cat("SE of the first fixed effect, held fit over floor fit (min, median, max):",
      signif(min(ratio_se, na.rm = TRUE), 3), signif(stats::median(ratio_se, na.rm = TRUE), 3),
      signif(max(ratio_se, na.rm = TRUE), 3), "\n")
  suppressMessages(library(illume))
  for (nm in names(rc <- reported())) {
    d <- rc[[nm]]
    f <- suppressMessages(suppressWarnings(ilm_model(y ~ 1, data = d, family = "gaussian",
                                                      ar = ilm_ar1(~ t | g), verbose = FALSE)))
    cat(sprintf("reported case %-11s sigma/sd(y) %.3g, held: %s\n", nm,
                exp(f$opt$par[["logdisp"]]) / stats::sd(d$y),
                paste(c(f$hessian_held, f$boundary_terms), collapse = ",")))
  }
}
