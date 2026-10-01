## Areal spatial terms: the BYM form, validated (item 115).
##
## PRE-REGISTERED DESIGN, committed before any code and before any run.
##
## ## The question
## The disease-surveillance case study fits a BYM-form areal term: an
##   intrinsic CAR field and an unstructured region effect, written in illume
##   as s(region, bs = "mrf", xt = list(nb = nb)) + (1 | region), the first
##   an mgcv Markov random field smooth (Wood 2017, 5.8.1), the second the
##   package's own random intercept. Nothing in illume's evidence covers it.
##   This study asks, for data generated from that model:
##   Q1 do the fixed effects' 95% intervals cover (the covariate's slope);
##   Q2 do the region-level relative risks' 95% intervals cover, the quantity
##      a disease map reports;
##   Q3 does illume agree with the reference fit of the same model;
##   Q4 how the fit's checks grade these models -- one count per region is
##      the rule in disease mapping, and obs_per_level[region] grades the
##      unstructured term FAIL there (seen in the design probe), so how many
##      fits are graded usable, and what coverage is among the ones that
##      are not;
##   Q5 whether ilm_variogram(coords = ) finds the spatial structure a model
##      without the term leaves behind, names the term that removes it, and
##      stops flagging once the term is in.
##
## ## The reference
## mgcv::gam() with the same Markov random field basis and the unstructured
##   effect as s(region, bs = "re"), method = "REML" (Laplace-approximate REML
##   for a count family). mgcv is a recommended package, already in the pinned
##   library lib-dfs (1.9-1); nothing is installed. It is the same model and
##   the same basis, so it tests illume's fitting, not the choice of prior.
## A disease-mapping reference proper -- an MCMC BYM fit, CARBayes::S.CARbym
##   on 100 replicates of two S1 cells, comparing the slope and the relative-
##   risk map -- is PENDING Craig's ruling (item 200). It needs a new package
##   and its dependencies in a pinned library, so nothing is installed until
##   he rules and the install is cleared. If it is added, it is specified in
##   an addendum committed before it runs; the arms below do not depend on
##   it.
##
## ## The maps
## Rook-adjacency lattices, 6 x 6 (R = 36) and 10 x 10 (R = 100), and the
##   centroids as coordinates. The structured field phi is drawn from the
##   intrinsic CAR with precision D - W on its non-null space and scaled to
##   its target SD; the unstructured theta ~ N(0, sd). A new field every
##   replicate. The covariate x ~ N(0, 1) per region, independent of phi, so
##   spatial confounding is not what is measured; slope 0.3.
##
## ## Arms and cells
## | Arm | Data | Cells |
## |---|---|---|
## | S1 counts | one Poisson count per region, y ~ Poisson(E exp(-0.2 + 0.3 x + phi + theta)), offset log(E), E ~ Gamma with mean Ebar | R in {36, 100} x Ebar in {5, 50} x total SD in {0.3, 0.6} x structured share rho in {0.2, 0.8} = 16 |
## | S2 gaussian | m observations per region, y = 1 + 0.3 x + phi + theta + N(0, 1), x per observation | R in {36, 100} x m in {3, 10} x rho in {0.2, 0.8}, total SD 0.6 = 8 |
## | S3 variogram | S1's counts at Ebar 50, rho 0.8: the model without the spatial term, y ~ x + offset(log(E)) + (1 | region), then the BYM model; ilm_variogram(coords = centroids, B = 39) on each | R in {36, 100} x total SD in {0.3, 0.6} = 4 |
##
## The total SD is sqrt(var(phi) + var(theta)) on the log (S1) or response
##   (S2) scale; rho = var(phi) / total. 500 replicates per S1 and S2 cell
##   (12,000 fits); the mgcv reference on the first 100 of each (2,400
##   fits); S3 at 100 replicates per cell, the BYM refits in its variogram
##   only at R = 36 (the R = 100 refits cost about 5 s each, times 39).
##
## ## What is recorded, per fit
## - the slope's estimate, SE, and 95% interval (z for S1; Satterthwaite t for
##   S2, the df tables' default) and whether it covers 0.3;
## - each region's relative risk (S1: exp of the linear predictor without the
##   offset; S2: the region mean) with its 95% interval from predict(interval
##   = "confidence", nsim = 500), and the share of regions covered;
## - the SDs of phi and theta as estimated, whether either was held at its
##   boundary, the checks' statuses (obs_per_level, latent_budget, hessian,
##   variance_boundary, optimizer, gradient) and ok;
## - each unheld variance parameter's SE on its log scale, from vcov(fit,
##   full = TRUE), and the fit's variance flag: any of them above 0.75, the
##   cut-off Arm B chose. The region intervals come from joint draws, and
##   this build predates the remedy for flagged fits (item 170); Arm B found
##   joint draws explode on them, so region coverage is reported by flag;
## - on the reference replicates, mgcv's slope, SE and the two SDs;
## - S3: the variogram's verdict, the bins flagged, and the advice line;
## - time for each fit.
##
## ## What counts, fixed now
## - "Calibrated": within 2 Monte Carlo SEs of 0.95. At 500 replicates the
##   slope's band is 0.931 to 0.969. The region-level coverage is a share of
##   R correlated regions per fit, so its MC SE is the SD of the per-fit
##   shares over replicates divided by sqrt(500), and its band is computed
##   from that.
## - Coverage is reported over ALL fits (the primary figure: a user sees a
##   fit whether or not a check fails) and separately over the fits graded
##   ok and not ok, held and unheld, and variance-flagged and not.
## - V1: the slope calibrated in every S1 cell with Ebar = 50 and every S2
##   cell. V2: the slope calibrated in every S1 cell with Ebar = 5 (reported
##   either way; low counts are where a Laplace fit is least sure).
## - V3: the region relative risks calibrated in every S1 and S2 cell.
## - V4: agreement with mgcv in every cell: the slope's median |illume -
##   mgcv| / mgcv's SE at most 0.1, and the median ratio of SEs within 0.9
##   to 1.1. Gaussian cells, where illume's REML and mgcv's are the same
##   criterion, fitted with reml = TRUE, are held to 1e-3 on the slope.
## - V5: S3 flags the model without the spatial term in at least 80% of
##   fits at total SD 0.6; V6: after the BYM term, flags in at most 10% (the
##   variogram's own verdict, at its default min_effect); V7: every flagged
##   fit's advice names the Markov random field term.
## - Q4: the share graded ok per cell, and the coverage in each grade. The
##   not-ok fits "cover as well as the ok ones" when, in every S1 cell where
##   both grades hold at least 50 fits, the slope's coverage and the region
##   coverage differ between the grades by less than 2 Monte Carlo SEs of
##   the difference (each grade's MC SE from its own count, combined in
##   quadrature). In a cell with fewer than 50 fits graded ok -- every fit
##   there, if the check fails them all, as the probe suggests -- the not-ok
##   fits are compared with nominal instead: they cover as well when they
##   are calibrated by V1's and V3's rules.
## - If V1 or V3 fails in a cell with Ebar = 50, or the not-ok fits cover as
##   well as the ok ones in that sense while the checks call them unusable,
##   that goes to Craig with a proposal before the case study relies on the
##   term or the check is changed.
##
## ## Fresh seeds
## S1 and S2 again at 250 replicates per cell (band 0.923 to 0.977) with a
##   seed offset, no reference, on the same build; the verdicts must hold
##   there too, and any that does not is reported as such.
##
## ## Build
## S1 and S2 need nothing beyond the df tables' build, and run on its pinned
##   library lib-dfs (illume acf93c2 with mgcv 1.9-1), counts and gaussian
##   alike -- the gaussian cells' t intervals are the df tables' default,
##   which that build carries. S3 tests the variogram's advice, which this
##   item rewrites to name the areal term (today it says illume fits no
##   spatial covariance), so it runs on a later pinned library built from
##   the branch with that change, with its own manifest.
##
## ## Size
## S1 and S2: 12,000 illume fits, about 1 s at R = 36 and 5 s at R = 100,
##   so about 10 core-hours with the predictions; 2,400 mgcv fits, about 3;
##   fresh seeds about 5; S3 about 3. Some 21 core-hours in all.
##
## Choices the design leaves to the code, fixed here before any run:
## - phi is the intrinsic CAR draw on the non-null space of D - W, scaled by
##   the field's typical marginal SD (the root mean diagonal of the
##   generalised inverse of D - W), so its SD varies from map to map as a
##   field's would, around the target; theta ~ N(0, target).
## - E ~ Gamma(shape 2, mean Ebar), a coefficient of variation of 0.71.
## - S2's x is per observation; its region "risk" is the region mean at x = 0,
##   1 + phi + theta.
## - The region intervals are predict(groups = "fitted", interval =
##   "confidence", nsim = 500, seed = 1) -- each region's own effects, not
##   the typical region's -- at x (S1: each region's own x, exposure 1; S2: x = 0), on the
##   response scale.
## - A fit's variance flag: any unheld variance parameter -- the MRF term's,
##   the region intercept's, and S2's residual SD -- with an SE on its log
##   scale above 0.75.
## - S3's verdict on a fit is the variogram's own: flagged when any bin is
##   WARN or FAIL (what its report counts). Its coordinates are the regions'
##   lattice centroids; B = 39 refits and seed 1 for every variogram.
##
## Usage: Rscript spatial_bym.R <arm S1|S2|S3> <nrep> <ncore> <outdir> [offset] [cells]
##        Rscript spatial_bym.R summarise <outdir>
## ---------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
SUMMARISE <- length(args) >= 1 && identical(args[1], "summarise")
if (!SUMMARISE) {
  ARM    <- args[1]
  NREP   <- as.integer(args[2])
  NCORE  <- as.integer(args[3])
  sp     <- args[4]
  OFFSET <- if (length(args) >= 5) as.integer(args[5]) else 0L
  ONLY   <- if (length(args) >= 6) as.integer(strsplit(args[6], ",")[[1]]) else NULL
  stopifnot(ARM %in% c("S1", "S2", "S3"))
} else sp <- args[2]
dir.create(sp, showWarnings = FALSE, recursive = TRUE)
suppressPackageStartupMessages(library(illume))
`%||%` <- function(a, b) if (is.null(a) || !length(a)) b else a

## ---- the cells ----------------------------------------------------------------
g <- function(...) expand.grid(..., stringsAsFactors = FALSE)
cells <- rbind(
  data.frame(arm = "S1", g(R = c(36, 100), Ebar = c(5, 50), tot = c(0.3, 0.6), rho = c(0.2, 0.8)), m = NA),
  data.frame(arm = "S2", g(R = c(36, 100), m = c(3, 10), rho = c(0.2, 0.8)), Ebar = NA, tot = 0.6),
  data.frame(arm = "S3", g(R = c(36, 100), tot = c(0.3, 0.6)), Ebar = 50, rho = 0.8, m = NA))
cells <- cells[, c("arm", "R", "Ebar", "m", "tot", "rho")]
cells$cell <- seq_len(nrow(cells))
stopifnot(nrow(cells) == 28L)
SLOPE <- 0.3; B0 <- c(S1 = -0.2, S2 = 1); NREF <- 100L

## ---- the maps ------------------------------------------------------------------
make_map <- function(R) {
  k <- as.integer(round(sqrt(R))); stopifnot(k * k == R)
  xy <- expand.grid(i = seq_len(k), j = seq_len(k))
  lev <- sprintf("r%03d", seq_len(R))
  nb <- lapply(seq_len(R), function(r) which(abs(xy$i - xy$i[r]) + abs(xy$j - xy$j[r]) == 1))
  names(nb) <- lev
  W <- matrix(0, R, R); for (r in seq_len(R)) W[r, nb[[r]]] <- 1
  e <- eigen(diag(rowSums(W)) - W, symmetric = TRUE)
  pos <- e$values > 1e-8
  ## the generalised inverse's diagonal: the field's marginal variances
  mv <- rowSums(sweep(e$vectors[, pos]^2, 2, e$values[pos], "/"))
  list(R = R, lev = lev, nb = nb, xy = as.matrix(xy), U = e$vectors[, pos],
       lam = e$values[pos], scale = sqrt(mean(mv)))
}
MAPS <- list(`36` = make_map(36), `100` = make_map(100))

## ---- the data ------------------------------------------------------------------
gen <- function(ce, seed) {
  set.seed(seed)
  mp <- MAPS[[as.character(ce$R)]]; R <- ce$R
  phi <- drop(mp$U %*% (stats::rnorm(length(mp$lam)) / sqrt(mp$lam))) / mp$scale *
    sqrt(ce$rho) * ce$tot
  theta <- stats::rnorm(R, 0, sqrt(1 - ce$rho) * ce$tot)
  reg <- factor(mp$lev, levels = mp$lev)
  if (ce$arm %in% c("S1", "S3")) {
    x <- stats::rnorm(R); E <- stats::rgamma(R, shape = 2, rate = 2 / ce$Ebar)
    eta <- B0[["S1"]] + SLOPE * x + phi + theta
    d <- data.frame(region = reg, x = x, E = E, y = stats::rpois(R, E * exp(eta)))
    nd <- data.frame(region = reg, x = x, E = 1)
    truth <- exp(eta)
  } else {
    d <- data.frame(region = rep(reg, each = ce$m))
    d$x <- stats::rnorm(nrow(d))
    d$y <- B0[["S2"]] + SLOPE * d$x + (phi + theta)[as.integer(d$region)] + stats::rnorm(nrow(d))
    nd <- data.frame(region = reg, x = 0)
    truth <- B0[["S2"]] + phi + theta
  }
  list(d = d, nd = nd, truth = truth, nb = mp$nb, xy = mp$xy)
}
fml <- list(
  S1 = y ~ x + offset(log(E)) + s(region, bs = "mrf", xt = list(nb = nb)) + (1 | region),
  S2 = y ~ x + s(region, bs = "mrf", xt = list(nb = nb)) + (1 | region))
ref_fml <- list(
  S1 = y ~ x + s(region, bs = "mrf", xt = list(nb = nb)) + s(region, bs = "re"),
  S2 = y ~ x + s(region, bs = "mrf", xt = list(nb = nb)) + s(region, bs = "re"))

## ---- one fit -------------------------------------------------------------------
status_of <- function(ck, prefix) {
  s <- ck$status[startsWith(ck$check, prefix)]
  if (!length(s)) NA_character_ else paste(unique(s), collapse = "/")
}
one <- function(job) {
  ce <- cells[cells$cell == job$cell, ]
  seed <- 7919L * ce$cell + job$rep + OFFSET
  sim <- gen(ce, seed)
  nb <- sim$nb
  key <- data.frame(cell = ce$cell, arm = ce$arm, R = ce$R, Ebar = ce$Ebar, m = ce$m,
                    tot = ce$tot, rho = ce$rho, rep = job$rep, seed = seed)
  fam <- if (ce$arm == "S1") "poisson" else "gaussian"
  ## the formulas look up the neighbour list here, in this fit's own
  ## environment (as.formula() leaves a formula's environment as it was)
  fo <- fml[[ce$arm]]; environment(fo) <- environment()
  fo_ref <- ref_fml[[ce$arm]]; environment(fo_ref) <- environment()
  t0 <- proc.time()[["elapsed"]]
  f <- tryCatch(suppressMessages(suppressWarnings(ilm_model(
         fo, data = sim$d, family = fam,
         reml = ce$arm == "S2", verbose = FALSE))), error = function(e) conditionMessage(e))
  key$t_fit <- proc.time()[["elapsed"]] - t0
  if (is.character(f)) { key$ok <- NA; key$err <- substr(f, 1, 200); return(key) }
  key$err <- NA_character_
  ck <- f$checks
  key$ok <- isTRUE(f$ok)
  key$held <- paste(f$hessian_held, collapse = "/")
  for (nm in c("obs_per_level", "latent_budget", "hessian", "variance_boundary",
               "optimizer", "gradient"))
    key[[paste0("ck_", nm)]] <- status_of(ck, nm)
  ## the slope
  ct <- tryCatch(suppressWarnings(ilm_coef_table(f)), error = function(e) NULL)
  if (!is.null(ct) && "x" %in% rownames(ct)) {
    key$est <- ct["x", "Estimate"]; key$se <- ct["x", "Std. Error"]
    key$df <- if ("df" %in% names(ct)) ct["x", "df"] else Inf
    q <- if (is.finite(key$df)) stats::qt(0.975, key$df) else stats::qnorm(0.975)
    key$cover <- abs(key$est - SLOPE) <= q * key$se
  }
  ## the variance parameters: SDs, and the SE of each unheld one's log scale
  pe <- f$opt$par; th <- pe[names(pe) == "theta"]
  key$lsd_mrf <- th[1]; key$lsd_region <- th[2]
  if (ce$arm == "S2") key$lsd_resid <- pe[["logdisp"]]
  V <- tryCatch(suppressWarnings(vcov(f, full = TRUE)), error = function(e) NULL)
  vr <- c(mrf = "s(region):L[1,1]", region = "region:L[1,1]")
  hl <- c(mrf = "s(region)", region = "region")
  sev <- c()
  if (!is.null(V)) {
    for (k in names(vr)) {
      se_k <- if (vr[[k]] %in% rownames(V)) sqrt(V[vr[[k]], vr[[k]]]) else NA_real_
      if (hl[[k]] %in% f$hessian_held) se_k <- NA_real_
      key[[paste0("s_", k)]] <- se_k; sev <- c(sev, se_k)
    }
    if (ce$arm == "S2") {
      rd <- setdiff(rownames(V)[grepl("sigma|disp", rownames(V))], character(0))[1]
      se_r <- if (!is.na(rd)) sqrt(V[rd, rd]) else NA_real_
      if ("dispersion" %in% f$hessian_held) se_r <- NA_real_
      key$s_resid <- se_r; sev <- c(sev, se_r)
    }
  }
  key$vflag <- any(sev > 0.75, na.rm = TRUE)
  ## the region relative risks (S1) or means (S2)
  t1 <- proc.time()[["elapsed"]]
  p <- tryCatch(suppressWarnings(stats::predict(f, newdata = sim$nd, groups = "fitted", interval = "confidence",
                                                nsim = 500, seed = 1)),
                error = function(e) conditionMessage(e))
  key$t_pred <- proc.time()[["elapsed"]] - t1
  if (is.list(p)) {
    lo <- as.numeric(p$lower); hi <- as.numeric(p$upper)
    key$region_cover <- mean(lo <= sim$truth & sim$truth <= hi)
    key$region_width <- stats::median(hi - lo)
    key$region_finite <- mean(is.finite(lo) & is.finite(hi))
  } else key$pred_err <- substr(p, 1, 200)
  ## the reference, on the first NREF replicates
  if (job$rep <= NREF && OFFSET == 0L) {
    t2 <- proc.time()[["elapsed"]]
    gm <- tryCatch(suppressWarnings(mgcv::gam(fo_ref,
            data = sim$d, family = if (ce$arm == "S1") stats::poisson() else stats::gaussian(),
            offset = if (ce$arm == "S1") log(sim$d$E) else NULL, method = "REML")),
            error = function(e) NULL)
    key$t_ref <- proc.time()[["elapsed"]] - t2
    if (!is.null(gm)) {
      key$ref_est <- stats::coef(gm)[["x"]]
      key$ref_se <- sqrt(stats::vcov(gm)["x", "x"])
      vc <- tryCatch({ utils::capture.output(v0 <- suppressWarnings(mgcv::gam.vcomp(gm, rescale = FALSE))); v0 },
                     error = function(e) NULL)
      if (is.matrix(vc)) {
        key$ref_lsd_mrf <- log(vc[1, 1]); key$ref_lsd_region <- log(vc[2, 1])
      } else if (is.list(vc) && !is.null(vc$vc)) {
        key$ref_lsd_mrf <- log(vc$vc[1]); key$ref_lsd_region <- log(vc$vc[2])
      }
    }
  }
  key
}

## ---- S3: the variogram on a model without the term, and with it ---------------
one_s3 <- function(job) {
  ce <- cells[cells$cell == job$cell, ]
  seed <- 7919L * ce$cell + job$rep + OFFSET
  sim <- gen(ce, seed); nb <- sim$nb
  key <- data.frame(cell = ce$cell, arm = ce$arm, R = ce$R, tot = ce$tot, rep = job$rep, seed = seed)
  xy <- sim$xy[as.integer(sim$d$region), , drop = FALSE]
  colnames(xy) <- c("cx", "cy"); xy <- as.data.frame(xy)
  vg <- function(f) tryCatch(suppressMessages(suppressWarnings(ilm_variogram(f, coords = xy,
          B = 39L, seed = 1L, plot = FALSE, verbose = FALSE, progress = FALSE))),
          error = function(e) conditionMessage(e))
  said <- function(v, pre) {
    if (is.character(v)) { key[[paste0(pre, "_err")]] <<- substr(v, 1, 200); return(invisible()) }
    st <- v$table$status
    key[[paste0(pre, "_flag")]] <<- any(st %in% c("WARN", "FAIL"))
    key[[paste0(pre, "_bins")]] <<- sum(st %in% c("WARN", "FAIL"))
    a <- illume:::ilm_variogram_advice(v)
    key[[paste0(pre, "_mrf")]] <<- grepl("bs = \"mrf\"", a, fixed = TRUE)
    key[[paste0(pre, "_advice")]] <<- substr(a, 1, 300)
  }
  fo0 <- y ~ x + offset(log(E)) + (1 | region); environment(fo0) <- environment()
  t0 <- proc.time()[["elapsed"]]
  f0 <- tryCatch(suppressMessages(suppressWarnings(ilm_model(fo0, data = sim$d, family = "poisson",
          verbose = FALSE))), error = function(e) conditionMessage(e))
  if (is.character(f0)) key$m0_err <- substr(f0, 1, 200) else said(vg(f0), "m0")
  key$t_m0 <- proc.time()[["elapsed"]] - t0
  ## the BYM model's variogram, at R = 36 only (each of its 39 refits takes
  ## about 5 s at R = 100)
  if (ce$R == 36) {
    fo1 <- fml$S1; environment(fo1) <- environment()
    t1 <- proc.time()[["elapsed"]]
    f1 <- tryCatch(suppressMessages(suppressWarnings(ilm_model(fo1, data = sim$d, family = "poisson",
            verbose = FALSE))), error = function(e) conditionMessage(e))
    if (is.character(f1)) key$m1_err <- substr(f1, 1, 200) else said(vg(f1), "m1")
    key$t_m1 <- proc.time()[["elapsed"]] - t1
  }
  key
}
summarise_s3 <- function(x) {
  do.call(rbind, lapply(sort(unique(x$cell)), function(cl) {
    z <- x[x$cell == cl, ]; r <- function(v) mean(v, na.rm = TRUE)
    data.frame(z[1, c("cell", "R", "tot")], fits = nrow(z),
               flag_without = r(z$m0_flag), mrf_named_when_flagged = r(z$m0_mrf[z$m0_flag %in% TRUE]),
               flag_with_bym = if (!is.null(z$m1_flag)) r(z$m1_flag) else NA,
               errors = sum(!is.na(z$m0_err %||% NA)) + sum(!is.na(z$m1_err %||% NA)),
               time_without = stats::median(z$t_m0), time_with = if (!is.null(z$t_m1)) stats::median(z$t_m1, na.rm = TRUE) else NA)
  }))
}
verdicts_s3 <- function(s) {
  all <- function(v) { v <- v[!is.na(v)]; if (length(v)) base::all(v) else NA }
  data.frame(verdict = c("V5 flags the model without the spatial term in at least 80% of fits at total SD 0.6",
                         "V6 flags the BYM model in at most 10% of fits",
                         "V7 every flagged fit's advice names the Markov random field term"),
             holds = c(all(s$flag_without[s$tot == 0.6] >= 0.8), all(s$flag_with_bym <= 0.10),
                       all(s$mrf_named_when_flagged == 1)))
}

## ---- summaries -----------------------------------------------------------------
## Coverage within 2 Monte Carlo SEs of 0.95: the slope's MC SE is binomial,
## the region coverage's the SD of the per-fit shares over sqrt(n).
summarise_fits <- function(x) {
  x <- x[!is.na(x$ok), ]
  ## columns a run did not produce (no prediction error, no reference on
  ## fresh seeds) are NA throughout
  for (n in c("pred_err", "ref_est", "ref_se", "ref_lsd_mrf", "ref_lsd_region"))
    if (is.null(x[[n]])) x[[n]] <- NA
  band <- function(v, share = FALSE) {
    v <- v[!is.na(v)]; n <- length(v)
    if (!n) return(c(n = 0, est = NA, mcse = NA, cal = NA))
    est <- mean(v); se <- if (share) stats::sd(v) / sqrt(n) else sqrt(0.95 * 0.05 / n)
    c(n = n, est = est, mcse = se, cal = abs(est - 0.95) <= 2 * se)
  }
  rows <- list()
  for (cl in sort(unique(x$cell))) {
    z <- x[x$cell == cl, ]
    for (grp in c("all", "ok", "not_ok", "held", "unheld", "flagged", "unflagged")) {
      w <- switch(grp, all = rep(TRUE, nrow(z)), ok = z$ok, not_ok = !z$ok,
                  held = nzchar(z$held), unheld = !nzchar(z$held),
                  flagged = z$vflag, unflagged = !z$vflag)
      w <- w %in% TRUE
      sl <- band(z$cover[w]); rg <- band(z$region_cover[w], share = TRUE)
      rows[[length(rows) + 1L]] <- data.frame(z[1, c("cell", "arm", "R", "Ebar", "m", "tot", "rho")],
        group = grp, fits = sum(w),
        slope_n = sl[["n"]], slope_cover = sl[["est"]], slope_mcse = sl[["mcse"]], slope_cal = as.logical(sl[["cal"]]),
        region_n = rg[["n"]], region_cover = rg[["est"]], region_mcse = rg[["mcse"]], region_cal = as.logical(rg[["cal"]]),
        region_width = stats::median(z$region_width[w], na.rm = TRUE),
        region_nonfinite = sum(z$region_finite[w] < 1, na.rm = TRUE),
        pred_errors = sum(!is.na(z$pred_err[w])),
        ref_n = sum(!is.na(z$ref_est[w])),
        ref_dslope = stats::median(abs(z$est[w] - z$ref_est[w]) / z$ref_se[w], na.rm = TRUE),
        ref_seratio = stats::median(z$se[w] / z$ref_se[w], na.rm = TRUE),
        ref_dslope_max = suppressWarnings(max(abs(z$est[w] - z$ref_est[w]), na.rm = TRUE)),
        ref_dlsd_mrf = stats::median(abs(z$lsd_mrf[w] - z$ref_lsd_mrf[w]), na.rm = TRUE),
        ref_dlsd_region = stats::median(abs(z$lsd_region[w] - z$ref_lsd_region[w]), na.rm = TRUE),
        time_fit = stats::median(z$t_fit[w]), time_pred = stats::median(z$t_pred[w], na.rm = TRUE),
        stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, rows)
}
## the pre-registered verdicts; NA where no fit bears on one
verdicts <- function(s, fresh = FALSE) {
  all <- function(v) if (length(v)) base::all(v) else NA
  a <- s[s$group == "all", ]
  s1 <- a[a$arm == "S1", ]; s2 <- a[a$arm == "S2", ]
  q4 <- vapply(unique(s1$cell), function(cl) {
    ok <- s[s$cell == cl & s$group == "ok", ]; no <- s[s$cell == cl & s$group == "not_ok", ]
    if (ok$fits >= 50 && no$fits >= 50) {
      d1 <- abs(ok$slope_cover - no$slope_cover) < 2 * sqrt(ok$slope_mcse^2 + no$slope_mcse^2)
      d2 <- abs(ok$region_cover - no$region_cover) < 2 * sqrt(ok$region_mcse^2 + no$region_mcse^2)
      d1 && d2
    } else isTRUE(no$slope_cal) && isTRUE(no$region_cal)
  }, TRUE)
  v <- data.frame(
    verdict = c("V1 slope calibrated in every S1 cell with Ebar = 50 and every S2 cell",
                "V2 slope calibrated in every S1 cell with Ebar = 5",
                "V3 region relative risks calibrated in every S1 and S2 cell",
                "V4 agreement with mgcv in every cell (slope within 0.1 SE, SE ratio 0.9 to 1.1; gaussian slope to 1e-3)",
                "Q4 the fits graded not ok cover as well as the ok ones, in every S1 cell"),
    holds = c(all(c(s1$slope_cal[s1$Ebar %in% 50], s2$slope_cal)),
              all(s1$slope_cal[s1$Ebar %in% 5]),
              all(a$region_cal),
              if (fresh) NA else all(c(a$ref_dslope <= 0.1, a$ref_seratio >= 0.9 & a$ref_seratio <= 1.1,
                                       s2$ref_dslope_max <= 1e-3)),
              all(q4)))
  v
}

if (SUMMARISE) {
  for (tag in c("main", "fresh")) {
    fs <- list.files(sp, pattern = paste0("^spatial_bym_S[12]_", tag, "[.]csv$"), full.names = TRUE)
    if (!length(fs)) next
    x <- lapply(fs, utils::read.csv, stringsAsFactors = FALSE)
    nm <- unique(unlist(lapply(x, names)))
    x <- do.call(rbind, lapply(x, function(d) { for (n in setdiff(nm, names(d))) d[[n]] <- NA; d[nm] }))
    s <- summarise_fits(x)
    utils::write.csv(s, file.path(sp, paste0("spatial_bym_", tag, "_summary.csv")), row.names = FALSE)
    v <- verdicts(s, fresh = tag == "fresh")
    utils::write.csv(v, file.path(sp, paste0("spatial_bym_", tag, "_verdicts.csv")), row.names = FALSE)
    cat("==", tag, "==\n"); print(v, row.names = FALSE)
  }
  fs3 <- list.files(sp, pattern = "^spatial_bym_S3_main[.]csv$", full.names = TRUE)
  if (length(fs3)) {
    s3 <- summarise_s3(utils::read.csv(fs3, stringsAsFactors = FALSE))
    utils::write.csv(s3, file.path(sp, "spatial_bym_S3_summary.csv"), row.names = FALSE)
    v3 <- verdicts_s3(s3)
    utils::write.csv(v3, file.path(sp, "spatial_bym_S3_verdicts.csv"), row.names = FALSE)
    cat("== S3 ==
"); print(s3, row.names = FALSE); print(v3, row.names = FALSE)
  }
  quit(save = "no")
}

## ---- run -----------------------------------------------------------------------
t_all <- Sys.time()
use <- cells[cells$arm == ARM, ]
if (!is.null(ONLY)) use <- use[use$cell %in% ONLY, ]
jobs <- merge(use[, "cell", drop = FALSE], data.frame(rep = seq_len(NREP)))
## the largest maps first, so the long fits do not all end the run
jobs <- jobs[order(-cells$R[match(jobs$cell, cells$cell)], jobs$cell, jobs$rep), ]
jl <- split(jobs, seq_len(nrow(jobs)))
tag <- if (OFFSET) "fresh" else if (!is.null(ONLY)) "smoke" else "main"
## CHECKPOINTS (an I/O change, no effect on results): each fit's result is
## saved as it finishes, and a run started again reads the fits already done.
## Every fit sets its own seed from its cell and replicate.
ck <- file.path(sp, paste0("checkpoints_", ARM, "_", tag))
dir.create(ck, showWarnings = FALSE, recursive = TRUE)
one_ck <- function(job) {
  fp <- file.path(ck, sprintf("cell%02d_rep%04d.rds", job$cell, job$rep))
  if (file.exists(fp)) return(readRDS(fp))
  r <- if (ARM == "S3") one_s3(job) else one(job)
  saveRDS(r, paste0(fp, ".part")); file.rename(paste0(fp, ".part"), fp)
  r
}
res <- if (NCORE > 1L) {
  cl <- parallel::makeCluster(NCORE)
  invisible(parallel::clusterEvalQ(cl, suppressPackageStartupMessages(library(illume))))
  parallel::clusterExport(cl, setdiff(ls(globalenv()), c("cl", "jobs", "jl")), envir = globalenv())
  r <- parallel::parLapplyLB(cl, jl, one_ck, chunk.size = 1L)
  parallel::stopCluster(cl)
  r
} else lapply(jl, one_ck)
nm <- unique(unlist(lapply(res, names)))
out <- do.call(rbind, lapply(res, function(d) { for (n in setdiff(nm, names(d))) d[[n]] <- NA; d[nm] }))
utils::write.csv(out, file.path(sp, paste0("spatial_bym_", ARM, "_", tag, ".csv")), row.names = FALSE)
cat("illume", format(utils::packageVersion("illume")), "from", find.package("illume"), "\n")
cat("fits:", nrow(jobs), " minutes:", round(as.numeric(difftime(Sys.time(), t_all, units = "mins")), 1), "\n")
