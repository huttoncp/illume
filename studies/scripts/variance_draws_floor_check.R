## Arm B after the AR(1) correlation floor: which main-run fits it changes,
## and whether the flagged fits' forecast comparison (joint draws against
## draws given theta) depends on them. Run after the main run, not part of
## its registered design; the caveat in the findings reports the result.
##
## Three steps, each from Arm B's own seeds (variance_draws.R's cells and
## simulate(), and its one() for the forecasts), run from the directory that
## holds this repository's worktree as wt-vd:
##   Rscript variance_draws_floor_check.R ident <out.rds> <ncore>
##       every B1 and B2 fit on the build in R_LIBS; run once on the main
##       run's pinned build and once on the build with the floor
##   Rscript variance_draws_floor_check.R rerun <ids.csv> <prefix> <ncore>
##       the whole per-fit pipeline for the fits the floor changes, on the
##       build with the floor
##   Rscript variance_draws_floor_check.R compare <before.rds> <after.rds> <prefix> <out.csv>
##       the fits changed, and the flagged-fit forecast table before and after
## The main run's pinned build is reproduced by the first step: its AR SD
## estimates agree with the main run's to 5e-5 (the csv's rounding), and one()
## on it gives the main run's forecasts to 5e-6.
## ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
mode <- args[1]
src <- readLines("wt-vd/studies/scripts/variance_draws.R")
OFFSET <- 0L; NREP <- 1L; NCORE <- 1L; sp <- tempdir(); ONLY <- NULL
eval(parse(text = src[grep("^## ---- the cells", src)[1]:(grep("^## ---- run", src)[1] - 1)]),
     envir = globalenv())
bind <- function(x) {
  x <- Filter(Negate(is.null), x); if (!length(x)) return(NULL)
  nm <- unique(unlist(lapply(x, names)))
  do.call(rbind, lapply(x, function(d) { for (n in setdiff(nm, names(d))) d[[n]] <- NA; d[nm] }))
}

if (identical(mode, "ident")) {
  out_file <- args[2]; ncore <- as.integer(args[3])
  jobs <- merge(cells[cells$arm %in% c("B1", "B2"), c("cell", "arm", "family")], data.frame(rep = 1:100))
  one_id <- function(i) {
    suppressPackageStartupMessages(library(illume))
    j <- jobs[i, ]; ce <- cells[cells$cell == j$cell, ]; ce$rep <- j$rep
    set.seed(9973L * ce$cell + j$rep + OFFSET); s <- simulate(ce)
    f <- tryCatch(suppressWarnings(suppressMessages(ilm_model(s$fml, data = s$d, family = ce$family,
           ar = if (s$ar) ilm_ar1(~ t | g) else NULL, verbose = FALSE))), error = function(e) NULL)
    if (is.null(f)) return(data.frame(j, ok = NA, rho_raw = NA, obj = NA))
    pe <- f$opt$par
    data.frame(j, ok = f$ok, rho_raw = unname(pe[names(pe) == "rho_raw"]), obj = f$opt$objective,
               pars = paste(signif(pe, 10), collapse = ","))
  }
  cl <- parallel::makeCluster(ncore)
  parallel::clusterExport(cl, setdiff(ls(globalenv()), "cl"), envir = globalenv())
  r <- bind(parallel::parLapplyLB(cl, seq_len(nrow(jobs)), one_id, chunk.size = 1L))
  parallel::stopCluster(cl)
  saveRDS(r, out_file)
}

if (identical(mode, "rerun")) {
  ids <- utils::read.csv(args[2]); pre <- args[3]; ncore <- as.integer(args[4])
  jl <- lapply(seq_len(nrow(ids)), function(i) { ce <- cells[cells$cell == ids$cell[i], ]; ce$rep <- ids$rep[i]; ce })
  cl <- parallel::makeCluster(ncore)
  parallel::clusterExport(cl, setdiff(ls(globalenv()), "cl"), envir = globalenv())
  res <- parallel::parLapplyLB(cl, jl, function(ce) one(ce), chunk.size = 1L)
  parallel::stopCluster(cl)
  for (k in c("par", "fc")) utils::write.csv(bind(lapply(res, `[[`, k)), paste0(pre, "_", k, ".csv"), row.names = FALSE)
}

if (identical(mode, "compare")) {
  b <- readRDS(args[2]); a <- readRDS(args[3]); pre <- args[4]; out_csv <- args[5]
  m <- merge(b, a, by = c("cell", "arm", "family", "rep"), suffixes = c(".b", ".a"))
  m$changed <- abs(m$obj.a - m$obj.b) > 1e-8 | abs(m$rho_raw.a - m$rho_raw.b) > 1e-6
  m$edge <- abs(tanh(m$rho_raw.b)) > 0.99; m$past <- abs(m$rho_raw.b) > 8
  cat("fits", nrow(m), " |rho| > 0.99", sum(m$edge, na.rm = TRUE), " past the bound", sum(m$past, na.rm = TRUE),
      " changed", sum(m$changed, na.rm = TRUE), " changed with |rho| <= 0.99", sum(m$changed & !m$edge, na.rm = TRUE), "
")
  print(table(arm = m$arm, changed = m$changed))
  ids <- m[m$changed %in% TRUE, c("cell", "rep")]
  rd <- function(f) utils::read.csv(gzfile(f))
  dir <- "wt-vd/studies/runs/0.0.8.9003/variance_draws/"
  par0 <- rd(paste0(dir, "variance_draws_main_par.csv.gz")); fc0 <- rd(paste0(dir, "variance_draws_main_fc.csv.gz"))
  parR <- utils::read.csv(paste0(pre, "_par.csv")); fcR <- utils::read.csv(paste0(pre, "_fc.csv"))
  key <- function(d) paste(d$cell, d$rep)
  swap <- function(a, r) { r <- r[, intersect(names(a), names(r))]; for (n in setdiff(names(a), names(r))) r[[n]] <- NA
    rbind(a[!key(a) %in% key(ids), ], r[names(a)]) }
  par1 <- swap(par0, parR); fc1 <- swap(fc0, fcR)
  fails <- c("unbounded", "wald_too_wide", "profile_failed")
  fx <- function(x) { x[!is.finite(x) | x < 0] <- Inf; x }
  tab <- function(par, fc, tag) {
    u <- par[!par$held & par$label %in% c("ok", fails) & is.finite(par$s), ]; u$fl <- u$s > 0.75
    ff <- stats::aggregate(fl ~ cell + rep, u, any)
    fc <- merge(fc, ff, by = c("cell", "rep"), all.x = TRUE); fc$fl[is.na(fc$fl)] <- FALSE
    out <- NULL
    for (rem in c("joint", "given_theta")) for (h in c(1, 3, 6)) {
      z <- fc[fc$remedy == rem & fc$h == h & fc$fl, ]; w <- z$n_rows; cm <- fx(z$crps_mean); fin <- is.finite(cm)
      out <- rbind(out, data.frame(run = tag, remedy = rem, h = h, fits = nrow(z),
        cov80 = round(weighted.mean(z$cov80, w), 3), cov95 = round(weighted.mean(z$cov95, w), 3),
        crps_mean_finite = signif(weighted.mean(cm[fin], w[fin]), 3), fits_exploded = sum(!fin),
        true_crps = signif(weighted.mean(z$ref_mean, w), 3), explosive = round(sum(z$n_explosive) / sum(w), 4)))
    }
    list(tab = out, fc = fc)
  }
  b <- tab(par0, fc0, "before"); a <- tab(par1, fc1, "after")
  options(width = 160)
  print(rbind(b$tab, a$tab), row.names = FALSE)
  ## the exploded joint forecasts before: how many are fits the floor changes?
  ex <- b$fc[b$fc$remedy == "joint" & b$fc$fl & !is.finite(fx(b$fc$crps_mean)), ]
  exf <- unique(key(ex))
  cat("\nfits whose joint forecast exploded before:", length(exf), "; of them changed by the floor:", sum(exf %in% key(ids)), "\n")
  exa <- a$fc[a$fc$remedy == "joint" & a$fc$fl & !is.finite(fx(a$fc$crps_mean)), ]
  cat("after:", length(unique(key(exa))), "\n")
  cat("changed fits flagged before:", sum(unique(key(b$fc[b$fc$fl, ])) %in% key(ids)),
      " after:", sum(unique(key(a$fc[a$fc$fl, ])) %in% key(ids)), "\n")
  utils::write.csv(rbind(b$tab, a$tab), out_csv, row.names = FALSE)
}
