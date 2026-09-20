## ---------------------------------------------------------------------------
## Monte Carlo coverage of illume's Wald confidence intervals for fixed effects.
##
## Design note.  The cells deliberately span the LATENT BUDGET -- observations
## per latent value, which governs how accurate the Laplace approximation is.
## Two multinomial cells sit below the ~3 obs/latent working rule, where
## coverage is expected to degrade.  That is deliberate: a coverage study that
## only samples the comfortable region cannot fail, and so proves nothing.
##
## The gaussian fixed-effects cell is the opposite check.  There the intervals
## are exact t intervals, so anything other than ~95% would indicate a bug
## rather than an approximation error.
##
## Raw per-replicate estimates are written to disk per cell, so a long run can
## be inspected while it is still going, and any summary can be recomputed
## later without refitting.
##
## The cells cover every model type the package fits, not just the ones it
## started with.  A coverage study that only exercises the families it was
## written for stops saying anything the moment a new one is added, and the
## new ones are exactly where a mistake is most likely.
##
## Cells are grouped by `kind`, which decides how the response is generated and
## what extra arguments the fit needs:
##
##   glm     the original families, response from the linear predictor
##   ar      an AR(1) latent field over evenly spaced time
##   car     an OU latent process over irregular time
##   censor  gaussian with a ceiling, so a Tobit fit
##   aft     time to an event, with right-censored follow-up
##   rp      the same, from a baseline no parametric family can bend to
##   disp    gaussian whose spread differs by a factor
##   smooth  gaussian with a penalised smooth alongside the parametric part
##   slope   a random slope as well as a random intercept
##
## A cell may name the coefficients it checks, in `check`.  For a flexible
## parametric baseline the spline coefficients are nuisance, and for a smooth
## the basis columns are; in both the question is whether the PARAMETRIC
## coefficients keep their coverage while the rest of the model absorbs shape.
##
## Usage: Rscript coverage_study.R <nrep> <ncore> <outdir> [cell1,cell2,...]
## ---------------------------------------------------------------------------

suppressPackageStartupMessages(library(parallel))

args   <- commandArgs(trailingOnly = TRUE)
NREP   <- if (length(args) >= 1) as.integer(args[1]) else 10L
NCORE  <- if (length(args) >= 2) as.integer(args[2]) else 6L
OUTDIR <- if (length(args) >= 3) args[3] else "."
ONLY   <- if (length(args) >= 4) strsplit(args[4], ",")[[1]] else NULL

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

LEVEL <- 0.95

## ---- cell definitions ------------------------------------------------------
## nlat/obs_per_lat are recorded so the results table can be read against the
## latent budget without recomputing it.
cells <- list(
  list(name = "gauss_fixed", family = "gaussian",   J = NA, re = FALSE,
       ncl = NA, per = NA, N = 400, sd_re = NA),
  list(name = "gauss_mm",    family = "gaussian",   J = NA, re = TRUE,
       ncl = 30, per = 20, sd_re = 0.6),
  list(name = "binom_mm",    family = "binomial",   J = NA, re = TRUE,
       ncl = 30, per = 20, sd_re = 0.6),
  list(name = "pois_mm",     family = "poisson",    J = NA, re = TRUE,
       ncl = 30, per = 20, sd_re = 0.6),
  list(name = "mn_J3_rich",  family = "multinomial", J = 3, re = TRUE,
       ncl = 30, per = 20, sd_re = 0.6),
  list(name = "mn_J3_thin",  family = "multinomial", J = 3, re = TRUE,
       ncl = 60, per = 5,  sd_re = 0.6),
  ## A latent-budget dose-response series at J = 5 (C = 4).  Observations per
  ## latent value are 1.5, 3.0, 5.0 and 10.0, straddling the ~3 working rule, so
  ## the study can say where the approximation starts to fail rather than only
  ## that it does somewhere.
  list(name = "mn_J5_thin",  family = "multinomial", J = 5, re = TRUE,
       ncl = 40, per = 6,  sd_re = 0.6),
  list(name = "mn_J5_mid",   family = "multinomial", J = 5, re = TRUE,
       ncl = 30, per = 12, sd_re = 0.6),
  list(name = "mn_J5_rich",  family = "multinomial", J = 5, re = TRUE,
       ncl = 30, per = 20, sd_re = 0.6),
  list(name = "mn_J5_vrich", family = "multinomial", J = 5, re = TRUE,
       ncl = 30, per = 40, sd_re = 0.6),
  list(name = "mn_J3_ar",    family = "multinomial", J = 3, re = TRUE,
       ncl = 25, Tt = 10, rep_per = 8, sd_re = 0.5, rho = 0.5, ar = TRUE),

  ## ---- families and structures added after the first run -------------------
  ## The negative binomial was supported from the start and never covered here.
  list(name = "nbinom_mm",  kind = "glm", family = "nbinom", re = TRUE,
       ncl = 30, per = 20, sd_re = 0.6, k = 2),
  ## AR(1) on a gaussian response: the original AR cell is multinomial, where
  ## the latent budget is the binding constraint. With a continuous response
  ## the Laplace approximation is exact, so this cell asks a different
  ## question -- whether the correlation machinery itself is right.
  list(name = "gauss_ar1",  kind = "ar", family = "gaussian", re = TRUE,
       ncl = 40, Tt = 10, rep_per = 3, sd_re = 0.6, rho = 0.6),
  ## CAR(1): one latent per observation, which is what continuous time always
  ## produces and what the latent-budget check used to fail.
  list(name = "gauss_car1", kind = "car", family = "gaussian", re = TRUE,
       ncl = 60, per = 6, tmax = 30, sd_re = 0.6, rho = 0.85, sd_e = 0.6),
  ## A ceiling at the 80th percentile, so a fifth of the data is censored.
  list(name = "tobit_ceil", kind = "censor", family = "gaussian", re = FALSE,
       N = 600, q = 0.8),
  ## The three accelerated failure time families, each at 35% censoring.
  list(name = "aft_weib",   kind = "aft", family = "weibull",     re = FALSE,
       N = 600, scale = 0.7, cr = 0.35),
  list(name = "aft_lnorm",  kind = "aft", family = "lognormal",   re = FALSE,
       N = 600, scale = 0.7, cr = 0.35),
  list(name = "aft_llogis", kind = "aft", family = "loglogistic", re = FALSE,
       N = 600, scale = 0.7, cr = 0.35),
  ## A frailty Weibull: survival with a random effect, which survreg cannot fit
  ## and so has no outside reference.
  list(name = "aft_frailty", kind = "aft", family = "weibull", re = TRUE,
       ncl = 60, per = 10, scale = 0.6, cr = 0.3, sd_re = 0.5),
  ## Royston-Parmar against a piecewise-constant hazard. The baseline is one no
  ## parametric family here can draw, so this asks whether the coefficients
  ## keep their coverage when the shape has to be absorbed by the spline.
  list(name = "rp_flex",    kind = "rp", family = "rp", re = FALSE,
       N = 800, rp_df = 5, cut = 3, h1 = 0.15, h2 = 0.9),
  ## A dispersion model: spreads of 0.5 and 2.0 by stratum.
  list(name = "disp_group", kind = "disp", family = "gaussian", re = FALSE,
       N = 600),
  ## A penalised smooth alongside the parametric terms. The smooth is a random
  ## effect, so this is also a latent-budget cell of a different shape.
  list(name = "gauss_smooth", kind = "smooth", family = "gaussian", re = FALSE,
       N = 600),
  ## A random slope as well as an intercept, which no earlier cell had.
  list(name = "gauss_slope", kind = "slope", family = "gaussian", re = TRUE,
       ncl = 40, per = 15, sd_re = 0.6, sd_sl = 0.4)
)
names(cells) <- vapply(cells, function(z) z$name, "")
if (!is.null(ONLY)) cells <- cells[ONLY]

## ---- true fixed effects ----------------------------------------------------
## Fixed once per cell, shared by every replicate, so "truth" is a constant and
## coverage is a property of the intervals rather than of a moving target.
truth_of <- function(cell) {
  C <- if (identical(cell$family, "multinomial")) cell$J - 1L else 1L
  p <- 4L                                   # (Intercept), x1, grpb, grpc
  set.seed(99L)
  matrix(round(runif(p * C, -0.7, 0.7), 3), p, C)
}

## ---- the covariate frame every cell shares ---------------------------------
## Same predictors everywhere, so a coefficient means the same thing from one
## cell to the next and the table can be read down a column.
covars <- function(N) {
  data.frame(x1 = rnorm(N),
             grp = factor(sample(c("a", "b", "c"), N, TRUE)))
}

## ---- generators for the model types added later ----------------------------
## Each returns a data frame carrying whatever the fit needs as attributes:
## "ar" a correlation structure, "censor" a censoring specification, "truth" a
## NAMED vector when the coefficients checked are a subset of coef().
gen_extra <- function(cell, seed) {
  Bt <- truth_of(cell)            # before the replicate seed; see gen()
  set.seed(seed)
  b <- as.vector(Bt)

  if (cell$kind == "ar") {
    ncl <- cell$ncl; Tt <- cell$Tt; rp <- cell$rep_per
    N <- ncl * Tt * rp
    dd <- data.frame(g = factor(rep(seq_len(ncl), each = Tt * rp)),
                     time = rep(rep(seq_len(Tt), each = rp), ncl))
    dd <- cbind(dd, covars(N))
    eta <- model.matrix(~ x1 + grp, dd) %*% Bt
    Ba <- numeric(ncl * Tt)
    for (g in seq_len(ncl)) {
      i0 <- (g - 1L) * Tt
      Ba[i0 + 1L] <- rnorm(1, 0, cell$sd_re)
      for (t in 2:Tt)
        Ba[i0 + t] <- cell$rho * Ba[i0 + t - 1L] +
          rnorm(1, 0, cell$sd_re * sqrt(1 - cell$rho^2))
    }
    idx <- (as.integer(dd$g) - 1L) * Tt + dd$time
    dd$y <- as.numeric(eta) + Ba[idx] + rnorm(N, 0, 1)
    attr(dd, "ar") <- illume::ilm_ar1(dd$time, dd$g, verbose = FALSE)
    return(dd)
  }

  if (cell$kind == "car") {
    ncl <- cell$ncl; per <- cell$per; N <- ncl * per
    tl <- lapply(seq_len(ncl), function(i) sort(sample.int(cell$tmax, per)))
    dd <- data.frame(g = factor(rep(seq_len(ncl), each = per)),
                     time = unlist(tl))
    dd <- cbind(dd, covars(N))
    eta <- model.matrix(~ x1 + grp, dd) %*% Bt
    ## an Ornstein-Uhlenbeck path at the observed times
    u <- unlist(lapply(tl, function(tt) {
      z <- numeric(length(tt)); z[1] <- rnorm(1)
      for (k in seq_along(tt)[-1]) {
        ph <- cell$rho^(tt[k] - tt[k - 1])
        z[k] <- ph * z[k - 1] + rnorm(1, 0, sqrt(1 - ph^2))
      }
      z
    }))
    dd$y <- as.numeric(eta) + cell$sd_re * u + rnorm(N, 0, cell$sd_e)
    attr(dd, "ar") <- illume::ilm_car1(dd$time, dd$g, verbose = FALSE)
    return(dd)
  }

  if (cell$kind == "censor") {
    dd <- covars(cell$N)
    eta <- model.matrix(~ x1 + grp, dd) %*% Bt
    ys <- as.numeric(eta) + rnorm(cell$N, 0, 1)
    up <- unname(quantile(ys, cell$q))
    dd$y <- pmin(ys, up)
    attr(dd, "censor") <- illume::ilm_censor(dd$y, upper = up)
    return(dd)
  }

  if (cell$kind == "disp") {
    dd <- covars(cell$N)
    dd$s <- factor(sample(c("lo", "hi"), cell$N, TRUE))
    eta <- model.matrix(~ x1 + grp, dd) %*% Bt
    dd$y <- as.numeric(eta) + rnorm(cell$N, 0, ifelse(dd$s == "hi", 2, 0.5))
    return(dd)
  }

  if (cell$kind == "smooth") {
    dd <- covars(cell$N)
    dd$z <- runif(cell$N, -3, 3)
    eta <- model.matrix(~ x1 + grp, dd) %*% Bt
    dd$y <- as.numeric(eta) + sin(dd$z) * 1.5 + rnorm(cell$N, 0, 1)
    return(dd)
  }

  if (cell$kind == "slope") {
    ncl <- cell$ncl; per <- cell$per; N <- ncl * per
    dd <- data.frame(g = factor(rep(seq_len(ncl), each = per)))
    dd <- cbind(dd, covars(N))
    eta <- model.matrix(~ x1 + grp, dd) %*% Bt
    b0 <- rnorm(ncl, 0, cell$sd_re); b1 <- rnorm(ncl, 0, cell$sd_sl)
    gi <- as.integer(dd$g)
    dd$y <- as.numeric(eta) + b0[gi] + b1[gi] * dd$x1 + rnorm(N, 0, 1)
    return(dd)
  }

  if (cell$kind == "aft") {
    N <- if (isTRUE(cell$re)) cell$ncl * cell$per else cell$N
    dd <- covars(N)
    if (isTRUE(cell$re)) dd$g <- factor(rep(seq_len(cell$ncl), each = cell$per))
    eta <- as.numeric(model.matrix(~ x1 + grp, dd) %*% Bt)
    if (isTRUE(cell$re))
      eta <- eta + rnorm(cell$ncl, 0, cell$sd_re)[as.integer(dd$g)]
    w <- switch(cell$family, weibull = log(rexp(N)), lognormal = rnorm(N),
                loglogistic = rlogis(N))
    tt <- exp(eta + cell$scale * w)
    q <- unname(quantile(tt, 1 - cell$cr))
    ct <- pmin(q, rexp(N, rate = 1 / (3 * q)))
    dd$y <- pmin(tt, ct)
    attr(dd, "censor") <- illume::ilm_surv(dd$y, as.integer(tt <= ct))
    return(dd)
  }

  if (cell$kind == "rp") {
    dd <- covars(cell$N)
    ## the covariate part acts on the log cumulative hazard, so these ARE the
    ## coefficients a proportional-hazards fit estimates
    lp <- as.numeric(model.matrix(~ x1 + grp, dd) %*% Bt) - Bt[1L]
    e <- -log(runif(cell$N))
    H1 <- cell$h1 * cell$cut
    t1 <- e / (cell$h1 * exp(lp))
    tt <- pmax(ifelse(t1 <= cell$cut, t1,
                      cell$cut + (e - H1 * exp(lp)) / (cell$h2 * exp(lp))), 1e-4)
    ct <- rexp(cell$N, rate = 1 / (2 * median(tt)))
    dd$y <- pmin(tt, ct)
    attr(dd, "censor") <- illume::ilm_surv(dd$y, as.integer(tt <= ct))
    ## the intercept and the spline are nuisance here: the baseline shape is
    ## not a coefficient anyone reports
    attr(dd, "truth") <- setNames(b[-1L], c("x1", "grpb", "grpc"))
    return(dd)
  }
  stop("unknown cell kind: ", cell$kind)
}

## ---- one simulated data set ------------------------------------------------
gen <- function(cell, seed) {
  if (!is.null(cell$kind) && cell$kind != "glm") return(gen_extra(cell, seed))
  gen_glm(cell, seed)
}

gen_glm <- function(cell, seed) {
  ## truth_of() sets its own seed so that truth is constant across replicates.
  ## It must therefore be called BEFORE the replicate seed, or it resets the
  ## stream and every replicate simulates identical data.
  C  <- if (identical(cell$family, "multinomial")) cell$J - 1L else 1L
  Bt <- truth_of(cell)
  set.seed(seed)

  if (isTRUE(cell$ar)) {
    ## group x time x replicate: several rows share each (group, time) cell, so
    ## the AR latent values are actually identified
    ncl <- cell$ncl; Tt <- cell$Tt; rp <- cell$rep_per
    N <- ncl * Tt * rp
    dd <- data.frame(
      g    = factor(rep(seq_len(ncl), each = Tt * rp)),
      time = rep(rep(seq_len(Tt), each = rp), ncl))
    dd$x1  <- rnorm(N)
    dd$grp <- factor(sample(c("a", "b", "c"), N, TRUE))
    X <- model.matrix(~ x1 + grp, dd)
    ## AR(1) latent field, one row per (group, time), C columns
    nb <- ncl * Tt
    Bar <- matrix(0, nb, C)
    for (g in seq_len(ncl)) {
      i0 <- (g - 1L) * Tt
      Bar[i0 + 1L, ] <- rnorm(C, 0, cell$sd_re)
      for (t in 2:Tt)
        Bar[i0 + t, ] <- cell$rho * Bar[i0 + t - 1L, ] +
          rnorm(C, 0, cell$sd_re * sqrt(1 - cell$rho^2))
    }
    idx <- (as.integer(dd$g) - 1L) * Tt + dd$time
    eta <- X %*% Bt + Bar[idx, , drop = FALSE]
    attr(dd, "ar") <- list(idx = as.integer(idx), n_group = ncl, Tt = Tt)
  } else if (isTRUE(cell$re)) {
    ncl <- cell$ncl; per <- cell$per; N <- ncl * per
    dd <- data.frame(g = factor(rep(seq_len(ncl), each = per)))
    dd$x1  <- rnorm(N)
    dd$grp <- factor(sample(c("a", "b", "c"), N, TRUE))
    X <- model.matrix(~ x1 + grp, dd)
    b <- matrix(rnorm(ncl * C, 0, cell$sd_re), ncl, C)
    eta <- X %*% Bt + b[as.integer(dd$g), , drop = FALSE]
  } else {
    N <- cell$N
    dd <- data.frame(x1 = rnorm(N),
                     grp = factor(sample(c("a", "b", "c"), N, TRUE)))
    X <- model.matrix(~ x1 + grp, dd)
    eta <- X %*% Bt
  }

  dd$y <- switch(cell$family,
    gaussian = as.numeric(eta) + rnorm(nrow(dd), 0, 1),
    poisson  = rpois(nrow(dd), exp(pmin(as.numeric(eta), 5))),
    nbinom   = rnbinom(nrow(dd), size = cell$k,
                       mu = exp(pmin(as.numeric(eta), 5))),
    binomial = rbinom(nrow(dd), 1, 1 / (1 + exp(-as.numeric(eta)))),
    multinomial = {
      Tc <- contr.sum(cell$J)
      P  <- exp(eta %*% t(Tc)); P <- P / rowSums(P)
      labs <- paste0("c", seq_len(cell$J))
      factor(labs[apply(P, 1, function(pr) sample.int(cell$J, 1L, prob = pr))],
             levels = labs)
    })
  dd
}

## ---- one replicate ---------------------------------------------------------
run_rep <- function(i, cell) {
  dd <- gen(cell, seed = 10000L + i)
  kind <- if (is.null(cell$kind)) "glm" else cell$kind
  ## the formula each kind needs; everything else is shared
  fm <- switch(kind,
    smooth = y ~ x1 + grp + t2(z),
    slope  = y ~ x1 + grp + (1 + x1 | g),
    if (isTRUE(cell$re)) y ~ x1 + grp + (1 | g) else y ~ x1 + grp)
  args <- list(formula = fm, data = dd, family = cell$family,
               ar = attr(dd, "ar"), censor = attr(dd, "censor"),
               verbose = FALSE)
  if (kind == "disp") args$dispformula <- ~ s
  if (kind == "rp")   args$rp_df <- cell$rp_df
  f <- tryCatch(suppressWarnings(do.call(illume::ilm_model, args)),
       error = function(e) structure(list(msg = conditionMessage(e)),
                                     class = "failed"))
  if (inherits(f, "failed"))
    return(list(ok = FALSE, err = f$msg, b = NULL, s = NULL))

  b <- coef(f)
  s <- suppressWarnings(sqrt(diag(vcov(f))))
  ok <- isTRUE(f$opt$convergence == 0L) && isTRUE(f$sdr$pdHess) &&
        length(s) == length(b) && all(is.finite(s)) && all(s > 0)
  crit <- if (isTRUE(f$exact_df)) qt(1 - (1 - LEVEL) / 2, f$resid_df)
          else qnorm(1 - (1 - LEVEL) / 2)

  ## A cell may check a subset: the spline coefficients of a flexible baseline
  ## and the basis columns of a smooth are nuisance, and reporting coverage for
  ## them would say nothing about the parametric part, which is the question.
  want <- names(attr(dd, "truth"))
  if (!is.null(want)) {
    ix <- match(want, names(b))
    if (anyNA(ix))
      return(list(ok = FALSE, err = "checked coefficients not in the fit",
                  b = NULL, s = NULL))
    b <- b[ix]; s <- s[ix]
  }
  list(ok = ok, err = NA_character_, b = unname(b), s = unname(s),
       crit = crit, nm = names(b))
}

## ---- run each cell ---------------------------------------------------------
cl <- makePSOCKcluster(NCORE)
on.exit(stopCluster(cl), add = TRUE)
invisible(clusterEvalQ(cl, suppressPackageStartupMessages(library(illume))))
clusterExport(cl, c("gen", "gen_glm", "gen_extra", "covars", "truth_of",
                    "run_rep", "LEVEL"), envir = environment())

for (cell in cells) {
  t0 <- Sys.time()
  clusterExport(cl, "cell", envir = environment())
  res <- parLapply(cl, seq_len(NREP), function(i) run_rep(i, cell))
  el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  ## the truth a cell checks against, which is the whole coefficient vector
  ## unless the generator named a subset
  tv0 <- attr(gen(cell, seed = 10000L + 1L), "truth")
  tv  <- if (is.null(tv0)) as.vector(truth_of(cell)) else unname(tv0)
  okv <- vapply(res, function(z) isTRUE(z$ok), TRUE)
  nb  <- length(tv)
  B <- t(vapply(res, function(z)
         if (is.null(z$b)) rep(NA_real_, nb) else z$b[seq_len(nb)], numeric(nb)))
  S <- t(vapply(res, function(z)
         if (is.null(z$s)) rep(NA_real_, nb) else z$s[seq_len(nb)], numeric(nb)))
  cr <- vapply(res, function(z) if (is.null(z$crit)) NA_real_ else z$crit, 1)
  nm <- res[[which(!vapply(res, function(z) is.null(z$nm), TRUE))[1]]]$nm
  nm <- if (is.null(nm)) paste0("b", seq_len(nb)) else nm[seq_len(nb)]

  keep <- okv & stats::complete.cases(B) & stats::complete.cases(S)
  Bk <- B[keep, , drop = FALSE]; Sk <- S[keep, , drop = FALSE]; ck <- cr[keep]
  lo <- Bk - ck * Sk; hi <- Bk + ck * Sk
  cov_i <- sweep(lo, 2, tv, "<=") & sweep(hi, 2, tv, ">=")

  out <- data.frame(
    cell      = cell$name,
    family    = cell$family,
    J         = if (is.null(cell$J) || is.na(cell$J)) NA_integer_ else cell$J,
    coef      = nm,
    truth     = tv,
    n_attempt = NREP,
    n_used    = sum(keep),
    coverage  = colMeans(cov_i),
    bias      = colMeans(Bk) - tv,
    emp_sd    = apply(Bk, 2, sd),
    mean_se   = colMeans(Sk),
    stringsAsFactors = FALSE)
  out$se_ratio <- out$mean_se / out$emp_sd
  out$mc_se    <- sqrt(out$coverage * (1 - out$coverage) / out$n_used)
  out$secs     <- el
  ## the latent budget this cell was built to sit at, carried into the results
  ## so the table can be read against it directly
  Cc <- if (identical(cell$family, "multinomial")) cell$J - 1L else 1L
  kind <- if (is.null(cell$kind)) "glm" else cell$kind
  if (isTRUE(cell$ar) || kind == "ar") {
    out$N <- cell$ncl * cell$Tt * cell$rep_per
    out$n_latent <- cell$ncl * cell$Tt * Cc
  } else if (kind == "car") {
    ## continuous time gives one latent per distinct observation time, which
    ## for this design is one per observation
    out$N <- cell$ncl * cell$per
    out$n_latent <- cell$ncl * (cell$per + 1L)
  } else if (kind == "slope") {
    out$N <- cell$ncl * cell$per
    out$n_latent <- cell$ncl * 2L
  } else if (isTRUE(cell$re)) {
    out$N <- cell$ncl * cell$per
    out$n_latent <- cell$ncl * Cc
  } else {
    out$N <- cell$N; out$n_latent <- 0L
  }
  out$kind <- kind
  out$obs_per_latent <- ifelse(out$n_latent > 0, out$N / out$n_latent, NA_real_)
  out$conv_rate <- sum(keep) / NREP

  saveRDS(list(cell = cell, raw = res, summary = out),
          file.path(OUTDIR, paste0("cov_", cell$name, ".rds")))
  write.csv(out, file.path(OUTDIR, paste0("cov_", cell$name, ".csv")),
            row.names = FALSE)
  ## Why replicates were dropped matters as much as how many.  An outright
  ## error and a non-positive-definite Hessian are different failures needing
  ## different explanations, and coverage below is conditional on survival.
  errs <- unlist(lapply(res, function(z) if (!is.na(z$err)) z$err else NULL))
  n_err <- length(errs); n_nonpd <- sum(!okv) - n_err
  cat(sprintf("%-12s  n=%3d/%3d  cover %.3f  se_ratio %.3f  %.0fs\n",
              cell$name, sum(keep), NREP, mean(out$coverage),
              mean(out$se_ratio), el))
  if (n_err || n_nonpd)
    cat(sprintf("               dropped: %d error(s), %d non-pd/non-finite%s\n",
                n_err, n_nonpd,
                if (n_err) paste0("  [", substr(names(sort(table(errs),
                  decreasing = TRUE))[1], 1, 70), "]") else ""))
  flush.console()
}
cat("done\n")
