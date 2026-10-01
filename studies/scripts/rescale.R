## Rescaling the fixed-effect columns inside the fit (item 6).
##
## PRE-REGISTERED DESIGN, committed before the change is measured.
##
## ## The problem, measured before anything was changed
## A zero-inflated negative binomial with a random intercept and income in
##   dollars (mean 33,000, SD 9,000) in both parts: 12 of 12 fits failed,
##   the Hessian not positive definite and three standard errors NaN, while
##   the checks advised simplifying the covariance. The same data with income
##   in thousands: 11 of 12 fits graded ok, every SE finite. A zero-inflated
##   Poisson with income only in the count part survived in dollars, but
##   only through the Hessian's recomputation, and its income SE was 0.7%
##   off the thousands fit's. Columns on very different scales make the
##   optimiser's steps and the Hessian's differencing badly conditioned; the
##   model is the same model either way.
##
## ## The change (Craig's ruling, item 6)
## Inside ilm_fit(), each non-constant column j of the fixed design (and of
##   the zero part's and the dispersion model's) is divided by s_j, its SD --
##   and, in the centring variant, has its mean m_j taken off first. The
##   optimisation and the Hessian run on those columns. At the end the
##   estimates and their covariance are converted back EXACTLY, by the linear
##   map beta = A beta_s (+ the intercept's shift when centring):
##     scale only    beta_j = beta_s,j / s_j
##     centring      beta_j = beta_s,j / s_j for j > 1, and the intercept
##                   beta_1 = beta_s,1 - sum_j m_j beta_s,j / s_j
##   V = A V_s A'. The joint precision's fixed block is mapped the same way.
##   Everything downstream sees the user's scale, as a REML fit's already
##   does after its own reshaping: opt$par and sdr$cov.fixed are rewritten
##   there, and a tape on the user's columns is built for code that
##   evaluates the objective at user-scale parameters (as obj_ml is for
##   REML). The remedies name rescaling first when a fit fails on columns
##   whose scales differ by more than a factor of 1e3.
## A model without an intercept is not centred (there is no intercept to
##   absorb the shift); nor is a column that is an indicator (0/1), whose
##   scale is already fine.
##
## ## Scope, stated (and said in the help)
## - Offsets are never rescaled: they enter the linear predictor with a
##   fixed coefficient of 1.
## - A smooth's penalised part is untouched. Its unpenalised (null-space)
##   columns are ordinary fixed columns, rescaled like any other; a test
##   shows the smooth's fitted values and edf unchanged against HAND.
## - Random-slope covariates are out of scope for this pass: their scale
##   sits in the random-effect covariance, not in a fixed column. A badly
##   scaled one keeps the remedy that names rescaling by hand.
## - The separation check's flat catch ("beyond 8 on the link scale") runs
##   on the internal standardised coefficients, where 8 means the same for
##   every column; a test shows the separation cases caught identically
##   before and after the change.
##
## ## The study: scale only against scale and centre
## | Arm | Model | Cells |
## |---|---|---|
## | R1 | the problem case: ZI negative binomial, (1 \| g), income in dollars in both parts | n in {300, 600} |
## | R2 | Poisson, (1 \| g), one covariate in dollars and one of SD 1e-4 | n in {300, 600} |
## | R3 | binomial, (1 \| g), income in dollars | n in {300, 600} |
## | R4 | well-scaled control: each family above with covariates of SD about 1 | 3 families |
## | R5 | a smooth beside a dollar-scale covariate: Poisson, y ~ income + s(x) + (1 \| g) | n = 600 |
## | R6 | out of scope, recorded: Poisson with a random slope on income in dollars, (1 + income \| g) | n = 600 |
##
## 11 cells, 100 replicates each. R5 checks exactness through the smooth
##   machinery (its fitted values and edf as well as the coefficients). R6
##   records the out-of-scope behaviour: the remedy naming rescaling by hand
##   fires, and nothing is silently wrong -- a fit graded ok there has finite
##   SEs that agree with HAND to E1's tolerance. Every data set is fitted four ways: the
##   current build on the raw columns (BEFORE), the current build on columns
##   the user rescaled by hand to SD about 1 (HAND), and the new build on the
##   raw columns in each variant (SCALE, CENTRE).
##
## What is recorded: whether the fit is graded ok; the count of non-finite
##   SEs; each coefficient and SE on the user's scale; the objective; the
##   checks that are not OK; time.
##
## ## What counts, fixed now
## - E1, exactness: wherever HAND and a new variant both converge (ok), their
##   user-scale coefficients and SEs agree to 1e-6 relative, the intercept
##   included under centring, and their objectives to 1e-8. This is the
##   stated test of the conversion back.
## - E2, the failures are fixed: in R1 to R3, each new variant is graded ok
##   in at least as many fits as HAND, less 2 (a variant that loses more
##   than 2 fits of 100 against the hand rescaling does not pass).
## - The noise floor, measured first: in R4, BEFORE refitted from a jittered
##   start (every parameter moved by N(0, 0.1)), against BEFORE from the
##   default start. The floor is the 99th percentile over fits of the
##   largest relative difference in a coefficient or SE.
## - E3, nothing well-scaled changes: in R4, each variant's coefficients and
##   SEs agree with BEFORE within twice that floor, and its grade is the same
##   in every fit. (E1 compares nearly identical problems and keeps 1e-6.)
## - R6, recorded: in every fit where the random slope's covariate is badly
##   scaled and the fit is not ok, the remedies name rescaling; no fit graded
##   ok differs from HAND beyond E1's tolerance.
## - The default is scale only if it passes E1 to E3, since it leaves the
##   intercept exactly where the user's model puts it; centring becomes the
##   default only if scale only fails E2 and centring passes all three. If
##   neither passes, that goes to Craig with the numbers before anything
##   ships.
##
## ## Build
## BEFORE and HAND on the pinned library of the branch the change is made on,
##   before the change (its parent commit); SCALE and CENTRE on a pinned
##   library of the change, the variant chosen by an internal option for the
##   study only. Manifests for both.
##
## ## Build and running
## One pinned library built from the branch's change commit, with its
##   manifest, and one of its parent for BEFORE and HAND; the runs as a
##   sequential chain on one core after the queued full checks, with each
##   fit's check lines recorded.
##
## ## Size
## 1,100 data sets, four fits each and the jittered refits of R4, about a
##   second a fit: about an hour and a half on one core.
##
## ## Amendment A1 (committed before any main-run fit)
## The smoke run (one replicate a cell) showed E1 as written compared two
##   different optimisation problems: HAND on the old build rescales two
##   columns by hand and leaves the rest raw, where the new build
##   standardises every column, so the two agreed only to the optimiser's
##   tolerance (gaps of 1e-5 to 1e-4, and more on parameters near zero) and
##   said nothing about the conversion back. So:
## (i) a fifth way, HANDNEW -- the hand-rescaled data fitted on the NEW
##     build, in each variant. Standardising X and standardising X/c give
##     the same internal columns, so SCALE against HANDNEW (and CENTRE
##     against HANDNEW under centring) is the same problem, and E1 isolates
##     the conversion back. E1 keeps 1e-6 and 1e-8.
## (ii) E1's measure: a coefficient's gap relative to the larger of its size
##     and its SE, so a parameter near zero does not inflate it; an SE's gap
##     relative to itself; the objective's gap absolute.
## (iii) No cell uses REML -- every one is fitted by maximum likelihood --
##     so the objective, invariant to rescaling the covariates under ML, is
##     compared as it is. (Under REML the restricted likelihood moves by
##     log|det A|; the fit puts that Jacobian back itself, and a REML cell
##     would be compared after it.)
## (iv) HAND on the old build stays: E2 is unchanged against it, and its gap
##     to SCALE is reported, not graded. E3 is unchanged.
## (v) A code error found in the same smoke run is fixed: R4's HAND
##     conversion applied the dollar factors to data already in hand units.
##   In the smoke run under A1, every cell met E1 but R5 (the smooth), at
##   1.07e-6 against 1e-6; the threshold stands as approved and the run
##   decides.
##
## ## Amendment A2 (committed before any main-run fit)
## The build under test is 944383a, not 7d77e10 as the design and A1 name
##   it. The full test suite found, before any main-run fit, that 7d77e10's
##   conversion back spread a held term's NA through the whole covariance (a
##   whole-matrix product, where 0 * NA is NA), so a held fit's variance
##   components had no usable covariance and its Satterthwaite df fell back
##   to z. 944383a maps only the coefficient blocks' rows and columns. The
##   change library (lib-rsa) is rebuilt from it; nothing else changes.
##
## ## Amendment A3 (committed before the rerun): the gate
## The full check of 944383a failed 7 tests, all from one cause: dividing
##   every column by its SD moved fits whose columns were already well
##   scaled (factors like 0.97), so the optimiser stopped elsewhere within
##   its tolerance, and a pathological beta landed in a worse basin. Settled
##   within the intent of Craig's ruling on item 6 (2026-09-29): the build
##   under test is now c0377f3, which rescales a
##   column only when its SD is outside [1e-2, 1e2]. So:
## (i) The AFTER arm is rerun on a new pinned library (lib-rsc, from c0377f3,
##     with its manifest), same seeds, same cells. BEFORE and HAND are the
##     old build's and are not rerun: their rows from the first run are used
##     as they are (the same code, library and seeds).
## (ii) HANDNEW changes so that it is still the same internal problem as the
##     variant it is compared with. Under the gate, income / 1e4 (SD 0.9) is
##     inside the band and left alone, while the fit on dollars runs on
##     income / sd(income). So HANDNEW_SCALE divides each rescaled column by
##     its own SD (income, and R2's tiny covariate), and HANDNEW_CENTRE
##     centres it too; both are inside the band, so the new build leaves them
##     exactly as given. Their estimates are put back on the user's scale by
##     the same linear map the fit uses -- beta_x = beta_h / s, and under
##     centring the intercept's beta_0 = beta_h0 - beta_h m / s -- with the SEs
##     from the mapped covariance of each coefficient block (the zero part's
##     likewise).
## (iii) E1 and E2 are regraded on the columns actually rescaled: every cell
##     of R1 to R3, R5 and R6 has a column outside the band (income, SD
##     9,000; R2's tiny, SD 1e-4); R4 has none, so there SCALE and CENTRE
##     are the old code path and E3 must hold with a gap of exactly 0 and
##     the same grade in every fit. R5's smooth: its null-space column (x,
##     SD 0.29) is inside the band, so only income moves. Thresholds as
##     registered: E1 1e-6 and 1e-8, E2 at most 2 fits fewer than HAND.
## (iv) The default rule stands: scale only if it passes E1 to E3; otherwise
##     the numbers go to Craig. The flatness test's fragility found while
##     diagnosing the failures (the beta hold flipping with the path) is
##     reported beside the results, not graded here.

## Choices the design leaves to the code, fixed here before any run:
## - HAND divides income by 1e4 (SD about 0.9) and multiplies the tiny
##   covariate by 1e4 (SD about 1); its coefficients and SEs are put back on
##   the user's scale by those factors before any comparison.
## - Comparisons are over the whole fixed parameter vector -- the fixed
##   effects, the zero part's coefficients, the variance and dispersion
##   parameters -- and their SEs from the fit's covariance.
## - The noise floor re-optimises BEFORE's own objective from a start moved
##   by N(0, 0.1) on every fixed parameter (seed 1), twice, and takes both
##   sets of SEs from RTMB's sdreport at the two optima, so reference and
##   jittered are on the same footing.
## - R4 uses income in units of 1e4 and the tiny covariate times 1e4.
## - Seeds: 1e4 * cell + rep.
##
## Usage: Rscript rescale.R before|after <nrep> <outdir>
##        Rscript rescale.R compare <outdir>
## ---------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
MODE <- args[1]
if (MODE == "compare") sp <- args[2] else { NREP <- as.integer(args[2]); sp <- args[3] }
dir.create(sp, showWarnings = FALSE, recursive = TRUE)
suppressPackageStartupMessages(library(illume))
`%||%` <- function(a, b) if (is.null(a)) b else a

cells <- data.frame(
  cell = 1:11,
  arm  = c("R1", "R1", "R2", "R2", "R3", "R3", "R4", "R4", "R4", "R5", "R6"),
  fam  = c("zinb", "zinb", "poisson", "poisson", "binomial", "binomial",
           "zinb", "poisson", "binomial", "poisson", "poisson"),
  n    = c(300, 600, 300, 600, 300, 600, 500, 500, 500, 600, 600),
  stringsAsFactors = FALSE)

gen <- function(ce, rep) {
  set.seed(1e4 * ce$cell + rep); n <- ce$n; G <- 25
  d <- data.frame(g = factor(sample(G, n, TRUE)), income = round(stats::rnorm(n, 33000, 9000)),
                  tiny = stats::rnorm(n, 0, 1e-4), x = stats::runif(n))
  u <- stats::rnorm(G, 0, 0.4)[d$g]
  lin <- 0.00003 * (d$income - 33000) + 2000 * d$tiny
  if (ce$arm == "R6") lin <- lin + stats::rnorm(G, 0, 0.00001)[d$g] * (d$income - 33000)
  if (ce$arm == "R5") lin <- lin + sin(3 * d$x)
  if (ce$fam == "zinb") {
    pz <- stats::plogis(-1 + 0.00004 * (d$income - 33000))
    d$y <- ifelse(stats::runif(n) < pz, 0, stats::rnbinom(n, mu = exp(0.4 + lin + u), size = 2))
  } else if (ce$fam == "poisson") d$y <- stats::rpois(n, exp(0.4 + lin + u))
  else d$y <- stats::rbinom(n, 1, stats::plogis(-0.3 + lin + u))
  if (ce$arm == "R4") { d$income <- d$income / 1e4; d$tiny <- d$tiny * 1e4 }
  d
}
spec <- function(ce, hand = FALSE) {
  inc <- if (hand) "inc_h" else "income"; tin <- if (hand) "tiny_h" else "tiny"
  rhs <- switch(ce$arm,
    R1 = c(inc, "(1 | g)"), R2 = c(inc, tin, "(1 | g)"), R3 = c(inc, "(1 | g)"),
    R4 = c(inc, tin, "(1 | g)"), R5 = c(inc, "s(x)", "(1 | g)"),
    R6 = c(inc, sprintf("(1 + %s | g)", inc)))
  list(f = stats::reformulate(rhs, "y"),
       zi = if (ce$fam == "zinb") stats::reformulate(inc) else NULL,
       family = if (ce$fam == "zinb") "nbinom" else ce$fam)
}
## the factor that puts a HAND coefficient back on the user's scale, per
## fixed parameter: the income columns by 1e-4, the tiny one by 1e4
hand_factor <- function(f, ce) {
  nm <- names(f$opt$par); k <- rep(1, length(nm))
  if (ce$arm == "R4") return(k)                 # R4's data are in hand units already
  ib <- which(nm == "beta"); xn <- colnames(f$X)
  k[ib[xn == "inc_h"]] <- 1e-4; k[ib[xn == "tiny_h"]] <- 1e4
  if (!is.null(f$Zzi)) { iz <- which(nm == "gzi"); k[iz[colnames(f$Zzi) == "inc_h"]] <- 1e-4 }
  k
}
## Amendment A3: the hand columns h = (x - m) / s, and the map back to the
## user's scale over the fixed parameters: beta_x = beta_h / s and, when
## centred, beta_0 = beta_h0 - beta_h m / s, in the fixed and zero parts.
## `tr` holds m and s for inc_h and tiny_h.
hand_cols <- function(ce, d, how) {
  tr <- list(inc_h = c(m = 0, s = 1), tiny_h = c(m = 0, s = 1))
  if (ce$arm != "R4" && how == "old") tr <- list(inc_h = c(m = 0, s = 1e4), tiny_h = c(m = 0, s = 1e-4))
  if (ce$arm != "R4" && how %in% c("sd", "sdc")) {
    tr$inc_h <- c(m = if (how == "sdc") mean(d$income) else 0, s = stats::sd(d$income))
    tr$tiny_h <- c(m = if (how == "sdc") mean(d$tiny) else 0, s = stats::sd(d$tiny))
  }
  d$inc_h <- (d$income - tr$inc_h[["m"]]) / tr$inc_h[["s"]]
  d$tiny_h <- (d$tiny - tr$tiny_h[["m"]]) / tr$tiny_h[["s"]]
  list(d = d, tr = tr)
}
hand_map <- function(f, tr) {
  nm <- names(f$opt$par); A <- diag(length(nm))
  blk <- function(pname, cn) {
    ib <- which(nm == pname); if (!length(ib) || is.null(cn)) return(invisible())
    i0 <- ib[cn == "(Intercept)"]
    for (h in names(tr)) {
      j <- ib[cn == h]; if (!length(j)) next
      A[j, j] <<- 1 / tr[[h]][["s"]]
      if (length(i0) && tr[[h]][["m"]] != 0) A[i0, j] <<- -tr[[h]][["m"]] / tr[[h]][["s"]]
    }
  }
  blk("beta", colnames(f$X)); blk("gzi", colnames(f$Zzi))
  A
}
one_fit <- function(ce, d, hand = FALSE, how = "old") {
  tr <- NULL
  if (hand) { hc <- hand_cols(ce, d, how); d <- hc$d; tr <- hc$tr }
  s <- spec(ce, hand)
  t0 <- proc.time()[["elapsed"]]
  f <- tryCatch(suppressWarnings(suppressMessages(ilm_model(s$f, data = d, family = s$family,
         ziformula = s$zi, verbose = FALSE))), error = function(e) conditionMessage(e))
  tt <- proc.time()[["elapsed"]] - t0
  if (is.character(f)) return(list(row = data.frame(ok = NA, err = substr(f, 1, 150), time = tt)))
  V <- as.matrix(f$sdr$cov.fixed)
  if (hand && how != "old") {
    ## Amendment A3: the linear map, with each coefficient block's covariance
    ## mapped on its own rows and columns (a held term's NA stays put)
    A <- hand_map(f, tr); par <- drop(A %*% f$opt$par)
    se <- sqrt(pmax(diag(V), 0)); se[!is.finite(diag(V))] <- NA
    I <- which(rowSums(A != diag(nrow(A))) > 0 | colSums(A != diag(nrow(A))) > 0)
    I <- sort(unique(c(I, which(names(f$opt$par) %in% c("beta", "gzi")))))
    Vb <- V[I, I, drop = FALSE]
    if (all(is.finite(Vb))) se[I] <- sqrt(pmax(diag(A[I, I, drop = FALSE] %*% Vb %*% t(A[I, I, drop = FALSE])), 0))
  } else {
    kf <- if (hand) hand_factor(f, ce) else rep(1, length(f$opt$par))
    par <- f$opt$par * kf
    se <- sqrt(pmax(diag(V), 0)) * kf
    se[!is.finite(diag(V))] <- NA
  }
  ck <- f$checks
  row <- data.frame(ok = isTRUE(f$ok), err = NA_character_, time = tt,
                    n_nonfinite_se = sum(!is.finite(se)), objective = f$opt$objective,
                    not_ok = paste(ck$check[ck$status != "OK"], collapse = " "),
                    rescale_remedy = if (!isTRUE(f$ok)) any(grepl("rescale .* by hand",
                      tryCatch(ilm_remedies(f)$remedy, error = function(e) ""))) else NA)
  extra <- list(par = unname(par), se = unname(se))
  if (ce$arm == "R5") {
    extra$edf <- unname(unlist(f$edf))
    extra$pred <- as.numeric(stats::predict(f, newdata = d[1:10, ], type = "link"))
  }
  list(row = row, extra = extra, fit = f)
}
## the noise floor: BEFORE's own objective re-optimised from a jittered start
jitter_floor <- function(f) {
  sdrep <- get("sdreport", asNamespace("RTMB"))
  obj <- f$obj; p0 <- f$opt$par
  set.seed(1); st <- p0 + stats::rnorm(length(p0), 0, 0.1)
  o <- tryCatch({ o1 <- stats::nlminb(st, obj$fn, obj$gr); stats::nlminb(o1$par, obj$fn, obj$gr) },
                error = function(e) NULL)
  if (is.null(o)) return(NA_real_)
  s0 <- tryCatch(sdrep(obj, par.fixed = p0), error = function(e) NULL)
  s1 <- tryCatch(sdrep(obj, par.fixed = o$par), error = function(e) NULL)
  invisible(tryCatch(obj$fn(p0), error = function(e) NULL))
  if (is.null(s0) || is.null(s1)) return(NA_real_)
  se0 <- sqrt(diag(s0$cov.fixed)); se1 <- sqrt(diag(s1$cov.fixed))
  max(abs(o$par - p0) / pmax(abs(p0), 1e-8), abs(se1 - se0) / pmax(se0, 1e-12), na.rm = TRUE)
}

if (MODE %in% c("before", "after")) {
  ways <- if (MODE == "before") c("before", "hand")
          else c("scale", "centre", "handnew_scale", "handnew_centre")
  rows <- list(); extras <- list()
  for (ci in seq_len(nrow(cells))) for (rep in seq_len(NREP)) {
    ce <- cells[ci, ]; d <- gen(ce, rep)
    for (w in ways) {
      if (w != "before" && w != "hand") options(illume.rescale = sub("^handnew_", "", w))
      r <- one_fit(ce, d, hand = w %in% c("hand", "handnew_scale", "handnew_centre"),
                   how = switch(w, handnew_scale = "sd", handnew_centre = "sdc", "old"))
      rr <- cbind(data.frame(cell = ce$cell, arm = ce$arm, fam = ce$fam, n = ce$n, rep = rep, way = w), r$row)
      if (w == "before" && ce$arm == "R4" && !is.null(r$fit)) rr$floor <- jitter_floor(r$fit)
      rows[[length(rows) + 1L]] <- rr
      extras[[paste(ce$cell, rep, w)]] <- r$extra
    }
  }
  nm <- unique(unlist(lapply(rows, names)))
  out <- do.call(rbind, lapply(rows, function(x) { for (k in setdiff(nm, names(x))) x[[k]] <- NA; x[nm] }))
  utils::write.csv(out, file.path(sp, paste0("rescale_", MODE, ".csv")), row.names = FALSE)
  saveRDS(extras, file.path(sp, paste0("rescale_", MODE, "_extras.rds")))
  cat("illume", format(utils::packageVersion("illume")), "from", find.package("illume"), "\n")
  cat(MODE, "fits:", nrow(out), "\n")
}

if (MODE == "compare") {
  rd <- lapply(c("before", "after"), function(m) utils::read.csv(file.path(sp, paste0("rescale_", m, ".csv"))))
  nm <- unique(unlist(lapply(rd, names)))
  a <- do.call(rbind, lapply(rd, function(x) { for (k in setdiff(nm, names(x))) x[[k]] <- NA; x[nm] }))
  ex <- c(readRDS(file.path(sp, "rescale_before_extras.rds")), readRDS(file.path(sp, "rescale_after_extras.rds")))
  ## a coefficient's gap relative to the larger of its size and its SE, so a
  ## parameter near zero does not inflate it; an SE's relative to itself
  gap <- function(x, y) {
    if (is.null(x) || is.null(y) || length(x$par) != length(y$par)) return(NA_real_)
    v <- c(abs(x$par - y$par) / pmax(abs(y$par), y$se, 1e-12), abs(x$se - y$se) / pmax(y$se, 1e-12))
    max(v, na.rm = TRUE)
  }
  key <- unique(a[, c("cell", "arm", "rep")])
  res <- do.call(rbind, lapply(seq_len(nrow(key)), function(i) {
    k <- key[i, ]; g <- function(w) a[a$cell == k$cell & a$rep == k$rep & a$way == w, ]
    e <- function(w) ex[[paste(k$cell, k$rep, w)]]
    h <- g("hand"); b <- g("before"); hs <- g("handnew_scale"); hc <- g("handnew_centre")
    data.frame(k, ok_before = b$ok, ok_hand = h$ok, ok_scale = g("scale")$ok, ok_centre = g("centre")$ok,
      floor = b$floor %||% NA,
      ## E1 (amendment A1): against the same data rescaled by hand on the NEW
      ## build, the same internal problem, so only the conversion differs
      gap_scale_hand = if (isTRUE(hs$ok) && isTRUE(g("scale")$ok)) gap(e("scale"), e("handnew_scale")) else NA,
      gap_centre_hand = if (isTRUE(hc$ok) && isTRUE(g("centre")$ok)) gap(e("centre"), e("handnew_centre")) else NA,
      dobj_scale_hand = if (isTRUE(hs$ok) && isTRUE(g("scale")$ok)) abs(g("scale")$objective - hs$objective) else NA,
      dobj_centre_hand = if (isTRUE(hc$ok) && isTRUE(g("centre")$ok)) abs(g("centre")$objective - hc$objective) else NA,
      ## and against HAND on the old build, reported, not graded
      gap_scale_oldhand = if (isTRUE(h$ok) && isTRUE(g("scale")$ok)) gap(e("scale"), e("hand")) else NA,
      gap_scale_before = gap(e("scale"), e("before")), gap_centre_before = gap(e("centre"), e("before")),
      grade_same_scale = identical(g("scale")$ok, b$ok), grade_same_centre = identical(g("centre")$ok, b$ok),
      remedy_scale = g("scale")$rescale_remedy, ok_r6_scale = g("scale")$ok,
      edf_gap_scale = if (!is.null(e("scale")$edf) && !is.null(e("handnew_scale")$edf))
        max(abs(e("scale")$edf - e("handnew_scale")$edf) / pmax(abs(e("handnew_scale")$edf), 1e-8)) else NA,
      pred_gap_scale = if (!is.null(e("scale")$pred) && !is.null(e("handnew_scale")$pred))
        max(abs(e("scale")$pred - e("handnew_scale")$pred)) else NA)
  }))
  utils::write.csv(res, file.path(sp, "rescale_per_dataset.csv"), row.names = FALSE)
  floor <- stats::quantile(res$floor[res$arm == "R4"], 0.99, na.rm = TRUE)
  by_cell <- do.call(rbind, lapply(split(res, res$cell), function(z) data.frame(
    cell = z$cell[1], arm = z$arm[1], n = nrow(z),
    ok_before = sum(z$ok_before, na.rm = TRUE), ok_hand = sum(z$ok_hand, na.rm = TRUE),
    ok_scale = sum(z$ok_scale, na.rm = TRUE), ok_centre = sum(z$ok_centre, na.rm = TRUE),
    e1_scale_max = suppressWarnings(max(z$gap_scale_hand, na.rm = TRUE)),
    e1_centre_max = suppressWarnings(max(z$gap_centre_hand, na.rm = TRUE)),
    e1_obj_scale_max = suppressWarnings(max(z$dobj_scale_hand, na.rm = TRUE)),
    e1_obj_centre_max = suppressWarnings(max(z$dobj_centre_hand, na.rm = TRUE)),
    oldhand_gap_median = suppressWarnings(stats::median(z$gap_scale_oldhand, na.rm = TRUE)),
    e3_scale_max = if (z$arm[1] == "R4") max(z$gap_scale_before, na.rm = TRUE) else NA,
    e3_centre_max = if (z$arm[1] == "R4") max(z$gap_centre_before, na.rm = TRUE) else NA,
    e3_grades_scale = if (z$arm[1] == "R4") all(z$grade_same_scale) else NA,
    e3_grades_centre = if (z$arm[1] == "R4") all(z$grade_same_centre) else NA,
    r5_edf_max = suppressWarnings(max(z$edf_gap_scale, na.rm = TRUE)),
    r5_pred_max = suppressWarnings(max(z$pred_gap_scale, na.rm = TRUE)),
    r6_not_ok = sum(!z$ok_r6_scale, na.rm = TRUE),
    r6_remedy_named = sum(z$remedy_scale %in% TRUE))))
  utils::write.csv(by_cell, file.path(sp, "rescale_by_cell.csv"), row.names = FALSE)
  e1 <- function(v) all(by_cell[[paste0("e1_", v, "_max")]][is.finite(by_cell[[paste0("e1_", v, "_max")]])] <= 1e-6) &&
    all(by_cell[[paste0("e1_obj_", v, "_max")]][is.finite(by_cell[[paste0("e1_obj_", v, "_max")]])] <= 1e-8)
  e2 <- function(v) { z <- by_cell[by_cell$arm %in% c("R1", "R2", "R3"), ]; all(z[[paste0("ok_", v)]] >= z$ok_hand - 2) }
  e3 <- function(v) { z <- by_cell[by_cell$arm == "R4", ]
    all(z[[paste0("e3_", v, "_max")]] <= 2 * floor) && all(z[[paste0("e3_grades_", v)]]) }
  r6 <- by_cell[by_cell$arm == "R6", ]
  v <- data.frame(verdict = c("E1 scale", "E1 centre", "E2 scale", "E2 centre", "E3 scale", "E3 centre",
                              "R6 remedy named in every failed fit", "noise floor (99th pct)"),
                  value = c(e1("scale"), e1("centre"), e2("scale"), e2("centre"), e3("scale"), e3("centre"),
                            r6$r6_remedy_named == r6$r6_not_ok, signif(floor, 3)))
  pass_s <- e1("scale") && e2("scale") && e3("scale"); pass_c <- e1("centre") && e2("centre") && e3("centre")
  default <- if (pass_s) "scale" else if (pass_c) "centre" else "neither: to Craig"
  utils::write.csv(v, file.path(sp, "rescale_verdicts.csv"), row.names = FALSE)
  print(by_cell, row.names = FALSE); print(v, row.names = FALSE); cat("default by the rule:", default, "\n")
}
