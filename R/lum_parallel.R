## illume: cross-platform parallelism.
##
## WHERE IT HELPS, AND WHERE IT DOES NOT.
## A single lum_fit() does NOT parallelise.  The inner loop is a Laplace
## approximation -- sparse Cholesky plus an inner Newton solve -- and RTMB builds
## its tape in R, so TMB's parallel_accumulator (a C++ template feature) is
## unavailable.  Measured on this machine: TMB::openmp() reports max threads = 1,
## and a J=6 / N=6000 fit took 31.7 s, 29.2 s and 28.2 s under openmp(1), (2) and
## (4) -- noise, not speedup.  So there is deliberately no `ncores` argument on
## lum_fit(): an argument that silently does nothing is worse than none.
##
## What DOES parallelise is every procedure built on independent refits:
##   lum_consistency()   B refits of simulated data
##   parametric bootstrap LRTs
##   LRT-based lum_anova()   one refit per dropped term
##   simulation envelopes for the residual diagnostics
## These share the pool below.

#' Start a worker pool
#'
#' Uses a PSOCK cluster rather than forking. Forking (`parallel::mclapply()`) is
#' faster where available but runs **silently sequentially on Windows**, which
#' would make a `ncores` argument a no-op on one of the three target platforms
#' without saying so. PSOCK behaves identically everywhere.
#'
#' @section What does and does not parallelise:
#' A single model fit does not. Its inner loop is a Laplace approximation, and
#' `RTMB` builds its computation graph in R, so TMB's parallel accumulation is
#' unavailable; measurement confirmed no speedup from additional threads. There
#' is deliberately no `ncores` argument on [lum_fit()], because an argument that
#' silently does nothing is worse than no argument.
#'
#' What does parallelise is everything built from independent refits:
#' [lum_consistency()], [lum_pb_lrt()], likelihood-ratio [lum_anova()], and
#' the simulated envelopes in the diagnostics.
#'
#' @param ncores Integer. Workers; 1 or fewer returns `NULL` (run sequentially).
#' @return A cluster object, or `NULL`.
#' @keywords internal
#' @noRd
lum_pool <- function(ncores) {
  ncores <- max(1L, as.integer(ncores))
  if (ncores <= 1L) return(NULL)
  avail <- tryCatch(parallel::detectCores(), error = function(e) 1L)
  if (is.na(avail)) avail <- 1L
  ncores <- min(ncores, avail)
  if (ncores <= 1L) return(NULL)
  cl <- parallel::makeCluster(ncores)
  ok <- tryCatch({
    parallel::clusterEvalQ(cl, suppressPackageStartupMessages(library(RTMB)))
    ## Prototype: ship the functions from the global environment.  A package
    ## would simply library(illume) on each worker instead.
    fns <- Filter(function(n) {
      v <- get0(n, envir = .GlobalEnv, inherits = FALSE); is.function(v)
    }, ls(.GlobalEnv))
    if (length(fns)) parallel::clusterExport(cl, fns, envir = .GlobalEnv)
    ## One TMB thread per worker.  The fit does not thread anyway (see above),
    ## so oversubscribing would only cause contention.
    parallel::clusterEvalQ(cl, try(TMB::openmp(1L), silent = TRUE))
    TRUE
  }, error = function(e) FALSE)
  if (!ok) { try(parallel::stopCluster(cl), silent = TRUE); return(NULL) }
  cl
}

#' Apply a function over a list, in parallel when a pool is supplied
#'
#' Load-balanced, because refits take uneven amounts of time -- some diverge
#' quickly, others grind -- and fixed chunking would leave workers idle.
#'
#' @param cl A cluster, or `NULL` for sequential.
#' @param X A list or vector.
#' @param FUN Function to apply.
#' @param ... Passed to `FUN`.
#' @return A list of results.
#' @keywords internal
#' @noRd
lum_lapply <- function(cl, X, FUN, ...) {
  if (is.null(cl)) lapply(X, FUN, ...) else parallel::parLapplyLB(cl, X, FUN, ...)
}

#' Refit a model to many simulated response vectors
#'
#' Only plain R objects are sent to workers. The fitted object holds external
#' pointers into TMB's C++ memory, which cannot meaningfully be serialised, and
#' only small summary vectors are returned, so whole models never cross a process
#' boundary.
#'
#' @param fit A fitted `"lum_model"` object.
#' @param ys Integer matrix of simulated responses, one column per replicate.
#' @param ncores Integer. Worker processes.
#' @param restarts Integer. Optimiser restarts.
#' @param verbose Logical. Report the pool size.
#' @return A list of [lum_summarise_fit()] vectors, `NULL` where a refit failed.
#' @keywords internal
#' @noRd
lum_refit_many <- function(fit, ys, ncores = 1L, restarts = 2L, verbose = FALSE) {
  ## Pull out plain R objects only.  fit$obj holds external pointers into TMB
  ## and must never be serialised to a worker.
  X <- fit$X; J <- fit$J; rl <- lum_re_list_of(fit)
  rs <- fit$re_struct; arr <- fit$ar; yl <- fit$ylevels
  fm <- fit$family; wt <- fit$weights
  cl <- lum_pool(ncores)
  on.exit(if (!is.null(cl)) try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
  if (verbose)
    cat(sprintf("refitting %d replicates on %d core(s)\n", ncol(ys),
                if (is.null(cl)) 1L else length(cl)))
  one <- function(b) {
    f <- try(lum_fit(X, ys[, b], J, rl, rs, arr, ylevels = yl,
                      weights = wt, family = fm,
                      verbose = FALSE, restarts = restarts), silent = TRUE)
    if (inherits(f, "try-error")) return(NULL)
    if (f$opt$convergence != 0 || !isTRUE(f$sdr$pdHess)) return(NULL)
    lum_summarise_fit(f)          # small named numeric; the fit itself stays put
  }
  lum_lapply(cl, seq_len(ncol(ys)), one)
}
