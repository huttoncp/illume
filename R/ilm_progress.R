## ---------------------------------------------------------------------------
## Progress reporting.
##
## Measured before being added: against a realistic loop -- 2000 bootstrap
## replicates on 20,000 rows -- a bar updated every 1% ran 2.48s against 2.53s
## without one, and updating every single iteration cost 2%. Both are inside the
## noise. The overhead only shows against a body so fast that a bar would be
## pointless anyway: 200,000 iterations of `r + 1` went 0.11s to 0.15s.
##
## So the bar is free wherever it is worth having, and the 1% step is there to
## avoid writing to a slow console rather than to save time.
##
## `utils` is already imported, so nothing new is needed for it.
##
## The default is interactive(): visible when a person is watching, silent in
## scripts, tests and knitr. A bar written into a vignette or a test log is
## noise, and would break every expect_silent() in the suite.
##
## The help page for the `progress` argument, ilm_progress_arg, is illumex's;
## ilm_progress() is a copy of illumex's, compared by
## tests/testthat/test-shared-helpers.R, so change the two together.
## ---------------------------------------------------------------------------

## A bar, or a silent stand-in with the same shape so callers need no branch.
#' @keywords internal
#' @noRd
ilm_progress <- function(n, progress = NULL, label = NULL) {
  on <- isTRUE(progress) ||
    (is.null(progress) && interactive() && n > 1L)
  if (!on || !is.finite(n) || n < 1L)
    return(list(tick = function(i) invisible(NULL),
                done = function() invisible(NULL)))
  if (!is.null(label)) message(label)
  pb <- utils::txtProgressBar(min = 0, max = n, style = 3)
  step <- max(1L, as.integer(n) %/% 100L)
  list(
    tick = function(i) {
      if (i %% step == 0L || i == n) utils::setTxtProgressBar(pb, i)
      invisible(NULL)
    },
    done = function() { utils::setTxtProgressBar(pb, n); close(pb); invisible(NULL) })
}

## lapply with a bar, over one process or many.
##
## parLapplyLB() hands out the work and returns only when all of it is done, so
## a per-element bar is not available across processes. Splitting the indices
## into chunks and dispatching them one chunk at a time gives a bar that moves,
## at the cost of a synchronisation point per chunk -- which is why there are
## only as many chunks as there are ticks worth showing.
#' @keywords internal
#' @noRd
ilm_lapply_progress <- function(cl, X, FUN, ..., progress = NULL) {
  n <- length(X)
  pb <- ilm_progress(n, progress)
  on.exit(pb$done(), add = TRUE)
  if (is.null(cl)) {
    out <- vector("list", n)
    for (i in seq_len(n)) { out[[i]] <- FUN(X[[i]], ...); pb$tick(i) }
    return(out)
  }
  nch <- min(n, max(1L, length(cl) * 4L))
  idx <- split(seq_len(n), cut(seq_len(n), nch, labels = FALSE))
  out <- vector("list", n); done <- 0L
  for (k in idx) {
    out[k] <- parallel::parLapplyLB(cl, X[k], FUN, ...)
    done <- done + length(k); pb$tick(done)
  }
  out
}
