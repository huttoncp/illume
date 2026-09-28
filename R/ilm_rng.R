## ---------------------------------------------------------------------------
## The random-number stream is the user's.
##
## A function given a seed sets it for its own draws -- so its result is the
## same every time -- and used to leave the stream there, so that a user's
## next rnorm() after ilm_contrast(), predict() or any check came out the same
## whatever came before it. That is a side effect nobody asked for, and one
## that makes a simulation of the user's own quietly less random.
##
## So a seeded function puts the stream back as it leaves: as it was, or
## absent if it was absent (a fresh session has no .Random.seed until the
## first draw, and should not be given one). The function's own draws are
## untouched -- set.seed() stays where it was -- so every seeded result is
## what it was. With seed = NULL nothing is set and nothing is restored: the
## function draws from the user's stream, as any R code does.
##
## illumex copies this helper; tests/testthat/test-shared-helpers.R keeps the
## two the same once it has.
## ---------------------------------------------------------------------------

## Register, in the CALLER's frame, a restore of the global random stream to
## run when the caller exits -- normally or by an error. The deferred call is
## base R's on.exit(), reached from the caller's frame the way withr::defer()
## reaches it, with add = TRUE so a handler the caller registers later cannot
## replace it, as long as that one adds too.
#' @keywords internal
#' @noRd
ilm_rng_restore <- function(seed, envir = parent.frame()) {
  if (is.null(seed)) return(invisible(FALSE))
  genv <- globalenv()
  had <- exists(".Random.seed", envir = genv, inherits = FALSE)
  old <- if (had) get(".Random.seed", envir = genv, inherits = FALSE)
  restore <- function() {
    if (had) assign(".Random.seed", old, envir = genv)
    else if (exists(".Random.seed", envir = genv, inherits = FALSE))
      rm(".Random.seed", envir = genv)
  }
  thunk <- as.call(list(function() restore()))
  do.call(base::on.exit, list(thunk, TRUE, FALSE), envir = envir)
  invisible(TRUE)
}
