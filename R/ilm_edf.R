## ---------------------------------------------------------------------------
## Effective degrees of freedom of each smooth.
##
## A smooth's variance is on the scale of its penalised basis, which a
## reparameterisation of the basis moves while the fitted curve stays put, so
## it says little about how wiggly the curve is. Its effective degrees of
## freedom say that, and no reparameterisation moves them.
## ---------------------------------------------------------------------------

## mgcv's definition (Wood 2017, section 6.1.2): the trace of the smooth's
## block of F = (X'WX + S)^-1 X'WX, taken over every coefficient of the fit at
## the estimated variance parameters. X is the whole design, fixed and random
## columns; S the penalty, which here is the random effects' prior precision
## P; and X'WX + S the Hessian H of the joint negative log-likelihood in the
## fixed and random effects at the modes. F = I - H^-1 S, so an unpenalised
## column -- one of the smooth's null space in the fixed effects -- counts
## exactly 1, and the penalised block counts its size less tr(H^-1 P) over it.
## For a multinomial the smooth has a column per category and these are its
## total over them.
##
## H in the random effects is TMB's inner Hessian. Under REML the fixed
## effects are among TMB's random parameters, so that is all of H; by ML their
## rows are taken from the gradient's tape, in the fixed effects only.
#' @keywords internal
#' @noRd
ilm_smooth_edf <- function(object) {
  sm <- object$smooths
  if (!length(sm)) return(NULL)
  C <- object$C
  out <- stats::setNames(rep(NA_real_, length(sm)), names(sm))
  ## the most it can be: every coefficient unpenalised, the null space's
  ## columns (one per category) and every penalised one
  mx <- vapply(names(sm), function(lab) {
    k <- match(ilm_smooth_terms(object, lab), names(object$re))
    ilm_null_dim(sm[[lab]]) * C +
      sum(vapply(k[!is.na(k)], function(j) length(object$b_idx[[j]]), 0))
  }, 0)
  attr(out, "max") <- mx
  H <- tryCatch(ilm_edf_hessian(object), error = function(e) NULL)
  if (is.null(H)) return(out)
  bv <- which(names(object$obj$env$par) == "bvec")
  ch <- tryCatch(Matrix::Cholesky(Matrix::forceSymmetric(H$H), LDL = FALSE,
                                  perm = TRUE), error = function(e) NULL)
  if (is.null(ch)) return(out)
  for (lab in names(sm)) {
    pen <- 0
    for (nm in ilm_smooth_terms(object, lab)) {
      k <- match(nm, names(object$re))
      if (is.na(k)) next
      pos <- match(bv[object$b_idx[[k]]], H$idx)
      P <- ilm_edf_prior(object, k)
      if (is.null(P)) next                  # a variance at zero: no freedom
      E <- Matrix::sparseMatrix(i = pos, j = seq_along(pos), x = 1,
                                dims = c(nrow(H$H), length(pos)))
      Hi <- as.matrix(Matrix::solve(ch, E))[pos, , drop = FALSE]
      pen <- pen + length(pos) - sum(Hi * P)
    }
    out[[lab]] <- ilm_null_dim(sm[[lab]]) * C + pen
  }
  out
}

## A smooth's random terms, named as ilm_model() names them: the label for
## one penalised block, the label and a number for each of several.
#' @keywords internal
#' @noRd
ilm_smooth_terms <- function(object, lab) {
  nr <- length(object$smooths[[lab]]$rand)
  if (nr == 1L) lab else paste0(lab, ".", seq_len(nr))
}

## The number of a smooth's unpenalised columns, its null space in the
## fixed effects.
#' @keywords internal
#' @noRd
ilm_null_dim <- function(sob) if (is.null(sob$Xf)) 0L else ncol(sob$Xf)

## What a smooth's edf is called wherever it is printed -- summary() and
## ilm_varcorr() both take theirs from here.
## `nm` is a random term; NULL unless it is the first (or only) block of a
## smooth.
#' @keywords internal
#' @noRd
ilm_edf_words <- function(object, nm) {
  e <- object$edf
  if (!length(e)) return(NULL)
  for (lab in names(e)) {
    tn <- ilm_smooth_terms(object, lab)
    if (!identical(nm, tn[1L])) next
    val <- if (is.na(e[[lab]])) "edf not available"
           else sprintf("edf %.2f of %d", e[[lab]], as.integer(attr(e, "max")[[lab]]))
    if (object$C > 1L) val <- paste0(val, ", total over the ", object$C,
                                     " categories")
    if (length(tn) > 1L) val <- paste0(val, " for ", lab)
    return(list(smooth = lab, label = val))
  }
  NULL
}

## The note that goes with a smooth's edf wherever it is printed.
#' @keywords internal
#' @noRd
ilm_edf_note <- function()
  paste("the SD is on the basis's scale and not comparable across bases; the",
        "edf is the effective number of parameters the curve uses, at most its",
        "null space plus its basis")

## The Hessian of the joint negative log-likelihood in the fixed effects and
## every random effect, sparse, laid out along `idx` (positions in the fit's
## full parameter vector), at the estimates and modes.
#' @keywords internal
#' @noRd
ilm_edf_hessian <- function(object) {
  env <- object$obj$env
  full <- ilm_full_par(object)
  rn <- names(full); ri <- env$random
  Hr <- env$spHess(full, random = TRUE)
  if (isTRUE(object$reml)) return(list(H = Hr, idx = ri))
  ib <- which(rn == "beta")
  ## the fixed effects' ROWS of the Hessian, one reverse sweep each through
  ## the gradient's tape: their columns cost a sweep per parameter instead,
  ## 100 times slower on a model with 400 random effects
  g <- RTMB::GetTape(object$obj)$jacfun()
  J <- t(RTMB::MakeTape(function(x) g(x)[ib], full)$jacobian(full))
  A <- J[ib, , drop = FALSE]; B <- J[ri, , drop = FALSE]
  H <- rbind(cbind(Matrix::Matrix((A + t(A)) / 2, sparse = TRUE),
                   Matrix::Matrix(t(B), sparse = TRUE)),
             cbind(Matrix::Matrix(B, sparse = TRUE), Hr))
  list(H = H, idx = c(ib, ri))
}

## The prior precision of random term k's block of bvec, in its layout: the
## block is the nl x C matrix of effects by category, column by column, with
## covariance Sigma (x) I -- or, for a reduced-rank term, the nl x r matrix of
## standardised scores, with covariance I. NULL for a variance at zero.
#' @keywords internal
#' @noRd
ilm_edf_prior <- function(object, k) {
  nl <- object$nlk[k]
  if (identical(object$ty[k], "rr")) return(base::diag(nl * object$rk[k]))
  S <- as.matrix(object$Sigma[[k]])
  Si <- tryCatch(solve(S), error = function(e) NULL)
  if (is.null(Si) || any(!is.finite(Si))) return(NULL)
  kronecker(Si, base::diag(nl))
}
