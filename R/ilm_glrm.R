## ---------------------------------------------------------------------------
## Generalized low rank models.
##
## PCA, and everything built on it, minimises squared error. On a numeric
## column that is the right thing. On a binary one it is not: squared loss on
## a 0/1 indicator says the distance from 0 to 1 is the same kind of quantity
## as the distance from 1.4 to 2.4, and a reconstruction can land at -0.3 or
## 1.7, which is not a category.
##
## FAMD -- what ilm_reduce() does by default -- handles mixed data by one-hot
## encoding the categories, scaling the indicators and running PCA on the
## result. It is fast, it has a closed form, and its treatment of a category is
## a Gaussian approximation to something that is not Gaussian.
##
## A generalized low rank model keeps the low-rank structure and changes the
## loss, one per column:
##
##     minimise  sum_j  L_j( (X Y)_j + mu_j ,  A_j )  +  lambda ( |X|^2 + |Y|^2 )
##
## with L quadratic for a numeric column, logistic for a binary one,
## multinomial across the levels of a categorical one, and Poisson for a count.
## Each column is then reconstructed on its own scale: a category comes back as
## a category.
##
## It costs an iterative fit where FAMD has a decomposition, and it is an
## OPTION here rather than the default for that reason -- on all-numeric data
## with quadratic loss the two are the same model, and the test suite checks
## that they agree with the SVD to 1e-8.
##
## Reference:
##   Udell, M., Horn, C., Zadeh, R. and Boyd, S. (2016). Generalized low rank
##     models. Foundations and Trends in Machine Learning 9, 1-118.
## ---------------------------------------------------------------------------

#' Which loss belongs to a column
#'
#' @keywords internal
#' @noRd
ilm_glrm_loss_of <- function(x) {
  if (is.logical(x)) return("logistic")
  if (is.factor(x) || is.character(x)) {
    nl <- nlevels(as.factor(x))
    return(if (nl <= 2L) "logistic" else "multinomial")
  }
  if (is.numeric(x)) return("quadratic")
  NA_character_
}

#' Encode the data into blocks of columns, one block per variable
#'
#' A numeric column is one column, standardised. A binary one is one column of
#' 0/1. A categorical one is a block with a column per level, carrying the
#' level index rather than indicators, because a multinomial loss needs to know
#' which level was seen rather than to score each indicator separately.
#'
#' @keywords internal
#' @noRd
ilm_glrm_encode <- function(sub, loss) {
  n <- nrow(sub)
  blocks <- list(); pos <- 0L
  for (v in names(sub)) {
    x <- sub[[v]]; lo <- loss[[v]]
    if (lo %in% c("quadratic", "poisson")) {
      a <- as.numeric(x)
      ctr <- if (lo == "quadratic") mean(a, na.rm = TRUE) else 0
      sc <- if (lo == "quadratic") {
        s <- stats::sd(a, na.rm = TRUE); if (!is.finite(s) || s <= 0) 1 else s
      } else 1
      blocks[[v]] <- list(loss = lo, cols = pos + 1L, target = (a - ctr) / sc,
                          centre = ctr, scale = sc, levels = NULL)
      pos <- pos + 1L
    } else if (lo == "logistic") {
      f <- as.factor(x)
      blocks[[v]] <- list(loss = lo, cols = pos + 1L,
                          target = as.numeric(as.integer(f) - 1L),
                          centre = 0, scale = 1, levels = levels(f))
      pos <- pos + 1L
    } else {
      f <- as.factor(x); L <- nlevels(f)
      blocks[[v]] <- list(loss = lo, cols = pos + seq_len(L),
                          target = as.integer(f), centre = 0, scale = 1,
                          levels = levels(f))
      pos <- pos + L
    }
  }
  list(blocks = blocks, d = pos, n = n)
}

#' Objective and gradient with respect to the linear predictor
#'
#' Returns the loss summed over observed cells and the matrix dL/dU, which is
#' all the alternating fit needs -- the chain rule to X and Y is the same for
#' every loss.
#'
#' @keywords internal
#' @noRd
ilm_glrm_grad <- function(U, enc, w = NULL) {
  G <- matrix(0, nrow(U), ncol(U))
  obj <- 0
  if (is.null(w)) w <- rep(1, nrow(U))
  for (b in enc$blocks) {
    j <- b$cols
    if (b$loss == "quadratic") {
      a <- b$target; ok <- !is.na(a)
      r <- U[, j] - a
      r[!ok] <- 0
      obj <- obj + 0.5 * sum(w * r^2)
      G[, j] <- w * r
    } else if (b$loss == "poisson") {
      a <- b$target; ok <- !is.na(a)
      u <- pmin(U[, j], 30)                    # exp(30) is already 1e13
      e <- exp(u)
      obj <- obj + sum((w * (e - a * u))[ok])
      g <- w * (e - a); g[!ok] <- 0
      G[, j] <- g
    } else if (b$loss == "logistic") {
      a <- b$target; ok <- !is.na(a)
      u <- U[, j]
      obj <- obj + sum((w * (logspace_add0(u) - a * u))[ok])
      g <- w * (stats::plogis(u) - a); g[!ok] <- 0
      G[, j] <- g
    } else {
      yi <- b$target; ok <- !is.na(yi)
      Ub <- U[, j, drop = FALSE]
      m <- apply(Ub, 1L, max)
      P <- exp(Ub - m)
      sm <- rowSums(P)
      P <- P / sm
      ## a row whose category was not observed still needs a valid subscript;
      ## it is masked out of both the objective and the gradient below, but
      ## an NA here is an error rather than a zero
      yy <- yi; yy[!ok] <- 1L
      pick <- cbind(seq_len(nrow(Ub)), yy)
      obj <- obj + sum((w * (log(sm) + m - Ub[pick]))[ok])
      Gb <- P
      Gb[pick] <- Gb[pick] - 1
      Gb <- Gb * w
      Gb[!ok, ] <- 0
      G[, j] <- Gb
    }
  }
  list(objective = obj, grad = G)
}

#' log(1 + exp(u)) without overflowing
#' @keywords internal
#' @noRd
logspace_add0 <- function(u) ifelse(u > 30, u, log1p(exp(pmin(u, 30))))

#' A low-rank model with a loss chosen per column
#'
#' The same idea as principal components -- describe every row with a handful
#' of numbers -- but with a loss appropriate to each column's type rather than
#' squared error everywhere. A binary column gets logistic loss, a categorical
#' one multinomial across its levels, a count Poisson, and a numeric one
#' squared error, which is where this reduces to ordinary PCA.
#'
#' @section When this is worth the iterative fit:
#'
#' On all-numeric data it is not: with quadratic loss everywhere it *is* PCA,
#' and [ilm_reduce()] gets there in closed form. It earns its cost when
#' categorical columns matter, because FAMD reaches them by one-hot encoding
#' and applying squared loss to the indicators, which is a Gaussian
#' approximation to something that is not Gaussian, and can reconstruct a
#' category as -0.3.
#'
#' @section Missing values are not a special case:
#'
#' A cell that was not observed contributes nothing to the loss, so the fit
#' uses whatever is there and the reconstruction fills the rest in. That is why
#' the same machinery serves [ilm_impute()].
#'
#' @param data A data frame.
#' @param cols Columns to use; see [ilm_selection].
#' @param rank Number of dimensions.
#' @param loss Optional named character vector overriding the automatic choice
#'   per column: `"quadratic"`, `"logistic"`, `"multinomial"` or `"poisson"`.
#'   An **ordered** factor defaults to `"multinomial"`, which ignores the
#'   ordering but never invents a spacing; `"quadratic"` on its integer codes
#'   uses the ordering and does invent one, and is available if that is the
#'   trade you want.
#' @param lambda Ridge penalty on the factors, or `NULL` (the default) to
#'   choose it by holding out observed cells. It is not a detail: a row with
#'   few observed cells has as many scores as observations and fits them
#'   exactly, so without enough penalty the reconstruction of everything else
#'   in that row goes wherever the algebra sends it. On a 400 by 4 matrix with
#'   20% missing, held-out error ran 0.82 at `lambda = 0.1`, 0.68 at 2, and
#'   1.16 at 25 -- against 0.78 for an iterative SVD and 1.18 for column means.
#'   The best value depends on the size, the width and the missingness, which
#'   is why it is chosen rather than fixed.
#' @param weights Optional non-negative row weights. A row counted twice
#'   contributes twice to the loss, which is how [ilm_impute()] draws a
#'   bootstrap replicate without resampling the rows themselves and losing
#'   the ones it has to reconstruct.
#' @param maxit Maximum alternating passes.
#' @param tol Relative change in the objective at which to stop.
#' @param seed Random seed for the starting point.
#' @param progress Show a progress bar; see [ilm_progress_arg].
#' @return An object of class `"ilm_glrm"` with `scores` (one row per
#'   observation), `archetypes`, the per-column `loss`, the `objective` trace
#'   and a `fitted` data frame reconstructing each column on its own scale.
#' @seealso [ilm_reduce()] for the closed-form FAMD route and
#'   [ilm_reduce(method = "glrm")][ilm_reduce] to use this one inside it.
#' @references Udell, M., Horn, C., Zadeh, R. and Boyd, S. (2016). Generalized
#'   low rank models. *Foundations and Trends in Machine Learning* 9, 1-118.
#' @examples
#' set.seed(1); n <- 200
#' f <- rnorm(n)
#' d <- data.frame(a = f + rnorm(n, 0, .4), b = -f + rnorm(n, 0, .4),
#'                 g = factor(ifelse(f + rnorm(n, 0, .5) > 0, "hi", "lo")))
#' g <- ilm_glrm(d, rank = 1)
#' head(g$scores)
#' @export
ilm_glrm <- function(data, cols = NULL, rank = 2L, loss = NULL, lambda = NULL,
                     weights = NULL, maxit = 300L, tol = 1e-7, seed = 1L,
                     progress = NULL) {
  if (!is.data.frame(data))
    stop("`data` must be a data frame; it is ", class(data)[1], call. = FALSE)
  keep <- ilm_resolve_cols(data, cols)
  sub <- data[keep]
  auto <- vapply(sub, ilm_glrm_loss_of, "")
  if (anyNA(auto)) {
    bad <- names(auto)[is.na(auto)]
    message("ilm_glrm(): dropping column(s) with no loss to give them: ",
            paste(bad, collapse = ", "))
    sub <- sub[!is.na(auto)]; auto <- auto[!is.na(auto)]
  }
  if (!ncol(sub)) stop("no usable columns", call. = FALSE)
  if (!is.null(loss)) {
    ok <- c("quadratic", "logistic", "multinomial", "poisson")
    bad <- setdiff(loss, ok)
    if (length(bad))
      stop("unknown loss: ", paste(unique(bad), collapse = ", "),
           ". Available: ", paste(ok, collapse = ", "), ".", call. = FALSE)
    unknown <- setdiff(names(loss), names(sub))
    if (length(unknown))
      stop("`loss` names column(s) that are not being used: ",
           paste(unknown, collapse = ", "), call. = FALSE)
    auto[names(loss)] <- unname(loss)
  }
  n <- nrow(sub)
  rank <- max(1L, as.integer(rank))
  if (is.null(weights)) weights <- rep(1, n)
  if (length(weights) != n || any(!is.finite(weights)) || any(weights < 0))
    stop("`weights` must be one non-negative finite number per row",
         call. = FALSE)
  if (is.null(lambda)) {
    lambda <- ilm_glrm_lambda(sub, as.list(auto), rank, maxit, tol, seed)
  } else if (!is.numeric(lambda) || length(lambda) != 1L || lambda < 0)
    stop("`lambda` must be a single non-negative number, or NULL to choose ",
         "one by cross-validation.", call. = FALSE)
  enc <- ilm_glrm_encode(sub, as.list(auto))
  d <- enc$d
  if (rank >= min(n, d))
    stop("`rank` is ", rank, " but the encoded data has ", d, " columns over ",
         n, " rows; a low-rank model needs a rank below both.", call. = FALSE)

  ## Start from the SVD of a crude numeric encoding. It costs one
  ## decomposition and lands the alternating fit in the right basin; from a
  ## random start the multinomial blocks routinely settle in a worse one.
  set.seed(seed)
  A0 <- matrix(0, n, d)
  for (b in enc$blocks) {
    if (b$loss %in% c("quadratic", "poisson")) {
      v <- b$target; v[is.na(v)] <- 0; A0[, b$cols] <- v
    } else if (b$loss == "logistic") {
      v <- b$target; v[is.na(v)] <- mean(v, na.rm = TRUE)
      A0[, b$cols] <- 2 * v - 1
    } else {
      yi <- b$target
      M <- matrix(0, n, length(b$cols))
      okr <- !is.na(yi)
      M[cbind(which(okr), yi[okr])] <- 1
      A0[, b$cols] <- sweep(M, 2L, colMeans(M), "-")
    }
  }
  sv <- svd(A0, nu = rank, nv = rank)
  X <- sv$u[, seq_len(rank), drop = FALSE] *
    rep(sqrt(sv$d[seq_len(rank)]), each = n)
  Y <- t(sv$v[, seq_len(rank), drop = FALSE] *
           rep(sqrt(sv$d[seq_len(rank)]), each = d))
  mu <- rep(0, d)

  ## Alternating descent with an adaptive step. A step is PROPOSED, evaluated,
  ## and kept only if the objective fell; otherwise the step is halved and the
  ## point is left alone.
  ##
  ## Convergence is declared only after an ACCEPTED step whose decrease was
  ## small. That is not fussiness: a first version tested the objective at the
  ## top of each pass, so the pass following a rejected step saw the same
  ## point, computed the same objective, and read "no change" as convergence.
  ## The very first rejection ended the fit, after five iterations, well short
  ## of the optimum.
  step <- 1 / max(1, sqrt(n))
  pen <- function(Xm, Ym) 0.5 * lambda * (sum(Xm^2) + sum(Ym^2))
  g <- ilm_glrm_grad(X %*% Y + rep(mu, each = n), enc, weights)
  cur <- g$objective + pen(X, Y)
  obj <- cur
  pb <- ilm_progress(maxit, progress, "fitting the low-rank model")
  for (it in seq_len(maxit)) {
    Xn <- X - step * (g$grad %*% t(Y) + lambda * X)
    gx <- ilm_glrm_grad(Xn %*% Y + rep(mu, each = n), enc, weights)
    Yn <- Y - step * (t(Xn) %*% gx$grad + lambda * Y)
    mun <- mu - step * colSums(gx$grad) / max(sum(weights), 1)
    gn <- ilm_glrm_grad(Xn %*% Yn + rep(mun, each = n), enc, weights)
    new <- gn$objective + pen(Xn, Yn)
    if (is.finite(new) && new < cur) {
      delta <- cur - new
      X <- Xn; Y <- Yn; mu <- mun; g <- gn; cur <- new
      obj <- c(obj, new)
      step <- step * 1.1
      if (delta < tol * max(1, abs(cur))) break
    } else {
      step <- step / 2
      if (step < 1e-14) break
    }
    pb$tick(it)
  }
  pb$done()

  U <- X %*% Y + rep(mu, each = n)
  fitted <- ilm_glrm_reconstruct(U, enc, names(sub))
  colnames(X) <- paste0("dim", seq_len(rank))
  rownames(Y) <- colnames(X)
  ## the residual scale of each quadratic block, on the ORIGINAL scale, so a
  ## draw around the reconstruction can be made without refitting anything
  sg <- vapply(enc$blocks, function(b) {
    if (b$loss != "quadratic") return(NA_real_)
    a <- b$target; ok <- !is.na(a)
    sqrt(max(sum((U[ok, b$cols] - a[ok])^2) / max(sum(ok) - rank, 1), 0)) *
      b$scale
  }, 0)
  structure(list(scores = as.data.frame(X), archetypes = Y, offset = mu,
                 linear_predictor = U, sigma = sg,
                 loss = auto, objective = obj, rank = rank, lambda = lambda,
                 fitted = fitted, encoding = enc, columns = names(sub),
                 iterations = length(obj) - 1L,
                 converged = length(obj) - 1L < maxit),
            class = "ilm_glrm")
}

#' Turn the fitted linear predictor back into columns on their own scales
#'
#' @keywords internal
#' @noRd
ilm_glrm_reconstruct <- function(U, enc, nms) {
  out <- vector("list", length(enc$blocks))
  names(out) <- names(enc$blocks)
  for (v in names(enc$blocks)) {
    b <- enc$blocks[[v]]; j <- b$cols
    out[[v]] <- if (b$loss == "quadratic") U[, j] * b$scale + b$centre
    else if (b$loss == "poisson") exp(pmin(U[, j], 30))
    else if (b$loss == "logistic")
      factor(b$levels[(stats::plogis(U[, j]) > 0.5) + 1L], levels = b$levels)
    else factor(b$levels[max.col(U[, j, drop = FALSE], ties.method = "first")],
                levels = b$levels)
  }
  as.data.frame(out, stringsAsFactors = FALSE)
}

#' @export
print.ilm_glrm <- function(x, ...) {
  cat(sprintf("Generalized low rank model: rank %d over %d columns\n",
              x$rank, length(x$columns)))
  tb <- table(x$loss)
  cat("  losses: ",
      paste(sprintf("%s x%d", names(tb), as.integer(tb)), collapse = ", "),
      "\n", sep = "")
  cat(sprintf("  %d iterations, objective %.4f%s\n", x$iterations,
              utils::tail(x$objective, 1L),
              if (x$converged) "" else "  [did NOT converge]"))
  if (!x$converged)
    cat("  The objective was still moving when maxit was reached. Raise",
        " maxit,\n  or lambda if it is wandering.\n", sep = "")
  invisible(x)
}

#' Choose the ridge penalty by holding out observed cells
#'
#' The same idea the rank selection in ilm_impute() uses, on a parameter whose
#' effect is smooth -- which matters, because that one lands on an erratic
#' curve and this one does not: across the grid below the held-out error
#' changes by a few per cent between neighbouring values and has a single
#' minimum, so being one grid point out costs almost nothing.
#'
#' @keywords internal
#' @noRd
ilm_glrm_lambda <- function(sub, loss, rank, maxit, tol, seed,
                            grid = c(0.1, 0.5, 1, 2, 5, 10), folds = 3L) {
  n <- nrow(sub)
  obs <- which(!is.na(as.matrix(
    data.frame(lapply(sub, function(z) as.numeric(as.factor(z)))))))
  if (length(obs) < 60L) return(1)
  set.seed(seed)
  fold <- sample(rep_len(seq_len(folds), length(obs)))
  err <- vapply(grid, function(lam) {
    e <- 0; m <- 0L
    for (f in seq_len(folds)) {
      hold <- obs[fold == f]
      sh <- sub
      rr <- ((hold - 1L) %% n) + 1L
      cc <- ((hold - 1L) %/% n) + 1L
      for (j in unique(cc)) sh[rr[cc == j], j] <- NA
      fit <- tryCatch(
        ilm_glrm(sh, rank = rank, loss = unlist(loss), lambda = lam,
                 maxit = min(maxit, 200L), tol = tol, seed = seed,
                 progress = FALSE),
        error = function(z) NULL)
      if (is.null(fit)) next
      ## scored on a common scale: squared error for a numeric column after
      ## standardising, and a miss counted as 1 for a category, so no single
      ## column's units decide the winner
      for (j in seq_along(sub)) {
        i <- rr[cc == j]
        if (!length(i)) next
        tv <- sub[[j]][i]; pv <- fit$fitted[[j]][i]
        if (is.numeric(tv)) {
          sc <- stats::sd(sub[[j]], na.rm = TRUE)
          if (!is.finite(sc) || sc <= 0) sc <- 1
          e <- e + sum(((as.numeric(pv) - tv) / sc)^2, na.rm = TRUE)
        } else {
          e <- e + sum(as.character(pv) != as.character(tv), na.rm = TRUE)
        }
        m <- m + length(i)
      }
    }
    if (m) e / m else Inf
  }, 0)
  if (all(!is.finite(err))) return(1)
  grid[which.min(err)]
}

#' Draw imputed values from a fitted generalized low rank model
#'
#' The reconstruction is a fitted value; an imputation has to be a DRAW around
#' it, or every one of the m datasets is identical and the pooled variance is
#' the variance of a single fit. A numeric cell gets Gaussian noise at the
#' block's residual scale; a category is SAMPLED from its fitted probabilities
#' rather than set to the most likely level, which would make the imputations
#' agree with each other far more than the data support.
#'
#' @keywords internal
#' @noRd
ilm_glrm_draw <- function(fit, data, which_cols) {
  out <- data
  U <- fit$linear_predictor
  for (v in which_cols) {
    b <- fit$encoding$blocks[[v]]
    if (is.null(b)) next
    mi <- is.na(data[[v]])
    if (!any(mi)) next
    if (b$loss == "quadratic") {
      pred <- U[mi, b$cols] * b$scale + b$centre
      sg <- fit$sigma[[v]]
      out[[v]][mi] <- pred + stats::rnorm(sum(mi), 0, if (is.na(sg)) 0 else sg)
    } else if (b$loss == "poisson") {
      out[[v]][mi] <- stats::rpois(sum(mi), exp(pmin(U[mi, b$cols], 30)))
    } else if (b$loss == "logistic") {
      p <- stats::plogis(U[mi, b$cols])
      out[[v]][mi] <- b$levels[stats::rbinom(sum(mi), 1L, p) + 1L]
    } else {
      Ub <- U[mi, b$cols, drop = FALSE]
      P <- exp(Ub - apply(Ub, 1L, max)); P <- P / rowSums(P)
      pick <- apply(P, 1L, function(pr) sample.int(length(pr), 1L, prob = pr))
      out[[v]][mi] <- b$levels[pick]
    }
  }
  out
}
