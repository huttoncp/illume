## ---------------------------------------------------------------------------
## Joint draws of every parameter of a fitted model.
##
## Everything a fit estimated -- fixed effects, variance and dispersion
## parameters, each group's random effects, each cell of a correlation over
## time -- is drawn together from the joint precision, so the draws carry the
## correlations between them: a smooth's penalised part against its null
## space, a group's effect against the intercept, a forecast state against
## the walk's variance. Laid out by name, and with the variance parameters
## also returned on their natural scale through the one transform the fit
## itself uses, so nothing built on the draws re-derives illume's
## parameterisation.
## ---------------------------------------------------------------------------

#' Joint draws of every parameter of a fitted model
#'
#' Draws the fixed effects, the variance and dispersion parameters, the random
#' effects and the cells of a correlation over time together, from the
#' Gaussian approximation to their joint distribution that the fit's standard
#' errors come from.
#'
#' @details
#' **The distribution.** The draws are centred on the estimates and the
#' conditional modes, with covariance the inverse of the joint precision --
#' the Laplace approximation's curvature in every parameter at once. For the
#' fixed parameters that is exactly `vcov(fit, full = TRUE)`; for the random
#' effects it adds the uncertainty the fixed parameters bring to the
#' conditional SDs of [ilm_ranef()].
#'
#' **At a boundary.** When a random term's covariance sits at the edge of its
#' range, the fit holds the direction the data cannot resolve at its estimate
#' (see `fit$hessian_held`), and so do the draws: exactly, by conditioning the
#' joint distribution on those directions. `held` says how many.
#'
#' **`given = "theta"`** holds the variance parameters of the random terms and
#' of a correlation over time at their estimates, and draws everything else
#' from its distribution given them -- the approach of the `merTools`
#' package. The dispersion and the family's other parameters are still drawn.
#'
#' **`given = "parameters"`** holds every parameter at its estimate -- the
#' fixed effects, the variance, dispersion and correlation parameters, and the
#' family's others -- and draws only the random effects and the cells of a
#' correlation over time, from their distribution given the parameters: the
#' curvature of the fit's inner optimisation at the conditional modes, which
#' for a gaussian model is that distribution exactly. Their SDs are
#' [ilm_ranef()]'s `sd`, except under REML, where `ilm_ranef()` integrates
#' over the fixed effects and these hold them, so these are smaller. The
#' variance components in `natural` are then those at the estimate. A model
#' without random effects has nothing left to draw: given its parameters the
#' distribution is a point mass, and every draw sits at the estimate.
#'
#' **`blocks`** returns only the named blocks of the parameter vector -- say
#' `c("beta", "bvec")` -- which saves memory; the draw is joint either way.
#'
#' A fit made without `joint = TRUE` has no joint precision stored, and it is
#' formed here from the fit's compiled objective; a fit read back from disk no
#' longer has that, and has to be refitted. A fit without random effects
#' needs neither: with nothing integrated out, its draws come from
#' `vcov(fit, full = TRUE)`.
#'
#' @param object A fitted `"ilm_model"`.
#' @param nsim Number of draws.
#' @param seed Optional random seed.
#' @param given `"none"` to draw everything, `"theta"` to hold the variance
#'   parameters at their estimates, or `"parameters"` to hold every parameter
#'   and draw only the random effects and the cells.
#' @param blocks Optional character vector of the blocks to return, named as in
#'   `names(object$obj$env$par)`: `"beta"`, `"theta"`, `"bvec"`, `"logdisp"`,
#'   `"B_ar"`, `"lchol_ar"`, `"rho_raw"` and so on.
#' @param natural Logical. Also transform each draw's variance parameters to
#'   the natural scale, as [ilm_varcorr()] reports them.
#' @return An object of class `"ilm_draws"`, a list with
#'   \describe{
#'     \item{`draws`}{a matrix, one row per parameter and one column per draw,
#'       on the fit's internal scale (log standard deviations and the like),
#'       with rows named by block.}
#'     \item{`mode`}{the centre: the estimates and conditional modes.}
#'     \item{`map`}{a data frame saying what each row is: `row`, `block`,
#'       `term` (a coefficient's or parameter's name, or a random term's), and
#'       for a random effect its `level` (the group's label), `dim` and, for a
#'       correlation over time, `cell` (the row of [ilm_cells()]).}
#'     \item{`natural`}{with `natural = TRUE`, each draw's variance
#'       parameters on their natural scale: `re`, one array per random term of
#'       its covariance by draw; `latent`, a correlation over time's
#'       covariance by draw with its `rho` (and `range`) or, for a random walk,
#'       `var_per_time`; and `dispersion`, one row per dispersion parameter.
#'       At the estimate these are exactly [ilm_varcorr()]'s.}
#'     \item{`held`}{the number of directions held at a boundary, and the
#'       terms they belong to; none under `given = "theta"` or
#'       `"parameters"`, which hold the variance parameters whole.}
#'     \item{`given`}{as asked.}
#'   }
#' @seealso [ilm_ranef()], whose `row` indexes these draws' full layout;
#'   [ilm_varcorr()]; [ilm_cells()].
#' @examples
#' set.seed(1)
#' d <- data.frame(id = factor(rep(1:12, each = 6)), x = rnorm(72))
#' d$y <- 1 + 0.5 * d$x + rnorm(12)[d$id] + rnorm(72)
#' fit <- ilm_model(y ~ x + (1 | id), data = d, family = "gaussian",
#'                  verbose = FALSE)
#' dr <- ilm_draws(fit, nsim = 200, seed = 1)
#' dim(dr$draws)
#' head(dr$map)
#' @export
ilm_draws <- function(object, nsim = 1000L, seed = NULL,
                      given = c("none", "theta", "parameters"), blocks = NULL,
                      natural = TRUE) {
  if (!inherits(object, "ilm_model"))
    stop("`object` must be a fitted ilm_model, not ", class(object)[1],
         call. = FALSE)
  given <- match.arg(given)
  nsim <- as.integer(nsim)
  if (length(nsim) != 1L || is.na(nsim) || nsim < 1L)
    stop("`nsim` must be a positive whole number", call. = FALSE)
  Q <- ilm_joint_prec(object)
  mu <- ilm_full_par(object)
  rn <- names(mu)
  if (!identical(rownames(Q), rn))
    stop("internal: the joint precision is not laid out along the fit's ",
         "parameter vector", call. = FALSE)
  n <- length(mu)
  if (!is.null(blocks)) {
    bad <- setdiff(blocks, unique(rn))
    if (length(bad))
      stop("no block called ", paste(sQuote(bad), collapse = ", "),
           " in this fit; its blocks are ", paste(unique(rn), collapse = ", "),
           call. = FALSE)
  }
  keep <- if (is.null(blocks)) seq_len(n) else which(rn %in% blocks)

  ## what is held at the estimate: whole coordinates (the variance parameters
  ## under given = "theta", or the held terms' parameters for a fit made
  ## before the directions were stored), or directions within them
  hold <- rep(FALSE, n)
  if (given == "theta") hold <- rn %in% c("theta", "lchol_ar", "rho_raw")
  ## every parameter at its estimate, and the random effects and cells drawn
  ## given them: their block of the precision is the inner Hessian at the
  ## modes. Under REML the fixed effects sit among TMB's random parameters,
  ## and are held all the same.
  if (given == "parameters") hold <- !rn %in% c("bvec", "B_ar")
  dirs <- NULL; n_held <- 0L
  if (given == "none" && length(object$hessian_held)) {
    if (!is.null(object$hessian_dirs)) {
      fp <- setdiff(seq_len(n), object$obj$env$random)
      dirs <- matrix(0, n, ncol(object$hessian_dirs))
      dirs[fp, ] <- object$hessian_dirs
      n_held <- ncol(dirs)
    } else {
      hc <- ilm_held_coords(object, rn)
      hold <- hold | hc; n_held <- sum(hc)
    }
  }
  free <- which(!hold)
  ch <- NULL
  if (length(free)) {
    Qf <- Matrix::forceSymmetric(Q[free, free, drop = FALSE])
    ch <- tryCatch(Matrix::Cholesky(Qf, LDL = FALSE, perm = TRUE),
                   error = function(e) NULL)
    if (is.null(ch))
      stop("the joint precision is not positive definite, so there is no ",
           "distribution to draw from; see fit$checks", call. = FALSE)
  }
  ## a held direction: condition the draw on it, x - S A' (A S A')^-1 A x
  ## with S the covariance of the free coordinates
  if (!is.null(dirs)) {
    A <- t(dirs[free, , drop = FALSE])
    SA <- as.matrix(Matrix::solve(ch, t(A), system = "A"))
    ASA <- A %*% SA
  }

  if (!is.null(seed)) set.seed(seed)
  out <- base::matrix(mu[keep], length(keep), nsim,
                      dimnames = list(rn[keep], NULL))
  nat <- list()
  step <- 250L
  for (s0 in seq(1L, nsim, by = step)) {
    cols <- s0:min(nsim, s0 + step - 1L)
    full <- base::matrix(mu, n, length(cols))
    if (length(free)) {
      Z <- base::matrix(stats::rnorm(length(free) * length(cols)),
                        length(free))
      X <- Matrix::solve(ch, Z, system = "Lt")
      X <- as.matrix(Matrix::solve(ch, X, system = "Pt"))
      if (!is.null(dirs)) X <- X - SA %*% solve(ASA, A %*% X)
      full[free, ] <- full[free, ] + X
    }
    out[, cols] <- full[keep, , drop = FALSE]
    if (natural) nat[[length(nat) + 1L]] <- ilm_natural_draws(object, full)
  }

  structure(list(draws = out, mode = mu[keep],
                 map = ilm_draws_map(object, keep),
                 natural = if (natural) ilm_natural_bind(nat) else NULL,
                 held = list(n = n_held, terms = object$hessian_held),
                 given = given),
            class = "ilm_draws")
}

## The joint precision: stored with the fit, or formed from its compiled
## objective and the fixed-parameter Hessian its covariance came from.
#' @keywords internal
#' @noRd
ilm_joint_prec <- function(object) {
  Q <- object$jointPrecision
  if (is.null(Q)) Q <- object$sdr$jointPrecision
  if (!is.null(Q)) return(Q)
  ## With nothing integrated out -- a model without random effects, fitted by
  ## maximum likelihood -- the whole parameter vector is fixed, and its
  ## distribution is the one its standard errors come from: vcov(), with the
  ## n / (n - p) that makes a gaussian fit's agree with lm(). sdreport()
  ## forms a joint precision only when there is a random part, and its
  ## absence used to be taken for a fit read back from disk, so no such
  ## model could be drawn from.
  if (!length(object$obj$env$random)) {
    rn <- names(object$obj$env$par)
    V <- suppressWarnings(stats::vcov(object, full = TRUE))
    Q <- if (!is.null(V) && all(is.finite(V)))
      tryCatch(solve(V), error = function(e) NULL)
    if (is.null(Q))
      stop("the fit's Hessian is not positive definite, so there is no ",
           "distribution to draw from; see fit$checks", call. = FALSE)
    if (!identical(dim(Q), rep(length(rn), 2L)))
      stop("internal: the fit's covariance does not match its parameter ",
           "vector", call. = FALSE)
    dimnames(Q) <- list(rn, rn)
    return(Q)
  }
  Hf <- object$hessian_fixed
  if (is.null(Hf)) {
    if (!isTRUE(object$sdr$pdHess))
      stop("the fit's Hessian is not positive definite, so there is no joint ",
           "distribution to draw from; see fit$checks", call. = FALSE)
    Hf <- solve(object$sdr$cov.fixed)
  }
  err <- NULL
  s2 <- tryCatch(suppressWarnings(sdreport(object$obj, par.fixed = object$opt$par,
                                           hessian.fixed = Hf,
                                           getJointPrecision = TRUE)),
                 error = function(e) { err <<- conditionMessage(e); NULL })
  if (is.null(s2) || is.null(s2$jointPrecision))
    stop("the joint precision could not be formed from this fit's compiled ",
         "objective", if (!is.null(err)) paste0(" (", err, ")"), ". A fit ",
         "read back from disk has lost that objective and has to be ",
         "refitted; with joint = TRUE the precision is kept with the fit.",
         call. = FALSE)
  s2$jointPrecision
}

## For a fit made before the held directions were stored: hold the held
## terms' whole parameter blocks instead, which holds at least as much.
#' @keywords internal
#' @noRd
ilm_held_coords <- function(object, rn) {
  h <- rep(FALSE, length(rn))
  it <- which(rn == "theta")
  for (nm in object$hessian_held) {
    if (identical(nm, "ar")) { h[rn %in% c("lchol_ar", "rho_raw")] <- TRUE; next }
    k <- match(nm, names(object$re))
    if (is.na(k)) next
    h[it[(object$toff[k] + 1L):object$toff[k + 1L]]] <- TRUE
  }
  h
}

## A full-vector draw put in coef(object, full = TRUE) order, the order
## ilm_rebuild() reads: each block by its type label, which under REML also
## finds the fixed effects in the random block.
#' @keywords internal
#' @noRd
ilm_coef_order <- function(object, full) {
  tl <- names(object$opt$par); rn <- names(object$obj$env$par)
  v <- numeric(length(tl))
  for (t in unique(tl)) v[tl == t] <- full[rn == t]
  v
}

## What each drawn row is.
#' @keywords internal
#' @noRd
ilm_draws_map <- function(object, keep) {
  rn <- names(object$obj$env$par)
  tl <- names(object$opt$par); pn <- object$pnames
  term <- rep(NA_character_, length(rn))
  for (t in unique(rn[rn %in% tl])) term[rn == t] <- pn[tl == t]
  level <- rep(NA_character_, length(rn)); dim <- level
  cell <- rep(NA_integer_, length(rn))
  r <- ilm_ranef_table(object, sd = FALSE)
  ok <- !is.na(r$row)
  term[r$row[ok]] <- r$term[ok]
  level[r$row[ok]] <- as.character(r$level[ok])
  dim[r$row[ok]] <- r$dim[ok]
  cell[r$row[ok]] <- r$cell[ok]
  data.frame(row = seq_along(keep), block = rn[keep], term = term[keep],
             level = level[keep], dim = dim[keep], cell = cell[keep],
             stringsAsFactors = FALSE)
}

## The variance parameters of a block of draws -- the columns of `full`, laid
## out along object$obj$env$par -- on their natural scale, every draw at once.
## Each covariance factor is filled from the draws' rows where the fit's own
## numeric builder puts each parameter, and the correlation over time and the
## dispersion go through the transforms ilm_rebuild() uses, so a draw at the
## mode gives exactly ilm_varcorr(). The labels are ilm_natural()'s at the
## estimate, laid out as the draws have always had them.
#' @keywords internal
#' @noRd
ilm_natural_draws <- function(object, full) {
  rn <- names(object$obj$env$par); m <- ncol(full); C <- object$C
  nat0 <- ilm_natural(object)
  lab <- function(V) c(dimnames(as.matrix(V[, ])), list(NULL))
  th <- full[rn == "theta", , drop = FALSE]
  re <- lapply(seq_along(object$re), function(k) {
    P <- th[object$toff[k] + seq_len(object$npc[k]), , drop = FALSE]
    S <- switch(object$ty[k],
      rr = ilm_prod_draws(P, function(v) ilm_mkLam_num(v, C, object$rk[k])),
      diag = ilm_prod_draws(P, function(v) ilm_mkD_num(v, C)),
      ilm_prod_draws(P, function(v) ilm_mkL_num(v, C)))
    d <- object$dk[k]
    if (C == 1L && d > 1L) {
      Pd <- th[object$toff[k] + object$npc[k] + seq_len(object$npd[k]), ,
               drop = FALSE]
      S <- ilm_prod_draws(Pd, function(v) ilm_mkLd_num(v, d, object$dcor[k])) *
        rep(S[1L, 1L, ], each = d * d)
    }
    dimnames(S) <- lab(nat0$re[[k]])
    S
  })
  names(re) <- names(nat0$re)

  latent <- NULL
  if (!is.null(object$ar)) {
    lt <- nat0$latent; rw <- identical(lt$type, "rw1")
    Sa <- ilm_prod_draws(full[rn == "lchol_ar", , drop = FALSE],
                         function(v) ilm_mkL_num(v, C))
    dimnames(Sa) <- lab(if (rw) lt$var_per_time else lt$Sigma)
    latent <- list(type = lt$type, meaning = lt$meaning)
    if (rw) {
      latent$var_per_time <- Sa
    } else {
      latent$Sigma <- Sa
      rr <- full[rn == "rho_raw", , drop = FALSE]
      tr <- if (nrow(rr)) ilm_rho_from_raw(lt$type, rr[1L, ])
            else list(rho = rep(lt$rho, m), range = rep(lt$range, m))
      latent$rho <- tr$rho
      if (identical(lt$type, "car1")) latent$range <- tr$range
    }
  }

  disp <- NULL
  if (!is.null(nat0$dispersion)) {
    v0 <- nat0$dispersion$value
    has_dm <- !is.null(object$Zd) || isTRUE(object$disp_mu)
    D <- if (object$family$n_disp > 0L && !has_dm && any(rn == "logdisp"))
      ilm_disp_scale(object, full[rn == "logdisp", , drop = FALSE])
    else if (has_dm)
      ## a modelled dispersion is summarised by its median over the rows,
      ## which takes the whole rebuilt fit, one draw at a time
      vapply(seq_len(m), function(j) unname(
        ilm_rebuild(object, ilm_coef_order(object, full[, j]))$dispersion),
        numeric(length(v0)))
    else rep(v0, m)
    disp <- base::matrix(D, length(v0), m, dimnames = list(names(v0), NULL))
  }
  list(re = re, latent = latent, dispersion = disp)
}

## L L' for each column of P, with L = build(P[, s]). Where each parameter
## sits in L is read off the builder once -- a parameter that moves an entry
## from 1 to e is there on the log scale, one that moves it from 0 to 1 is
## taken as it is, and an entry no parameter moves is fixed -- and every
## product is then formed across the draws at once.
#' @keywords internal
#' @noRd
ilm_prod_draws <- function(P, build) {
  np <- nrow(P); m <- ncol(P)
  L0 <- build(numeric(np)); a <- nrow(L0); b <- ncol(L0)
  Lm <- base::matrix(as.vector(L0), a * b, m)
  for (k in seq_len(np)) {
    v <- numeric(np); v[k] <- 1
    w <- which(build(v) != L0)
    if (length(w) != 1L)
      stop("internal: a covariance parameter moves ", length(w),
           " entries of its factor", call. = FALSE)
    Lm[w, ] <- if (L0[w] == 1) exp(P[k, ]) else P[k, ]
  }
  if (m && !isTRUE(all.equal(as.vector(build(P[, 1L])), Lm[, 1L],
                             tolerance = 1e-12)))
    stop("internal: the covariance factor of the draws is not the builder's",
         call. = FALSE)
  V <- base::matrix(0, a * a, m)
  cj <- (seq_len(b) - 1L) * a
  for (i in seq_len(a)) for (j in seq_len(i)) {
    s <- colSums(Lm[i + cj, , drop = FALSE] * Lm[j + cj, , drop = FALSE])
    V[i + (j - 1L) * a, ] <- s
    V[j + (i - 1L) * a, ] <- s
  }
  array(V, c(a, a, m))
}

## Chunks of ilm_natural_draws() joined along the draws.
#' @keywords internal
#' @noRd
ilm_natural_bind <- function(chunks) {
  if (length(chunks) == 1L) return(chunks[[1L]])
  first <- chunks[[1L]]
  along <- function(get) {
    a <- get(first)
    out <- array(unlist(lapply(chunks, function(x) as.vector(get(x)))),
                 c(dim(a)[1:2], sum(vapply(chunks, function(x) dim(get(x))[3L],
                                           0L))))
    dimnames(out) <- dimnames(a)
    out
  }
  re <- lapply(names(first$re), function(nm) along(function(x) x$re[[nm]]))
  names(re) <- names(first$re)
  latent <- first$latent
  for (f in intersect(c("Sigma", "var_per_time"), names(latent)))
    latent[[f]] <- along(function(x) x$latent[[f]])
  for (f in intersect(c("rho", "range"), names(latent)))
    latent[[f]] <- unlist(lapply(chunks, function(x) x$latent[[f]]))
  disp <- if (!is.null(first$dispersion))
    do.call(cbind, lapply(chunks, `[[`, "dispersion"))
  list(re = re, latent = latent, dispersion = disp)
}

#' @export
print.ilm_draws <- function(x, ...) {
  cat("<ilm_draws>", ncol(x$draws), "joint draws of", nrow(x$draws),
      "parameters", switch(x$given,
        theta = "(variance parameters held at their estimates)",
        parameters = "(every parameter held at its estimate)",
        ""), "\n")
  tab <- table(factor(x$map$block, levels = unique(x$map$block)))
  cat(paste0("  ", names(tab), ": ", as.integer(tab), collapse = "\n"), "\n")
  if (x$held$n > 0L)
    cat("  held at a boundary:", x$held$n, "direction(s) of",
        paste(x$held$terms, collapse = ", "), "\n")
  invisible(x)
}
