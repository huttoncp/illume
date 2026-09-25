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
#' **`blocks`** returns only the named blocks of the parameter vector -- say
#' `c("beta", "bvec")` -- which saves memory; the draw is joint either way.
#'
#' A fit made without `joint = TRUE` has no joint precision stored, and it is
#' formed here from the fit's compiled objective; a fit read back from disk no
#' longer has that, and has to be refitted.
#'
#' @param object A fitted `"ilm_model"`.
#' @param nsim Number of draws.
#' @param seed Optional random seed.
#' @param given `"none"` to draw everything, or `"theta"` to hold the variance
#'   parameters at their estimates.
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
#'       terms they belong to.}
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
                      given = c("none", "theta"), blocks = NULL,
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
  Qf <- Matrix::forceSymmetric(Q[free, free, drop = FALSE])
  ch <- tryCatch(Matrix::Cholesky(Qf, LDL = FALSE, perm = TRUE),
                 error = function(e) NULL)
  if (is.null(ch))
    stop("the joint precision is not positive definite, so there is no ",
         "distribution to draw from; see fit$checks", call. = FALSE)
  ## a held direction: condition the draw on it, x - S A' (A S A')^-1 A x
  ## with S the covariance of the free coordinates
  if (!is.null(dirs)) {
    A <- t(dirs[free, , drop = FALSE])
    SA <- as.matrix(Matrix::solve(ch, t(A), system = "A"))
    ASA <- A %*% SA
  }

  if (!is.null(seed)) set.seed(seed)
  out <- matrix(mu[keep], length(keep), nsim,
                dimnames = list(rn[keep], NULL))
  nat <- if (natural) vector("list", nsim) else NULL
  step <- 250L
  for (s0 in seq(1L, nsim, by = step)) {
    cols <- s0:min(nsim, s0 + step - 1L)
    Z <- matrix(stats::rnorm(length(free) * length(cols)), length(free))
    X <- Matrix::solve(ch, Z, system = "Lt")
    X <- as.matrix(Matrix::solve(ch, X, system = "Pt"))
    if (!is.null(dirs)) X <- X - SA %*% solve(ASA, A %*% X)
    full <- matrix(mu, n, length(cols))
    full[free, ] <- full[free, ] + X
    out[, cols] <- full[keep, , drop = FALSE]
    if (natural)
      for (j in seq_along(cols))
        nat[[cols[j]]] <- ilm_natural(object, ilm_coef_order(object, full[, j]))
  }

  structure(list(draws = out, mode = mu[keep],
                 map = ilm_draws_map(object, keep),
                 natural = if (natural) ilm_natural_stack(nat) else NULL,
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
  Hf <- object$hessian_fixed
  if (is.null(Hf)) {
    if (!isTRUE(object$sdr$pdHess))
      stop("the fit's Hessian is not positive definite, so there is no joint ",
           "distribution to draw from; see fit$checks", call. = FALSE)
    Hf <- solve(object$sdr$cov.fixed)
  }
  s2 <- tryCatch(suppressWarnings(sdreport(object$obj, par.fixed = object$opt$par,
                                           hessian.fixed = Hf,
                                           getJointPrecision = TRUE)),
                 error = function(e) NULL)
  if (is.null(s2) || is.null(s2$jointPrecision))
    stop("this fit's compiled objective is gone -- it was read back from ",
         "disk? -- so its joint precision cannot be formed. Refit it, with ",
         "joint = TRUE to keep the precision with the fit.", call. = FALSE)
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
    if (identical(nm, "dispersion")) { h[rn == "logdisp"] <- TRUE; next }
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

## ilm_natural() of each draw, stacked into arrays by draw.
#' @keywords internal
#' @noRd
ilm_natural_stack <- function(nat) {
  first <- nat[[1L]]; ns <- length(nat)
  stack <- function(get) {
    m <- get(first)
    a <- array(vapply(nat, function(x) as.numeric(get(x)), numeric(length(m))),
               c(dim(as.matrix(m)), ns))
    dimnames(a) <- c(dimnames(as.matrix(m)), list(NULL))
    a
  }
  re <- lapply(names(first$re), function(nm) stack(function(x) x$re[[nm]][, ]))
  names(re) <- names(first$re)
  latent <- NULL
  if (!is.null(first$latent)) {
    lt <- first$latent
    latent <- list(type = lt$type, meaning = lt$meaning)
    if (identical(lt$type, "rw1")) {
      latent$var_per_time <- stack(function(x) x$latent$var_per_time[, ])
    } else {
      latent$Sigma <- stack(function(x) x$latent$Sigma[, ])
      latent$rho <- vapply(nat, function(x) x$latent$rho, 0)
      if (identical(lt$type, "car1"))
        latent$range <- vapply(nat, function(x) x$latent$range, 0)
    }
  }
  disp <- if (!is.null(first$dispersion))
    matrix(vapply(nat, function(x) x$dispersion$value,
                  numeric(length(first$dispersion$value))),
           nrow = length(first$dispersion$value),
           dimnames = list(names(first$dispersion$value), NULL))
  list(re = re, latent = latent, dispersion = disp)
}

#' @export
print.ilm_draws <- function(x, ...) {
  cat("<ilm_draws>", ncol(x$draws), "joint draws of", nrow(x$draws),
      "parameters", if (identical(x$given, "theta"))
        "(variance parameters held at their estimates)" else "", "\n")
  tab <- table(factor(x$map$block, levels = unique(x$map$block)))
  cat(paste0("  ", names(tab), ": ", as.integer(tab), collapse = "\n"), "\n")
  if (x$held$n > 0L)
    cat("  held at a boundary:", x$held$n, "direction(s) of",
        paste(x$held$terms, collapse = ", "), "\n")
  invisible(x)
}
