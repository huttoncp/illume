## illume: simulation and the empirical Laplace-consistency check.
##
## TMB::checkConsistency() cannot be used here.  It relies on obj$simulate(),
## which redraws quantities produced by distribution calls inside the tape.
## ilm_model()'s likelihood has none: the response enters as a fixed indicator matrix
## and the random-effect densities are hand-written quadratic forms, so nothing
## in the tape is redrawable and obj$simulate() errors outright.  The check below
## does the same job directly -- simulate from the fit, refit, compare -- and
## reuses ilm_simulate(), which the package needs anyway for population-averaged
## marginal means, parametric-bootstrap LRTs and the residual diagnostics.

#' Symmetric square root of a covariance matrix
#'
#' Used instead of a Cholesky factor because it tolerates singular matrices,
#' which arise legitimately here -- a reduced-rank covariance is singular by
#' construction, and a smooth component penalised to zero produces one too.
#' `chol()` would fail on both.
#'
#' @param S A symmetric positive semi-definite matrix.
#' @return A symmetric matrix `A` with `A %*% A` equal to `S`.
#' @keywords internal
#' @noRd
ilm_msqrt <- function(S) {
  e <- eigen(S, symmetric = TRUE); v <- pmax(e$values, 0)
  e$vectors %*% diag(sqrt(v), nrow(S)) %*% t(e$vectors)
}

## One draw of the latent AR(1), CAR(1) or random-walk process, at the
## observed rows.
##
## Simulating must walk the same chain the likelihood scores, or every
## envelope built on top of it is calibrated against the wrong process. Shared
## by ilm_simulate(), which walks the fitted rows, and ilm_power(), which walks
## a simulated study's -- so the two cannot drift apart.
#' @keywords internal
#' @noRd
ilm_ar_draw <- function(ar, rho, Sar, C) {
  La <- ilm_msqrt(Sar)
  Ba <- matrix(0, ar$n_cell, C)
  if (identical(ar$type, "rw1")) {
    ## zero at each group's first cell, as the fit holds it, then independent
    ## steps with covariance gap * Sigma, in cell order like CAR(1)
    for (k in seq_along(ar$rest))
      Ba[ar$rest[k], ] <- Ba[ar$prev[k], ] +
        sqrt(ar$gap[k]) * (rnorm(C) %*% La)
  } else if (identical(ar$type, "car1")) {
    Ba[ar$first, ] <- matrix(rnorm(length(ar$first) * C), ncol = C) %*% La
    phi <- rho ^ ar$gap
    ## transitions are in cell order, so the predecessor is always already
    ## filled by the time its successor is reached
    for (k in seq_along(ar$rest))
      Ba[ar$rest[k], ] <- phi[k] * Ba[ar$prev[k], ] +
        sqrt(1 - phi[k]^2) * (rnorm(C) %*% La)
  } else {
    for (g in seq_len(ar$n_group)) {
      r0 <- (g - 1L) * ar$Tt
      Ba[r0 + 1L, ] <- rnorm(C) %*% La
      for (tt in 2:ar$Tt)
        Ba[r0 + tt, ] <- rho * Ba[r0 + tt - 1L, ] + sqrt(1 - rho^2) * (rnorm(C) %*% La)
    }
  }
  Ba[ar$idx, , drop = FALSE]
}

#' Rebuild the random-effects specification from a fitted model
#' @param fit A fitted `"ilm_model"` object.
#' @return A list suitable for passing back to [ilm_fit()].
#' @keywords internal
#' @noRd
ilm_re_list_of <- function(fit) {
  rl <- lapply(fit$re, function(e) {
    if (e$kind == "basis") list(basis = e$basis)
    else if (e$d > 1L)     list(group = e$group, Z = e$Z)
    else                   e$group
  })
  names(rl) <- names(fit$re); rl
}

#' Simulate new outcomes from a fitted model
#'
#' Draws fresh random effects from their estimated distributions, forms the
#' linear predictor, and samples a new response for each observation from the
#' model's family -- a category for a multinomial or ordinal outcome.
#'
#' @section Why this is central:
#' Four separate tools depend on it: population-averaged predictions,
#' parametric-bootstrap tests, the simulated envelopes in every residual
#' diagnostic, and [ilm_consistency()]. Simulation is how this package builds
#' reference distributions, rather than relying on theoretical ones that may not
#' hold.
#'
#' `TMB::checkConsistency()` cannot be used in its place. It relies on the
#' objective function being able to simulate its own data, but this model's
#' likelihood contains no distribution calls that TMB can invert -- the response
#' enters as a fixed matrix and the random-effect densities are written out by
#' hand -- so there is nothing in the computation graph for TMB to redraw.
#'
#' @param fit A fitted `"ilm_model"` object.
#' @param nsim Integer. Number of simulated datasets.
#' @param seed Integer. Random seed, for reproducibility.
#' @return An integer matrix with one column per simulated dataset, each holding
#'   category codes.
#' @seealso [ilm_consistency()], [ilm_pb_lrt()], [ilm_appraise()].
#' @export
ilm_simulate <- function(fit, nsim = 1L, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  C <- fit$C; J <- fit$J; Tc <- contr.sum(J); X <- fit$X; N <- nrow(X)
  K <- length(fit$re)
  fam <- fit$family
  mn  <- is.null(fam) || identical(fam$name, "multinomial")
  ## a category index is an integer; every other family returns a measurement,
  ## and storing those in an integer matrix would silently truncate them
  out <- if (mn) matrix(0L, N, nsim) else matrix(0, N, nsim)
  wts <- if (is.null(fit$weights)) rep(1, N) else fit$weights
  dsp <- if (!is.null(fit$dispersion)) log(unname(fit$dispersion)) else numeric(0)
  nms <- names(fit$re)
  for (s in seq_len(nsim)) {
    eta <- X %*% fit$beta
    for (k in seq_len(K)) {
      e <- fit$re[[k]]; nl <- fit$nlk[k]; d <- fit$dk[k]; w <- fit$wk[k]
      isrr <- identical(fit$re_struct[[k]]$type, "rr")
      Ad <- if (d > 1L) t(chol(fit$Sigma_d[[nms[k]]])) else matrix(1, 1, 1)
      Bc <- if (isrr) diag(1, w) else ilm_msqrt(fit$Sigma[[k]])
      B <- matrix(0, nl * d, w)
      for (g in seq_len(nl)) {
        Bg <- Ad %*% matrix(rnorm(d * w), d, w) %*% Bc      # d x w
        for (i in seq_len(d)) B[(i - 1L) * nl + g, ] <- Bg[i, ]
      }
      if (e$kind == "basis") {
        ctb <- e$basis %*% B
        if (isrr) ctb <- ctb %*% t(fit$Lambda[[k]])
        eta <- eta + ctb
      } else {
        for (i in seq_len(d)) {
          Bi <- B[((i - 1L) * nl + 1L):(i * nl), , drop = FALSE]
          ctb <- Bi[e$group, , drop = FALSE]
          if (isrr) ctb <- ctb %*% t(fit$Lambda[[k]])
          eta <- eta + e$Z[, i] * ctb
        }
      }
    }
    if (!is.null(fit$ar)) {
      eta <- eta + ilm_ar_draw(fit$ar, fit$rho, fit$Sigma[["ar"]], C)
    }
    lsig_sim <- if (!is.null(fit$Zd) || isTRUE(fit$disp_mu))
      log(ilm_disp_vec(fit)) else NULL
    rp_off <- if (!is.null(fit$rp)) as.vector(eta[, 1]) -
      as.vector(ilm_rcs(log(as.numeric(fit$y)), fit$rp$knots) %*%
                fit$beta[fit$rp$cols, 1L]) else NULL
    if (isTRUE(fit$ordinal)) {
      ## Draw the latent variable and see which band it falls in, which is the
      ## model's own account of where a category comes from -- and keeps the
      ## random effects and any AR structure already in eta.
      z <- as.numeric(eta[, 1]) + fit$family$qfun(stats::runif(N))
      out[, s] <- as.integer(rowSums(outer(z, as.numeric(fit$zeta), `>`))) + 1L
    } else if (mn) {
      P <- exp(eta %*% t(Tc)); P <- P / rowSums(P)
      out[, s] <- apply(P, 1, function(pr) sample.int(J, 1L, prob = pr))
    } else {
      yy <- if (!is.null(rp_off)) ilm_rp_sim(fit, off = rp_off)
            else as.numeric(fam$sim(eta, wts, dsp, logsig = lsig_sim))
      ## A zero part is not a post-hoc thinning of an ordinary count draw
      ## under a hurdle: the positives come from a distribution that cannot
      ## produce a zero, so they are redrawn rather than filtered.
      if (!is.null(fit$Zzi)) {
        mu <- as.numeric(fam$linkinv(eta[, 1]))
        dv <- if (!is.null(lsig_sim)) exp(lsig_sim) else
              if (length(dsp)) rep(exp(dsp[1]), length(mu)) else rep(1, length(mu))
        yy <- ilm_zi_rng(fit, mu, dv, ilm_zi_p(fit), yy)
      }
      out[, s] <- ilm_censor_apply(fit$censor, yy)
    }
  }
  out
}

#' Interpretable summary of a fit, for comparison across refits
#'
#' Reduces a fitted model to fixed effects plus per-term category standard
#' deviations. Small enough to return from a parallel worker, and on a scale that
#' can be compared directly across replicates.
#'
#' @param fit A fitted `"ilm_model"` object.
#' @return A named numeric vector.
#' @keywords internal
#' @noRd
ilm_summarise_fit <- function(fit) {
  v <- as.vector(fit$beta)
  nm <- paste0("beta", seq_along(v))
  for (k in seq_along(fit$Sigma)) {
    sdv <- sqrt(diag(fit$Sigma[[k]]))
    v <- c(v, sdv); nm <- c(nm, paste0("sd_", names(fit$Sigma)[k], seq_along(sdv)))
  }
  if (!is.na(fit$rho)) { v <- c(v, fit$rho); nm <- c(nm, "rho") }
  setNames(v, nm)
}

#' Check whether the model can recover itself
#'
#' Simulates datasets from the fitted model, refits each one, and asks whether
#' the re-estimates centre on the original values. A systematic offset indicates
#' that the Laplace approximation is biased for this model and this amount of
#' data.
#'
#' @section How to read the result:
#' The **refit rate comes first and matters most**. If many replicates fail to
#' refit, the model cannot reliably recover itself from data it generated, and
#' the replicates that did converge are a success-conditioned sample -- their
#' apparent lack of bias understates the problem rather than excusing it. A low
#' refit rate is a finding, not a technical hiccup.
#'
#' Parameters at the **zero boundary** are reported separately from bias. A
#' variance or loading fitted at essentially zero cannot recentre on itself:
#' re-estimates are bounded below by zero, so their mean must exceed it. That is
#' arithmetic, not approximation error, and the usual symmetric test does not
#' apply. Such a parameter is flagged `"BOUNDARY"`, which tells you the component
#' is unsupported by the data rather than that the fitting method failed.
#'
#' With many parameters tested, some will exceed a two-sigma threshold by chance;
#' the printed output states how many to expect.
#'
#' @param fit A fitted `"ilm_model"` object.
#' @param B Integer. Simulated datasets. Detectable bias scales roughly as
#'   `1/sqrt(B)`, so `B = 30` resolves offsets of about 0.37 standard deviations.
#' @param seed Integer. Random seed.
#' @param ncores Integer. Worker processes. Every replicate is simulated up front
#'   on the main process, so results do not depend on this.
#' @param verbose Logical. Print the table.
#' @return Invisibly, a list with the comparison table, the raw draws, the number
#'   of usable refits and the refit rate.
#' @references
#' Joe, H. (2008). Accuracy of Laplace approximation for discrete response mixed
#' models. *Computational Statistics & Data Analysis*, 52(12), 5066--5074.
#'
#' Self, S. G., & Liang, K.-Y. (1987). Asymptotic properties of maximum
#' likelihood estimators and likelihood ratio tests under nonstandard
#' conditions. *Journal of the American Statistical Association*, 82(398),
#' 605--610. (On why boundary parameters need separate treatment.)
#' @export
ilm_consistency <- function(fit, B = 50L, seed = 1L, ncores = 1L, verbose = TRUE) {
  t0 <- proc.time()[3]
  ## Simulated up front on the master, so the parallel section is deterministic
  ## given its input: results do not depend on ncores.
  ys <- ilm_simulate(fit, B, seed)
  p0 <- ilm_summarise_fit(fit)
  TH <- matrix(NA_real_, length(p0), B, dimnames = list(names(p0), NULL))
  rr <- ilm_refit_many(fit, ys, ncores = ncores, restarts = 2L, verbose = verbose)
  for (b in seq_len(B)) if (!is.null(rr[[b]])) TH[, b] <- rr[[b]]
  nok <- sum(!is.na(TH[1, ])); conv_rate <- nok / B
  mu <- rowMeans(TH, na.rm = TRUE); sdv <- apply(TH, 1, sd, na.rm = TRUE)
  se <- sdv / sqrt(nok); bias <- mu - p0; z <- bias / se
  ## NB: ifelse() takes its length from the TEST, so a scalar test would collapse
  ## the whole vector to length 1 and silently recycle it.
  st <- if (nok < 10) rep("INCONCLUSIVE", length(z))
        else ifelse(abs(z) >= 4, "FAIL", ifelse(abs(z) >= 2, "WARN", "OK"))
  ## BOUNDARY, not bias.  The sd_* entries are non-negative magnitudes (for an
  ## rr term they are |Lambda_j|).  One fitted at essentially zero cannot be
  ## recovered without upward bias -- resimulated estimates are bounded below by
  ## zero, so their mean must exceed it -- and the symmetric z-test is invalid
  ## there.  Observed at J = 10: sd_site9 fitted 0.0095, refit mean 0.0369,
  ## z = 5.44, which is the boundary and not the Laplace approximation.
  is_sd <- grepl("^sd_", names(p0))
  if (any(is_sd)) {
    big <- max(p0[is_sd], na.rm = TRUE)
    ## Small beside the largest, or below the 1e-3 at which the fit itself
    ## calls a variance zero. Beside itself alone, a model's only standard
    ## deviation is never small: a random intercept fitted at 0.0001 was
    ## reported as biased (z = 5.00, FAIL) and the Laplace approximation
    ## blamed, in a gaussian model, which has none.
    near0 <- is_sd & (p0 < 0.05 * big | p0 < 1e-3)
    st[near0 & st %in% c("WARN", "FAIL")] <- "BOUNDARY"
  }
  res <- data.frame(parameter = names(p0), fitted = round(p0, 4),
                    mean_refit = round(mu, 4), bias = round(bias, 4),
                    rel_bias = round(bias / abs(p0), 3), z = round(z, 2),
                    status = st, row.names = NULL)
  if (verbose) {
    cat(sprintf("\nconsistency check: %d of %d replicates refitted successfully (%.0f%%), %.0f s\n",
                nok, B, 100 * conv_rate, proc.time()[3] - t0))
    print(res, row.names = FALSE)
    nf <- sum(st == "FAIL"); nw <- sum(st == "WARN"); np <- length(z)
    ## The refit RATE is itself a diagnostic, and it comes first: if most
    ## replicates fail, the surviving ones are a success-conditioned sample and
    ## their bias estimates understate the problem rather than exonerating it.
    if (conv_rate < 0.8)
      cat(sprintf("\n>> UNSTABLE: only %d of %d datasets simulated FROM this fit could be\n   refitted. The model cannot reliably recover itself from its own data.\n   The rows above come only from replicates that converged, so they are\n   conditioned on success and understate any bias -- do not read them as\n   reassurance. Check the latent_budget diagnostic.\n", nok, B))
    nbd <- sum(st == "BOUNDARY")
    if (nbd)
      cat(sprintf("\n>> %d parameter(s) sit at the ZERO BOUNDARY (a variance or loading fitted\n   at essentially nothing). Their apparent upward bias is forced: a\n   non-negative statistic cannot recentre on zero. This is not evidence of a\n   Laplace problem -- but it does mean those components are unsupported by the\n   data and the term could be simplified.\n", nbd))
    if (nok < 10)
      cat("\n>> INCONCLUSIVE: too few replicates converged to judge bias at all.\n")
    else if (nf > 0)
      cat(sprintf("\n>> %d parameter(s) show systematic bias under resimulation.\n   The Laplace approximation is not reliable for this model/data combination;\n   check the latent_budget diagnostic, which governs exactly this.\n", nf))
    else if (nw > 0)
      cat(sprintf("\n>> %d parameter(s) borderline at |z| >= 2. With %d parameters tested,\n   about %.1f such flags are expected by chance even with no bias.\n   Raise B to sharpen the test.\n", nw, np, 0.0455 * np))
    else if (conv_rate >= 0.8)
      cat("\n>> no detectable bias: estimates recentre on the fitted values.\n")
  }
  invisible(list(table = res, draws = TH, n_ok = nok, conv_rate = conv_rate))
}
