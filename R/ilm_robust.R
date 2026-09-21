## ---------------------------------------------------------------------------
## Sandwich standard errors.
##
## The model-based covariance is (-H)^-1, and it is right only if the model is
## right. A sandwich replaces the middle of it with what the data actually did:
##
##     V  =  (-H)^-1  [ sum_g s_g s_g' ]  (-H)^-1
##
## where s_g is the sum of the score contributions in cluster g. If the mean
## structure is correct but the variance structure is not -- unequal spread,
## correlation within schools or firms or subjects -- the point estimates are
## still consistent and only the standard errors are wrong, and this fixes
## them without having to say what the right variance structure would be.
##
## Two things about it are routinely got wrong and both are handled here.
##
## THE REFERENCE DISTRIBUTION. The sandwich is consistent as the number of
## CLUSTERS grows, not as the number of observations does. With 15 clusters and
## a normal reference, a nominal 95% interval covers well below 95% however
## many rows each cluster holds. The default here is a t distribution on the
## Bell-McCaffrey degrees of freedom, which is usually well below G - 1.
##
## THE SMALL-SAMPLE CORRECTION. CR0 is the raw sandwich and is biased downwards.
## CR1 multiplies by a constant, which helps a little. CR2 adjusts each
## cluster's residuals by its own leverage, which is what actually works, and
## is the default.
##
## For a model with random effects the within-cluster correlation is already
## in the likelihood, and the per-cluster score is not the sum of per-row
## scores, so this does not apply -- see the error ilm_estfun() raises.
##
## References:
##   Liang, K.-Y. and Zeger, S. L. (1986). Longitudinal data analysis using
##     generalized linear models. Biometrika 73, 13-22.
##   Bell, R. M. and McCaffrey, D. F. (2002). Bias reduction in standard errors
##     for linear regression with multi-stage samples. Survey Methodology 28,
##     169-181.
##   Cameron, A. C. and Miller, D. L. (2015). A practitioner's guide to
##     cluster-robust inference. Journal of Human Resources 50, 317-372.
## ---------------------------------------------------------------------------

#' Per-observation score contributions
#'
#' For every family here the log-likelihood reaches a row only through its
#' linear predictor, so the score factors as `x_i * u_i` with `u_i` a single
#' number per row -- the derivative of that row's log-likelihood with respect
#' to its own linear predictor. That is the same quantity `sandwich::estfun()`
#' builds from the working residuals, which is why the two agree exactly.
#'
#' @param object A fitted [ilm_model()].
#' @return An `N x p` matrix, one row per observation and one column per
#'   fixed-effect coefficient.
#' @keywords internal
#' @noRd
ilm_estfun <- function(object) {
  if (length(object$re))
    stop("this model has random effects, so the within-cluster correlation is ",
         "already in the likelihood and a cluster's score is not the sum of ",
         "its rows' scores. A sandwich is for when you do NOT want to model ",
         "that correlation: fit the same mean structure without the random ",
         "effect and cluster on the same grouping, or keep this model and ",
         "trust it. ilm_check_variance() and ilm_variogram() say whether the ",
         "variance structure is holding up.", call. = FALSE)
  if (object$C > 1L)
    stop("a multinomial fit has one score per category dimension and no single ",
         "linear predictor to differentiate along; cluster-robust standard ",
         "errors are not available for it here.", call. = FALSE)
  if (!is.null(object$Zzi))
    stop("a zero part gives each row two linear predictors, and the score for ",
         "one is not separable from the other; cluster-robust standard errors ",
         "are not available for a zero-inflated or hurdle fit here.",
         call. = FALSE)

  fam <- if (!is.null(object$family)) object$family$name else "gaussian"
  X <- object$X
  y <- as.numeric(object$y)
  N <- nrow(X)
  w <- if (is.null(object$weights)) rep(1, N) else object$weights
  eta <- as.numeric(ilm_eta_hat(object, FALSE)[, 1])
  d <- ilm_disp_vec(object)

  ## u = d(log-likelihood of one row) / d(that row's linear predictor)
  u <- if (isTRUE(object$ordinal)) {
    zt <- as.numeric(object$zeta); J <- object$J
    pf <- object$family$pfun
    dens <- switch(object$family$link,
      logit   = function(z) stats::dlogis(z),
      probit  = function(z) stats::dnorm(z),
      cloglog = function(z) exp(z - exp(z)))
    yi <- as.integer(object$y)
    hi <- ifelse(yi < J, zt[pmin(yi, J - 1L)] - eta, Inf)
    lo <- ifelse(yi > 1L, zt[pmax(yi - 1L, 1L)] - eta, -Inf)
    p <- ifelse(yi < J, pf(hi), 1) - ifelse(yi > 1L, pf(lo), 0)
    fhi <- ifelse(yi < J, dens(hi), 0)
    flo <- ifelse(yi > 1L, dens(lo), 0)
    ## P = F(hi) - F(lo) and both cuts carry -eta, so dP/deta = f(lo) - f(hi)
    (flo - fhi) / pmax(p, .Machine$double.eps)
  } else switch(fam,
    ## the working residual, family by family; each is (y - mu) times whatever
    ## the link and the variance function contribute
    gaussian = (y - eta) / pmax(d, .Machine$double.eps)^2,
    poisson  = y - exp(eta),
    binomial = y - stats::plogis(eta),
    nbinom   = {
      mu <- exp(eta); k <- if (is.null(d)) 1 else d
      (y - mu) * k / (k + mu)
    },
    beta     = ilm_beta_u(y, eta, if (is.null(d)) 1 else d),
    stop("cluster-robust standard errors are not implemented for the ", fam,
         " family. They are available for gaussian, binomial, poisson, ",
         "nbinom, beta and the ordinal families.", call. = FALSE))
  X * (w * u)
}

#' Which rows belong to which cluster
#'
#' @keywords internal
#' @noRd
ilm_cluster_id <- function(object, cluster) {
  N <- nrow(object$X)
  if (inherits(cluster, "formula")) {
    v <- all.vars(cluster)
    ## The clustering variable is usually NOT a predictor -- that is rather
    ## the point of clustering on it -- so the model frame will not have it.
    ## Fall back to the data the model was fitted from, and drop the same rows
    ## the fit dropped, which is the step that otherwise goes wrong silently:
    ## a cluster vector one row longer than the design still recycles.
    dat <- object$model
    if (!all(v %in% names(dat))) {
      raw <- tryCatch(eval(object$call$data, environment(object$formula)),
                      error = function(e) NULL)
      if (is.null(raw) || !all(v %in% names(raw)))
        stop("`cluster` names ", paste(setdiff(v, names(dat)), collapse = ", "),
             ", which is neither in the model frame nor in the data the model ",
             "was fitted from. Pass the values directly as a vector, one per ",
             "row of the data.", call. = FALSE)
      om <- object$na.action
      raw <- raw[v]
      if (!is.null(om) && length(om)) raw <- raw[-as.integer(om), , drop = FALSE]
      if (nrow(raw) != N)
        stop("`cluster` was found in the original data but has ", nrow(raw),
             " rows against the model's ", N, ". Pass the values directly as ",
             "a vector to be sure they line up.", call. = FALSE)
      dat <- raw
    }
    f <- interaction(dat[v], drop = TRUE, sep = ":")
  } else {
    if (length(cluster) != N)
      stop("`cluster` has ", length(cluster), " values but the model was fitted ",
           "to ", N, " rows. If some rows were dropped as missing, pass a ",
           "formula instead -- ", "`cluster = ~ id` -- so the same rows are ",
           "used.", call. = FALSE)
    f <- factor(cluster)
  }
  if (anyNA(f))
    stop("`cluster` has missing values; every row has to belong to some ",
         "cluster for the sum over clusters to mean anything.", call. = FALSE)
  droplevels(f)
}

#' Cluster-robust and heteroskedasticity-robust covariance
#'
#' A sandwich covariance for the fixed effects: the model's own covariance on
#' the outside, and what the data actually did on the inside. Use it when the
#' mean structure is credible but the variance structure is not -- spread that
#' changes with a predictor, or observations that are correlated within a
#' school, firm, household or subject.
#'
#' @section Which correction, and why it matters:
#'
#' A sandwich is consistent in the number of **clusters**, not the number of
#' rows, and with few clusters the raw version is biased downwards:
#'
#' * `"CR0"` is the raw sandwich, with no correction.
#' * `"CR1"` multiplies it by `G/(G-1) * (N-1)/(N-p)`, the constant Stata
#'   applies. It is better than nothing and does not depend on the design.
#' * `"CR2"` (the default) adjusts each cluster's residuals by that cluster's
#'   own leverage, after Bell and McCaffrey. It is the one that works when the
#'   clusters are unbalanced, which is when the problem is worst.
#'
#' With `cluster = NULL` each row is its own cluster and these become the
#' familiar heteroskedasticity-robust HC0, HC1 and HC2.
#'
#' @section How few is too few:
#'
#' There is no threshold that makes the problem go away, but below roughly 40
#' clusters the correction and the reference distribution both start to matter
#' a great deal, and below about 15 no adjustment reliably rescues the interval
#' -- the estimator is being asked to learn a covariance from a dozen numbers.
#' [ilm_robust()] reports the count and says so. The honest alternative there
#' is a design with more clusters, or a model that says what the correlation is
#' rather than working around it.
#'
#' @param object A fitted [ilm_model()] with no random effects.
#' @param cluster A one-sided formula naming variables in the model frame
#'   (`~ school`, or `~ school + year` for their interaction), or a vector with
#'   one entry per row. `NULL` treats every row as its own cluster.
#' @param type `"CR2"`, `"CR1"` or `"CR0"`; see above.
#' @return A `p x p` covariance matrix, with the number of clusters and the
#'   type in its attributes.
#' @seealso [ilm_robust()] for a coefficient table built on it.
#' @references Bell, R. M. and McCaffrey, D. F. (2002). Bias reduction in
#'   standard errors for linear regression with multi-stage samples.
#'   *Survey Methodology* 28, 169-181.
#'
#'   Cameron, A. C. and Miller, D. L. (2015). A practitioner's guide to
#'   cluster-robust inference. *Journal of Human Resources* 50, 317-372.
#' @examples
#' set.seed(1); G <- 40; n <- 10
#' d <- data.frame(g = factor(rep(seq_len(G), each = n)), x = rnorm(G * n))
#' d$y <- 1 + 0.5 * d$x + rep(rnorm(G, 0, 1.5), each = n) + rnorm(G * n)
#' f <- ilm_model(y ~ x, data = d, family = "gaussian", verbose = FALSE)
#' sqrt(diag(vcov(f)))                          # model-based, too small here
#' sqrt(diag(ilm_vcov_cluster(f, ~ g)))         # cluster-robust
#' @export
ilm_vcov_cluster <- function(object, cluster = NULL,
                             type = c("CR2", "CR1", "CR0")) {
  type <- match.arg(type)
  if (!inherits(object, "ilm_model"))
    stop("`object` must be a fitted ilm_model, not ", class(object)[1],
         call. = FALSE)
  S <- ilm_estfun(object)
  X <- object$X
  N <- nrow(X); p <- ncol(X)
  ## The bread comes from vcov(), not from sdr$cov.fixed directly, because
  ## the two differ by n / (n - p) for a gaussian fit and ilm_estfun() divides
  ## by the SAME dispersion vcov() is scaled to. A sandwich cancels the
  ## dispersion exactly when both halves use one value for it and silently
  ## does not when they use two -- which cost about 2% on every standard error
  ## before the scales were made to match.
  bread <- as.matrix(suppressWarnings(stats::vcov(object)))
  bread <- bread[seq_len(p), seq_len(p), drop = FALSE]

  g <- if (is.null(cluster)) factor(seq_len(N)) else
    ilm_cluster_id(object, cluster)
  G <- nlevels(g)
  if (G < 2L)
    stop("a sandwich needs at least two clusters; this one has ", G, ".",
         call. = FALSE)

  if (type == "CR2") {
    ## Bell-McCaffrey: within each cluster, pre-multiply the scores by
    ## (I - H_gg)^(-1/2), where H_gg is that cluster's own block of the hat
    ## matrix. A cluster with high leverage has its residuals shrunk towards
    ## zero by the fit, and this is what undoes that -- which is exactly the
    ## bias CR0 suffers, and why a constant multiplier cannot fix it.
    W <- ilm_estfun_weight(object)
    Xw <- X * sqrt(W)
    B <- bread
    for (lv in levels(g)) {
      i <- which(g == lv)
      Hgg <- Xw[i, , drop = FALSE] %*% B %*% t(Xw[i, , drop = FALSE])
      A <- ilm_imh_half(Hgg)
      S[i, ] <- A %*% S[i, , drop = FALSE]
    }
  }

  ## the meat: one summed score vector per cluster, squared up
  Sg <- if (is.null(cluster)) S else
    do.call(rbind, lapply(split(seq_len(N), g), function(i)
      colSums(S[i, , drop = FALSE])))
  meat <- crossprod(Sg)
  adj <- if (type == "CR1") (G / (G - 1)) * ((N - 1) / (N - p)) else 1
  V <- adj * (bread %*% meat %*% bread)
  dimnames(V) <- list(colnames(X), colnames(X))
  structure(V, n_clusters = G, type = type,
            cluster_sizes = as.integer(table(g)))
}

#' The GLM weight each row carries, for the hat matrix
#'
#' @keywords internal
#' @noRd
ilm_estfun_weight <- function(object) {
  fam <- if (!is.null(object$family)) object$family$name else "gaussian"
  eta <- as.numeric(ilm_eta_hat(object, FALSE)[, 1])
  N <- length(eta)
  w <- if (is.null(object$weights)) rep(1, N) else object$weights
  d <- ilm_disp_vec(object)
  v <- if (isTRUE(object$ordinal)) rep(1, N) else switch(fam,
    gaussian = 1 / pmax(d, .Machine$double.eps)^2,
    poisson  = exp(eta),
    binomial = {p <- stats::plogis(eta); p * (1 - p)},
    nbinom   = {mu <- exp(eta); k <- if (is.null(d)) 1 else d
                mu * k / (k + mu)},
    beta     = ilm_beta_w(eta, if (is.null(d)) 1 else d),
    rep(1, N))
  w * v
}

#' (I - H)^(-1/2) for a small symmetric block
#'
#' Through the eigen-decomposition rather than a solve, because a cluster whose
#' leverage reaches one leaves `I - H` singular and the inverse square root
#' undefined. Clamping the eigenvalues away from zero there keeps the
#' adjustment finite instead of returning infinities for a cluster that the
#' design cannot separate from the fit.
#'
#' @keywords internal
#' @noRd
ilm_imh_half <- function(H) {
  M <- diag(nrow(H)) - (H + t(H)) / 2
  e <- eigen(M, symmetric = TRUE)
  lam <- pmax(e$values, 1e-8)
  e$vectors %*% diag(1 / sqrt(lam), nrow(H)) %*% t(e$vectors)
}

#' Fixed effects with sandwich standard errors
#'
#' The same coefficients, re-tested against a covariance that does not assume
#' the variance structure is right. See [ilm_vcov_cluster()] for what the
#' corrections do.
#'
#' @section The degrees of freedom are the point:
#'
#' A sandwich is consistent in the number of clusters, so with few of them the
#' normal reference is badly optimistic. The default here is a t distribution
#' on the Bell-McCaffrey degrees of freedom, computed per coefficient from the
#' design rather than fixed at `G - 1`. Those can be far smaller than `G - 1`
#' when one cluster dominates a predictor -- a treatment assigned to three
#' schools out of thirty does not have twenty-nine degrees of freedom behind
#' it, and the number this reports is the warning.
#'
#' @param object A fitted [ilm_model()] with no random effects.
#' @param cluster Clustering, as in [ilm_vcov_cluster()]. `NULL` gives
#'   heteroskedasticity-robust standard errors.
#' @param type `"CR2"`, `"CR1"` or `"CR0"`.
#' @param df `"bm"` for Bell-McCaffrey (the default), `"G-1"`, `"normal"`, or a
#'   single number to use for every coefficient.
#' @param level Confidence level.
#' @return A data frame of `term`, `estimate`, `se`, `df`, `t`, `p`, `lower`,
#'   `upper`, with the model-based standard error alongside for comparison.
#' @seealso [ilm_vcov_cluster()].
#' @examples
#' set.seed(1); G <- 30; n <- 12
#' d <- data.frame(g = factor(rep(seq_len(G), each = n)), x = rnorm(G * n))
#' d$y <- 1 + 0.4 * d$x + rep(rnorm(G, 0, 1.5), each = n) + rnorm(G * n)
#' f <- ilm_model(y ~ x, data = d, family = "gaussian", verbose = FALSE)
#' ilm_robust(f, ~ g)
#' @export
ilm_robust <- function(object, cluster = NULL, type = c("CR2", "CR1", "CR0"),
                       df = c("bm", "G-1", "normal"), level = 0.95) {
  type <- match.arg(type)
  if (is.character(df)) df <- match.arg(df)
  V <- ilm_vcov_cluster(object, cluster, type)
  G <- attr(V, "n_clusters")
  b <- stats::coef(object)
  se <- sqrt(pmax(diag(V), 0))
  dfv <- if (is.numeric(df)) rep(df[1], length(b))
         else switch(df,
           normal = rep(Inf, length(b)),
           `G-1`  = rep(G - 1, length(b)),
           bm     = ilm_bm_df(object, cluster, G))
  crit <- ifelse(is.finite(dfv), stats::qt(1 - (1 - level) / 2, dfv),
                 stats::qnorm(1 - (1 - level) / 2))
  tv <- b / se
  pv <- ifelse(is.finite(dfv), 2 * stats::pt(-abs(tv), dfv),
               2 * stats::pnorm(-abs(tv)))
  out <- data.frame(term = names(b), estimate = unname(b),
                    se = unname(se),
                    se_model = unname(sqrt(diag(as.matrix(stats::vcov(object))))),
                    df = unname(dfv), t = unname(tv), p = unname(pv),
                    lower = unname(b - crit * se), upper = unname(b + crit * se),
                    row.names = NULL, stringsAsFactors = FALSE)
  structure(out, class = c("ilm_robust", "data.frame"),
            n_clusters = G, type = type, df_kind = df,
            clustered = !is.null(cluster))
}

#' Bell-McCaffrey degrees of freedom, one per coefficient
#'
#' The Satterthwaite degrees of freedom of the cluster-robust variance of a
#' single coefficient, treating the errors as if the null model held. It falls
#' out of the eigenvalues of the matrix whose quadratic form the variance is.
#'
#' @keywords internal
#' @noRd
ilm_bm_df <- function(object, cluster, G) {
  X <- object$X; N <- nrow(X); p <- ncol(X)
  g <- if (is.null(cluster)) factor(seq_len(N)) else
    ilm_cluster_id(object, cluster)
  W <- ilm_estfun_weight(object)
  Xw <- X * sqrt(W)
  XtXinv <- tryCatch(solve(crossprod(Xw)), error = function(e) NULL)
  if (is.null(XtXinv)) return(rep(G - 1, p))
  idx <- split(seq_len(N), g)
  ## A_g = (I - H_gg)^(-1/2) applied to cluster g's rows of the weighted design
  Ag <- lapply(idx, function(i) {
    Xi <- Xw[i, , drop = FALSE]
    ilm_imh_half(Xi %*% XtXinv %*% t(Xi)) %*% Xi
  })
  vapply(seq_len(p), function(j) {
    l <- XtXinv[, j]
    ## the variance of l'beta is a quadratic form in the errors; its
    ## Satterthwaite df is (sum lambda)^2 / sum lambda^2 for the eigenvalues
    ## of that form, and stacking the per-cluster pieces gives them
    Gm <- do.call(rbind, lapply(seq_along(idx), function(k) {
      i <- idx[[k]]
      v <- rep(0, N); v[i] <- Ag[[k]] %*% l
      v
    }))
    M <- tcrossprod(Gm)
    ev <- eigen(M, symmetric = TRUE, only.values = TRUE)$values
    ev <- ev[ev > max(ev) * 1e-10]
    if (!length(ev)) return(G - 1)
    max(sum(ev)^2 / sum(ev^2), 1)
  }, 0)
}

#' @export
print.ilm_robust <- function(x, ...) {
  G <- attr(x, "n_clusters")
  cat(sprintf("Fixed effects with %s standard errors (%s)\n",
              if (attr(x, "clustered")) "cluster-robust" else "heteroskedasticity-robust",
              attr(x, "type")))
  if (attr(x, "clustered")) cat(sprintf("  %d clusters\n", G))
  d <- as.data.frame(x)
  d$estimate <- signif(d$estimate, 5); d$se <- signif(d$se, 4)
  d$se_model <- signif(d$se_model, 4); d$df <- round(d$df, 1)
  d$t <- round(d$t, 3); d$p <- signif(d$p, 3)
  d$lower <- signif(d$lower, 4); d$upper <- signif(d$upper, 4)
  print(d, row.names = FALSE)
  ratio <- x$se / x$se_model
  cat(sprintf("\n  se / model-based se: %.2f to %.2f\n",
              min(ratio), max(ratio)))
  if (attr(x, "clustered") && G < 15L)
    cat("\n  ", G, " clusters is very few. No correction reliably rescues an\n",
        "  interval down here -- the estimator is learning a covariance from\n",
        "  ", G, " numbers. Treat these as optimistic.\n", sep = "")
  else if (attr(x, "clustered") && G < 40L)
    cat("\n  Under 40 clusters, the correction and the reference distribution\n",
        "  both matter. CR2 with Bell-McCaffrey degrees of freedom is the\n",
        "  least bad combination and is what was used here.\n", sep = "")
  invisible(x)
}
