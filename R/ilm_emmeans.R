## ---------------------------------------------------------------------------
## Estimated marginal means, and contrasts among them.
##
## Almost every fitted model in applied work is followed by "...and now compare
## the groups". ilm_anova() says whether a term matters and ilm_ame() says what
## a unit of a predictor is worth; neither answers "what is the predicted mean
## for group A, and is it different from group B".
##
## A marginal mean is a prediction on a REFERENCE GRID: every combination of the
## factors asked for, with the other factors averaged over and the numeric
## predictors held at their means. What "averaged over" means is a choice, and
## it is the choice people get wrong -- equal weights answer "what would the
## groups look like in a balanced design", proportional weights answer "what do
## they look like in this population". Both are defensible and they differ
## whenever the design is unbalanced, so the argument is explicit and the
## printed output says which was used.
##
## Every quantity here is a linear combination of the coefficients, L %*% beta,
## so its variance is L V L' exactly -- no simulation, no delta method. That
## holds on the LINK scale, which is why the means are reported there by
## default and back-transformation is asked for rather than assumed.
## ---------------------------------------------------------------------------

## The reference grid: one row per combination of every factor in the model,
## numeric predictors at their means unless `at` says otherwise.
#' @keywords internal
#' @noRd
ilm_ref_grid <- function(object, at = NULL) {
  mf <- object$model
  if (is.null(mf)) stop("the fit did not keep its model frame", call. = FALSE)
  tl <- attr(stats::terms(object), "term.labels")
  vars <- intersect(all.vars(stats::delete.response(stats::terms(object))),
                    names(mf))
  vals <- list()
  for (v in vars) {
    z <- mf[[v]]
    vals[[v]] <- if (!is.null(at) && v %in% names(at)) at[[v]]
                 else if (is.numeric(z)) mean(z, na.rm = TRUE)
                 else factor(levels(factor(z)), levels = levels(factor(z)))
  }
  g <- expand.grid(vals, stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)
  for (v in names(g)) if (is.factor(mf[[v]]) || is.character(mf[[v]]))
    g[[v]] <- factor(g[[v]], levels = levels(factor(mf[[v]])))
  g
}

#' Cell weights for averaging over the grid
#'
#' The difference between the two frequency schemes is the whole reason for
#' offering both: a product of MARGINAL counts gives every level of `specs` the
#' same mix of everything else, so a contrast is the effect alone, while the
#' JOINT count gives each level its own mix, so a contrast also carries
#' whatever the composition does.
#'
#' @param g A reference grid.
#' @param mf The fit's model frame, which supplies the observed frequencies.
#' @param weights One of `"equal"`, `"cells"`, `"proportional"`.
#' @return A numeric weight per grid row.
#' @keywords internal
#' @noRd
ilm_emm_cellw <- function(g, mf, weights) {
  fv <- names(g)[vapply(g, is.factor, TRUE)]
  w <- if (weights == "equal" || !length(fv)) rep(1, nrow(g)) else
    if (weights == "cells") {
      tb <- table(interaction(mf[fv], drop = FALSE, sep = "\r"))
      key <- as.character(interaction(g[fv], drop = FALSE, sep = "\r"))
      as.numeric(tb[match(key, names(tb))])
    } else {
      ## outer product of the one-way margins
      Reduce(`*`, lapply(fv, function(v) {
        tb <- table(mf[[v]])
        as.numeric(tb[match(as.character(g[[v]]), names(tb))])
      }))
    }
  w[!is.finite(w)] <- 0
  w
}

#' Average model-matrix rows within each level of `specs`
#'
#' Shared by [ilm_emmeans()] and [ilm_trends()]. The only thing that differs
#' between them is what the matrix holds: the grid's own model matrix for a
#' mean, its derivative with respect to a covariate for a slope. The averaging
#' is identical, and so is everything downstream of it.
#'
#' @param mm A matrix with one row per grid row.
#' @param g The reference grid.
#' @param specs Variables whose levels the rows are grouped by.
#' @param w Cell weights from [ilm_emm_cellw()].
#' @return A list with `L` and the level labels `lv`.
#' @keywords internal
#' @noRd
ilm_emm_avg <- function(mm, g, specs, w) {
  keyspec <- if (length(specs)) interaction(g[specs], drop = TRUE, sep = "\r")
             else factor(rep("", nrow(g)))
  lv <- levels(keyspec)
  L <- matrix(0, length(lv), ncol(mm), dimnames = list(NULL, colnames(mm)))
  for (i in seq_along(lv)) {
    k <- which(keyspec == lv[i])
    ww <- w[k]
    if (sum(ww) <= 0) ww <- rep(1, length(k))
    L[i, ] <- colSums(mm[k, , drop = FALSE] * (ww / sum(ww)))
  }
  list(L = L, lv = lv)
}

#' Estimated marginal means
#'
#' The model's predicted mean for each level of the variables asked for, with
#' the other predictors averaged over rather than left at whatever value
#' happens to be first.
#'
#' @section What is averaged, and how:
#'
#' A marginal mean is a prediction on a grid of every factor combination, with
#' numeric predictors held at their means. The factors *not* named in `specs`
#' have to be averaged over, and how they are averaged changes the answer
#' whenever the design is unbalanced:
#'
#' * `weights = "equal"` treats every cell alike, answering "what would these
#'   groups look like in a balanced design". This is the usual default and
#'   matches what most software calls an estimated marginal mean.
#' * `weights = "proportional"` weights each cell by the product of the
#'   *marginal* frequencies of the variables being averaged over, answering
#'   "what do these groups look like in a population with this sample's
#'   composition".
#' * `weights = "cells"` weights by the *joint* frequency actually observed,
#'   so each group is averaged over its own covariate distribution.
#'
#' The first two keep a comparison clean: the averaged-over part is identical
#' for every level of `specs`, so a difference between two marginal means is
#' the model's effect and nothing else. `"cells"` does not, because each group
#' gets its own mix. When the averaged-over variables are associated with
#' `specs`, a difference of cell-weighted means carries that difference in
#' composition as well as the effect, and [ilm_contrast()] says so. That is
#' sometimes the quantity you want -- it describes the groups as they are --
#' but it is not an adjusted comparison.
#'
#' @section An ordered response:
#'
#' For an ordinal fit the link-scale values are marginal means of the **latent
#' scale** -- the linear predictor the thresholds cut up -- and not of the
#' categories, which have no mean to take. A contrast between two of them is a
#' difference in log odds of being in a higher category, constant across cuts
#' by the same assumption [ilm_check_proportional()] tests.
#'
#' `type = "response"` gives each category's PROBABILITY instead, computed in
#' every cell of the grid and averaged over the cells, with a delta-method
#' standard error that carries the thresholds' uncertainty as well as the
#' slopes'. [ilm_contrast()] then compares groups within each category as
#' differences in probability. This agrees with `emmeans` on a `MASS::polr()`
#' fit with `mode = "prob"`.
#'
#' @section A multinomial response:
#'
#' Every category gets its own row for every level of `specs`.
#'
#' On the **link** scale the value is the category's centred log-odds: its
#' log-probability less the average log-probability over all the categories,
#' which is what the sum-to-zero coefficients describe. These are exact linear
#' combinations, as for any other family, and a contrast between two groups
#' within a category is a difference of log-odds against the same average.
#'
#' On the **response** scale the value is the category's PROBABILITY, computed
#' in every cell of the grid and then averaged over the cells with the chosen
#' weights -- the probabilities of each group sum to 1 -- with a delta-method
#' standard error and an interval formed on the logit scale, so it stays
#' inside 0 and 1. [ilm_contrast()] then compares groups within each category
#' as differences in probability. This is what `emmeans` computes for an
#' `nnet::multinom()` fit with `mode = "prob"`, and the two agree; the
#' coefficients differ between the packages, because `nnet` codes against a
#' baseline category, but the probabilities do not.
#'
#' @section Which scale:
#'
#' Marginal means are computed on the **link** scale, where they are exact
#' linear combinations of the coefficients and their variance is `L V L'` with
#' no approximation. `type = "response"` back-transforms the endpoints, which
#' keeps the interval's coverage but means the reported centre is a median
#' rather than a mean on that scale. For a model with a random effect that
#' centre is also conditional on the group rather than population-averaged --
#' [ilm_ame()] is the population-averaged quantity.
#'
#' @param object A fitted [ilm_model()].
#' @param specs Variables to keep, as a character vector. Everything else is
#'   averaged over or held at its mean.
#' @param at Named list fixing particular values, for instance
#'   `list(age = c(40, 60))`.
#' @param weights `"equal"`, `"proportional"` or `"cells"`; see above.
#' @param type `"link"` or `"response"`.
#' @param level Confidence level.
#' @return An object of class `"ilm_emm"`: a data frame of the grid with
#'   `estimate`, `se`, `lower`, `upper`, plus the contrast machinery it carries.
#' @seealso [ilm_contrast()] to compare them, [ilm_ame()] for the average
#'   marginal effect of a predictor.
#' @examples
#' set.seed(1); n <- 200
#' d <- data.frame(g = factor(sample(c("a", "b", "c"), n, TRUE)), x = rnorm(n))
#' d$y <- 1 + 0.5 * (d$g == "b") + 0.2 * d$x + rnorm(n)
#' fit <- ilm_model(y ~ g + x, data = d, family = "gaussian", verbose = FALSE)
#' ilm_emmeans(fit, "g")
#' @export
ilm_emmeans <- function(object, specs, at = NULL,
                        weights = c("equal", "proportional", "cells"),
                        type = c("link", "response"), level = 0.95) {
  weights <- match.arg(weights); type <- match.arg(type)
  if (!inherits(object, "ilm_model"))
    stop("`object` must be a fitted ilm_model, not ", class(object)[1],
         call. = FALSE)
  mn <- is.null(object$family) || identical(object$family$name, "multinomial")
  mf <- object$model
  specs <- as.character(specs)
  miss <- setdiff(specs, names(mf))
  if (length(miss))
    stop("variable(s) not in the model: ", paste(miss, collapse = ", "),
         ". The model has: ",
         paste(setdiff(names(mf), names(mf)[1]), collapse = ", "),
         call. = FALSE)
  ## The response is a column of the model frame too, so it passed the check
  ## above and then failed inside the grid on "undefined columns selected".
  if (names(mf)[1L] %in% specs)
    stop("`", names(mf)[1L], "` is the response. Marginal means are taken ",
         "over the levels of predictors: ",
         paste(setdiff(names(mf), names(mf)[1L]), collapse = ", "),
         if (identical(object$family$name, "multinomial") || isTRUE(object$ordinal))
           ". Every category of the response gets its own row whatever is asked for."
         else ".", call. = FALSE)

  g <- ilm_ref_grid(object, at)
  ## the model matrix of the grid, built with the FIT's terms and contrasts so
  ## the columns line up with the coefficients
  tt <- stats::delete.response(stats::terms(object))
  mmg <- ilm_drop_intercept(
    stats::model.matrix(tt, data = g, contrasts.arg = object$contrasts), object)
  b <- stats::coef(object)
  ## a multinomial fit has one coefficient per column PER CATEGORY
  nC <- if (mn) object$C else 1L
  if (ncol(mmg) * nC != length(b))
    stop("the reference grid does not match the fitted coefficients; a term ",
         "here is not a plain variable (a smooth or a matrix column), and a ",
         "marginal mean is not defined for it", call. = FALSE)
  V <- suppressWarnings(as.matrix(stats::vcov(object)))

  w <- ilm_emm_cellw(g, mf, weights)
  av <- ilm_emm_avg(mmg, g, specs, w)
  if (mn)
    return(ilm_emm_multinom(object, g, mmg, av, specs, w, b, V, type, level,
                            weights))
  if (isTRUE(object$ordinal) && identical(type, "response"))
    return(ilm_emm_ordinal(object, g, mmg, av, specs, w, level, weights))
  L <- av$L; lv <- av$lv
  est <- as.numeric(L %*% b)
  Vem <- L %*% V %*% t(L)
  se <- sqrt(pmax(diag(Vem), 0))
  crit <- if (isTRUE(object$exact_df))
            stats::qt(1 - (1 - level) / 2, object$resid_df)
          else stats::qnorm(1 - (1 - level) / 2)

  out <- if (length(specs))
    as.data.frame(do.call(rbind, strsplit(lv, "\r", fixed = TRUE)),
                  stringsAsFactors = FALSE)
  else data.frame(.all = "", stringsAsFactors = FALSE)
  names(out) <- if (length(specs)) specs else ".all"
  out$estimate <- est; out$se <- se
  out$lower <- est - crit * se; out$upper <- est + crit * se

  fam <- if (!is.null(object$family)) object$family$name else "gaussian"
  if (type == "response") {
    inv <- ilm_emm_linkinv(object)
    out$estimate <- inv(est); out$lower <- inv(out$lower)
    out$upper <- inv(out$upper); out$se <- NA_real_
  }
  rownames(out) <- NULL
  structure(out, class = c("ilm_emm", "data.frame"), L = L, V = Vem,
            specs = specs, weights = weights, type = type, level = level,
            family = fam, object = object)
}

## The multinomial case: a row for every category in every level of `specs`.
##
## On the link scale each row is a category's centred log-odds, T_c[j, ] B'
## L_i -- exact and linear, so its gradient in vec(B) is kron(T_c[j, ], L_i).
## On the response scale it is the category's probability in each grid cell,
## averaged over the cells with the chosen weights, and the gradient is the
## same average of each cell's softmax Jacobian:
##
##   dp/dvec(B) = ((diag(p) - p p') T_c) (x) x_cell'
##
## Contrasts are formed on whichever scale was asked for, from these gradients.
#' @keywords internal
#' @noRd
ilm_emm_multinom <- function(object, g, mmg, av, specs, w, b, V, type, level,
                             weights) {
  J <- object$J; C <- object$C; p <- ncol(mmg)
  Tc <- stats::contr.sum(J); cats <- object$ylevels
  B <- matrix(b, p, C)
  z <- stats::qnorm(1 - (1 - level) / 2)
  lv <- av$lv; nl <- length(lv)
  keyspec <- if (length(specs)) interaction(g[specs], drop = TRUE, sep = "\r")
             else factor(rep("", nrow(g)))
  est <- numeric(nl * J); G <- matrix(0, nl * J, p * C)
  for (i in seq_len(nl)) {
    rows <- (i - 1L) * J + seq_len(J)
    if (type == "link") {
      Li <- av$L[i, , drop = FALSE]
      est[rows] <- as.numeric(Li %*% B %*% t(Tc))
      for (j in seq_len(J)) G[rows[j], ] <- kronecker(Tc[j, , drop = FALSE], Li)
    } else {
      k <- which(keyspec == lv[i]); ww <- w[k]
      if (sum(ww) <= 0) ww <- rep(1, length(k))
      ww <- ww / sum(ww)
      for (m in seq_along(k)) {
        x <- mmg[k[m], , drop = FALSE]
        pr <- exp(as.numeric(x %*% B %*% t(Tc))); pr <- pr / sum(pr)
        est[rows] <- est[rows] + ww[m] * pr
        G[rows, ] <- G[rows, ] +
          ww[m] * kronecker((diag(pr) - pr %o% pr) %*% Tc, x)
      }
    }
  }
  Vem <- G %*% V %*% t(G)
  se <- sqrt(pmax(diag(Vem), 0))
  if (type == "link") {
    lo <- est - z * se; hi <- est + z * se
  } else {
    ## on the logit scale and back, so the interval cannot leave [0, 1]
    pe <- pmin(pmax(est, 1e-12), 1 - 1e-12)
    sl <- se / (pe * (1 - pe))
    lo <- stats::plogis(stats::qlogis(pe) - z * sl)
    hi <- stats::plogis(stats::qlogis(pe) + z * sl)
  }
  sp <- if (length(specs))
    as.data.frame(do.call(rbind, strsplit(lv, "\r", fixed = TRUE)),
                  stringsAsFactors = FALSE)
  else data.frame(.all = "", stringsAsFactors = FALSE)
  names(sp) <- if (length(specs)) specs else ".all"
  out <- sp[rep(seq_len(nl), each = J), , drop = FALSE]
  out$category <- factor(rep(cats, nl), levels = cats)
  out$estimate <- est; out$se <- se; out$lower <- lo; out$upper <- hi
  rownames(out) <- NULL
  structure(out, class = c("ilm_emm", "data.frame"), L = G, V = Vem,
            est_c = est, category = as.character(out$category),
            specs = specs, weights = weights, type = type, level = level,
            family = "multinomial", object = object)
}

## An ordinal fit on the response scale: each category's probability,
##
##   P(Y = j) = F(theta_j - eta) - F(theta_{j-1} - eta),
##
## in every grid cell, averaged over the cells as for a multinomial. The
## thresholds are estimated too, so the gradient runs over them as well as the
## slopes -- through the first-value-and-log-gaps form they are fitted in --
## and the variance comes from the fit's full covariance.
#' @keywords internal
#' @noRd
ilm_emm_ordinal <- function(object, g, mmg, av, specs, w, level, weights) {
  fam <- object$family; zeta <- as.numeric(object$zeta)
  K <- length(zeta); J <- K + 1L; cats <- object$ylevels
  beta <- stats::coef(object); p <- length(beta)
  tl <- names(object$opt$par)
  ib <- seq_len(p); iz <- which(tl == "zeta_raw")
  zr <- unname(object$opt$par[iz])
  Vf <- suppressWarnings(as.matrix(stats::vcov(object, full = TRUE)))
  ## d theta_k / d zeta_raw: theta_k = r_1 + sum_{m = 2..k} exp(r_m)
  Dz <- matrix(0, K, K); Dz[, 1L] <- 1
  if (K > 1L) for (k in 2:K) Dz[k, 2:k] <- exp(zr[2:k])
  dens <- switch(fam$link, logit = stats::dlogis, probit = stats::dnorm,
                 cloglog = function(z) exp(z - exp(z)))
  z <- stats::qnorm(1 - (1 - level) / 2)
  lv <- av$lv; nl <- length(lv)
  keyspec <- if (length(specs)) interaction(g[specs], drop = TRUE, sep = "\r")
             else factor(rep("", nrow(g)))
  est <- numeric(nl * J); G <- matrix(0, nl * J, ncol(Vf))
  for (i in seq_len(nl)) {
    rows <- (i - 1L) * J + seq_len(J)
    k <- which(keyspec == lv[i]); ww <- w[k]
    if (sum(ww) <= 0) ww <- rep(1, length(k))
    ww <- ww / sum(ww)
    for (m in seq_along(k)) {
      x <- mmg[k[m], ]
      eta <- sum(x * beta)
      pr <- diff(c(0, fam$pfun(zeta - eta), 1))
      fz <- dens(zeta - eta)
      ## dP_j/dbeta = -(f_j - f_{j-1}) x, with f_0 = f_J = 0
      fe <- c(0, fz, 0)
      gb <- -(fe[-1L] - fe[-(J + 1L)]) %o% x
      ## dP_j/dtheta_k = f_k [j = k] - f_k [j = k + 1]
      At <- matrix(0, J, K)
      At[cbind(seq_len(K), seq_len(K))] <- fz
      At[cbind(seq_len(K) + 1L, seq_len(K))] <- -fz
      est[rows] <- est[rows] + ww[m] * pr
      G[rows, ib] <- G[rows, ib] + ww[m] * gb
      G[rows, iz] <- G[rows, iz] + ww[m] * (At %*% Dz)
    }
  }
  Vem <- G %*% Vf %*% t(G)
  se <- sqrt(pmax(diag(Vem), 0))
  pe <- pmin(pmax(est, 1e-12), 1 - 1e-12)
  sl <- se / (pe * (1 - pe))
  sp <- if (length(specs))
    as.data.frame(do.call(rbind, strsplit(lv, "\r", fixed = TRUE)),
                  stringsAsFactors = FALSE)
  else data.frame(.all = "", stringsAsFactors = FALSE)
  names(sp) <- if (length(specs)) specs else ".all"
  out <- sp[rep(seq_len(nl), each = J), , drop = FALSE]
  out$category <- factor(rep(cats, nl), levels = cats)
  out$estimate <- est; out$se <- se
  out$lower <- stats::plogis(stats::qlogis(pe) - z * sl)
  out$upper <- stats::plogis(stats::qlogis(pe) + z * sl)
  rownames(out) <- NULL
  structure(out, class = c("ilm_emm", "data.frame"), L = G, V = Vem,
            est_c = est, category = as.character(out$category),
            specs = specs, weights = weights, type = "response",
            level = level, family = fam$name, object = object)
}

#' @keywords internal
#' @noRd
## The inverse link, from the family object, which knows its own. A switch on
## the family's NAME used to stand here, and fell through to the identity for
## any family it did not list: a beta model's "response" means were its
## link-scale values, 0.37 where the proportion was 0.59.
ilm_emm_linkinv <- function(object) {
  f <- object$family
  if (is.null(f) || !is.function(f$linkinv)) return(function(z) z)
  f$linkinv
}

#' @export
print.ilm_emm <- function(x, ...) {
  sp <- attr(x, "specs")
  cat("<ilm_emm> marginal means of", paste(sp, collapse = " x "),
      sprintf("(%s scale, %s weights)\n\n", attr(x, "type"),
              attr(x, "weights")))
  print(as.data.frame(x), row.names = FALSE, digits = 4)
  ord <- startsWith(attr(x, "family"), "ordinal")
  catg <- identical(attr(x, "family"), "multinomial") || ord
  if (catg && attr(x, "type") == "response") {
    cat("\n  Category probabilities, averaged over the grid; each group's sum\n",
        "  to 1. Intervals are formed on the logit scale.",
        if (length(attr(x, "object")$re))
          "\n  They are for a typical group, with the random effects at zero;\n  ilm_ame(groups = \"population\") gives the comparison averaged over\n  the groups."
        else "", "\n", sep = "")
  } else if (identical(attr(x, "family"), "multinomial")) {
    cat("\n  Each category's CENTRED log-odds: its log-probability less the\n",
        "  average over all the categories. A single value is not a\n",
        "  probability; differences between groups within a category are\n",
        "  comparable. type = \"response\" gives the probabilities.\n", sep = "")
  } else if (ord)
    cat("\n  On the LATENT scale, whose origin the thresholds set, so a single\n",
        "  mean is not interpretable on its own -- the differences are.\n",
        sep = "")
  else if (attr(x, "type") == "link" && attr(x, "family") != "gaussian")
    cat("\n  On the link scale. type = \"response\" back-transforms.\n")
  cat("  Compare them with ilm_contrast().\n")
  invisible(x)
}

## ---- comparing them --------------------------------------------------------

## Build the contrast matrix for a named scheme.
#' @keywords internal
#' @noRd
ilm_contrast_matrix <- function(lab, method, ref = NULL) {
  k <- length(lab)
  if (k < 2L)
    stop("at least 2 marginal means are needed to contrast them", call. = FALSE)
  if (method == "pairwise") {
    ij <- utils::combn(k, 2L)
    C <- matrix(0, ncol(ij), k)
    for (i in seq_len(ncol(ij))) { C[i, ij[1, i]] <- -1; C[i, ij[2, i]] <- 1 }
    rownames(C) <- paste(lab[ij[2, ]], "-", lab[ij[1, ]])
    return(C)
  }
  if (method == "trt.vs.ctrl") {
    r <- if (is.null(ref)) 1L else match(as.character(ref), lab)
    if (is.na(r))
      stop("`ref` (", ref, ") is not one of: ", paste(lab, collapse = ", "),
           call. = FALSE)
    o <- setdiff(seq_len(k), r)
    C <- matrix(0, length(o), k)
    for (i in seq_along(o)) { C[i, r] <- -1; C[i, o[i]] <- 1 }
    rownames(C) <- paste(lab[o], "-", lab[r])
    return(C)
  }
  ## orthogonal polynomial trends, which assume the levels are ordered AND
  ## equally spaced -- true for a dose series, false for arbitrary groups
  P <- stats::contr.poly(k)
  C <- t(P)
  rownames(C) <- colnames(P)
  C
}

#' Compare estimated marginal means
#'
#' Differences between the means from [ilm_emmeans()], with intervals and a
#' multiplicity adjustment.
#'
#' @section The adjustment:
#'
#' Comparing every pair of five groups is ten tests, and ten intervals each
#' nominally 95% do not jointly cover at 95%. The default is the single-step
#' studentized maximum: the contrasts' joint covariance is known exactly here,
#' `C L V L' C'`, so the reference distribution is simulated from it directly
#' rather than bootstrapped. This is the same construction
#' [illumex::ilm_boot_diff()] uses on raw data, where it reproduced `TukeyHSD()` to
#' 0.006 on the design Tukey is exact for.
#'
#' `"bonferroni"` is the conservative fallback and `"none"` is there for
#' comparisons chosen in advance. With a single contrast there is nothing to
#' adjust and the column reads `"none"` whatever was asked.
#'
#' @param object An [ilm_emmeans()] result.
#' @param method `"pairwise"`, `"trt.vs.ctrl"` or `"poly"`. The last assumes
#'   the levels are ordered and equally spaced.
#' @param ref Reference level for `"trt.vs.ctrl"`.
#' @param adjust `"max_t"`, `"bonferroni"` or `"none"`.
#' @param level Confidence level for the family.
#' @param nsim Draws used to find the studentized-maximum critical value.
#' @param seed Random seed for that simulation.
#' @return A data frame with `contrast`, `estimate`, `se`, `lower`, `upper`,
#'   `p_value`, `p_adj` and `adjust`.
#' @seealso [ilm_emmeans()], [illumex::ilm_boot_diff()] for the same comparison made
#'   without a model.
#' @examples
#' set.seed(1); n <- 200
#' d <- data.frame(g = factor(sample(c("a", "b", "c"), n, TRUE)), x = rnorm(n))
#' d$y <- 1 + 0.5 * (d$g == "b") + 0.2 * d$x + rnorm(n)
#' fit <- ilm_model(y ~ g + x, data = d, family = "gaussian", verbose = FALSE)
#' ilm_contrast(ilm_emmeans(fit, "g"))
#' @export
ilm_contrast <- function(object, method = c("pairwise", "trt.vs.ctrl", "poly"),
                         ref = NULL,
                         adjust = c("max_t", "bonferroni", "none"),
                         level = 0.95, nsim = 20000L, seed = 1L) {
  method <- match.arg(method); adjust <- match.arg(adjust)
  if (!inherits(object, "ilm_emm"))
    stop("`object` must be an ilm_emmeans() result, not ", class(object)[1],
         call. = FALSE)
  sp <- attr(object, "specs")
  ## cell weights average each group over its own covariate mix, so a
  ## difference between two of them is not an adjusted comparison
  if (identical(attr(object, "weights"), "cells"))
    warning("these marginal means were averaged with weights = \"cells\", so ",
            "each group carries its own mix of the other variables. A ",
            "difference between them is the effect plus whatever that ",
            "difference in composition contributes, not an adjusted ",
            "comparison. Use weights = \"equal\" or \"proportional\" for one.",
            call. = FALSE)
  lab <- if (length(sp))
    apply(as.data.frame(object)[sp], 1L, paste, collapse = " ") else
      rep("overall", nrow(object))
  grp <- attr(object, "category")
  C <- if (is.null(grp)) ilm_contrast_matrix(lab, method, ref) else {
    ## a multinomial: the same comparisons among the groups, within each
    ## category, and never across categories
    do.call(rbind, lapply(unique(grp), function(k) {
      r <- which(grp == k)
      Ck <- ilm_contrast_matrix(lab[r], method, ref)
      M <- matrix(0, nrow(Ck), length(grp))
      M[, r] <- Ck
      rownames(M) <- paste0(k, ": ", rownames(Ck))
      M
    }))
  }
  Vem <- attr(object, "V")
  ## The marginal means on the LINK scale, whatever the print showed -- except
  ## for a multinomial, whose means carry the scale they were asked for: a
  ## difference in probability is what a comparison of categories is usually
  ## wanted as, and its delta-method variance is exact to first order.
  est0 <- attr(object, "est_c")
  if (is.null(est0))
    est0 <- as.numeric(attr(object, "L") %*% stats::coef(attr(object, "object")))
  est <- as.numeric(C %*% est0)
  Vc <- C %*% Vem %*% t(C)
  se <- sqrt(pmax(diag(Vc), 0))
  m <- length(est)
  fit <- attr(object, "object")
  df <- if (isTRUE(fit$exact_df)) fit$resid_df else Inf
  if (m == 1L) adjust <- "none"

  tstat <- ifelse(se > 0, est / se, 0)
  praw <- 2 * stats::pt(-abs(tstat), df)
  crit <- stats::qt(1 - (1 - level) / 2, df)
  padj <- praw

  if (adjust == "max_t") {
    ## The contrasts' joint distribution is known here, so the critical value
    ## comes from simulating that distribution rather than resampling data.
    R <- stats::cov2cor(Vc + diag(1e-12, m))
    ev <- eigen((R + t(R)) / 2, symmetric = TRUE)
    A <- ev$vectors %*% diag(sqrt(pmax(ev$values, 0)), m, m)
    set.seed(seed)
    Zs <- A %*% matrix(stats::rnorm(m * nsim), m, nsim)
    ## multivariate t, not normal: the contrasts are divided by an ESTIMATED
    ## standard error, so the reference carries that estimate's uncertainty
    ## too. The two agree once the degrees of freedom are large, and diverge
    ## exactly where it matters -- a small sample, where a normal reference
    ## gives a critical value that is too small and intervals that are too
    ## narrow.
    if (is.finite(df) && df > 0)
      Zs <- Zs / rep(sqrt(stats::rchisq(nsim, df) / df), each = m)
    mx <- apply(abs(Zs), 2L, max)
    crit <- stats::quantile(mx, level, names = FALSE)
    padj <- vapply(abs(tstat), function(t0) mean(mx >= t0), 1)
  } else if (adjust == "bonferroni") {
    crit <- stats::qt(1 - (1 - level) / (2 * m), df)
    padj <- pmin(1, m * praw)
  }

  out <- data.frame(contrast = rownames(C), estimate = est, se = se,
                    lower = est - crit * se, upper = est + crit * se,
                    p_value = praw, p_adj = padj, adjust = adjust,
                    n_contrasts = m, stringsAsFactors = FALSE)
  rownames(out) <- NULL
  structure(out, class = c("ilm_contrast", "data.frame"),
            level = level, type = attr(object, "type"),
            family = attr(object, "family"), weights = attr(object, "weights"))
}

#' @export
print.ilm_contrast <- function(x, ...) {
  cat("<ilm_contrast>", nrow(x), "comparison(s), adjust =", x$adjust[1], "\n\n")
  d <- as.data.frame(x)[, c("contrast", "estimate", "se", "lower", "upper",
                            "p_adj")]
  print(d, row.names = FALSE, digits = 4)
  if (x$adjust[1] == "max_t")
    cat("\n  Intervals hold jointly at",
        sprintf("%.0f%%", 100 * attr(x, "level")),
        "across all", nrow(x), "comparisons.\n")
  invisible(x)
}
