## illume::ilm_model() -- multinomial linear mixed model, matrix-normal core.
## Prototype.  Arbitrary J; crossed random intercepts AND random slopes, each
## with a per-term category covariance structure ("us", "diag", "rr");
## optional AR(1).

#' Number of free parameters in a covariance structure
#'
#' An unstructured (symmetric) C-by-C covariance matrix has C(C+1)/2 distinct
#' entries, because the upper and lower triangles mirror each other. `ilm_nrr()`
#' gives the count for a reduced-rank structure instead.
#'
#' @param C Integer. Number of category dimensions (`J - 1` for `J` categories).
#' @param r Integer. Rank of the reduced-rank structure, smaller than C.
#' @return Integer parameter count.
#' @keywords internal
#' @noRd
ilm_ncov <- function(C) as.integer(C * (C + 1L) / 2L)              # unstructured

#' Number of free parameters in a covariance structure
#'
#' An unstructured (symmetric) C-by-C covariance matrix has C(C+1)/2 distinct
#' entries, because the upper and lower triangles mirror each other. `ilm_nrr()`
#' gives the count for a reduced-rank structure instead.
#'
#' @param C Integer. Number of category dimensions (`J - 1` for `J` categories).
#' @param r Integer. Rank of the reduced-rank structure, smaller than C.
#' @return Integer parameter count.
#' @keywords internal
#' @noRd
ilm_nrr  <- function(C, r) as.integer(C * r - r * (r - 1L) / 2L)   # reduced rank

#' Build a covariance factor from unconstrained parameters
#'
#' Optimisers work best when parameters can take any real value, but a
#' covariance matrix must be positive definite. The standard solution is to
#' estimate the **Cholesky factor** `L` (a lower-triangular matrix) and form
#' `Sigma = L L'`, which is positive definite for any `L` with a positive
#' diagonal. Taking `exp()` of the diagonal guarantees that positivity, so the
#' optimiser sees an unconstrained problem.
#'
#' `ilm_mkL()` builds a full unstructured factor, `ilm_mkD()` a diagonal one (no
#' correlation between categories), and `ilm_mkLam()` a C-by-r loadings matrix for a
#' reduced-rank structure in which `Sigma = Lambda Lambda'`.
#'
#' @param v Numeric vector of unconstrained parameters.
#' @param C Integer. Number of category dimensions.
#' @param r Integer. Rank; `ilm_mkLam()` only.
#' @return A matrix, automatic-differentiation aware for use inside the
#'   likelihood.
#' @references
#' Pinheiro, J. C., & Bates, D. M. (1996). Unconstrained parametrizations for
#' variance-covariance matrices. *Statistics and Computing*, 6, 289--296.
#' @keywords internal
#' @noRd
ilm_mkL <- function(v, C) {                       # unconstrained -> lower Cholesky
  Lv <- AD(numeric(C * C)); k <- 1L
  for (j in 1:C) for (i in j:C) {
    Lv[(j - 1L) * C + i] <- if (i == j) exp(v[k]) else v[k]; k <- k + 1L
  }
  matrix(Lv, C, C)
}

#' Build a covariance factor from unconstrained parameters
#'
#' Optimisers work best when parameters can take any real value, but a
#' covariance matrix must be positive definite. The standard solution is to
#' estimate the **Cholesky factor** `L` (a lower-triangular matrix) and form
#' `Sigma = L L'`, which is positive definite for any `L` with a positive
#' diagonal. Taking `exp()` of the diagonal guarantees that positivity, so the
#' optimiser sees an unconstrained problem.
#'
#' `ilm_mkL()` builds a full unstructured factor, `ilm_mkD()` a diagonal one (no
#' correlation between categories), and `ilm_mkLam()` a C-by-r loadings matrix for a
#' reduced-rank structure in which `Sigma = Lambda Lambda'`.
#'
#' @param v Numeric vector of unconstrained parameters.
#' @param C Integer. Number of category dimensions.
#' @param r Integer. Rank; `ilm_mkLam()` only.
#' @return A matrix, automatic-differentiation aware for use inside the
#'   likelihood.
#' @references
#' Pinheiro, J. C., & Bates, D. M. (1996). Unconstrained parametrizations for
#' variance-covariance matrices. *Statistics and Computing*, 6, 289--296.
#' @keywords internal
#' @noRd
ilm_mkD <- function(v, C) {                       # diagonal "Cholesky"
  Lv <- AD(numeric(C * C))
  for (i in 1:C) Lv[(i - 1L) * C + i] <- exp(v[i])
  matrix(Lv, C, C)
}

#' Build a covariance factor from unconstrained parameters
#'
#' Optimisers work best when parameters can take any real value, but a
#' covariance matrix must be positive definite. The standard solution is to
#' estimate the **Cholesky factor** `L` (a lower-triangular matrix) and form
#' `Sigma = L L'`, which is positive definite for any `L` with a positive
#' diagonal. Taking `exp()` of the diagonal guarantees that positivity, so the
#' optimiser sees an unconstrained problem.
#'
#' `ilm_mkL()` builds a full unstructured factor, `ilm_mkD()` a diagonal one (no
#' correlation between categories), and `ilm_mkLam()` a C-by-r loadings matrix for a
#' reduced-rank structure in which `Sigma = Lambda Lambda'`.
#'
#' @param v Numeric vector of unconstrained parameters.
#' @param C Integer. Number of category dimensions.
#' @param r Integer. Rank; `ilm_mkLam()` only.
#' @return A matrix, automatic-differentiation aware for use inside the
#'   likelihood.
#' @references
#' Pinheiro, J. C., & Bates, D. M. (1996). Unconstrained parametrizations for
#' variance-covariance matrices. *Statistics and Computing*, 6, 289--296.
#' @keywords internal
#' @noRd
ilm_mkLam <- function(v, C, r) {                  # C x r loadings, lower-trapezoidal
  Lv <- AD(numeric(C * r)); k <- 1L
  for (j in 1:r) for (i in j:C) {
    Lv[(j - 1L) * C + i] <- if (i == j) exp(v[k]) else v[k]; k <- k + 1L
  }
  matrix(Lv, C, r)
}

#' Within-group covariance factor, with its scale fixed
#'
#' For a random-slope term the full random-effect covariance is a **Kronecker
#' product**, `Sigma_cat` times `Sigma_d`: one part describes how effects vary
#' across outcome categories, the other how they vary across the within-group
#' dimensions (intercept, slope, and so on).
#'
#' A Kronecker product carries a scale indeterminacy: multiplying one factor by
#' a constant and dividing the other by the same constant leaves the model
#' unchanged. Left alone this appears as a parameter correlation of exactly 1
#' and the fit cannot separate the two pieces. Fixing `Sigma_d[1,1] = 1` removes
#' it. `Sigma_cat` then carries the scale and is interpretable as the covariance
#' of the *intercept* effects across categories, while `Sigma_d[i,i]` for
#' `i > 1` is that dimension's variance **relative to** the intercept.
#'
#' @param v Numeric vector of unconstrained parameters.
#' @param d Integer. Number of within-group dimensions; 2 for intercept plus
#'   slope.
#' @param d_cor Logical. `TRUE` estimates an intercept-slope correlation;
#'   `FALSE` forces them uncorrelated, equivalent to lme4's
#'   `(1 | g) + (0 + time | g)`. Use `FALSE` when the correlation runs to the
#'   -1 or +1 boundary.
#' @return A d-by-d lower-triangular matrix with entry `[1,1]` fixed at 1.
#' @references
#' Dawid, A. P. (1981). Some matrix-variate distribution theory: notational
#' considerations and a Bayesian application. *Biometrika*, 68(1), 265--274.
#' @keywords internal
#' @noRd
ilm_mkLd <- function(v, d, d_cor = TRUE) {
  Lv <- AD(numeric(d * d)); Lv[1] <- 1; k <- 1L
  if (!d_cor) {
    for (i in 2:d) { Lv[(i - 1L) * d + i] <- exp(v[k]); k <- k + 1L }
  } else {
    for (j in 1:d) for (i in j:d) {
      if (i == 1L && j == 1L) next
      Lv[(j - 1L) * d + i] <- if (i == j) exp(v[k]) else v[k]; k <- k + 1L
    }
  }
  matrix(Lv, d, d)
}

#' Parameter count for a within-group covariance
#'
#' One fewer than the unconstrained count, because `Sigma_d[1,1]` is held at 1
#' to resolve the Kronecker scale indeterminacy described in [ilm_mkLd()].
#'
#' @param d Integer. Number of within-group dimensions.
#' @param d_cor Logical. Whether an intercept-slope correlation is estimated.
#' @return Integer parameter count.
#' @keywords internal
#' @noRd
ilm_npar_d <- function(d, d_cor = TRUE)
  if (d > 1L) { if (d_cor) ilm_ncov(d) - 1L else d - 1L } else 0L

#' Plain-numeric versions of the covariance factor builders
#'
#' Behave identically to [ilm_mkL()] and its relatives but operate on ordinary
#' numbers rather than automatic-differentiation types. Used after fitting, to
#' rebuild covariance matrices from the estimated parameters.
#'
#' @param v Numeric vector of unconstrained parameters.
#' @param C Integer. Number of category dimensions.
#' @param d Integer. Number of within-group dimensions.
#' @param r Integer. Rank.
#' @param d_cor Logical; `ilm_mkLd_num()` only.
#' @return A numeric matrix.
#' @keywords internal
#' @noRd
ilm_mkL_num <- function(v, C) {
  L <- matrix(0, C, C); k <- 1L
  for (j in 1:C) for (i in j:C) { L[i, j] <- if (i == j) exp(v[k]) else v[k]; k <- k + 1L }
  L
}

#' Plain-numeric versions of the covariance factor builders
#'
#' Behave identically to [ilm_mkL()] and its relatives but operate on ordinary
#' numbers rather than automatic-differentiation types. Used after fitting, to
#' rebuild covariance matrices from the estimated parameters.
#'
#' @param v Numeric vector of unconstrained parameters.
#' @param C Integer. Number of category dimensions.
#' @param d Integer. Number of within-group dimensions.
#' @param r Integer. Rank.
#' @param d_cor Logical; `ilm_mkLd_num()` only.
#' @return A numeric matrix.
#' @keywords internal
#' @noRd
ilm_mkLd_num <- function(v, d, d_cor = TRUE) {
  L <- matrix(0, d, d); L[1, 1] <- 1; k <- 1L
  if (!d_cor) { for (i in 2:d) { L[i, i] <- exp(v[k]); k <- k + 1L }; return(L) }
  for (j in 1:d) for (i in j:d) {
    if (i == 1L && j == 1L) next
    L[i, j] <- if (i == j) exp(v[k]) else v[k]; k <- k + 1L
  }
  L
}

#' Plain-numeric versions of the covariance factor builders
#'
#' Behave identically to [ilm_mkL()] and its relatives but operate on ordinary
#' numbers rather than automatic-differentiation types. Used after fitting, to
#' rebuild covariance matrices from the estimated parameters.
#'
#' @param v Numeric vector of unconstrained parameters.
#' @param C Integer. Number of category dimensions.
#' @param d Integer. Number of within-group dimensions.
#' @param r Integer. Rank.
#' @param d_cor Logical; `ilm_mkLd_num()` only.
#' @return A numeric matrix.
#' @keywords internal
#' @noRd
ilm_mkD_num   <- function(v, C) diag(exp(v), C, C)

#' Plain-numeric versions of the covariance factor builders
#'
#' Behave identically to [ilm_mkL()] and its relatives but operate on ordinary
#' numbers rather than automatic-differentiation types. Used after fitting, to
#' rebuild covariance matrices from the estimated parameters.
#'
#' @param v Numeric vector of unconstrained parameters.
#' @param C Integer. Number of category dimensions.
#' @param d Integer. Number of within-group dimensions.
#' @param r Integer. Rank.
#' @param d_cor Logical; `ilm_mkLd_num()` only.
#' @return A numeric matrix.
#' @keywords internal
#' @noRd
ilm_mkLam_num <- function(v, C, r) {
  L <- matrix(0, C, r); k <- 1L
  for (j in 1:r) for (i in j:C) { L[i, j] <- if (i == j) exp(v[k]) else v[k]; k <- k + 1L }
  L
}

#' Parameter names, in the order the optimiser stores them
#'
#' Downstream packages (`car`, `marginaleffects`, `emmeans`) all require that
#' `names(coef())` line up exactly with `colnames(vcov())`, so these names are
#' generated once and reused everywhere.
#'
#' Their fill order **must** mirror the matching matrix builders ([ilm_mkL()],
#' [ilm_mkD()], [ilm_mkLam()], [ilm_mkLd()]). Were the two ever to drift apart, every
#' parameter would be silently mislabelled, which is worse than having no names
#' at all. [ilm_fit()] therefore checks that the number of names matches the
#' number of parameters and stops immediately if it does not.
#'
#' @param pre Character. Term name used as a prefix, for example `"subj"`.
#' @param C Integer. Number of category dimensions.
#' @param d Integer. Number of within-group dimensions.
#' @param r Integer. Rank.
#' @param d_cor Logical; `ilm_nm_ld()` only.
#' @return A character vector of parameter names.
#' @keywords internal
#' @noRd
ilm_nm_tri  <- function(pre, C) { o <- character(0)
  for (j in 1:C) for (i in j:C) o <- c(o, sprintf("%s:L[%d,%d]", pre, i, j)); o }

#' Parameter names, in the order the optimiser stores them
#'
#' Downstream packages (`car`, `marginaleffects`, `emmeans`) all require that
#' `names(coef())` line up exactly with `colnames(vcov())`, so these names are
#' generated once and reused everywhere.
#'
#' Their fill order **must** mirror the matching matrix builders ([ilm_mkL()],
#' [ilm_mkD()], [ilm_mkLam()], [ilm_mkLd()]). Were the two ever to drift apart, every
#' parameter would be silently mislabelled, which is worse than having no names
#' at all. [ilm_fit()] therefore checks that the number of names matches the
#' number of parameters and stops immediately if it does not.
#'
#' @param pre Character. Term name used as a prefix, for example `"subj"`.
#' @param C Integer. Number of category dimensions.
#' @param d Integer. Number of within-group dimensions.
#' @param r Integer. Rank.
#' @param d_cor Logical; `ilm_nm_ld()` only.
#' @return A character vector of parameter names.
#' @keywords internal
#' @noRd
ilm_nm_diag <- function(pre, C) sprintf("%s:logsd[%d]", pre, seq_len(C))

#' Parameter names, in the order the optimiser stores them
#'
#' Downstream packages (`car`, `marginaleffects`, `emmeans`) all require that
#' `names(coef())` line up exactly with `colnames(vcov())`, so these names are
#' generated once and reused everywhere.
#'
#' Their fill order **must** mirror the matching matrix builders ([ilm_mkL()],
#' [ilm_mkD()], [ilm_mkLam()], [ilm_mkLd()]). Were the two ever to drift apart, every
#' parameter would be silently mislabelled, which is worse than having no names
#' at all. [ilm_fit()] therefore checks that the number of names matches the
#' number of parameters and stops immediately if it does not.
#'
#' @param pre Character. Term name used as a prefix, for example `"subj"`.
#' @param C Integer. Number of category dimensions.
#' @param d Integer. Number of within-group dimensions.
#' @param r Integer. Rank.
#' @param d_cor Logical; `ilm_nm_ld()` only.
#' @return A character vector of parameter names.
#' @keywords internal
#' @noRd
ilm_nm_lam  <- function(pre, C, r) { o <- character(0)
  for (j in 1:r) for (i in j:C) o <- c(o, sprintf("%s:Lambda[%d,%d]", pre, i, j)); o }

#' Parameter names, in the order the optimiser stores them
#'
#' Downstream packages (`car`, `marginaleffects`, `emmeans`) all require that
#' `names(coef())` line up exactly with `colnames(vcov())`, so these names are
#' generated once and reused everywhere.
#'
#' Their fill order **must** mirror the matching matrix builders ([ilm_mkL()],
#' [ilm_mkD()], [ilm_mkLam()], [ilm_mkLd()]). Were the two ever to drift apart, every
#' parameter would be silently mislabelled, which is worse than having no names
#' at all. [ilm_fit()] therefore checks that the number of names matches the
#' number of parameters and stops immediately if it does not.
#'
#' @param pre Character. Term name used as a prefix, for example `"subj"`.
#' @param C Integer. Number of category dimensions.
#' @param d Integer. Number of within-group dimensions.
#' @param r Integer. Rank.
#' @param d_cor Logical; `ilm_nm_ld()` only.
#' @return A character vector of parameter names.
#' @keywords internal
#' @noRd
ilm_nm_ld   <- function(pre, d, d_cor) { o <- character(0)
  if (!d_cor) for (i in 2:d) o <- c(o, sprintf("%s:Ld[%d,%d]", pre, i, i))
  else for (j in 1:d) for (i in j:d) {
    if (i == 1L && j == 1L) next
    o <- c(o, sprintf("%s:Ld[%d,%d]", pre, i, j)) }
  o }

#' Properties of a category covariance structure
#'
#' Report how many parameters a structure costs (`ilm_str_npar()`), how many latent
#' values it needs per group (`ilm_str_width()`), and a short printing label
#' (`ilm_str_label()`).
#'
#' The width matters as much as the parameter count. A `diag` structure and an
#' `rr(1)` structure can cost the same number of parameters, but `rr(1)` needs
#' only one latent value per group where `diag` needs `C`. Fewer latent values
#' means more data informing each one, which is what governs whether the Laplace
#' approximation is accurate.
#'
#' @param s A structure specification:
#'   `list(type = "us" | "diag" | "rr", rank = <integer, rr only>)`.
#' @param C Integer. Number of category dimensions.
#' @return An integer, or a character label for `ilm_str_label()`.
#' @keywords internal
#' @noRd
ilm_str_npar  <- function(s, C) switch(s$type, us = ilm_ncov(C), diag = as.integer(C), rr = ilm_nrr(C, s$rank))

#' Properties of a category covariance structure
#'
#' Report how many parameters a structure costs (`ilm_str_npar()`), how many latent
#' values it needs per group (`ilm_str_width()`), and a short printing label
#' (`ilm_str_label()`).
#'
#' The width matters as much as the parameter count. A `diag` structure and an
#' `rr(1)` structure can cost the same number of parameters, but `rr(1)` needs
#' only one latent value per group where `diag` needs `C`. Fewer latent values
#' means more data informing each one, which is what governs whether the Laplace
#' approximation is accurate.
#'
#' @param s A structure specification:
#'   `list(type = "us" | "diag" | "rr", rank = <integer, rr only>)`.
#' @param C Integer. Number of category dimensions.
#' @return An integer, or a character label for `ilm_str_label()`.
#' @keywords internal
#' @noRd
ilm_str_width <- function(s, C) switch(s$type, us = C, diag = C, rr = s$rank)

#' Properties of a category covariance structure
#'
#' Report how many parameters a structure costs (`ilm_str_npar()`), how many latent
#' values it needs per group (`ilm_str_width()`), and a short printing label
#' (`ilm_str_label()`).
#'
#' The width matters as much as the parameter count. A `diag` structure and an
#' `rr(1)` structure can cost the same number of parameters, but `rr(1)` needs
#' only one latent value per group where `diag` needs `C`. Fewer latent values
#' means more data informing each one, which is what governs whether the Laplace
#' approximation is accurate.
#'
#' @param s A structure specification:
#'   `list(type = "us" | "diag" | "rr", rank = <integer, rr only>)`.
#' @param C Integer. Number of category dimensions.
#' @return An integer, or a character label for `ilm_str_label()`.
#' @keywords internal
#' @noRd
ilm_str_label <- function(s) if (s$type == "rr") sprintf("rr(%d)", s$rank) else s$type

#' Suggest a reduced rank the data can support
#'
#' Given a grouping factor with `nl` levels, returns the largest rank whose
#' parameter count sits comfortably within the available information. This turns
#' a failing check into a specific, actionable suggestion -- "use `rr(3)`" --
#' rather than vague advice to simplify.
#'
#' The target is `nl/6` parameters rather than `nl/3`. The latter is roughly
#' where fits were observed to break down, so aiming at half of it leaves the
#' recommendation in safe territory instead of at the cliff edge.
#'
#' @param C Integer. Number of category dimensions.
#' @param nl Integer. Number of levels of the grouping factor.
#' @return Integer rank, or `0` if even rank 1 is more than the data support.
#' @keywords internal
#' @noRd
ilm_rec_rank <- function(C, nl) {
  r <- 0L
  for (rr in seq_len(max(1L, C - 1L))) if (ilm_nrr(C, rr) <= nl / 6) r <- rr
  r
}

#' Normalise the random-effects specification
#'
#' Accepts the three ways a random term may be supplied and converts them to one
#' internal form: a bare grouping vector (random intercept), a
#' `list(group =, Z =)` (random slopes), or a `list(basis =)` (a penalised
#' smooth, after reparameterisation by `ilm_smooth()`).
#'
#' @param re_list Named list of random-effect specifications.
#' @param N Integer. Number of rows of data.
#' @return A named list holding, for each term, `kind` (`"group"` or `"basis"`),
#'   `group`, `Z`, `basis`, `d` (within-group dimensions) and `nl` (number of
#'   levels, or of basis functions).
#' @keywords internal
#' @noRd
ilm_norm_re <- function(re_list, N) lapply(re_list, function(e) {
  if (is.list(e) && !is.null(e$basis)) {          # spline / basis term
    Bs <- as.matrix(e$basis)
    list(kind = "basis", basis = Bs, group = NULL, Z = NULL, d = 1L, nl = ncol(Bs))
  } else if (is.list(e)) {                        # random slopes
    Z <- as.matrix(e$Z); g <- factor(e$group)
    list(kind = "group", basis = NULL, group = as.integer(g), Z = Z,
         d = ncol(Z), nl = nlevels(g))
  } else {                                        # random intercept
    g <- factor(e)
    list(kind = "group", basis = NULL, group = as.integer(g),
         Z = matrix(1, N, 1), d = 1L, nl = nlevels(g))
  }
})

#' Convert an mgcv smooth into fixed and random parts
#'
#' A penalised spline can be rewritten as a small set of **fixed** columns -- the
#' unpenalised "null space", typically a linear trend -- plus a set of **random**
#' columns whose coefficients are independent and normally distributed. This is
#' the classical mixed-model representation of a smoother, and it is what allows
#' a spline to be fitted by exactly the same machinery as a random intercept: no
#' new covariance code is required, only a design matrix.
#'
#' `mgcv::smooth2random(type = 2)` performs the reparameterisation.
#'
#' @param spec A smooth specification, for example the result of
#'   `mgcv::s(x, k = 10)` or `mgcv::t2(x, z)`. Note that `mgcv::te()` is
#'   rejected by `smooth2random()`; use `t2()` for tensor products.
#' @param data A data frame containing the smooth's variables.
#' @return A list with `Xf` (fixed null-space columns), `rand` (a list of random
#'   basis blocks, more than one for a tensor product), `sm` (the original
#'   smooth object) and `re` (the full reparameterisation, retained because
#'   [predict.ilm_model()] needs it to rebuild the basis at new covariate values).
#' @references
#' Wood, S. N. (2017). *Generalized Additive Models: An Introduction with R*,
#' 2nd ed. Chapman & Hall/CRC. Section 5.4 develops the mixed-model
#' representation of penalised smoothers.
#'
#' Wood, S. N. (2003). Thin plate regression splines. *Journal of the Royal
#' Statistical Society, Series B*, 65(1), 95--114.
#' @keywords internal
#' @noRd
ilm_smooth <- function(spec, data) {
  ## A variable whose name needs backticks is given a syntactic stand-in for
  ## mgcv, which cannot take one (see ilm_add_standins()). The label keeps the
  ## name as written, and `name_map` goes with the smooth so that prediction
  ## can put the same stand-ins onto new data.
  vars <- ilm_spec_vars(spec)
  bad <- vars[ilm_bq(vars) != vars]
  map <- NULL
  if (length(bad)) {
    stand <- paste0(".ilm_sm", seq_along(bad), ".")
    map <- stats::setNames(bad, stand)
    spec <- ilm_swap_spec(spec, bad, stand, label = FALSE)
    data <- ilm_add_standins(data, map)
  }
  sm <- mgcv::smoothCon(spec, data = data, absorb.cons = TRUE, scale.penalty = TRUE)[[1]]
  re <- mgcv::smooth2random(sm, "", type = 2)
  list(Xf = re$Xf, rand = re$rand, sm = sm, re = re, name_map = map)
}

#' Build the table of model checks
#'
#' Each check records a status, what was measured, why it matters, and what to do
#' about it. Statuses are `"OK"`, `"WARN"`, `"FAIL"`, `"BOUNDARY"` (a parameter
#' pinned at zero, where the usual symmetric tests do not apply) and
#' `"INCONCLUSIVE"` (not enough information to judge). The last is deliberate:
#' reporting "assumption satisfied" when the honest answer is "cannot tell" is
#' worse than saying nothing.
#'
#' @param ck An existing checks data frame.
#' @param check Character. Name of the check.
#' @param status Character. One of the statuses above.
#' @param detail Character. What was measured.
#' @param cause Character. Why it matters.
#' @param suggestion Character. What to do about it.
#' @return A checks data frame.
#' @keywords internal
#' @noRd
ilm_new_checks <- function() data.frame(check=character(), status=character(),
                                    detail=character(), cause=character(),
                                    suggestion=character(), stringsAsFactors=FALSE)

#' Build the table of model checks
#'
#' Each check records a status, what was measured, why it matters, and what to do
#' about it. Statuses are `"OK"`, `"WARN"`, `"FAIL"`, `"BOUNDARY"` (a parameter
#' pinned at zero, where the usual symmetric tests do not apply) and
#' `"INCONCLUSIVE"` (not enough information to judge). The last is deliberate:
#' reporting "assumption satisfied" when the honest answer is "cannot tell" is
#' worse than saying nothing.
#'
#' @param ck An existing checks data frame.
#' @param check Character. Name of the check.
#' @param status Character. One of the statuses above.
#' @param detail Character. What was measured.
#' @param cause Character. Why it matters.
#' @param suggestion Character. What to do about it.
#' @return A checks data frame.
#' @keywords internal
#' @noRd
ilm_add_check <- function(ck, check, status, detail, cause = "", suggestion = "")
  rbind(ck, data.frame(check=check, status=status, detail=detail, cause=cause,
                       suggestion=suggestion, stringsAsFactors=FALSE))

#' Checks that can be made before fitting
#'
#' These examine the data and the requested structure only, so they cost nothing
#' and can warn you before a long fit begins. They cover sparse outcome
#' categories, whether each grouping factor has enough levels to support its
#' covariance structure, observations per group, the type of any weights
#' supplied, and the latent budget.
#'
#' @section The latent budget:
#' The most useful check of the set. A mixed model introduces one unobserved
#' ("latent") value per random effect per group. The Laplace approximation is
#' accurate when each latent value is informed by a reasonable amount of data and
#' degrades when it is not.
#'
#' Crucially this is a **global** constraint. A term can pass its own
#' level-count check comfortably and still break the fit, because the model as a
#' whole has more latent values than the data can support. Per-term checks cannot
#' see that; this one can.
#'
#' Thresholds were calibrated by simulation rather than assumed. Fits stayed well
#' behaved at roughly three or more observations per latent value and collapsed
#' below that.
#'
#' @param y Integer vector of category codes.
#' @param J Integer. Number of outcome categories.
#' @param re Normalised random-effects list from [ilm_norm_re()].
#' @param re_struct Named list of covariance structures.
#' @param ar Optional AR(1) specification.
#' @param weights Optional numeric vector of observation weights. When supplied,
#'   information is measured by `sum(weights)` rather than the row count.
#' @return A data frame of checks.
#' @references
#' Joe, H. (2008). Accuracy of Laplace approximation for discrete response mixed
#' models. *Computational Statistics & Data Analysis*, 52(12), 5066--5074.
#' @keywords internal
#' @noRd
ilm_precheck <- function(y, J, re, re_struct, ar = NULL, weights = NULL,
                     family = NULL) {
  C <- J - 1L
  ## With weights, information is carried by sum(w), not the row count: a
  ## thousand rows each worth one trial and ten rows each worth a hundred are
  ## not the same amount of data, and the latent budget must use the latter.
  if (is.null(weights)) weights <- rep(1, length(y))
  N <- sum(weights); ck <- ilm_new_checks()
  ## sparse-category checking only makes sense for the multinomial
  if (is.null(family) || identical(family$name, "multinomial")) {
  tabc <- vapply(seq_len(J), function(j) sum(weights[y == j]), 0); mn <- min(tabc)
  ck <- ilm_add_check(ck, "category_counts",
    if (mn < 10) "FAIL" else if (mn < 30) "WARN" else "OK",
    sprintf("smallest category n = %g of %g total (J = %d)", mn, N, J),
    if (mn < 30) "sparse outcome category" else "",
    if (mn < 30) "collapse rare categories, or expect unstable category-specific effects" else "")
  }

  ## Frequency weights ONLY.  Measured, with a NON-informative 10% sample and
  ## w = 10 (where weighting should be a harmless no-op):
  ##   fixed effects inflated 20-130% away from truth (0.9 -> 1.07, 0.2 -> 0.47)
  ##   standard errors understated ~3.5x (ratio 0.28 vs expected 1/sqrt(10)=0.32)
  ## The estimate damage is specific to MIXED models: the weight multiplies the
  ## conditional likelihood INSIDE the Laplace integral, so each cluster looks to
  ## carry w times its real information and the random effects are under-shrunk.
  ## Correct survey handling needs scaled level-specific weights (Pfeffermann et
  ## al. 1998; Rabe-Hesketh & Skrondal 2006) plus a sandwich variance.  Neither
  ## is implemented, so this refuses to stay quiet about it.
  noninteger <- any(abs(weights - round(weights)) > 1e-8)
  ## For a binomial response the weights ARE trial counts, so a constant 60
  ## trials a row is ordinary rather than suspicious; neither rule below
  ## applies to it.
  binom <- identical(family$name, "binomial")
  flat <- !binom && length(weights) > 1L && all(weights == weights[1]) &&
    weights[1] > 1
  ## Sampling weights are often rounded, so being whole numbers proves nothing.
  ## What does give them away is the total: replicate counts sum to the number
  ## of observations they stand for, while sampling weights sum to a
  ## POPULATION. Grouped binomial data is the one legitimate case where the
  ## weights are genuinely large -- they are trial counts -- so it is excluded
  ## rather than warned about every time.
  inflated <- !noninteger && !all(weights == 1) &&
    sum(weights) > 10 * length(weights) && !binom
  wst <- if (noninteger) "FAIL" else if (flat || inflated) "WARN" else "OK"
  mixed <- length(re) > 0L
  ck <- ilm_add_check(ck, "weights_type", wst,
    sprintf("%s weights; range %g to %g, total %g",
            if (all(weights == 1)) "none (all 1)" else if (noninteger) "non-integer" else "integer",
            min(weights), max(weights), sum(weights)),
    if (wst != "OK") paste0(
      if (noninteger) "non-integer weights cannot be replicate counts, so these look like sampling weights"
      else if (inflated) sprintf("the weights total %.0f across %d rows, which is a population rather than a count of replicates, so these look like sampling weights",
                                 sum(weights), length(weights))
      else "every row carries the same weight > 1, which is a scaling rather than replicate counts",
      "; ilm_model() treats weights as FREQUENCIES. The standard errors are then too small by at least sqrt(n / sum(w)), because the likelihood believes it saw sum(w) observations",
      if (mixed)
        ", and in a model with random effects the COEFFICIENTS move too: the weight multiplies the conditional likelihood inside the Laplace integral, so each cluster appears to carry w times its real information and the random effects are under-shrunk"
      else ". The COEFFICIENTS are fine in a fixed-effects fit -- a weighted likelihood is design-consistent for the population parameter -- which is what makes this easy to miss") else "",
    if (wst != "OK") paste0(
      "aggregate genuine replicate counts instead; for a complex sample pass design = ilm_design(weights = , ids = , strata = ) and read the standard errors from ilm_svy_coef()",
      if (mixed) ", which refuses a model with random effects for the reason above" else "") else "")

  lat <- integer(0)
  for (nm in names(re)) {
    e <- re[[nm]]; s <- re_struct[[nm]]; nl <- e$nl
    np <- ilm_str_npar(s, C) + ilm_npar_d(e$d, isTRUE(s$d_cor))
    w  <- ilm_str_width(s, C); lat[nm] <- nl * e$d * w
    ratio <- nl / np
    ## Thresholds calibrated on the J sweep: 4.0 levels/parameter still fit, 2.7
    ## collapsed to a rank-deficient covariance.  They apply to GROUPING factors,
    ## where each level is an independent replicate.  Basis functions are not
    ## replicates -- the 2-D t2() test fit cleanly at 3.0 per parameter -- so the
    ## rule is not applied there until it has its own calibration.
    st <- if (e$kind == "basis") "OK"
          else if (ratio < 3) "FAIL" else if (ratio < 6) "WARN" else "OK"
    sug <- ""
    if (st != "OK") {
      if (s$type %in% c("us", "diag")) {
        r <- ilm_rec_rank(C, nl)
        ## Prefer rr over diag: at equal parameter count rr also cuts the latent
        ## width from C to r, which ran ~2.6x faster over 3 seeds in
        ## test_diag_vs_rr.R.  Both converged there, so this is a cost argument,
        ## not a demonstrated stability one.
        sug <- if (r >= 1)
          sprintf("use rr(%d) for this term: %d parameters instead of %d, and latent width %d instead of %d",
                  r, ilm_nrr(C, r), ilm_str_npar(s, C), r, C)
        else sprintf("use rr(1): %d parameters, latent width 1 (diag costs the same but leaves the width at %d)", C, C)
      } else if (isTRUE(s$rank > 1L)) {
        sug <- sprintf("lower the rank below %d for this term, or drop this grouping factor", s$rank)
      } else sug <- "the rank is already 1; pool levels of this grouping factor, or drop it"
    }
    unit <- if (e$kind == "basis") "basis functions" else "levels"
    ck <- ilm_add_check(ck, paste0("re_levels[", nm, "]"), st,
      sprintf("%d %s for %d covariance parameters (%.1f per parameter); %s, C = %d%s",
              nl, unit, np, ratio, ilm_str_label(s), C,
              if (e$d > 1L) sprintf(", %d random-effect dimensions", e$d) else ""),
      if (st != "OK") sprintf("too few %s for this category covariance structure; below about 3 per parameter it fails outright", unit) else "",
      sug)
    if (e$kind != "basis") {          # a basis term has no grouping levels
      ope <- N / nl
      ck <- ilm_add_check(ck, paste0("obs_per_level[", nm, "]"),
        if (ope < 2) "FAIL" else if (ope < 5) "WARN" else "OK",
        sprintf("%.1f observations per level", ope),
        if (ope < 5) "little information per group" else "",
        if (ope < 5) "drop this grouping factor or pool levels" else "")
    }
  }
  if (!is.null(ar)) lat["ar"] <- ar$n_cell * C

  ## ---- GLOBAL latent budget -------------------------------------------------
  ## The binding constraint is total observations per latent VALUE, not any
  ## per-term level ratio.  Two independent experiments agree on the threshold:
  ##   AR experiment  -- 1 obs/latent failed, 4 recovered cleanly
  ##   J sweep        -- 3.13 obs/latent fit, 2.88 collapsed
  ## A term can pass its own level check and still break the fit this way.
  ## Nothing is integrated out in a fixed-effects-only model, so the latent
  ## budget does not apply and reporting it would be meaningless.
  tot <- sum(lat)
  if (tot == 0L) return(ck)
  ratio <- N / tot

  ## THE BUDGET DOES NOT APPLY TO A GAUSSIAN MODEL.
  ##
  ## Everything above measures how much the Laplace approximation is being
  ## asked to do. With a gaussian response and gaussian latents the model is
  ## linear-Gaussian, the integral over the latents is a gaussian integral, and
  ## the Laplace "approximation" evaluates it EXACTLY -- there is no
  ## approximation error to run out of observations for. Measured: illume's
  ## log-likelihood agrees with lme4 and nlme, which compute that integral in
  ## closed form, to four decimal places.
  ##
  ## This matters because CAR(1) puts one latent value under every observation
  ## by construction, so the thresholds below -- calibrated on MULTINOMIAL
  ## experiments, where one categorical observation says very little about its
  ## own latent -- would FAIL every correctly specified continuous-time model.
  ## Over 400 replicates at exactly 1.00 observations per latent, Wald coverage
  ## for a fixed effect came back 0.955, 0.937 and 0.949 against a nominal
  ## 0.95, with standard-error ratios of 0.98, 0.98 and 0.97.
  ##
  ## What can still go wrong for a gaussian model is identifiability rather
  ## than approximation -- a latent process and the residual describing the
  ## same variation -- and that surfaces as a non-positive-definite Hessian,
  ## which has its own check after the fit.
  gaus <- !is.null(family) &&
    identical(if (is.list(family)) family$name else family, "gaussian")
  st <- if (gaus) "OK" else
        if (ratio < 3) "FAIL" else if (ratio < 5) "WARN" else "OK"
  sug <- ""
  if (st != "OK") {
    big <- names(lat)[which.max(lat)]
    target <- N / 5                                  # latent budget for the OK band
    if (big %in% names(re)) {
      e <- re[[big]]; s <- re_struct[[big]]
      keep <- max(1L, as.integer(floor((lat[big] - (tot - target)) / (e$nl * e$d))))
      sug <- sprintf("total latent budget for %g observations is about %d; term '%s' contributes %d (%s, width %d). rr(%d) would cut it to %d",
                     N, as.integer(target), big, lat[big], ilm_str_label(s),
                     ilm_str_width(s, C), keep, as.integer(e$nl * e$d * keep))
    } else sug <- sprintf("total latent budget for %g observations is about %d; the AR term contributes %d -- coarsen its time grid",
                          N, as.integer(target), lat[["ar"]])
  }
  ck <- ilm_add_check(ck, "latent_budget", st,
    sprintf("%.2f observations per latent value (%g observations, %d latent values: %s)%s",
            ratio, N, tot, paste(sprintf("%s %d", names(lat), lat), collapse = ", "),
            if (gaus) "; not a constraint for a gaussian response, where the Laplace approximation is exact" else ""),
    if (st != "OK") "too few observations per latent value; the Laplace approximation attenuates the variance components and the covariance estimates go rank deficient, even when every individual term passes its own level check" else "",
    sug)

  if (!is.null(ar)) {
    nlat <- ar$n_cell; r2 <- N / nlat
    car <- identical(ar$type, "car1")
    ## For a gaussian response this ratio is reported but not failed on, for
    ## the reason given above: one latent per observation is what continuous
    ## time produces, and coverage there is nominal.
    thin <- !gaus && r2 < 4
    ck <- ilm_add_check(ck, "obs_per_ar_latent",
      if (gaus) "OK" else if (r2 < 2) "FAIL" else if (r2 < 4) "WARN" else "OK",
      sprintf("%.2f observations per %s latent time point (%d time points)%s",
              r2, if (car) "CAR(1)" else "AR", as.integer(nlat),
              if (gaus && r2 < 4) "; expected for continuous time and not a problem for a gaussian response" else ""),
      if (thin) "the latent process carries about one categorical observation per latent value; Laplace attenuates the variance components and rho is driven toward the boundary" else "",
      if (thin) {
        if (car) "coarsen the time passed to ilm_car1() so observations share a latent value, or replace the term with s(time) plus a random slope"
        else "coarsen the AR time grid, use a reduced-rank AR, or replace AR with s(time) plus a random slope"
      } else "")
  }
  ck
}

#' Checks that require a fitted model
#'
#' Covers optimiser convergence, the gradient at the solution, whether the
#' Hessian is positive definite, collapsed variances, rank deficiency of each
#' category covariance, the within-group intercept-slope covariance, an AR
#' parameter at its boundary, and near-aliased parameter pairs.
#'
#' Two of these are subtler than they appear. A category covariance can be
#' numerically **singular while every pairwise correlation remains modest**, so
#' rank is assessed through eigenvalues rather than correlations. And a
#' reduced-rank term is *supposed* to be singular, so there the test is applied
#' within the term's own declared subspace instead, which avoids a guaranteed
#' false alarm.
#'
#' Checks also cross-reference one another: a rank-deficient covariance reports a
#' different likely cause depending on whether a pre-fit check had already
#' flagged that same term.
#'
#' @param opt Optimiser output from [stats::nlminb()].
#' @param obj The `RTMB` objective object.
#' @param sdr Output of `TMB::sdreport()`.
#' @param C Integer. Number of category dimensions.
#' @param has_ar Logical. Whether an AR(1) or CAR(1) term was fitted.
#' @param gap For CAR(1), the gaps between consecutive observations, so the
#'   boundary check can judge the correlation across a typical gap rather than
#'   across one time unit.
#' @param pre The pre-fit checks, used for cross-referencing.
#' @param Sig Named list of fitted category covariance matrices.
#' @param Sigd Named list of fitted within-group covariance matrices.
#' @param re_struct Named list of covariance structures.
#' @param kinds Named list marking each term `"group"` or `"basis"`.
#' @param pnames Character vector of parameter names.
#' @return A data frame of checks.
#' @keywords internal
#' @noRd
ilm_postcheck <- function(opt, obj, sdr, C, has_ar, pre, Sig, Sigd, re_struct, kinds = NULL,
                      pnames = NULL, gap = NULL, hess = NULL,
                      boundary = character(0)) {
  how <- if (is.null(hess)) "tmb" else hess$how
  held <- if (is.null(hess)) character(0) else hess$held
  kind_of <- function(nm) if (!is.null(kinds) && nm %in% names(kinds)) kinds[[nm]] else "group"
  ck <- ilm_new_checks(); g <- max(abs(obj$gr(opt$par)))
  ## nlminb's code 8, "false convergence", means it could not verify a descent
  ## direction from where it stopped -- not that it is far from a solution. The
  ## first-order condition is the gradient, and where that is small the
  ## stopping code is a heuristic disagreeing with the arithmetic. Measured on
  ## 60 flexible parametric survival fits, the 21 reporting code 8 had a
  ## maximum gradient of 3.7e-03 against 3.1e-03 for those reporting success,
  ## and recovered the same coefficients to within Monte Carlo error: -0.542,
  ## 0.273, 0.699 against -0.539, 0.251, 0.696, for a truth of -0.541, 0.258,
  ## 0.690. Grading those FAIL discards a third of perfectly good fits.
  ck <- ilm_add_check(ck, "optimizer",
    if (opt$convergence == 0) "OK" else if (g <= 1e-2) "WARN" else "FAIL",
    sprintf("nlminb code %d (%s)", opt$convergence, opt$message),
    if (opt$convergence != 0)
      paste0("optimizer stopped without meeting its tolerance",
             if (g <= 1e-2)
               "; the gradient is small, so this is its stopping rule rather than the fit"
             else "") else "", "")
  ck <- ilm_add_check(ck, "gradient",
    if (g > 1e-2) "FAIL" else if (g > 1e-3) "WARN" else "OK",
    sprintf("max |gradient| = %.2e", g),
    if (g > 1e-3) "not at a stationary point" else "",
    if (g > 1e-3) "restart from the current estimates, or simplify the random structure" else "")
  pd <- isTRUE(sdr$pdHess)
  ## a held term's rows are NA by design, so only the rest is judged
  cfd <- diag(sdr$cov.fixed)
  judged <- if (length(held)) is.finite(cfd) else rep(TRUE, length(cfd))
  anyNaN <- any(!is.finite(cfd[judged])) || any(cfd[judged] <= 0)
  lb <- pre$status[pre$check == "latent_budget"]
  if (identical(how, "boundary") && !anyNaN) {
    hq <- paste(sQuote(held, FALSE), collapse = ", ")
    ck <- ilm_add_check(ck, "hessian", "BOUNDARY",
      sprintf("not positive definite along the covariance of %s, which is held at its estimate", hq),
      paste0("the covariance of ", hq, " sits at the edge of its range -- a ",
             "variance of zero, or a correlation of +/-1 -- so the likelihood ",
             "is flat in that direction and that covariance has no standard error"),
      paste0("the fixed effects, their standard errors and tests ARE usable: ",
             "they are computed with that covariance held at its estimate, as ",
             "lme4 does when the full Hessian fails. The covariance itself is ",
             "not. If the term is not needed, drop it: at a variance of zero ",
             "the fixed effects do not change"))
  } else
  ck <- ilm_add_check(ck, "hessian",
    if (!pd || anyNaN) "FAIL" else "OK",
    sprintf("positive definite: %s%s; non-finite or non-positive variances: %s",
            pd, if (identical(how, "recomputed") && pd)
              " (recomputed: the first, coarser Hessian was not)" else "",
            anyNaN),
    if (!pd || anyNaN) paste0("likelihood is flat in at least one direction (a ridge)",
      if (length(lb) && lb %in% c("WARN", "FAIL")) "; the latent budget check also failed, which is the likely driver" else "") else "",
    if (!pd || anyNaN) "simplify the covariance structure; standard errors are unusable until this is resolved" else "")
  ## A fixed-effects-only model reports nothing, and asking for an empty
  ## summary warns; treat it as no reported quantities rather than an error.
  rp <- tryCatch(suppressWarnings(summary(sdr, "report")),
                 error = function(e) matrix(numeric(0), 0, 4))
  est <- if (nrow(rp)) rp[, 1] else numeric(0)
  nmv <- if (nrow(rp)) rownames(rp) else character(0)
  ## Grouping terms only.  For a basis (spline) term a component shrinking to
  ## zero is how penalisation WORKS -- it is reported by smooth_shrinkage below.
  gsd <- unlist(lapply(names(Sig), function(nm)
    if (kind_of(nm) == "basis") NULL else sqrt(diag(Sig[[nm]]))))
  if (length(gsd)) {
    bad <- gsd < 1e-3
    ck <- ilm_add_check(ck, "variance_boundary", if (any(bad)) "WARN" else "OK",
      sprintf("smallest random-effect SD across grouping terms = %.2e", min(gsd)),
      if (any(bad)) "a random-effect variance has collapsed to zero" else "",
      if (any(bad)) "drop that term; the data do not support it" else "")
  }
  ## Penalised-out smooth components: informative, not a failure -- unless the
  ## whole smooth vanishes, which means the term has no effect at all.
  for (nm in names(Sig)) {
    if (kind_of(nm) != "basis") next
    sv <- sqrt(diag(Sig[[nm]])); z <- sv < 1e-3
    ck <- ilm_add_check(ck, paste0("smooth_shrinkage[", nm, "]"),
      if (all(z)) "FAIL" else "OK",
      sprintf("%d of %d category components penalised to zero (SDs: %s)",
              sum(z), length(sv), paste(sprintf("%.3g", sv), collapse = ", ")),
      if (all(z)) "the entire smooth was penalised out: these data show no smooth effect for this term" else "",
      if (all(z)) sprintf("drop term '%s' from the model", nm) else "")
  }
  ## Eigenvalue rank check.  For C > 2 a category covariance can be numerically
  ## singular while every pairwise correlation stays modest.  A reduced-rank term
  ## is SUPPOSED to be singular, so there we test conditioning within its own
  ## r-dimensional subspace instead.
  for (nm in names(Sig)) {
    if (kind_of(nm) == "basis") next      # rank deficiency here is penalisation
    s <- if (nm %in% names(re_struct)) re_struct[[nm]] else list(type = "us")
    ev <- sort(eigen(Sig[[nm]], only.values = TRUE, symmetric = TRUE)$values, decreasing = TRUE)
    keep <- if (s$type == "rr") seq_len(s$rank) else seq_along(ev)
    rt <- min(ev[keep]) / max(ev[keep])
    lv <- pre$status[pre$check == paste0("re_levels[", nm, "]")]
    ck <- ilm_add_check(ck, paste0("sigma_rank[", nm, "]"),
      if (rt < 1e-6) "FAIL" else if (rt < 1e-3) "WARN" else "OK",
      if (s$type == "rr") sprintf("eigenvalue ratio within the rank-%d subspace = %.1e", s$rank, rt)
      else sprintf("eigenvalue ratio min/max = %.1e", rt),
      if (rt < 1e-3) paste0(
        "estimated category covariance is numerically rank deficient (fewer effective dimensions than it was given)",
        if (length(lv) && lv %in% c("WARN", "FAIL"))
          "; the pre-fit level count for this term was already flagged"
        else if (length(lb) && lb %in% c("WARN", "FAIL"))
          "; this term's own level count was adequate, so the global latent budget is the likely cause" else "") else "",
      if (rt < 1e-3) {
        if (s$type == "rr") sprintf("lower the rank below %d for this term", s$rank)
        else "use a reduced-rank (rr) category covariance for this term"
      } else "")
  }
  ## Within-group (intercept/slope) covariance.  Distinct from sigma_rank, which
  ## looks at the CATEGORY dimension; this one catches an intercept-slope
  ## correlation running to +/-1, which makes Sigma_d singular.
  ## Sigma_d can be singular for two DIFFERENT reasons, which need different
  ## remedies: a slope variance collapsing to zero, or a correlation at +/-1.
  ## Report them separately -- an eigenvalue ratio alone cannot tell them apart.
  for (nm in names(Sigd)) {
    Sd <- Sigd[[nm]]; sv <- sqrt(diag(Sd)); rel <- sv[-1] / sv[1]
    Rd <- Sd / outer(sv, sv); mxc <- max(abs(Rd[upper.tri(Rd)]))
    if (any(rel < 1e-3)) {
      st <- "FAIL"
      dt <- sprintf("smallest slope SD, relative to the intercept, = %.2e (correlation %.3f)", min(rel), mxc)
      wy <- "a random-slope variance has collapsed to zero: these data do not support an independent slope for this term"
      tr <- sprintf("remove the random slope from term '%s', keeping its random intercept", nm)
    } else if (mxc > 0.99) {
      st <- if (mxc > 0.999) "FAIL" else "WARN"
      dt <- sprintf("largest |intercept-slope correlation| = %.4f", mxc)
      wy <- "intercept and slope are (near) perfectly correlated, so Sigma_d is singular and the slope adds no independent dimension"
      tr <- sprintf("set d_cor = FALSE for term '%s' to drop the intercept-slope correlation, or remove the random slope", nm)
    } else {
      st <- "OK"
      dt <- sprintf("slope/intercept SD ratio %.3f, |correlation| %.3f", min(rel), mxc)
      wy <- ""; tr <- ""
    }
    ck <- ilm_add_check(ck, paste0("sigma_within[", nm, "]"), st, dt, wy, tr)
  }
  if (has_ar) {
    rho <- est[nmv == "rho"]
    lat <- pre$status[pre$check == "obs_per_ar_latent"]
    ## Judge the correlation across a TYPICAL gap, not across one time unit.
    ## Under CAR(1) rho is per unit of time, so with fine units (days, seconds)
    ## it sits at 0.999-something for a process that is nowhere near a random
    ## walk over the window actually observed.
    car <- !is.null(gap) && length(gap)
    eff <- if (length(rho) && car) rho^stats::median(gap) else rho
    lbl <- if (car) sprintf("rho = %.4f per time unit, %.4f across the median gap of %s",
                            if (length(rho)) rho else NA_real_,
                            if (length(eff)) eff else NA_real_,
                            signif(stats::median(gap), 3))
           else sprintf("rho = %.4f", if (length(rho)) rho else NA_real_)
    ck <- ilm_add_check(ck, "rho_boundary",
      if (length(eff) && abs(eff) > 0.99) "FAIL" else if (length(eff) && abs(eff) > 0.95) "WARN" else "OK",
      lbl,
      if (length(eff) && abs(eff) > 0.95)
        paste0("the process is near a random walk",
               if (length(lat) && lat %in% c("WARN", "FAIL"))
                 "; the obs-per-AR-latent check also failed, which is the likely driver" else "") else "",
      if (length(eff) && abs(eff) > 0.95)
        paste0("coarsen the ", if (car) "CAR(1)" else "AR", " time grid or drop the term") else "")
  }
  cf <- sdr$cov.fixed
  pnm <- if (!is.null(pnames) && !is.null(cf) && length(pnames) == ncol(cf)) pnames
         else make.unique(names(obj$par))
  ## a term held at its boundary has no covariance to correlate; judge the rest
  if (!is.null(cf) && length(held)) {
    ok <- is.finite(diag(cf))
    cf <- cf[ok, ok, drop = FALSE]; pnm <- pnm[ok]
  }
  if (!is.null(cf) && length(cf) && all(is.finite(cf)) && all(diag(cf) > 0)) {
    cm <- cov2cor(cf); cm[!upper.tri(cm)] <- 0
    mx <- max(abs(cm)); ij <- which(abs(cm) == mx, arr.ind = TRUE)[1, ]
    pair <- sprintf("%s <-> %s", pnm[ij[1]], pnm[ij[2]])
    ck <- ilm_add_check(ck, "parameter_aliasing",
      if (mx > 0.995) "FAIL" else if (mx > 0.95) "WARN" else "OK",
      sprintf("largest |parameter correlation| = %.3f (%s)", mx, pair),
      if (mx > 0.95) sprintf("%s cannot be separated by these data", pair) else "",
      if (mx > 0.95) "remove one of the competing terms, or fix one of the two parameters" else "")
  } else {
    ck <- ilm_add_check(ck, "parameter_aliasing", "INCONCLUSIVE",
      "covariance of fixed parameters unavailable", "Hessian not usable", "")
  }
  ## A covariance at its boundary -- whether it had to be held for the Hessian
  ## or the Hessian was positive definite anyway -- is not a failure of the
  ## fit when the fixed effects are usable. The checks that say WHY it is
  ## there keep their explanation and remedy, but report BOUNDARY rather than
  ## FAIL, so the verdict agrees with the standard errors above it. Measured:
  ## a site covariance fitted at a correlation of -1, with a positive definite
  ## Hessian and standard errors within 4% of the model without the term,
  ## used to be graded FAIL, with "the standard errors are not usable".
  at <- union(held, boundary)
  if (!identical(how, "none") && length(at)) {
    tags <- c(sprintf("sigma_rank[%s]", at), sprintf("sigma_within[%s]", at),
              sprintf("smooth_shrinkage[%s]", at),
              if ("ar" %in% at) "rho_boundary")
    hit <- ck$check %in% tags & ck$status == "FAIL"
    ck$status[hit] <- "BOUNDARY"
  }
  ck
}

## ---- when the Hessian says no -------------------------------------------------
##
## TMB's sdreport() judges the Hessian of the Laplace objective by differencing
## its gradient ONCE, with a fixed step of 1e-3. That is coarse wherever the
## curvature is small -- near a variance boundary above all -- and it calls
## Hessians indefinite that are not. glmmTMB recomputes the Hessian more
## accurately before giving up (numDeriv::jacobian() in its finalizeTMB()),
## which is much of why it reports fewer failures on the same data. This is the
## same idea without the dependency: central differences of the EXACT gradient
## RTMB already provides, with one Richardson step.
##
## When a variance component really does sit at the edge of its range -- an SD
## of zero, a correlation of +/-1, an AR(1) correlation at its limit -- the
## likelihood is flat in that direction and no Hessian over EVERY parameter is
## positive definite, however accurately it is computed. The fixed effects are
## still identified, though, and their covariance holding that term's
## covariance at its estimate is well defined: the Hessian with those rows and
## columns removed. That is what lme4 falls back to for a GLMM whose full
## Hessian fails, and at a variance of zero it is the covariance of the model
## without the term -- which is the model the data are describing.
##
## Only covariance parameters are ever held, and only whole terms, the ones at
## the boundary first. If what is left is still singular -- separation,
## aliased columns, a dispersion that ran off -- nothing is held, and the fit
## stays a failure. glmmTMB accepts any Hessian whose smallest eigenvalue beats
## machine epsilon, which lets a boundary variance through with a standard
## error in the billions; a boundary is treated as one here instead.

## Hessian of `fn` from its exact gradient: central differences at steps s and
## s/2, combined so the leading error term cancels.
#' @keywords internal
#' @noRd
ilm_hessian <- function(gr, par, h = 1e-4) {
  n <- length(par); H <- matrix(NA_real_, n, n)
  for (i in seq_len(n)) {
    s <- h * max(abs(par[i]), 1)
    e <- numeric(n); e[i] <- 1
    d1 <- (gr(par + s * e) - gr(par - s * e)) / (2 * s)
    d2 <- (gr(par + s / 2 * e) - gr(par - s / 2 * e)) / s
    H[, i] <- (4 * d2 - d1) / 3
  }
  (H + t(H)) / 2
}

## The covariance parameters, by term, as positions in the optimised vector;
## which of those terms sit at the boundary; and which positions may never be
## held -- the fixed effects, the dispersion, thresholds and zero part, which
## are the model a user reads rather than its variance structure.
#' @keywords internal
#' @noRd
ilm_cov_blocks <- function(re, Sig, Sigd, ty, rk, dk, toff, ar, pe, pn) {
  it <- which(pn == "theta")
  blocks <- list(); flagged <- character(0)
  for (k in seq_along(re)) {
    nm <- names(re)[k]
    blocks[[nm]] <- it[(toff[k] + 1L):toff[k + 1L]]
    S <- Sig[[nm]]; sdv <- sqrt(pmax(diag(S), 0))
    ev <- sort(eigen(S, symmetric = TRUE, only.values = TRUE)$values,
               decreasing = TRUE)
    kk <- if (ty[k] == "rr") seq_len(rk[k]) else seq_along(ev)
    bad <- any(sdv < 1e-3) ||
      min(ev[kk]) < 1e-6 * max(ev[1], .Machine$double.eps)
    if (!bad && dk[k] > 1L && !is.null(Sigd[[nm]])) {
      Sd <- Sigd[[nm]]; sv <- sqrt(diag(Sd)); R <- Sd / outer(sv, sv)
      bad <- any(sv[-1] / sv[1] < 1e-3) || max(abs(R[upper.tri(R)])) > 0.999
    }
    if (bad) flagged <- c(flagged, nm)
  }
  if (!is.null(ar)) {
    blocks[["ar"]] <- which(pn %in% c("lchol_ar", "rho_raw"))
    rr <- pe[pn == "rho_raw"]
    car <- identical(ar$type, "car1")
    rho <- if (car) exp(-1 / exp(rr)) else tanh(rr)
    eff <- if (car) rho^stats::median(ar$gap) else rho
    Sa <- Sig[["ar"]]
    eva <- eigen(Sa, symmetric = TRUE, only.values = TRUE)$values
    if (abs(eff) > 0.99 || any(sqrt(pmax(diag(Sa), 0)) < 1e-3) ||
        min(eva) < 1e-6 * max(eva, .Machine$double.eps))
      flagged <- c(flagged, "ar")
  }
  list(blocks = blocks, flagged = flagged,
       keep = which(!pn %in% c("theta", "lchol_ar", "rho_raw")))
}

## The three outcomes: "recomputed" (the accurate Hessian is positive definite
## and nothing is at a boundary, so the fit is fully usable), "boundary" (the
## terms in `held` are held at their estimates and everything else is
## usable), or "none" (still a failure). A fit TMB was already happy with is
## "tmb". The sdreport comes back redone from the Hessian that was accepted, so
## everything derived from it -- thresholds, rho, the joint precision REML and
## smooths need -- is consistent with it.
#' @keywords internal
#' @noRd
ilm_hess_recover <- function(obj, opt, sdr, cb, joint) {
  if (isTRUE(sdr$pdHess)) return(list(sdr = sdr, how = "tmb", held = character(0)))
  out <- list(sdr = sdr, how = "none", held = character(0))
  if (!length(opt$par)) return(out)
  H <- tryCatch(ilm_hessian(function(p) as.numeric(obj$gr(p)), opt$par),
                error = function(e) NULL)
  ## those gradient calls moved the tape; put it back at the optimum
  invisible(tryCatch(obj$fn(opt$par), error = function(e) NULL))
  if (is.null(H) || !all(is.finite(H))) return(out)
  ## Positive definite on a scale-free test, not merely chol()-able. Finite
  ## differences leave an exactly singular Hessian slightly positive, so an
  ## aliased pair of columns (x2 = 2 * x1) passed chol() and was "rescued";
  ## glmmTMB's threshold -- the smallest eigenvalue above machine epsilon --
  ## lets the same case through. Rescaled to a unit diagonal, differencing
  ## noise sits near 1e-8 and a genuine dependence at the aliasing check's
  ## own FAIL line (|correlation| 0.995) sits near 5e-3; 1e-6 separates them.
  pd <- function(M) {
    dg <- diag(M)
    if (!length(dg) || any(!is.finite(dg)) || any(dg <= 0)) return(FALSE)
    S <- M / sqrt(outer(dg, dg))
    ev <- tryCatch(eigen(S, symmetric = TRUE, only.values = TRUE)$values,
                   error = function(e) NA_real_)
    all(is.finite(ev)) && min(ev) > 1e-6 &&
      !inherits(try(chol(M), silent = TRUE), "try-error")
  }
  redo <- function(Hm) tryCatch(suppressWarnings(
    sdreport(obj, par.fixed = opt$par, hessian.fixed = Hm,
             getJointPrecision = joint)), error = function(e) NULL)
  ## a boundary variance is handled as one whatever the Hessian says: an
  ## inverse that is merely positive definite gives it a standard error in the
  ## billions, which is a number rather than information
  if (!length(cb$flagged) && pd(H)) {
    s2 <- redo(H)
    if (!is.null(s2) && isTRUE(s2$pdHess))
      return(list(sdr = s2, how = "recomputed", held = character(0)))
  }
  for (on in unique(Filter(length, list(cb$flagged, names(cb$blocks))))) {
    d <- unlist(cb$blocks[on], use.names = FALSE)
    if (!length(d) || any(cb$keep %in% d)) next
    kp <- setdiff(seq_len(nrow(H)), d)
    if (!length(kp) || !pd(H[kp, kp, drop = FALSE])) next
    ## held: its rows and columns cut loose, with a curvature so large that
    ## no uncertainty is carried through from it
    Hm <- H; Hm[d, ] <- 0; Hm[, d] <- 0
    Hm[cbind(d, d)] <- 1e12 * max(1, abs(diag(H)))
    s2 <- redo(Hm)
    if (is.null(s2)) next
    cf <- s2$cov.fixed; cf[d, ] <- NA_real_; cf[, d] <- NA_real_
    s2$cov.fixed <- cf
    ## the model's own Hessian was not positive definite, and that stays on
    ## the record; what is usable is said by `how`
    s2$pdHess <- FALSE
    return(list(sdr = s2, how = "boundary", held = on))
  }
  out
}

## Whether the FIXED effects of a fit carry usable standard errors: a positive
## definite Hessian, or a boundary term held at its estimate.
#' @keywords internal
#' @noRd
ilm_fixed_usable <- function(object)
  isTRUE(object$sdr$pdHess) || length(object$hessian_held) > 0L

#' Print a table of checks
#'
#' @param ck A checks data frame.
#' @param title Character heading.
#' @return `ck`, invisibly.
#' @keywords internal
#' @noRd
ilm_print_checks <- function(ck, title) {
  sym <- c(OK = "  ok  ", WARN = " WARN ", FAIL = " FAIL ", INCONCLUSIVE = "  ??  ",
           BOUNDARY = "BOUND ")
  cat("\n", title, "\n", strrep("-", 78), "\n", sep = "")
  for (i in seq_len(nrow(ck))) {
    cat(sprintf("[%s] %-24s %s\n", sym[ck$status[i]], ck$check[i], ck$detail[i]))
    if (nzchar(ck$cause[i]))      cat("        why: ", ck$cause[i], "\n", sep = "")
    if (nzchar(ck$suggestion[i])) cat("        try: ", ck$suggestion[i], "\n", sep = "")
  }
  invisible(ck)
}

#' Fit a multinomial mixed model from a design matrix
#'
#' The computational engine. Most users should call [ilm_model()] with a formula
#' instead; this is the entry point when you already have a design matrix, and it
#' is what the formula interface calls internally.
#'
#' @section How the model is parameterised:
#' With `J` outcome categories there are `C = J - 1` free linear predictors.
#' Coefficients use **sum-to-zero** coding, so each one is that category's
#' deviation from the average across categories, not a contrast against a
#' baseline. This differs from [nnet::multinom()] and `brms`.
#'
#' Random effects are held as a matrix rather than a long vector, with a
#' *matrix-normal* prior. Writing it this way means the covariance over groups
#' and the covariance over categories never have to be combined into one large
#' Kronecker product, which keeps both memory and computation manageable and
#' lets AR(1), spatial and i.i.d. structures share a single code path.
#'
#' @section Choosing a category covariance structure:
#' Each random term gets its own structure through `re_struct`:
#' \describe{
#'   \item{`"us"`}{unstructured: every variance and correlation free. Costs
#'     `C(C+1)/2` parameters, which grows quickly -- 45 at `J = 10`.}
#'   \item{`"diag"`}{diagonal: categories uncorrelated. Costs `C`.}
#'   \item{`"rr"` with `rank = r`}{reduced rank: `Sigma = Lambda Lambda'` with
#'     `Lambda` of size C-by-r. Costs fewer parameters **and** fewer latent
#'     values, which is usually the better trade when a term is stretched.}
#' }
#' If a structure is too rich for the data the pre-fit checks will say so and
#' name a specific rank to try.
#'
#' @section Models with no random effects:
#' `re_list` may be empty. The latent parameter block is then omitted entirely
#' rather than declared with length zero, and for a gaussian response the fit
#' reports exact t and F inference instead of the large-sample approximations.
#' See [ilm_model()].
#'
#' @section Fitting method:
#' The random effects are integrated out by the **Laplace approximation**,
#' implemented through `RTMB`. Exact quadrature is not an option here: with
#' crossed random effects the integral does not factorise, so its dimension is
#' the total number of latent values, often in the thousands. Laplace is an
#' approximation, and its accuracy depends on having enough data per latent
#' value -- see the `latent_budget` entry in `fit$checks`, and
#' [ilm_consistency()] to test it directly on your own fit.
#'
#' @param X Numeric design matrix for the fixed effects, with `N` rows.
#' @param y Integer vector of length `N` giving the observed category, coded
#'   `1` to `J`.
#' @param J Integer. Number of outcome categories, for the multinomial family
#'   only; `NULL` for every other family.
#' @param family Response distribution: a name, or the object returned by
#'   [ilm_family()]. See [ilm_family()] for what each one assumes.
#' @param re_list Named list of random terms. Each element is a grouping vector
#'   (random intercept), a `list(group =, Z =)` (random slopes), or a
#'   `list(basis =)` (a smooth, from `ilm_smooth()`).
#' @param re_struct Optional named list of category covariance structures,
#'   parallel to `re_list`. Defaults to `"us"` for every term. Set `d_cor =
#'   FALSE` within an element to drop an intercept-slope correlation.
#' @param ar Optional AR(1) specification:
#'   `list(idx =, n_group =, Tt =)`.
#' @param censor Optional censoring specification from [ilm_censor()], marking
#'   observations known only as an interval -- at or below a floor, at or above
#'   a ceiling. Supported for the gaussian family, where it gives a Tobit model.
#' @param rp Internal. The flexible parametric baseline: knots, the derivative
#'   design and which columns of `X` hold the spline. Built by [ilm_model()]
#'   from `rp_df`.
#' @param Zd Optional design matrix for the dispersion model, one row per
#'   observation. Its columns become a linear predictor for the logarithm of
#'   the dispersion, so its intercept replaces the single dispersion parameter.
#'   Built by [ilm_model()] from `dispformula`.
#' @param disp_mu Logical. Add a term in the logarithm of the fitted mean to
#'   that predictor, making the dispersion a power of the mean.
#' @param ylevels Optional character vector of category labels, used in output.
#' @param weights Optional numeric vector of **frequency** weights: the number of
#'   replicate observations each row stands for, exactly as in a binomial `glm()`
#'   fitted to grouped data. A weighted fit is then numerically identical to
#'   expanding each row, standard errors included.
#'
#'   Valid only when the pooled observations share every random effect in the
#'   model -- you may aggregate within subject and covariate pattern, but not
#'   across subjects. These are **not** sampling weights; see the
#'   `weights_type` check and the note below.
#' @param Zzi Optional design matrix for the zero part of a count model, one
#'   row per observation, modelling the logit of an excess-zero probability.
#'   Built by [ilm_model()] from `ziformula`.
#' @param zi_type `"inflated"` for a mixture, where a zero may have come either
#'   from the zero process or from the count itself, or `"hurdle"`, where every
#'   zero comes from the zero process and the positives come from a count that
#'   cannot be zero. Ignored when `Zzi` is `NULL`.
#' @param verbose Logical. Print the checks while fitting.
#' @param restarts Integer. Number of optimiser restarts from the previous
#'   solution, which helps on difficult surfaces.
#' @param joint Logical. Also compute the joint precision over fixed and random
#'   parameters. Needed by [predict.ilm_model()] to propagate uncertainty in penalised
#'   smooth coefficients; [ilm_model()] switches it on automatically when the model
#'   contains smooths.
#' @param reml Logical. Integrate the fixed effects out along with the random
#'   ones, giving restricted maximum likelihood. Under a flat prior that
#'   integral IS the restricted likelihood, and for a linear-gaussian model the
#'   Laplace approximation to it is exact -- so this is REML rather than an
#'   approximation to it. Gaussian responses only; elsewhere the integral is
#'   still well defined but has none of REML's properties, so it is refused.
#'   See [ilm_model()] for when to switch it on.
#'
#' @return An object of class `"ilm_model"`: a list whose most useful elements are
#'   `checks` (the diagnostic table), `Sigma` (fitted category covariances),
#'   `beta` (fixed effects as a p-by-C matrix), `opt` (optimiser output), `sdr`
#'   (`TMB::sdreport()` output) and `ok` (whether every check passed).
#'
#' @section A warning about sampling weights:
#' Passing inverse-probability or survey weights here is not supported and will
#' mislead you. In a *mixed* model the weight multiplies the conditional
#' likelihood inside the Laplace integral, so each cluster appears to carry more
#' information than it does; both the standard errors and the coefficients are
#' affected. Proper handling needs scaled level-specific weights and a
#' design-based variance estimator, neither of which is implemented. The
#' `weights_type` check flags weights that look like sampling weights.
#'
#' @references
#' Kristensen, K., Nielsen, A., Berg, C. W., Skaug, H., & Bell, B. M. (2016).
#' TMB: Automatic differentiation and Laplace approximation.
#' *Journal of Statistical Software*, 70(5), 1--21.
#'
#' Agresti, A. (2013). *Categorical Data Analysis*, 3rd ed. Wiley.
#'
#' Pfeffermann, D., Skinner, C. J., Holmes, D. J., Goldstein, H., & Rasbash, J.
#' (1998). Weighting for unequal selection probabilities in multilevel models.
#' *Journal of the Royal Statistical Society, Series B*, 60(1), 23--40.
#'
#' Rabe-Hesketh, S., & Skrondal, A. (2006). Multilevel modelling of complex
#' survey data. *Journal of the Royal Statistical Society, Series A*, 169(4),
#' 805--827.
#'
#' @seealso [ilm_model()] for the formula interface, [summary.ilm_model()],
#'   [ilm_consistency()].
#' @export
ilm_fit <- function(X, y, J = NULL, re_list, re_struct = NULL, ar = NULL,
                     ylevels = NULL, weights = NULL, family = "multinomial",
                     verbose = TRUE, restarts = 3L, joint = FALSE,
                     censor = NULL, Zd = NULL, disp_mu = FALSE,
                     rp = NULL, Zzi = NULL, zi_type = c("inflated", "hurdle"),
                     reml = FALSE) {
  zi_type <- match.arg(zi_type)
  fam <- if (is.list(family)) family else ilm_family(family)
  ## Integrating the fixed effects out under a flat prior is available for any
  ## family, and glmmTMB does exactly this -- `if (REML) randomArg <-
  ## c(randomArg, "beta")`. What differs is how much it delivers. For a linear
  ## model the integral is exact and the result is textbook REML, unbiased with
  ## the n - p divisor. For anything else the Laplace approximation to that
  ## integral is itself an approximation, so what comes out is an
  ## approximately-restricted likelihood: it still removes most of the
  ## downward bias in the variance components, and it is not entitled to REML's
  ## exact properties. Recorded on the fit as `reml_exact` so nothing has to
  ## infer which of the two it is holding.
  reml_exact <- isTRUE(reml) && identical(fam$name, "gaussian")
  ## Censoring: derive the codes from THIS response, so a simulated replicate
  ## is censored by the same rule the data were rather than inheriting the
  ## observed pattern. See ilm_censor().
  if (!is.null(censor)) {
    if (!isTRUE(fam$censorable))
      stop("the ", fam$name, " family has no censored form. Censoring needs a ",
           "continuous response whose distribution function can be evaluated; ",
           "available: gaussian, and the survival families weibull, lognormal, ",
           "loglogistic and rp.", call. = FALSE)
    if (length(censor) != nrow(X))
      stop("`censor` has ", length(censor), " values but the model matrix has ",
           nrow(X), " rows", call. = FALSE)
  }
  cens <- ilm_censor_for(censor, y)
  ## An ordered response has thresholds where every other family has an
  ## intercept, and a column of ones would be confounded with all of them at
  ## once rather than with any one of them.
  has_ord <- isTRUE(fam$ordinal)
  if (has_ord) {
    if (is.null(J) || J < 3L)
      stop("an ordinal family needs at least 3 ordered categories; with 2 the ",
           "thresholds collapse to a single intercept and the model is a ",
           "binomial one -- use family = \"binomial\".", call. = FALSE)
    if ("(Intercept)" %in% colnames(X))
      stop("internal: an ordinal design still carries an intercept column",
           call. = FALSE)
  }
  ## A dispersion model replaces the single dispersion parameter with a linear
  ## predictor for its logarithm. It needs a dispersion to model in the first
  ## place, which a poisson or binomial response does not have -- their spread
  ## is fixed by their mean, and the remedy for extra spread there is a
  ## different family, not a formula.
  if (!is.null(Zd) || isTRUE(disp_mu)) {
    if (fam$n_disp < 1L)
      stop("the ", fam$name, " family has no dispersion parameter to model: ",
           "its variance is determined by its mean. For counts that are more ",
           "spread out than poisson allows, use family = \"nbinom\"; for a ",
           "binary response, extra spread has to come from a random effect.",
           call. = FALSE)
    if (!is.null(Zd) && nrow(Zd) != nrow(X))
      stop("the dispersion design has ", nrow(Zd), " rows but the model ",
           "matrix has ", nrow(X), call. = FALSE)
  }
  ## A zero part needs a response for which a zero is a distinct event that the
  ## distribution already gives a probability to, so that "more zeros than
  ## this" is a statement with content. That is counts.
  has_zi <- !is.null(Zzi)
  if (has_zi) {
    if (!isTRUE(fam$zi_ok))
      stop("`ziformula` applies to counts (poisson, nbinom), not to the ",
           fam$name, " family: a zero there is an ordinary value of a ",
           "continuous response and carries no separate probability to ",
           "inflate.", call. = FALSE)
    if (nrow(Zzi) != nrow(X))
      stop("the zero-part design has ", nrow(Zzi), " rows but the model ",
           "matrix has ", nrow(X), call. = FALSE)
    ## A mixture adds the probability that the response distribution produced
    ## a zero by itself. A continuous one never does -- its density at zero is
    ## a density, not a probability -- so "inflated" has nothing to add and
    ## every zero necessarily came from the zero process. That is a hurdle,
    ## and asking for the other thing is asking for something that does not
    ## exist rather than for a worse approximation of it.
    if (isTRUE(fam$continuous) && !identical(zi_type, "hurdle"))
      stop("the ", fam$name, " family is continuous, so a zero-inflated ",
           "MIXTURE is not defined for it: there is no probability that the ",
           "response produced a zero on its own to mix with. Every zero comes ",
           "from the zero process, which is zi_type = \"hurdle\".",
           call. = FALSE)
  }
  ## An accelerated failure time model works on log(t), so a non-positive time
  ## is not a hard case, it is a contradiction. Say so rather than returning
  ## NaN from inside the optimiser.
  if (isTRUE(fam$positive)) {
    bad <- !is.na(y) & y <= 0
    if (any(bad))
      stop("the ", fam$name, " family models log(time), so every response must ",
           "be strictly positive; ", sum(bad), " value",
           if (sum(bad) > 1L) "s are" else " is", " zero or negative",
           if (any(y[bad] == 0, na.rm = TRUE))
             ". A recorded time of exactly 0 usually means an event before the first assessment: give it the smallest time the study could have measured, or left-censor it with ilm_censor()." else "",
           call. = FALSE)
  }
  ## Accept a spec from ilm_ar1()/ilm_car1(), or the bare list the fitter took
  ## before those existed.
  ar <- ilm_as_cor(ar)
  if (!is.null(ar) && length(ar$idx) != nrow(X))
    stop("the correlation structure covers ", length(ar$idx),
         " observations but the model matrix has ", nrow(X),
         " rows. Build it from the same rows the model is fitted to.",
         call. = FALSE)
  ## J is the number of outcome categories, and is meaningful only for the
  ## multinomial; every other family uses a single linear predictor.
  if (is.null(J)) J <- 2L
  C <- fam$C_of(J); N <- nrow(X); p <- ncol(X)
  Tc <- if (fam$name == "multinomial") contr.sum(J) else matrix(1, 1, 1)
  re <- ilm_norm_re(re_list, N)
  if (is.null(re_struct)) re_struct <- lapply(re, function(z) list(type = "us"))
  re_struct <- re_struct[names(re)]
  for (nm in names(re_struct)) {
    if (is.null(re_struct[[nm]]$rank))  re_struct[[nm]]$rank  <- NA_integer_
    if (is.null(re_struct[[nm]]$d_cor)) re_struct[[nm]]$d_cor <- TRUE
    ty_nm <- re_struct[[nm]]$type
    if (is.null(ty_nm) || !ty_nm %in% c("us", "diag", "rr"))
      stop(sprintf("re_struct$%s$type must be one of \"us\", \"diag\" or \"rr\"", nm),
           call. = FALSE)
    ## rank = 2 is a double, which is what anyone would actually type, but the
    ## width and parameter counts below are gathered with vapply(..., 1L) and
    ## fail on a double.  Coerce once, here, having checked it is sensible.
    r <- re_struct[[nm]]$rank
    if (ty_nm == "rr") {
      if (is.na(r))
        stop(sprintf("re_struct$%s: a \"rr\" structure needs a rank, e.g. list(type = \"rr\", rank = 2)", nm),
             call. = FALSE)
      if (!is.numeric(r) || length(r) != 1L || !is.finite(r) || r != round(r) || r < 1)
        stop(sprintf("re_struct$%s$rank must be a single whole number >= 1", nm),
             call. = FALSE)
      if (r > C)
        stop(sprintf("re_struct$%s$rank is %d, but there are only %d category dimensions; use type = \"us\" for a full-rank covariance",
                     nm, as.integer(r), C), call. = FALSE)
      re_struct[[nm]]$rank <- as.integer(r)
    } else re_struct[[nm]]$rank <- NA_integer_
  }

  if (is.null(weights)) weights <- rep(1, N)
  if (length(weights) != N) stop("weights must have one entry per row of X")
  if (any(weights < 0)) stop("weights must be non-negative")
  pre <- ilm_precheck(y, J, re, re_struct, ar, weights, fam)
  if (verbose) ilm_print_checks(pre, "pre-fit data checks")
  if (any(pre$status == "FAIL") && verbose)
    cat("\n>> pre-fit checks FAILED; fitting anyway, but treat the result as unreliable.\n")

  ## For the multinomial the response becomes a COUNT matrix, with wrow scaling
  ## the log-normaliser.  Every other family takes the response as supplied.
  if (fam$name == "multinomial") {
    yobs <- matrix(0, N, J); yobs[cbind(seq_len(N), y)] <- weights
  } else {
    ilm_check_response(y, fam, has_zero_part = has_zi)
    yobs <- as.numeric(y)
  }
  ## For an ordered response the likelihood needs, per row, which threshold
  ## bounds it above and which below, and whether either of those is actually
  ## a limit rather than a cut. All of it comes from the response, so it is
  ## built once here rather than inside the tape.
  ord_idx <- if (has_ord) list(
    iu = pmin(as.integer(y), J - 1L), il = pmax(as.integer(y) - 1L, 1L),
    mu = as.numeric(as.integer(y) < J), ml = as.numeric(as.integer(y) > 1L))
    else NULL
  K   <- length(re)
  nlk <- vapply(re, `[[`, 1L, "nl"); dk <- vapply(re, `[[`, 1L, "d")
  wk  <- vapply(names(re), function(nm) ilm_str_width(re_struct[[nm]], C), 1L)
  npc <- vapply(names(re), function(nm) ilm_str_npar(re_struct[[nm]], C), 1L)
  dcor <- vapply(re_struct, function(s) isTRUE(s$d_cor), TRUE)
  npd  <- as.integer(mapply(ilm_npar_d, dk, dcor))
  ty  <- vapply(re_struct, `[[`, "", "type"); rk <- vapply(re_struct, `[[`, 1L, "rank")
  bl  <- as.integer(nlk * dk * wk); boff <- c(0L, cumsum(bl))
  tl  <- as.integer(npc + npd);     toff <- c(0L, cumsum(tl))

  ## names in obj$par order: beta (column-major over categories), theta, [ar]
  xn <- colnames(X); if (is.null(xn)) xn <- paste0("X", seq_len(p))
  if (is.null(ylevels))
    ylevels <- if (fam$name == "multinomial") paste0("cat", seq_len(J)) else ""
  ## a univariate family has one unnamed linear predictor, so its coefficient
  ## names are simply the predictor names
  cpre <- if (fam$name == "multinomial") paste0(ylevels[seq_len(C)], ":") else rep("", C)
  pnames <- c(
    unlist(lapply(seq_len(C), function(cc) paste0(cpre[cc], xn))),
    unlist(lapply(seq_len(K), function(k) {
      nm <- names(re)[k]
      a <- switch(ty[k], us = ilm_nm_tri(nm, C), diag = ilm_nm_diag(nm, C),
                  rr = ilm_nm_lam(nm, C, rk[k]))
      if (dk[k] > 1L) a <- c(a, ilm_nm_ld(nm, dk[k], dcor[k]))
      a
    })))
  has_dm <- !is.null(Zd) || isTRUE(disp_mu)
  if (fam$n_disp > 0L && !has_dm) pnames <- c(pnames, fam$disp_names)
  if (has_dm) {
    if (!is.null(Zd)) pnames <- c(pnames, paste0("disp:", colnames(Zd)))
    if (isTRUE(disp_mu)) pnames <- c(pnames, "disp:mu_power")
    ## a second dispersion parameter (nbinom has none, but keep the block
    ## general) is still a plain scalar
    if (fam$n_disp > 1L)
      pnames <- c(pnames, fam$disp_names[-1L])
  }
  if (has_ord) pnames <- c(pnames, paste0("zeta:", ilm_ord_cut_names(ylevels, J)))
  if (has_zi) pnames <- c(pnames, paste0("zi:", colnames(Zzi)))
  if (!is.null(ar)) pnames <- c(pnames, ilm_nm_tri("ar", C), "ar:rho_raw")

  dl <- list(X = X, yobs = yobs, wrow = weights, Tct = t(Tc),
             n_disp = fam$n_disp,
             ## the derivative of the linear predictor with respect to log
             ## time, which a flexible parametric density needs and nothing
             ## else does; zero columns everywhere but the spline
             Drp = if (is.null(rp)) matrix(0, nrow(X), 0L) else rp$D,
             has_rp = !is.null(rp),
             Zdisp = if (is.null(Zd)) matrix(0, nrow(X), 0L) else Zd,
             has_dm = has_dm, disp_mu = isTRUE(disp_mu),
             Zzi = if (is.null(Zzi)) matrix(0, nrow(X), 0L) else Zzi,
             has_ord = has_ord,
             ## rows of cumulative sums: row j adds the first j - 1 gaps, so
             ## the first threshold gets none and the ordering is structural
             ord_L = if (has_ord)
               outer(seq_len(J - 1L), seq_len(J - 2L), `>`) * 1 else
               matrix(0, 0, 0),
             ord_iu = if (has_ord) ord_idx$iu else integer(0),
             ord_il = if (has_ord) ord_idx$il else integer(0),
             ord_mu = if (has_ord) ord_idx$mu else numeric(0),
             ord_ml = if (has_ord) ord_idx$ml else numeric(0),
             has_zi = has_zi, zi_hurdle = identical(zi_type, "hurdle"),
             zi_cont = isTRUE(fam$continuous),
             ## split once, outside the likelihood: the zero rows and the
             ## positive rows take different terms, and indexing them keeps the
             ## hurdle's truncation factor from ever being evaluated on a row
             ## where it is undefined
             i_zero = if (has_zi) which(as.numeric(yobs) == 0) else integer(0),
             i_pos  = if (has_zi) which(as.numeric(yobs) >  0) else integer(0),
             yzero  = if (has_zi) rep(0, nrow(X)) else numeric(0),
             grp = lapply(re, `[[`, "group"), Zl = lapply(re, `[[`, "Z"),
             kind = vapply(re, `[[`, "", "kind"), bas = lapply(re, `[[`, "basis"),
             b_idx = lapply(seq_len(K), function(k) (boff[k] + 1L):boff[k + 1L]),
             t_idx = lapply(seq_len(K), function(k) (toff[k] + 1L):toff[k + 1L]),
             nlk = as.integer(nlk), dk = as.integer(dk), wk = as.integer(wk),
             npc = as.integer(npc), ty = as.character(ty), rk = as.integer(rk),
             dcor = as.logical(dcor), K = K, C = C, has_ar = !is.null(ar),
             is_car = isTRUE(ar$type == "car1"))
  if (!is.null(ar)) {
    if (identical(ar$type, "car1")) {
      ## cells are ordered by group then time, so the preceding observation of
      ## a group is always the cell before it -- the same fact the evenly
      ## spaced case relies on, and why the index arithmetic is shared
      dl$ar_idx <- ar$idx; dl$idx1 <- ar$first
      dl$idx_t <- ar$rest; dl$idx_lag <- ar$prev
      dl$ar_gap <- as.numeric(ar$gap)
      dl$n_g <- ar$n_group; dl$Tt <- NA_integer_
    } else {
      i1 <- (seq_len(ar$n_group) - 1L) * ar$Tt + 1L
      dl$ar_idx <- ar$idx; dl$idx1 <- i1
      dl$idx_t <- setdiff(seq_len(ar$n_cell), i1); dl$idx_lag <- dl$idx_t - 1L
      dl$ar_gap <- numeric(0)
      dl$n_g <- ar$n_group; dl$Tt <- ar$Tt
    }
    dl$nre_ar <- length(dl$idx_t)
  }

  fam_nll <- fam$nll                    # captured in the closure, not via dl
  fam_logden <- fam$logden
  fam_linkinv <- fam$linkinv
  f <- function(pars) {
    getAll(pars, dl)
    nll <- 0; eta <- X %*% beta; sdv <- AD(numeric(0))
    for (k in seq_len(K)) {
      d <- dk[k]; nl <- nlk[k]; w <- wk[k]
      blk <- matrix(bvec[b_idx[[k]]], nl * d, w)   # rows: dimension-major
      th  <- theta[t_idx[[k]]]
      rows <- function(i) blk[((i - 1L) * nl + 1L):(i * nl), , drop = FALSE]
      if (ty[k] == "rr") {                         # columns are already N(0,1)
        Lam <- ilm_mkLam(th[1:npc[k]], C, rk[k]); lcld <- 0
        M <- lapply(1:d, rows)
        S <- Lam %*% t(Lam)
      } else {
        L <- if (ty[k] == "diag") ilm_mkD(th[1:npc[k]], C) else ilm_mkL(th[1:npc[k]], C)
        Lci <- solve(t(L)); lcld <- sum(log(diag(L)))
        M <- lapply(1:d, function(i) rows(i) %*% Lci)
        S <- L %*% t(L)
      }
      ## quadratic form.  With row covariance I_nl (x) Sigma_d and column
      ## covariance Sigma_c, tr(Sc^-1 B' (I (x) Sd^-1) B) collapses to
      ## sum_ij Sd^-1[i,j] * sum(M_i * M_j) -- no Kronecker is ever formed.
      if (d == 1L) {
        nll <- nll + 0.5 * sum(M[[1]] * M[[1]])
      } else {
        Ld <- ilm_mkLd(th[(npc[k] + 1L):(npc[k] + ilm_npar_d(d, dcor[k]))], d, dcor[k])
        ## sum_ij Sd^-1[i,j] <M_i, M_j> = sum_k || sum_i A[k,i] M_i ||^2 with
        ## A = Ld^-1.  A triangular solve, so Sigma_d is never inverted -- which
        ## matters because the intercept-slope correlation can reach +/-1.
        A <- solve(Ld)
        for (kk in 1:d) {
          Nk <- A[kk, 1] * M[[1]]
          if (kk > 1L) for (i in 2:kk) Nk <- Nk + A[kk, i] * M[[i]]
          nll <- nll + 0.5 * sum(Nk * Nk)
        }
        nll <- nll + w * nl * sum(log(diag(Ld)))
      }
      nll <- nll + nl * d * lcld
      if (kind[k] == "basis") {
        ctb <- bas[[k]] %*% rows(1L)             # dense basis, no grouping index
        if (ty[k] == "rr") ctb <- ctb %*% t(Lam)
        eta <- eta + ctb
      } else {
        Zk <- Zl[[k]]; gk <- grp[[k]]
        for (i in 1:d) {
          ctb <- rows(i)[gk, , drop = FALSE]
          if (ty[k] == "rr") ctb <- ctb %*% t(Lam)
          eta <- eta + Zk[, i] * ctb
        }
      }
      sdv <- c(sdv, sqrt(diag(S)))
    }
    if (has_ar) {
      La <- ilm_mkL(lchol_ar, C)
      ## the first observation of each group, marginally
      W1 <- B_ar[idx1, , drop = FALSE] %*% solve(t(La))
      nll <- nll + 0.5 * sum(W1 * W1) + n_g * sum(log(diag(La)))
      ## Every transition after that is Markov.  The only difference between
      ## the two structures is whether the correlation from one observation to
      ## the next is a single number or depends on the gap.
      if (is_car) {
        ## CAR(1): phi_k = rho ^ d_k, innovation variance 1 - phi_k^2.  Exact
        ## for the Ornstein-Uhlenbeck process, so an irregular gap costs
        ## nothing in accuracy.
        ##
        ## Parameterised on log(range) rather than on rho.  rho is the
        ## correlation one TIME UNIT apart, so with fine units (seconds, days)
        ## it sits at 0.9999-something and the gradient in it is negligible;
        ## the range is in the units of time and just shifts when they change.
        rng <- exp(rho_raw)
        phi <- exp(-ar_gap / rng)
        rho <- exp(-1 / rng)
        rsd <- B_ar[idx_t, , drop = FALSE] - phi * B_ar[idx_lag, , drop = FALSE]
        W <- (rsd / sqrt(1 - phi^2)) %*% solve(t(La))
        nll <- nll + 0.5 * sum(W * W) + nre_ar * sum(log(diag(La))) +
               (C / 2) * sum(log(1 - phi^2))
      } else {
        rho <- tanh(rho_raw)
        rsd <- B_ar[idx_t, , drop = FALSE] - rho * B_ar[idx_lag, , drop = FALSE]
        W <- (rsd %*% solve(t(La))) / sqrt(1 - rho^2)
        nll <- nll + 0.5 * sum(W * W) +
               nre_ar * (sum(log(diag(La))) + (C / 2) * log(1 - rho^2))
      }
      eta <- eta + B_ar[ar_idx, , drop = FALSE]
      Sa <- La %*% t(La); sdv <- c(sdv, sqrt(diag(Sa)))
      ADREPORT(rho)
    }
    ## The only family-specific line in the whole likelihood.  Everything above
    ## builds eta and accumulates the random-effect prior, and none of it
    ## depends on the response distribution.
    ## logdisp exists only when there is no dispersion model, or when the
    ## family has a second dispersion parameter the model does not cover
    dsp <- if (n_disp > 0L && !has_dm) logdisp else numeric(0)
    ## A dispersion model gives every row its own log dispersion. The reserved
    ## `mu` term makes it a power of the fitted mean, which is the remedy for
    ## the pattern ilm_check_variance() calls a trend; everything else in the
    ## dispersion design is ordinary data.
    lsig <- NULL
    if (has_dm) {
      lsig <- if (ncol(Zdisp) > 0L) as.vector(Zdisp %*% gamma) else 0 * eta[, 1]
      if (disp_mu) {
        mu1 <- fam_linkinv(eta[, 1])
        lsig <- lsig + mu_pow * log(abs(mu1) + 1e-8)
      }
      dsp <- if (n_disp > 1L) c(0, logdisp) else numeric(0)
    }
    etad <- if (has_rp) as.vector(Drp %*% beta[, 1]) else NULL
    ## The thresholds are a first value and a set of log increments, so they
    ## come out ordered whatever the optimiser does with them.
    zt <- NULL
    if (has_ord) {
      zt <- zeta_raw[1] + ord_L %*% exp(zeta_raw[-1])
      ADREPORT(zt)
    }
    if (has_zi) {
      ## Second linear predictor, for the logit of the excess-zero
      ## probability. Everything below is on the log scale through
      ## logspace_add and logspace_sub, which compute log(exp(a) + exp(b)) and
      ## log(exp(a) - exp(b)) without forming either exponential -- the
      ## mixture weights and the truncation factor both underflow otherwise.
      lpz <- as.vector(Zzi %*% gzi)
      l1p <- logspace_add(0 * lpz, lpz)        # log(1 + exp(lpz))
      lden  <- fam_logden(eta, yobs,  dsp, logsig = lsig)
      lden0 <- fam_logden(eta, yzero, dsp, logsig = lsig)
      if (zi_hurdle) {
        ## a zero is a zero and nothing else, and a positive value comes from a
        ## distribution that cannot produce one. For a CONTINUOUS response
        ## that second part needs no rescaling: the density already puts no
        ## mass at zero, so 1 - f(0) is 1 and dividing by it would be dividing
        ## by one minus a density, which is not a probability at all.
        nll <- nll - sum(wrow[i_zero] * (lpz[i_zero] - l1p[i_zero]))
        nll <- nll - sum(wrow[i_pos] *
          (lden[i_pos] - l1p[i_pos] -
             (if (zi_cont) 0 else
                logspace_sub(0 * lden0[i_pos], lden0[i_pos]))))
      } else {
        ## a zero has two possible origins and the likelihood adds them
        nll <- nll - sum(wrow[i_zero] *
          (logspace_add(lpz[i_zero], lden0[i_zero]) - l1p[i_zero]))
        nll <- nll - sum(wrow[i_pos] * (lden[i_pos] - l1p[i_pos]))
      }
    } else
    nll <- nll + fam_nll(eta,
                         if (has_ord) structure(yobs, ord_idx = list(
                           iu = ord_iu, il = ord_il,
                           mu = ord_mu, ml = ord_ml)) else yobs,
                         wrow, dsp, Tct = Tct, cens = cens,
                         logsig = lsig, etad = etad, zeta = zt)
    sd_ <- sdv; ADREPORT(sd_)
    if (n_disp > 0L && !has_dm) { disp_ <- exp(logdisp); ADREPORT(disp_) }
    nll
  }

  pars <- list(beta = matrix(0, p, C))
  ## With no random or smooth terms there is nothing to integrate out, so the
  ## latent parameter block is omitted entirely rather than declared with length
  ## zero -- MakeADFun cannot make an empty vector random.
  if (sum(tl) > 0L) pars$theta <- rep(0, sum(tl))
  if (sum(bl) > 0L) pars$bvec  <- rep(0, sum(bl))
  if (fam$n_disp > 0L) {
    ## with a dispersion model the intercept of that model IS the log
    ## dispersion, so only the extra parameters beyond the first are kept
    nd <- if (has_dm) fam$n_disp - 1L else fam$n_disp
    if (nd > 0L) pars$logdisp <- rep(0, nd)
  }
  ## A flexible baseline starts with the log cumulative hazard rising at rate
  ## one in log time, which is the exponential special case. From beta = 0 the
  ## derivative is zero everywhere, which is the boundary of the valid region
  ## and a poor place to differentiate.
  if (!is.null(rp)) pars$beta[rp$cols[1L], 1L] <- 1
  if (has_dm) {
    if (!is.null(Zd)) pars$gamma <- rep(0, ncol(Zd))
    if (isTRUE(disp_mu)) pars$mu_pow <- 0

    ## WARM START. A dispersion model is a refinement of the ordinary fit, so
    ## start it from there. It matters most for the `mu` term: from beta = 0
    ## the fitted mean is zero for every row, log|mu| is the log of the
    ## numerical floor, and the power multiplying it can go anywhere. Measured
    ## on data with a true power of 1.0, a cold start converged falsely with
    ## the power at -115818 and the slope collapsed to zero.
    base <- try(suppressWarnings(
      ilm_fit(X, y, J, re_list, re_struct, ar, ylevels = ylevels,
              weights = weights, family = fam, verbose = FALSE,
              restarts = 1L, censor = censor)), silent = TRUE)
    if (!inherits(base, "try-error")) {
      bp <- base$opt$par; bn <- names(bp)
      pars$beta <- matrix(bp[bn == "beta"], p, C)
      if (sum(tl) > 0L && sum(bn == "theta") == sum(tl))
        pars$theta <- unname(bp[bn == "theta"])
      br <- base$sdr$par.random
      if (!is.null(br) && sum(bl) > 0L && sum(names(br) == "bvec") == sum(bl))
        pars$bvec <- unname(br[names(br) == "bvec"])
      ## the intercept of the dispersion model is the log dispersion the
      ## ordinary fit found
      if (!is.null(Zd) && "(Intercept)" %in% colnames(Zd) &&
          any(bn == "logdisp"))
        pars$gamma[match("(Intercept)", colnames(Zd))] <-
          unname(bp[bn == "logdisp"])[1L]
    }
    ## A power of the mean is a power of |mu|, which has no useful derivative
    ## where mu passes through zero.
    if (isTRUE(disp_mu)) {
      mu0 <- fam$linkinv(as.vector(X %*% pars$beta[, 1]))
      if (any(mu0 > 0) && any(mu0 < 0))
        warning("`dispformula = ~ mu` makes the spread a power of the fitted ",
                "mean, and the fitted mean changes sign across these data, so ",
                "the term is not well behaved near the crossing. Model the ",
                "spread on a covariate instead, or shift the response so the ",
                "mean stays on one side of zero.", call. = FALSE)
    }
  }
  if (has_ord) {
    ## Start the thresholds at the cuts of the OBSERVED category proportions,
    ## which is the exact answer for a model with no predictors, and so a good
    ## place to begin for one with them.
    tab <- tapply(weights, factor(y, levels = seq_len(J)), sum)
    tab[is.na(tab)] <- 0
    cp <- cumsum(tab / sum(tab))[-J]
    cp <- pmin(pmax(cp, 1 / (2 * N)), 1 - 1 / (2 * N))
    th <- fam$qfun(cp)
    d <- pmax(diff(th), 1e-3)
    pars$zeta_raw <- c(th[1], log(d))
  }
  if (has_zi) {
    pars$gzi <- rep(0, ncol(Zzi))
    ## Start the zero part where the data already say it is. For a hurdle with
    ## a constant probability the observed proportion of zeros IS the maximum
    ## likelihood estimate of that probability, so this starts at the answer.
    ## For a mixture it is an upper bound -- some of those zeros came from the
    ## count part -- so start at half of it rather than above the truth.
    if ("(Intercept)" %in% colnames(Zzi)) {
      p0 <- mean(as.numeric(yobs) == 0, na.rm = TRUE)
      if (!identical(zi_type, "hurdle")) p0 <- p0 / 2
      p0 <- min(max(p0, 0.01), 0.95)
      pars$gzi[match("(Intercept)", colnames(Zzi))] <- log(p0 / (1 - p0))
    }
    ## and start the count part from the ordinary fit, for the same reason the
    ## dispersion model does: from beta = 0 every fitted mean is 1, and with a
    ## zero part free to explain the zeros the two can trade against each other
    ## into a flat region of the surface
    bz <- try(suppressWarnings(
      ilm_fit(X, y, J, re_list, re_struct, ar, ylevels = ylevels,
              weights = weights, family = fam, verbose = FALSE,
              restarts = 1L, censor = censor, Zd = Zd,
              disp_mu = disp_mu)), silent = TRUE)
    if (!inherits(bz, "try-error")) {
      bp <- bz$opt$par; bn <- names(bp)
      pars$beta <- matrix(bp[bn == "beta"], p, C)
      if (sum(tl) > 0L && sum(bn == "theta") == sum(tl))
        pars$theta <- unname(bp[bn == "theta"])
      if (!is.null(pars$logdisp) && sum(bn == "logdisp") == length(pars$logdisp))
        pars$logdisp <- unname(bp[bn == "logdisp"])
      if (!is.null(pars$gamma) && sum(bn == "gamma") == length(pars$gamma))
        pars$gamma <- unname(bp[bn == "gamma"])
      br <- bz$sdr$par.random
      if (!is.null(br) && sum(bl) > 0L && sum(names(br) == "bvec") == sum(bl))
        pars$bvec <- unname(br[names(br) == "bvec"])
    }
  }
  rnd <- if (sum(bl) > 0L) "bvec" else character(0)
  if (!is.null(ar)) {
    pars$lchol_ar <- rep(0, ilm_ncov(C))
    ## for CAR(1) rho_raw is log(range); starting at the median gap puts the
    ## correlation between consecutive observations near exp(-1)
    pars$rho_raw <- if (identical(ar$type, "car1"))
      log(stats::median(ar$gap)) else 0.5
    pars$B_ar <- matrix(0, ar$n_cell, C); rnd <- c(rnd, "B_ar")
  }
  t0 <- proc.time()[3]
  ## REML integrates the FIXED effects out along with the random ones. Under a
  ## flat prior that integral is the restricted likelihood, and for a
  ## linear-gaussian model the Laplace approximation to it is exact -- so
  ## adding "beta" to `random` is not an approximation to REML, it IS REML.
  ## Checked against lme4: variance components agree to 6e-08.
  rnd_fit <- if (reml) c(rnd, "beta") else rnd
  obj <- MakeADFun(f, pars, random = if (length(rnd_fit)) rnd_fit else NULL,
                   silent = TRUE)
  ## guard: if a fill order in nm_* ever drifts from ilm_mkL/ilm_mkD/ilm_mkLam/ilm_mkLd, the
  ## names would silently mislabel every parameter.  Fail loudly instead.
  n_expected <- if (reml) length(pnames) - p * C else length(pnames)
  if (n_expected != length(obj$par))
    stop(sprintf("internal: %d parameter names for %d parameters",
                 n_expected, length(obj$par)))
  ctl <- list(iter.max = 3000, eval.max = 3000)
  opt <- nlminb(obj$par, obj$fn, obj$gr, control = ctl)
  for (k in seq_len(restarts))
    opt <- nlminb(opt$par, obj$fn, obj$gr, control = ctl)

  ## the fitted covariance structures, which depend on the estimates alone
  pn <- names(obj$par)
  structs <- function(pe) {
    thv <- pe[pn == "theta"]
    Lams <- vector("list", K); names(Lams) <- names(re)
    Sig <- lapply(seq_len(K), function(k) {
      th <- thv[(toff[k] + 1L):toff[k + 1L]][1:npc[k]]
      if (ty[k] == "rr") { L <- ilm_mkLam_num(th, C, rk[k]); Lams[[k]] <<- L; L %*% t(L) }
      else { L <- if (ty[k] == "diag") ilm_mkD_num(th, C) else ilm_mkL_num(th, C); L %*% t(L) }
    })
    names(Sig) <- names(re)
    Sigd <- list()
    for (k in seq_len(K)) if (dk[k] > 1L) {
      th <- thv[(toff[k] + 1L):toff[k + 1L]][(npc[k] + 1L):(npc[k] + ilm_npar_d(dk[k], dcor[k]))]
      Ld <- ilm_mkLd_num(th, dk[k], dcor[k]); Sigd[[names(re)[k]]] <- Ld %*% t(Ld)
    }
    if (!is.null(ar)) { La <- ilm_mkL_num(pe[pn == "lchol_ar"], C); Sig[["ar"]] <- La %*% t(La) }
    list(Sig = Sig, Sigd = Sigd, Lams = Lams)
  }
  s <- structs(opt$par)
  cb <- ilm_cov_blocks(re, s$Sig, s$Sigd, ty, rk, dk, toff, ar, opt$par, pn)

  ## AT A BOUNDARY the optimiser chases a log standard deviation towards minus
  ## infinity along a ridge the likelihood is flat on, and nlminb stops with
  ## "false" or "singular convergence" however well everything else is
  ## determined. On the messy-data regime with a true SD of 0.05 that was 29
  ## fits in 400, whose intervals covered at 0.957 all the same. lme4 avoids
  ## it by constraining a variance at zero; the equivalent here is to hold
  ## the boundary term where it got to -- equal lower and upper bounds, so the
  ## parameter vector keeps its shape -- and let the rest finish converging.
  ## The likelihood is flat in the held direction, so nothing of substance
  ## moves, and nothing is kept unless the objective is at least as good.
  if (opt$convergence != 0L && length(cb$flagged)) {
    d <- unlist(cb$blocks[cb$flagged], use.names = FALSE)
    if (length(d) && !any(cb$keep %in% d)) {
      lo <- rep(-Inf, length(opt$par)); up <- rep(Inf, length(opt$par))
      lo[d] <- up[d] <- opt$par[d]
      o2 <- tryCatch(nlminb(opt$par, obj$fn, obj$gr, lower = lo, upper = up,
                            control = ctl), error = function(e) NULL)
      if (!is.null(o2) && is.finite(o2$objective) &&
          o2$objective <= opt$objective + 1e-8) {
        o2$message <- paste0(o2$message, "; the covariance of ",
                             paste(cb$flagged, collapse = ", "),
                             " held at its boundary estimate")
        opt <- o2
        s <- structs(opt$par)
        cb <- ilm_cov_blocks(re, s$Sig, s$Sigd, ty, rk, dk, toff, ar,
                             opt$par, pn)
      }
      ## the tape last evaluated wherever nlminb stopped; put it at the optimum
      invisible(tryCatch(obj$fn(opt$par), error = function(e) NULL))
    }
  }
  Sig <- s$Sig; Sigd <- s$Sigd; Lams <- s$Lams
  pe <- opt$par; thv <- pe[pn == "theta"]

  ## joint = TRUE also returns the joint precision over (fixed, random), which
  ## is what lets predict() propagate uncertainty in the PENALISED SMOOTH
  ## coefficients instead of holding them at their conditional modes.
  ## Under REML it is not optional: beta now lives in the random block, so the
  ## only route to its covariance is through the joint precision.
  sdr <- suppressWarnings(sdreport(obj, par.fixed = opt$par,
                                   getJointPrecision = joint || reml))

  ## When the Hessian says no: recompute it more accurately, and hold a term
  ## that sits at its boundary rather than lose the fixed effects along with
  ## it (see ilm_hess_recover()). This has to come before the REML block
  ## below, which reads the joint precision out of sdr.
  hess <- ilm_hess_recover(obj, opt, sdr, cb, joint || reml)
  sdr <- hess$sdr

  ## Under REML beta is in the random block, so pull it and its covariance out
  ## of the joint precision. Nothing is reshaped yet: the checks below and the
  ## covariance derivations all run against obj, which still expects the native
  ## REML parameter vector.
  obj_ml <- NULL; reml_beta <- NULL; reml_Vb <- NULL
  if (reml) {
    jp <- sdr$jointPrecision
    rn <- rownames(jp)
    ibx <- which(rn == "beta")
    iux <- setdiff(which(rn %in% unique(rnd)), ibx)
    ## marginal covariance of beta: Schur-complement the other random effects
    ## out of the joint precision
    Hbb <- as.matrix(jp[ibx, ibx, drop = FALSE])
    reml_Vb <- if (length(iux)) {
      Hbu <- as.matrix(jp[ibx, iux, drop = FALSE])
      Huu <- as.matrix(jp[iux, iux, drop = FALSE])
      solve(Hbb - Hbu %*% solve(Huu, t(Hbu)))
    } else solve(Hbb)
    srr <- summary(sdr, "random")
    reml_beta <- srr[rownames(srr) == "beta", 1]
    ## An ML-SHAPED objective, built once and never optimised. V_beta(theta) is
    ## a function of theta and the data and does not care how theta was
    ## estimated, so evaluating this at the REML estimates gives exactly the
    ## derivatives ilm_denom_df() needs. Costs one extra tape, about 20 ms.
    obj_ml <- MakeADFun(f, pars, random = if (length(rnd)) rnd else NULL,
                        silent = TRUE)
  }
  sec <- proc.time()[3] - t0

  post <- ilm_postcheck(opt, obj, sdr, C, !is.null(ar), pre, Sig, Sigd, re_struct,
                        gap = if (identical(ar$type, "car1")) ar$gap else NULL,
                    kinds = as.list(vapply(re, `[[`, "", "kind")), pnames = pnames,
                    hess = hess, boundary = cb$flagged)
  if (verbose) ilm_print_checks(post, "post-fit convergence checks")
  st <- c(pre$status, post$status)
  if (verbose) {
    cat("\nstructure: ", paste(sprintf("%s = %s%s (%d par, %d latent)", names(re),
        vapply(re_struct, ilm_str_label, ""),
        ifelse(dk > 1L, sprintf(" x %dd", dk), ""), tl, bl), collapse = " | "), "\n", sep = "")
    if (any(st == "FAIL"))
      cat(">> MODEL FIT UNRELIABLE:", sum(st == "FAIL"), "check(s) failed. See the 'why' lines above.\n")
    else if (any(st == "BOUNDARY"))
      cat(">> fixed effects usable. ",
          {
            at <- setdiff(union(hess$held, cb$flagged),
                          names(re)[vapply(re, `[[`, "", "kind") == "basis"])
            if (length(at))
              paste0("The covariance of ", paste(at, collapse = ", "),
                     " sits at a boundary, so that estimate should not be ",
                     "interpreted. ")
            else ""
          },
          "See the BOUNDARY lines above.\n", sep = "")
    else if (any(st == "WARN")) cat(">> fit completed with", sum(st == "WARN"), "warning(s).\n")
    else cat(">> all checks passed.\n")
  }
  ## A gaussian model with nothing integrated out is an ordinary linear model,
  ## where exact t inference is available and strictly better than the
  ## large-sample normal approximation the Laplace machinery would otherwise
  ## give.  Recorded here; vcov(), ilm_coef_table() and ilm_anova() act on it.
  ## Exact t and F inference rests on the response being a gaussian linear
  ## model, where the residual sum of squares is chi-square and independent of
  ## the coefficients. A censored fit is not that model: its likelihood mixes
  ## densities with tail probabilities, the sampling distribution is only
  ## asymptotically normal, and the N/(N-p) correction below is derived for
  ## ordinary least squares and does not apply. Fall back to Wald, which is
  ## what survreg() and every other censored fitter reports.
  ## A dispersion model rules out exact inference for the same reason censoring
  ## does. Exact t and F rest on a constant variance, where the residual sum of
  ## squares is chi-square and independent of the coefficients. With the
  ## variance itself estimated from the data as a function of covariates, that
  ## independence goes and the t distribution is an approximation. gls()
  ## reports t here by convention; this package reserves t for the case where
  ## it is exact.
  exact_df <- fam$name == "gaussian" && sum(bl) == 0L && is.null(ar) &&
    (is.null(cens) || !any(cens != 0L)) && !has_dm

  ## ---- put a REML fit back into the shape the rest of the package expects --
  ## Everything downstream reads opt$par and sdr$cov.fixed and assumes beta sits
  ## at the front of both. Rebuilding that layout once, here, is far safer than
  ## teaching thirty-odd call sites about two parameter orderings -- but it has
  ## to happen AFTER the checks and the covariance derivations above, which run
  ## against obj and need the native REML vector. Reshaping earlier silently
  ## recycles a short logical index over a long vector and corrupts every
  ## variance component.
  if (reml) {
    nb <- length(reml_beta); nt <- length(opt$par)
    opt$par <- c(reml_beta, opt$par)
    ## TYPE labels, as MakeADFun produces them -- not pnames. Downstream code
    ## selects on `names(opt$par) == "theta"`, and coef() renames with pnames.
    names(opt$par) <- c(rep("beta", nb), names(obj$par))
    Vfull <- matrix(0, nb + nt, nb + nt)
    Vfull[seq_len(nb), seq_len(nb)] <- reml_Vb
    Vfull[nb + seq_len(nt), nb + seq_len(nt)] <- sdr$cov.fixed
    ## the off-diagonal stays zero on purpose: under REML the fixed effects and
    ## the variance components are independent, which is exactly what
    ## restricting the likelihood to contrasts orthogonal to X buys
    sdr$cov.fixed <- Vfull
  }
  structure(list(obj = obj, opt = opt, sdr = sdr, checks = rbind(pre, post),
                 ## the ML-shaped twin, present only under REML, used by
                 ## ilm_denom_df() to differentiate V_beta(theta)
                 obj_ml = obj_ml, reml = reml, reml_exact = reml_exact,
                 ## how the standard errors were obtained: "tmb" (the first
                 ## Hessian was fine), "recomputed" (a more accurate one was),
                 ## "boundary" (the covariance of the terms in hessian_held is
                 ## held at its estimate) or "none"; see ilm_hess_recover()
                 hessian_how = hess$how, hessian_held = hess$held,
                 ## terms whose covariance sits at the edge of its range,
                 ## held or not; see ilm_cov_blocks()
                 boundary_terms = cb$flagged,
                 exact_df = exact_df, resid_df = if (exact_df) N - p else NA_integer_,
                 Sigma = Sig, Sigma_d = Sigd, re_struct = re_struct, sec = sec,
                 censor = censor, n_censored = if (is.null(cens)) 0L else sum(cens != 0L),
                 J = J, C = C, n_latent = sum(bl), n_covpar = sum(tl),
                 ## how many scalars TMB actually integrated out. logLik() needs
                 ## it to put back the (q/2)log(2*pi) the hand-written gaussian
                 ## priors below leave out of the joint density.
                 ## under REML beta is in the random block too, and it is not
                 ## one of the latent values logLik() is correcting for
                 n_integrated = length(obj$env$random) - if (reml) p * C else 0L,
                 ## index bookkeeping kept so set_coef() can rebuild every
                 ## derived quantity from a perturbed parameter vector
                 npc = as.integer(npc), npd = as.integer(npd), toff = toff,
                 ty = as.character(ty), rk = as.integer(rk), dcor = as.logical(dcor),
                 b_idx = dl$b_idx, nlk = as.integer(nlk), dk = as.integer(dk),
                 wk = as.integer(wk), re = re, X = X, y = y, ar = ar, Lambda = Lams,
                 weights = if (all(weights == 1)) NULL else weights,
                 pnames = pnames, ylevels = ylevels, family = fam,
                 ## For an exact model the residual SD is reported on the
                 ## unbiased (n - p) scale, so it agrees with lm() and with
                 ## the rescaled standard errors.  Elsewhere it is the
                 ## maximum likelihood estimate.
                 dispersion = if (fam$n_disp > 0L && !has_dm) {
                   d <- exp(pe[pn == "logdisp"])
                   if (exact_df) d <- d * sqrt(N / (N - p))
                   stats::setNames(d, fam$disp_names)
                 } else if (has_dm) {
                   ## the dispersion is no longer one number; report the value
                   ## at the median row so anything expecting a scalar gets a
                   ## representative one, and carry the model beside it
                   stats::setNames(
                     stats::median(ilm_disp_rows(
                       Zd, unname(pe[pn == "gamma"]),
                       if (isTRUE(disp_mu)) unname(pe[pn == "mu_pow"]) else NA_real_,
                       fam, X, matrix(pe[pn == "beta"], p, C))),
                     fam$disp_names[1])
                 } else NULL,
                 disp_formula = if (has_dm) attr(Zd, "formula") else NULL,
                 disp_gamma = if (has_dm) unname(pe[pn == "gamma"]) else NULL,
                 disp_mu_pow = if (has_dm && isTRUE(disp_mu))
                   unname(pe[pn == "mu_pow"]) else NA_real_,
                 disp_coef = if (has_dm)
                   stats::setNames(pe[pn %in% c("gamma", "mu_pow")],
                                   pnames[substr(pnames, 1L, 5L) == "disp:"])
                   else NULL,
                 Zd = Zd, disp_mu = isTRUE(disp_mu), rp = rp,
                 ordinal = has_ord,
                 zeta = if (has_ord) {
                   zr <- unname(pe[pn == "zeta_raw"])
                   stats::setNames(zr[1] + c(0, cumsum(exp(zr[-1]))),
                                   ilm_ord_cut_names(ylevels, J))
                 } else NULL,
                 zeta_se = if (has_ord) {
                   ## the thresholds are a transform of the fitted parameters,
                   ## so their standard errors come from ADREPORT rather than
                   ## from the parameter covariance directly
                   sv <- summary(sdr, "report")
                   sv <- sv[rownames(sv) == "zt", , drop = FALSE]
                   stats::setNames(unname(sv[, 2]), ilm_ord_cut_names(ylevels, J))
                 } else NULL,
                 Zzi = Zzi, zi_type = if (has_zi) zi_type else NULL,
                 zi_formula = if (has_zi) attr(Zzi, "formula") else NULL,
                 zi_gamma = if (has_zi)
                   stats::setNames(unname(pe[pn == "gzi"]), colnames(Zzi))
                   else NULL,
                 zi_se = if (has_zi) {
                   ss <- sqrt(pmax(diag(as.matrix(sdr$cov.fixed)), 0))
                   stats::setNames(unname(ss[pn == "gzi"]), colnames(Zzi))
                 } else NULL,
                 jointPrecision = if (joint) sdr$jointPrecision else NULL,
                 beta = matrix(pe[pn == "beta"], p, C),
                 ## rho is the correlation ONE TIME UNIT apart under both
                 ## structures, so the two are directly comparable. For CAR(1)
                 ## the fitted parameter is log(range), and rho = exp(-1/range).
                 rho = unname(if (is.null(ar)) NA_real_
                       else if (identical(ar$type, "car1"))
                         exp(-1 / exp(pe[pn == "rho_raw"]))
                       else tanh(pe[pn == "rho_raw"])),
                 ar_range = unname(if (!is.null(ar) && identical(ar$type, "car1"))
                   exp(pe[pn == "rho_raw"]) else NA_real_),
                 ok = !any(st == "FAIL")), class = "ilm_model")
}

## ---- the dispersion model --------------------------------------------------

## Per-row dispersion implied by a fitted dispersion model. Needed in two
## places that cannot share code -- inside the likelihood, where everything is
## an AD type, and afterwards on plain numbers -- so this is the plain-number
## one.
#' @keywords internal
#' @noRd
ilm_disp_rows <- function(Zd, gamma, mu_pow, fam, X, beta) {
  ls <- rep(0, nrow(X))
  if (!is.null(Zd) && ncol(Zd) > 0L && length(gamma) == ncol(Zd))
    ls <- ls + as.vector(Zd %*% gamma)
  if (is.finite(mu_pow)) {
    mu <- fam$linkinv(as.vector(X %*% beta[, 1]))
    ls <- ls + mu_pow * log(abs(mu) + 1e-8)
  }
  exp(ls)
}

## ---- refitting a fitted model ----------------------------------------------

## Every internal refit -- the likelihood-ratio tests, the parametric
## bootstrap, every simulation envelope -- has to rebuild the SAME model with a
## different response or a reduced design. Each new structural argument is
## another thing those sites can forget, and forgetting one is silent: the
## family default gave multinomial fits to gaussian data once, and omitting the
## censoring made a reduced likelihood that ignored it, which put the
## likelihood-ratio test's rejection rate at 100% under the null.
##
## So there is one place that knows what "the same model" means.
#' @keywords internal
#' @noRd
ilm_refit_like <- function(object, X = NULL, y = NULL, keep = NULL,
                           restarts = 1L, verbose = FALSE) {
  if (is.null(X)) X <- object$X
  if (is.null(y)) y <- object$y
  rp <- object$rp
  if (!is.null(rp) && !is.null(keep)) {
    ## the derivative design has one column per column of X, so it is subset
    ## the same way; the spline columns are always kept, because their
    ## term map is NA, and they move to wherever they now sit
    rp$D <- rp$D[, keep, drop = FALSE]
    rp$cols <- match(object$rp$cols, which(keep))
    if (anyNA(rp$cols))
      stop("internal: a baseline spline column was dropped from a refit",
           call. = FALSE)
  }
  Zd <- object$Zd
  ## `reml` too: a refit of a REML fit by maximum likelihood is a different
  ## estimator of the same model, so a Type III recode or a consistency check
  ## built on one quietly reported ML's variance components for a REML fit
  ilm_fit(X, y, object$J, ilm_re_list_of(object), object$re_struct, object$ar,
          ylevels = object$ylevels, weights = object$weights,
          family = object$family, verbose = verbose, restarts = restarts,
          censor = object$censor, Zd = Zd, disp_mu = isTRUE(object$disp_mu),
          rp = rp, Zzi = object$Zzi,
          zi_type = if (is.null(object$zi_type)) "inflated" else object$zi_type,
          reml = isTRUE(object$reml))
}

## Everything ilm_refit_like() reads, as plain R objects and nothing else --
## fit$obj holds external pointers into TMB and cannot travel to a parallel
## worker. Both parallel refit paths used to list these fields themselves,
## and both lists went out of date the same way: first without the censoring,
## later without the zero part. A zero-inflated model refitted without it is
## a different model, and the likelihood-ratio test between the two put a
## predictor with no effect at chi-square 68.8 where the right answer was
## 0.02. So the stub is built here, next to the one function that reads it.
#' @keywords internal
#' @noRd
ilm_refit_stub <- function(object)
  list(X = object$X, y = object$y, J = object$J, C = object$C,
       assign = object$assign, term_labels = object$term_labels,
       re = object$re, re_struct = object$re_struct, ar = object$ar,
       ylevels = object$ylevels, family = object$family,
       weights = object$weights, censor = object$censor, Zd = object$Zd,
       disp_mu = isTRUE(object$disp_mu), rp = object$rp, Zzi = object$Zzi,
       zi_type = object$zi_type, reml = isTRUE(object$reml))

#' Fitted dispersion, one value per observation
#'
#' The dispersion a model implies for each row. A constant repeated when the
#' model has a single dispersion parameter, and the fitted dispersion model
#' when it has one.
#'
#' Every piece of machinery that needs a standard deviation -- the residuals,
#' the simulators, the survival curves -- asks for it here rather than reading
#' `object$dispersion`, which is one number and stops being the whole story as
#' soon as a dispersion formula is in play.
#'
#' @param object A fitted `"ilm_model"` object.
#' @return A numeric vector with one entry per observation, or `NULL` when the
#'   family has no dispersion parameter.
#' @keywords internal
#' @noRd
ilm_disp_vec <- function(object) {
  if (is.null(object$dispersion)) return(NULL)
  if (is.null(object$Zd) && !isTRUE(object$disp_mu))
    return(rep(unname(object$dispersion[[1]]), nrow(object$X)))
  ilm_disp_rows(object$Zd, object$disp_gamma,
                if (is.null(object$disp_mu_pow)) NA_real_ else object$disp_mu_pow,
                object$family, object$X, object$beta)
}
