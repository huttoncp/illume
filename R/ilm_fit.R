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
  sm <- mgcv::smoothCon(spec, data = data, absorb.cons = TRUE, scale.penalty = TRUE)[[1]]
  re <- mgcv::smooth2random(sm, "", type = 2)
  list(Xf = re$Xf, rand = re$rand, sm = sm, re = re)
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
  flat <- length(weights) > 1L && all(weights == weights[1]) && weights[1] > 1
  wst <- if (noninteger) "FAIL" else if (flat) "WARN" else "OK"
  ck <- ilm_add_check(ck, "weights_type", wst,
    sprintf("%s weights; range %g to %g, total %g",
            if (all(weights == 1)) "none (all 1)" else if (noninteger) "non-integer" else "integer",
            min(weights), max(weights), sum(weights)),
    if (wst != "OK") paste0(
      if (noninteger) "non-integer weights cannot be replicate counts, so these look like sampling weights"
      else "every row carries the same weight > 1, which is a scaling rather than replicate counts",
      "; ilm_model() treats weights as FREQUENCIES. Applied as sampling weights they bias the fixed effects AND understate the standard errors") else "",
    if (wst != "OK")
      "aggregate genuine replicate counts instead; for survey data use a design-based method -- scaled level-specific weights plus a sandwich variance are not implemented here" else "")

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
  if (!is.null(ar)) lat["ar"] <- ar$n_group * ar$Tt * C

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
  st <- if (ratio < 3) "FAIL" else if (ratio < 5) "WARN" else "OK"
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
    sprintf("%.2f observations per latent value (%g observations, %d latent values: %s)",
            ratio, N, tot, paste(sprintf("%s %d", names(lat), lat), collapse = ", ")),
    if (st != "OK") "too few observations per latent value; the Laplace approximation attenuates the variance components and the covariance estimates go rank deficient, even when every individual term passes its own level check" else "",
    sug)

  if (!is.null(ar)) {
    nlat <- ar$n_group * ar$Tt; r2 <- N / nlat
    ck <- ilm_add_check(ck, "obs_per_ar_latent",
      if (r2 < 2) "FAIL" else if (r2 < 4) "WARN" else "OK",
      sprintf("%.2f observations per AR latent time point (%d time points)", r2, as.integer(nlat)),
      if (r2 < 4) "latent AR carries about one categorical observation per latent value; Laplace attenuates the variance components and rho is driven toward the boundary" else "",
      if (r2 < 4) "coarsen the AR time grid, use a reduced-rank AR, or replace AR with s(time) plus a random slope" else "")
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
#' @param has_ar Logical. Whether an AR(1) term was fitted.
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
                      pnames = NULL) {
  kind_of <- function(nm) if (!is.null(kinds) && nm %in% names(kinds)) kinds[[nm]] else "group"
  ck <- ilm_new_checks(); g <- max(abs(obj$gr(opt$par)))
  ck <- ilm_add_check(ck, "optimizer",
    if (opt$convergence == 0) "OK" else "FAIL",
    sprintf("nlminb code %d (%s)", opt$convergence, opt$message),
    if (opt$convergence != 0) "optimizer stopped without meeting its tolerance" else "", "")
  ck <- ilm_add_check(ck, "gradient",
    if (g > 1e-2) "FAIL" else if (g > 1e-3) "WARN" else "OK",
    sprintf("max |gradient| = %.2e", g),
    if (g > 1e-3) "not at a stationary point" else "",
    if (g > 1e-3) "restart from the current estimates, or simplify the random structure" else "")
  pd <- isTRUE(sdr$pdHess)
  anyNaN <- any(!is.finite(diag(sdr$cov.fixed))) || any(diag(sdr$cov.fixed) <= 0)
  lb <- pre$status[pre$check == "latent_budget"]
  ck <- ilm_add_check(ck, "hessian",
    if (!pd || anyNaN) "FAIL" else "OK",
    sprintf("positive definite: %s; non-finite or non-positive variances: %s", pd, anyNaN),
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
    ck <- ilm_add_check(ck, "rho_boundary",
      if (length(rho) && abs(rho) > 0.99) "FAIL" else if (length(rho) && abs(rho) > 0.95) "WARN" else "OK",
      sprintf("rho = %.4f", if (length(rho)) rho else NA_real_),
      if (length(rho) && abs(rho) > 0.95)
        paste0("AR process is near a random walk",
               if (length(lat) && lat %in% c("WARN", "FAIL"))
                 "; the obs-per-AR-latent check also failed, which is the likely driver" else "") else "",
      if (length(rho) && abs(rho) > 0.95) "coarsen the AR grid or drop the AR term" else "")
  }
  cf <- sdr$cov.fixed
  if (!is.null(cf) && all(is.finite(cf)) && all(diag(cf) > 0)) {
    cm <- cov2cor(cf); cm[!upper.tri(cm)] <- 0
    mx <- max(abs(cm)); ij <- which(abs(cm) == mx, arr.ind = TRUE)[1, ]
    pnm <- if (!is.null(pnames) && length(pnames) == ncol(cf)) pnames
           else make.unique(names(obj$par))
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
  ck
}

#' Print a table of checks
#'
#' @param ck A checks data frame.
#' @param title Character heading.
#' @return `ck`, invisibly.
#' @keywords internal
#' @noRd
ilm_print_checks <- function(ck, title) {
  sym <- c(OK = "  ok  ", WARN = " WARN ", FAIL = " FAIL ", INCONCLUSIVE = "  ??  ")
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
#' @param verbose Logical. Print the checks while fitting.
#' @param restarts Integer. Number of optimiser restarts from the previous
#'   solution, which helps on difficult surfaces.
#' @param joint Logical. Also compute the joint precision over fixed and random
#'   parameters. Needed by [predict.ilm_model()] to propagate uncertainty in penalised
#'   smooth coefficients; [ilm_model()] switches it on automatically when the model
#'   contains smooths.
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
                     verbose = TRUE, restarts = 3L, joint = FALSE) {
  fam <- if (is.list(family)) family else ilm_family(family)
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
    ilm_check_response(y, fam)
    yobs <- as.numeric(y)
  }
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
  if (fam$n_disp > 0L) pnames <- c(pnames, fam$disp_names)
  if (!is.null(ar)) pnames <- c(pnames, ilm_nm_tri("ar", C), "ar:rho_raw")

  dl <- list(X = X, yobs = yobs, wrow = weights, Tct = t(Tc),
             n_disp = fam$n_disp,
             grp = lapply(re, `[[`, "group"), Zl = lapply(re, `[[`, "Z"),
             kind = vapply(re, `[[`, "", "kind"), bas = lapply(re, `[[`, "basis"),
             b_idx = lapply(seq_len(K), function(k) (boff[k] + 1L):boff[k + 1L]),
             t_idx = lapply(seq_len(K), function(k) (toff[k] + 1L):toff[k + 1L]),
             nlk = as.integer(nlk), dk = as.integer(dk), wk = as.integer(wk),
             npc = as.integer(npc), ty = as.character(ty), rk = as.integer(rk),
             dcor = as.logical(dcor), K = K, C = C, has_ar = !is.null(ar))
  if (!is.null(ar)) {
    i1 <- (seq_len(ar$n_group) - 1L) * ar$Tt + 1L
    dl$ar_idx <- ar$idx; dl$idx1 <- i1
    dl$idx_t <- setdiff(seq_len(ar$n_group * ar$Tt), i1); dl$idx_lag <- dl$idx_t - 1L
    dl$n_g <- ar$n_group; dl$Tt <- ar$Tt
  }

  fam_nll <- fam$nll                    # captured in the closure, not via dl
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
      La <- ilm_mkL(lchol_ar, C); rho <- tanh(rho_raw)
      W1 <- B_ar[idx1, , drop = FALSE] %*% solve(t(La))
      nll <- nll + 0.5 * sum(W1 * W1) + n_g * sum(log(diag(La)))
      rsd <- B_ar[idx_t, , drop = FALSE] - rho * B_ar[idx_lag, , drop = FALSE]
      W <- (rsd %*% solve(t(La))) / sqrt(1 - rho^2); nre <- n_g * (Tt - 1L)
      nll <- nll + 0.5 * sum(W * W) + nre * (sum(log(diag(La))) + (C / 2) * log(1 - rho^2))
      eta <- eta + B_ar[ar_idx, , drop = FALSE]
      Sa <- La %*% t(La); sdv <- c(sdv, sqrt(diag(Sa)))
      ADREPORT(rho)
    }
    ## The only family-specific line in the whole likelihood.  Everything above
    ## builds eta and accumulates the random-effect prior, and none of it
    ## depends on the response distribution.
    dsp <- if (n_disp > 0L) logdisp else numeric(0)
    nll <- nll + fam_nll(eta, yobs, wrow, dsp, Tct = Tct)
    sd_ <- sdv; ADREPORT(sd_)
    if (n_disp > 0L) { disp_ <- exp(logdisp); ADREPORT(disp_) }
    nll
  }

  pars <- list(beta = matrix(0, p, C))
  ## With no random or smooth terms there is nothing to integrate out, so the
  ## latent parameter block is omitted entirely rather than declared with length
  ## zero -- MakeADFun cannot make an empty vector random.
  if (sum(tl) > 0L) pars$theta <- rep(0, sum(tl))
  if (sum(bl) > 0L) pars$bvec  <- rep(0, sum(bl))
  if (fam$n_disp > 0L) pars$logdisp <- rep(0, fam$n_disp)
  rnd <- if (sum(bl) > 0L) "bvec" else character(0)
  if (!is.null(ar)) {
    pars$lchol_ar <- rep(0, ilm_ncov(C)); pars$rho_raw <- 0.5
    pars$B_ar <- matrix(0, ar$n_group * ar$Tt, C); rnd <- c(rnd, "B_ar")
  }
  t0 <- proc.time()[3]
  obj <- MakeADFun(f, pars, random = if (length(rnd)) rnd else NULL,
                   silent = TRUE)
  ## guard: if a fill order in nm_* ever drifts from ilm_mkL/ilm_mkD/ilm_mkLam/ilm_mkLd, the
  ## names would silently mislabel every parameter.  Fail loudly instead.
  if (length(pnames) != length(obj$par))
    stop(sprintf("internal: %d parameter names for %d parameters",
                 length(pnames), length(obj$par)))
  opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(iter.max = 3000, eval.max = 3000))
  for (k in seq_len(restarts))
    opt <- nlminb(opt$par, obj$fn, obj$gr, control = list(iter.max = 3000, eval.max = 3000))
  ## joint = TRUE also returns the joint precision over (fixed, random), which
  ## is what lets predict() propagate uncertainty in the PENALISED SMOOTH
  ## coefficients instead of holding them at their conditional modes.
  sdr <- suppressWarnings(sdreport(obj, getJointPrecision = joint))
  sec <- proc.time()[3] - t0

  pn <- names(obj$par); pe <- opt$par; thv <- pe[pn == "theta"]
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

  post <- ilm_postcheck(opt, obj, sdr, C, !is.null(ar), pre, Sig, Sigd, re_struct,
                    kinds = as.list(vapply(re, `[[`, "", "kind")), pnames = pnames)
  if (verbose) ilm_print_checks(post, "post-fit convergence checks")
  st <- c(pre$status, post$status)
  if (verbose) {
    cat("\nstructure: ", paste(sprintf("%s = %s%s (%d par, %d latent)", names(re),
        vapply(re_struct, ilm_str_label, ""),
        ifelse(dk > 1L, sprintf(" x %dd", dk), ""), tl, bl), collapse = " | "), "\n", sep = "")
    if (any(st == "FAIL"))
      cat(">> MODEL FIT UNRELIABLE:", sum(st == "FAIL"), "check(s) failed. See the 'why' lines above.\n")
    else if (any(st == "WARN")) cat(">> fit completed with", sum(st == "WARN"), "warning(s).\n")
    else cat(">> all checks passed.\n")
  }
  ## A gaussian model with nothing integrated out is an ordinary linear model,
  ## where exact t inference is available and strictly better than the
  ## large-sample normal approximation the Laplace machinery would otherwise
  ## give.  Recorded here; vcov(), ilm_coef_table() and ilm_anova() act on it.
  exact_df <- fam$name == "gaussian" && sum(bl) == 0L && is.null(ar)
  structure(list(obj = obj, opt = opt, sdr = sdr, checks = rbind(pre, post),
                 exact_df = exact_df, resid_df = if (exact_df) N - p else NA_integer_,
                 Sigma = Sig, Sigma_d = Sigd, re_struct = re_struct, sec = sec,
                 J = J, C = C, n_latent = sum(bl), n_covpar = sum(tl),
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
                 dispersion = if (fam$n_disp > 0L) {
                   d <- exp(pe[pn == "logdisp"])
                   if (exact_df) d <- d * sqrt(N / (N - p))
                   stats::setNames(d, fam$disp_names)
                 } else NULL,
                 jointPrecision = if (joint) sdr$jointPrecision else NULL,
                 beta = matrix(pe[pn == "beta"], p, C),
                 rho = if (!is.null(ar)) tanh(pe[pn == "rho_raw"]) else NA_real_,
                 ok = !any(st == "FAIL")), class = "ilm_model")
}
