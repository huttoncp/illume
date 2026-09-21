#' Describe which outcome a binomial fit is modelling
#'
#' The coefficients of a binomial model describe the probability of one of the
#' two categories, and for a factor or character response that is the **second**
#' level, matching [stats::glm()]. Nothing in the coefficient names records
#' which one, unlike a multinomial fit, so the summary states it.
#'
#' @param o A `"ilm_model"` object with a binomial family.
#' @return A one-line character description.
#' @keywords internal
#' @noRd
ilm_bin_target <- function(o) {
  rsp <- if (!is.null(o$formula)) deparse(o$formula[[2]]) else "y"
  if (!is.null(o$ylevels) && length(o$ylevels) == 2L) {
    ## a logical response has levels named TRUE and FALSE, which are not
    ## strings and should not be printed as though they were
    q <- if (identical(o$ylevels, c("FALSE", "TRUE"))) "" else "'"
    sprintf("modelling P(%s = %s%s%s), with %s%s%s as the reference",
            rsp, q, o$ylevels[2], q, q, o$ylevels[1], q)
  } else if (!is.null(o$y) && all(o$y %in% c(0, 1)))
    sprintf("modelling P(%s = 1)", rsp)
  else
    "response is a proportion; weights give the number of trials"
}

## illume: summary() and print().
##
## One deliberate design choice: summary() SURFACES THE CHECK VERDICTS.  The
## whole point of the diagnostic layer is that this model fails silently -- a
## non-positive-definite Hessian still yields plausible coefficients -- so the
## checks belong where everyone looks, not behind a separate call.

#' Format a covariance matrix as standard deviations and correlations
#'
#' @param S A covariance matrix.
#' @param labs Character vector of row labels.
#' @param digits Integer. Rounding.
#' @return A data frame for printing.
#' @keywords internal
#' @noRd
ilm_fmt_corr <- function(S, labs, digits = 3) {
  sdv <- sqrt(diag(S)); R <- S / outer(sdv, sdv)
  out <- format(round(R, digits), nsmall = digits)
  out[upper.tri(out, diag = TRUE)] <- ""
  keep <- out[, -ncol(out), drop = FALSE]
  colnames(keep) <- labs[-length(labs)]          # name the correlation columns
  cbind(data.frame(SD = format(round(sdv, digits), nsmall = digits),
                   row.names = labs, check.names = FALSE), keep)
}

#' Summarise a fitted multinomial mixed model
#'
#' Prints the model, fit statistics, random-effect covariances, a coefficient
#' table with Wald tests, and -- unusually for a `summary()` method -- the
#' **assumption checks**.
#'
#' @details
#' The checks appear here on purpose. This model class fails quietly: the
#' coefficient table can look completely ordinary, significance stars and all,
#' while the covariance matrix behind it is unusable and every standard error is
#' meaningless. Putting the verdicts behind a separate function would mean the
#' people most likely to be misled are the least likely to look. When any check
#' fails, `summary()` says so directly beneath the coefficients.
#'
#' The printed note about sum-to-zero contrasts is also deliberate: a
#' coefficient here is a deviation from the across-category average, not a
#' contrast against a baseline category, and readers used to
#' [nnet::multinom()] will otherwise assume the latter.
#'
#' Correlation matrices are suppressed above `max_corr_dim` categories, where
#' they become too large to read; use `fit$Sigma` for the full matrices.
#'
#' @param object,x A fitted `"ilm_model"` object.
#' @param digits Integer. Significant digits in the coefficient table.
#' @param max_corr_dim Integer. Largest `C` for which correlation matrices are
#'   printed in full.
#' @param ... Unused.
#' @return `summary()` returns an object of class `"summary.ilm_model"`; its print
#'   method returns it invisibly.
#' @seealso [ilm_coef_table()], [ilm_appraise()].
#' @rdname summary.ilm_model
#' @export
summary.ilm_model <- function(object, ...) {
  structure(list(object = object), class = "summary.ilm_model")
}

#' Summarise a fitted multinomial mixed model
#'
#' Prints the model, fit statistics, random-effect covariances, a coefficient
#' table with Wald tests, and -- unusually for a `summary()` method -- the
#' **assumption checks**.
#'
#' @details
#' The checks appear here on purpose. This model class fails quietly: the
#' coefficient table can look completely ordinary, significance stars and all,
#' while the covariance matrix behind it is unusable and every standard error is
#' meaningless. Putting the verdicts behind a separate function would mean the
#' people most likely to be misled are the least likely to look. When any check
#' fails, `summary()` says so directly beneath the coefficients.
#'
#' The printed note about sum-to-zero contrasts is also deliberate: a
#' coefficient here is a deviation from the across-category average, not a
#' contrast against a baseline category, and readers used to
#' [nnet::multinom()] will otherwise assume the latter.
#'
#' Correlation matrices are suppressed above `max_corr_dim` categories, where
#' they become too large to read; use `fit$Sigma` for the full matrices.
#'
#' @param object,x A fitted `"ilm_model"` object.
#' @param digits Integer. Significant digits in the coefficient table.
#' @param max_corr_dim Integer. Largest `C` for which correlation matrices are
#'   printed in full.
#' @param ... Unused.
#' @return `summary()` returns an object of class `"summary.ilm_model"`; its print
#'   method returns it invisibly.
#' @seealso [ilm_coef_table()], [ilm_appraise()].
#' @rdname summary.ilm_model
#' @export
print.summary.ilm_model <- function(x, digits = 4, max_corr_dim = 6L, ...) {
  o <- x$object; C <- o$C; J <- o$J
  fam <- if (!is.null(o$family)) o$family$name else "multinomial"
  ## "mixed" only if something is actually integrated out; a fixed-effects
  ## gaussian model is an ordinary linear model and is fitted exactly
  mixed <- length(o$re) > 0L || !is.null(o$ar)
  cat(if (isTRUE(o$exact_df))
        "Linear model fit by maximum likelihood (exact t and F inference)\n"
      else if (mixed)
        "Linear mixed model fit by maximum likelihood (Laplace approximation)\n"
      else
        "Generalized linear model fit by maximum likelihood\n")
  if (fam == "multinomial")
    cat(sprintf(" Family: multinomial (%d categories: %s)\n", J,
                paste(o$ylevels, collapse = ", ")))
  else
    cat(sprintf(" Family: %s (%s link)\n", fam, o$family$link))
  ## Which of the two categories the coefficients describe.  A multinomial fit
  ## carries the category in every coefficient name, so the question does not
  ## arise; a binomial fit does not, and with a factor response the direction
  ## is otherwise invisible at the console.
  if (fam == "binomial")
    writeLines(paste0("        ", ilm_bin_target(o)))
  if (!is.null(o$formula)) cat("Formula:", deparse(o$formula), "\n")
  ll <- logLik(o); dev <- -2 * as.numeric(ll)
  nn <- if (!is.null(o$weights)) sum(o$weights) else nrow(o$X)
  cat(sprintf("\n     AIC      BIC   logLik deviance df.resid\n%8.1f %8.1f %8.1f %8.1f %8d\n",
              suppressWarnings(AIC(o)), suppressWarnings(BIC(o)),
              as.numeric(ll), dev, as.integer(nn - attr(ll, "df"))))

  if (length(o$re)) cat("\nRandom effects:\n")
  for (k in seq_along(o$re)) {
    nm <- names(o$re)[k]; s <- o$re_struct[[k]]; e <- o$re[[k]]
    unit <- if (e$kind == "basis") "basis functions" else "levels"
    cat(sprintf(" %s  [%s%s]  %d %s\n", nm, ilm_str_label(s),
                if (e$d > 1L) sprintf(", %d-dim", e$d) else "", e$nl, unit))
    S <- o$Sigma[[k]]
    if (C <= max_corr_dim) {
      print(ilm_fmt_corr(S, o$ylevels[seq_len(C)], 3))
    } else {
      cat("   category SDs: ",
          paste(sprintf("%.3f", sqrt(diag(S))), collapse = " "), "\n", sep = "")
      cat("   (correlation matrix suppressed at C = ", C, "; see fit$Sigma)\n", sep = "")
    }
    if (!is.null(o$Sigma_d[[nm]])) {
      Sd <- o$Sigma_d[[nm]]; sv <- sqrt(diag(Sd))
      cat(sprintf("   within-group: SD ratios %s | correlation %.3f\n",
                  paste(sprintf("%.3f", sv), collapse = ", "),
                  Sd[2, 1] / (sv[1] * sv[2])))
      cat("   (the first within-group SD is fixed at 1: Sigma_cat carries the scale)\n")
    }
  }
  if (!is.null(o$ar)) cat(sprintf(" ar1 over time: rho = %.4f\n", o$rho))
  gl <- vapply(o$re, function(e) if (e$kind == "basis") NA_integer_ else e$nl, 1L)
  gl0 <- vapply(o$re, function(e) if (e$kind == "basis") NA_integer_ else e$nl, 1L)
  if (!length(gl0) || all(is.na(gl0)))
    cat(sprintf("\nNumber of obs: %d%s\n", nrow(o$X),
                if (!is.null(o$weights)) sprintf(" rows (%g weighted)", sum(o$weights)) else ""))
  else
  cat(sprintf("\nNumber of obs: %d%s; groups: %s\n", nrow(o$X),
              if (!is.null(o$weights)) sprintf(" rows (%g weighted)", sum(o$weights)) else "",
              paste(sprintf("%s %d", names(gl)[!is.na(gl)], gl[!is.na(gl)]), collapse = ", ")))
  if (isTRUE(o$exact_df))
    cat(sprintf("Residual degrees of freedom: %d\n", o$resid_df))

  ## dispersion is family-specific: gaussian has a residual SD, the negative
  ## binomial an overdispersion parameter, and binomial/Poisson/multinomial
  ## have none at all
  ## With a dispersion model the single number is only the value at a typical
  ## row, so the model behind it is printed too rather than leaving the reader
  ## to think the spread is constant.
  if (!is.null(o$disp_coef)) {
    cat("
Dispersion model: ", deparse(o$disp_formula), "
", sep = "")
    cf <- o$disp_coef
    for (nm in names(cf))
      cat(sprintf(" %-22s %8.4f
", sub("^disp:", "", nm), cf[[nm]]))
    cat(" (on the log scale; the value below is that at the median row)
")
  }
  if (!is.null(o$dispersion)) {
    cat("\nDispersion:\n")
    for (nm in names(o$dispersion)) {
      lab <- sub("^log_", "", nm)
      note <- if (lab == "sigma") "   (residual standard deviation)"
              else if (lab == "k") "   (smaller means more overdispersion)"
              else ""
      cat(sprintf(" %-8s %.4f%s\n", lab, o$dispersion[[nm]], note))
    }
  }
  cat("\nFixed effects:\n")
  ct <- ilm_coef_table(o)
  stats::printCoefmat(as.matrix(ct), digits = digits, signif.stars = TRUE,
                      has.Pvalue = TRUE, P.values = TRUE)
  ## A zero part is half the model and invisible in the table above, since
  ## those coefficients belong to the count. Printing it here, with what the
  ## probability actually refers to, is the difference between a reader seeing
  ## a zero-inflated fit and seeing an ordinary count model.
  if (!is.null(o$zi_gamma)) {
    hurdle <- identical(o$zi_type, "hurdle")
    cat(sprintf("\nZero part (%s): %s\n",
                if (hurdle) "hurdle" else "zero-inflated",
                deparse(o$zi_formula)))
    zt <- ilm_zi_coef(o)
    zm <- as.matrix(zt[, c("estimate", "se", "z", "p")])
    dimnames(zm) <- list(zt$term,
                         c("Estimate", "Std. Error", "z value", "Pr(>|z|)"))
    stats::printCoefmat(zm, digits = digits, signif.stars = TRUE,
                        has.Pvalue = TRUE, P.values = TRUE)
    pz <- ilm_zi_p(o)
    cat(sprintf(" on the logit scale; fitted P(excess zero) %.3f to %.3f, median %.3f\n",
                min(pz), max(pz), stats::median(pz)))
    cat(if (hurdle)
      " every zero comes from this part; the count above cannot produce one\n"
      else
      " these are STRUCTURAL zeros only -- the count above produces zeros too\n")
  }
  ## only meaningful when there is more than one category dimension
  if (fam == "multinomial") {
    cat("---\n")
    cat("Category contrasts are SUM-TO-ZERO: each coefficient is that category's\n")
    cat("deviation from the across-category average, NOT a contrast against a baseline.\n")
  }

  ## verdicts, in the place people actually look
  ck <- o$checks; bad <- ck[ck$status != "OK", , drop = FALSE]
  cat("\nModel checks: ")
  if (!nrow(bad)) cat("all passed.\n")
  else {
    nf <- sum(bad$status == "FAIL"); nw <- sum(bad$status == "WARN")
    ni <- sum(bad$status %in% c("INCONCLUSIVE", "BOUNDARY"))
    cat(sprintf("%d FAIL, %d WARN, %d inconclusive\n", nf, nw, ni))
    for (i in seq_len(min(4L, nrow(bad))))
      cat(sprintf("  [%s] %s: %s\n", bad$status[i], bad$check[i], bad$detail[i]))
    if (nrow(bad) > 4L) cat(sprintf("  ... and %d more; see fit$checks\n", nrow(bad) - 4L))
    if (nf) cat("  >> FAIL means the standard errors above are not usable.\n")
  }
  invisible(x)
}

#' Compact display of a fitted model
#'
#' A one-line description plus the log-likelihood and AIC. Flags
#' `[CHECKS FAILED]` when any assumption check did not pass, so a failure is
#' visible even from a bare `print()`.
#'
#' @param x A fitted `"ilm_model"` object.
#' @param ... Unused.
#' @return `x`, invisibly.
#' @export
print.ilm_model <- function(x, ...) {
  cat(sprintf("ilm_model fit: %d categories, %d obs, %d fixed + %d covariance parameters%s%s\n",
              x$J, nrow(x$X), ncol(x$X) * x$C, x$n_covpar,
              ## a zero part is not among the fixed effects counted above, and
              ## a one-line print that does not mention it reads as an
              ## ordinary count model
              if (is.null(x$zi_gamma)) "" else
                sprintf(" + %d %s", length(x$zi_gamma),
                        if (identical(x$zi_type, "hurdle")) "hurdle" else "zero-inflation"),
              if (x$ok) "" else "  [CHECKS FAILED]"))
  cat(sprintf("  logLik %.2f | AIC %.1f\n", -x$opt$objective, suppressWarnings(AIC(x))))
  invisible(x)
}
