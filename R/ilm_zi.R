## ---------------------------------------------------------------------------
## Zero-inflated and hurdle counts.
##
## Two models for the same symptom and they are not the same model.
##
## A ZERO-INFLATED count is a mixture. Some rows are structural zeros, drawn
## from a process that can only ever produce a zero; the rest come from an
## ordinary count distribution which is itself free to produce a zero. A zero
## in the data could have come from either, and the likelihood adds the two
## routes:
##
##     P(Y = 0) = p + (1 - p) f(0),      P(Y = y) = (1 - p) f(y),  y > 0
##
## A HURDLE count is two separate processes. Whether the response clears zero
## is one model, and how far past zero it goes is another, fitted to a count
## that CANNOT be zero. A zero has exactly one origin:
##
##     P(Y = 0) = p,   P(Y = y) = (1 - p) f(y) / (1 - f(0)),  y > 0
##
## Which one to use is a question about the subject, not about fit. Structural
## zeros mean a subpopulation that was never at risk -- someone who does not
## own a car cannot report a collision. A hurdle means everyone is at risk and
## the first event is governed by something other than the rate afterwards.
## The two give different answers about p: under a hurdle p is the whole
## probability of a zero and is identified by the zero count alone, while under
## zero inflation p is the structural share only, and the count model absorbs
## the rest.
##
## Both are fitted by ilm_model(ziformula = ), which adds a second linear
## predictor for the logit of p, alongside the one for the mean.
## ---------------------------------------------------------------------------

#' The excess-zero probability implied by a fit, one value per row
#'
#' @param object A fitted `"ilm_model"` with a `ziformula`.
#' @param Z Optional design matrix for the zero part; the fit's own by default.
#' @return Numeric vector of probabilities.
#' @keywords internal
#' @noRd
ilm_zi_p <- function(object, Z = NULL) {
  if (is.null(Z)) Z <- object$Zzi
  if (is.null(Z)) return(NULL)
  stats::plogis(as.vector(Z %*% object$zi_gamma))
}

#' Probability of a zero from the count part alone
#'
#' @keywords internal
#' @noRd
ilm_zi_f0 <- function(object, mu) {
  fam <- object$family
  if (is.null(fam$p0))
    stop("family ", fam$name, " has no zero probability", call. = FALSE)
  fam$p0(mu, ilm_zi_disp(object, length(mu)))
}

#' Dispersion aligned to an arbitrary number of prediction rows
#'
#' One value per row when the dispersion is constant, whatever the number of
#' rows asked for. A dispersion model gives a different value to every row of
#' the FIT, and there is no way to carry that onto rows it never saw, so that
#' combination is refused rather than silently recycled.
#'
#' @keywords internal
#' @noRd
ilm_zi_disp <- function(object, n) {
  d <- ilm_disp_vec(object)
  if (is.null(d)) return(rep(1, n))
  d <- unname(d)
  if (length(d) == n) return(d)
  if (length(unique(round(d, 12))) == 1L) return(rep(d[1], n))
  stop("this model has both a dispersion formula and a zero part, and the ",
       "dispersion is a different value for every row it was fitted to. ",
       "Predicting for ", n, " new rows would need their dispersions, which ",
       "are not carried here.", call. = FALSE)
}

#' Marginal mean of a zero-inflated or hurdle response
#'
#' The mean of the response, which is not the mean of the count part: the zero
#' process pulls it down under zero inflation, and under a hurdle it pulls the
#' count part up as well, because a count that cannot be zero has a larger mean
#' than one that can.
#'
#' @keywords internal
#' @noRd
ilm_zi_mean <- function(object, mu, Z = NULL) {
  pz <- ilm_zi_p(object, Z)
  if (is.null(pz)) return(mu)
  if (identical(object$zi_type, "hurdle")) {
    f0 <- ilm_zi_f0(object, mu)
    (1 - pz) * mu / pmax(1 - f0, .Machine$double.eps)
  } else (1 - pz) * mu
}

#' Distribution function of a zero-inflated or hurdle count
#'
#' Needed by the randomised quantile residuals, which take the CDF at the
#' observation and at the value below it and draw uniformly across the jump.
#' Getting this wrong is silent: the residuals still look like residuals, they
#' are simply no longer uniform under a correct model.
#'
#' @keywords internal
#' @noRd
ilm_zi_cdf <- function(object, q, mu, pz) {
  fam <- object$family
  d <- ilm_zi_disp(object, length(mu))
  Fc <- fam$pcdf(q, mu, d)
  out <- if (identical(object$zi_type, "hurdle")) {
    f0 <- fam$p0(mu, d)
    ## below zero there is no mass at all; at zero the hurdle's own p; above
    ## it, p plus the truncated count rescaled onto what is left
    ifelse(q < 0, 0,
           pz + (1 - pz) * pmax(Fc - f0, 0) / pmax(1 - f0, .Machine$double.eps))
  } else pz + (1 - pz) * Fc
  ifelse(q < 0, 0, pmin(pmax(out, 0), 1))
}

#' Draw from a zero-inflated or hurdle count
#'
#' @keywords internal
#' @noRd
ilm_zi_rng <- function(object, mu, disp, pz, ycount) {
  n <- length(mu)
  structural <- stats::runif(n) < pz
  if (identical(object$zi_type, "hurdle")) {
    ## every non-structural row clears the hurdle, so it is drawn from a count
    ## that cannot be zero -- not from the untruncated one and then rejected
    y <- object$family$rtrunc(n, mu, disp)
    y[structural] <- 0L
    as.numeric(y)
  } else {
    ycount[structural] <- 0
    as.numeric(ycount)
  }
}

#' Coefficients of the zero part of a zero-inflated or hurdle fit
#'
#' The second linear predictor a `ziformula` adds, on the logit scale: a
#' positive coefficient means more excess zeros. What "excess" means depends on
#' which model was fitted -- under `zi_type = "inflated"` these are the
#' structural zeros only, over and above the ones the count part produces by
#' itself, while under `"hurdle"` they are every zero there is.
#'
#' @param object A fitted [ilm_model()] with a `ziformula`.
#' @param level Confidence level for the interval.
#' @return A data frame of `term`, `estimate`, `se`, `lower`, `upper`, `z` and
#'   `p`, with an `odds_ratio` column since the scale is a logit.
#' @seealso [ilm_zi_prob()] for the fitted probabilities themselves,
#'   [ilm_check_zeros()] for whether a zero part is called for at all.
#' @examples
#' set.seed(1); n <- 400
#' d <- data.frame(x = rnorm(n), z = rnorm(n))
#' lam <- exp(0.6 + 0.4 * d$x)
#' d$y <- ifelse(runif(n) < plogis(-0.5 + 0.8 * d$z), 0, rpois(n, lam))
#' f <- ilm_model(y ~ x, data = d, family = "poisson", ziformula = ~ z,
#'                verbose = FALSE)
#' ilm_zi_coef(f)
#' @export
ilm_zi_coef <- function(object, level = 0.95) {
  if (!inherits(object, "ilm_model"))
    stop("`object` must be a fitted ilm_model, not ", class(object)[1],
         call. = FALSE)
  if (is.null(object$Zzi))
    stop("this model has no zero part. Fit one with ilm_model(ziformula = ~ 1) ",
         "for a constant excess-zero probability, or ziformula = ~ x to let it ",
         "depend on a predictor.", call. = FALSE)
  est <- object$zi_gamma
  se <- object$zi_se
  if (is.null(se)) se <- rep(NA_real_, length(est))
  crit <- stats::qnorm(1 - (1 - level) / 2)
  z <- est / se
  data.frame(term = colnames(object$Zzi), estimate = est, se = se,
             lower = est - crit * se, upper = est + crit * se,
             odds_ratio = exp(est), z = z,
             p = 2 * stats::pnorm(-abs(z)), row.names = NULL)
}

#' Fitted excess-zero probability, one value per row
#'
#' @param object A fitted [ilm_model()] with a `ziformula`.
#' @param newdata Optional data to predict for; the model frame by default.
#' @return Numeric vector of probabilities between 0 and 1.
#' @seealso [ilm_zi_coef()].
#' @examples
#' set.seed(1); n <- 300
#' d <- data.frame(x = rnorm(n))
#' d$y <- ifelse(runif(n) < 0.3, 0, rpois(n, exp(0.5 + 0.3 * d$x)))
#' f <- ilm_model(y ~ x, data = d, family = "poisson", ziformula = ~ 1,
#'                verbose = FALSE)
#' summary(ilm_zi_prob(f))
#' @export
ilm_zi_prob <- function(object, newdata = NULL) {
  if (!inherits(object, "ilm_model"))
    stop("`object` must be a fitted ilm_model, not ", class(object)[1],
         call. = FALSE)
  if (is.null(object$Zzi))
    stop("this model has no zero part; see ilm_model(ziformula = ).",
         call. = FALSE)
  Z <- if (is.null(newdata)) object$Zzi else
    ilm_zi_design(object$zi_formula, newdata, colnames(object$Zzi))
  ilm_zi_p(object, Z)
}

#' Build the zero-part design matrix for new data
#'
#' @keywords internal
#' @noRd
ilm_zi_design <- function(formula, data, want = NULL) {
  tl <- attr(stats::terms(formula), "term.labels")
  f2 <- if (length(tl)) stats::reformulate(tl, env = environment(formula))
        else stats::as.formula("~ 1", environment(formula))
  Z <- stats::model.matrix(f2, data)
  if (is.null(want)) return(Z)
  miss <- setdiff(want, colnames(Z))
  if (length(miss))
    stop("the new data cannot reproduce the zero part of the model: ",
         paste(miss, collapse = ", "), " missing. A factor level absent from ",
         "the new data is the usual cause.", call. = FALSE)
  Z[, want, drop = FALSE]
}
