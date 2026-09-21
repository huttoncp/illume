## ---------------------------------------------------------------------------
## Proportions on the unit interval.
##
## A proportion that is not a count of anything -- percent cover, share of
## time, a score already scaled to [0, 1] -- has no denominator for a binomial
## fit and does not behave like a gaussian one: its spread has to shrink as it
## approaches either end, and a normal model puts mass outside the interval.
##
## Beta regression is the standard answer. The endpoints are the awkward part,
## because the density has no mass at 0 or 1 and a value sitting exactly there
## has no likelihood. There are two honest responses and they are different
## claims about the data:
##
##   The boundary is a SEPARATE PROCESS. Some units were never at risk of
##   being anywhere but zero. That is a two-part model and the zero part
##   belongs in its own formula -- ilm_model(ziformula = , zi_type = "hurdle").
##
##   The boundary is ROUNDING. A proportion recorded as 0 was really some
##   small number the instrument could not resolve. Then shifting it inside is
##   defensible, and ilm_squeeze() does it by the amount Smithson and
##   Verkuilen propose, which shrinks every value towards 1/2 rather than
##   nudging only the offending ones and leaving a gap.
##
## Guessing between the two silently would be the worst option, so
## ilm_check_response() names both and fits neither.
## ---------------------------------------------------------------------------

#' Move proportions off 0 and 1
#'
#' The beta family is defined on the open interval, so a proportion recorded as
#' exactly 0 or 1 has no likelihood. When those values are a rounding of
#' interior ones rather than a separate process, the usual remedy is the
#' transformation of Smithson and Verkuilen (2006):
#'
#' \deqn{y' = \frac{y (n - 1) + 1/2}{n}}
#'
#' which compresses the whole scale towards one half by a factor of
#' `(n - 1) / n`. Every value moves, not only the offending ones: nudging the
#' endpoints alone would leave a gap just inside them and change the shape of
#' the distribution more than this does.
#'
#' It is a choice, not a repair. The function reports how far it moved things
#' so the choice is visible, and with `n` small the shift is not negligible --
#' at `n = 50` an observed 0 becomes 0.01.
#'
#' @param y A numeric vector of proportions in `[0, 1]`.
#' @param n Sample size to use in the transformation; the number of
#'   observations by default.
#' @param quiet Suppress the message reporting the shift.
#' @return A numeric vector strictly inside `(0, 1)`.
#' @references Smithson, M. and Verkuilen, J. (2006). A better lemon squeezer?
#'   Maximum-likelihood regression with beta-distributed dependent variables.
#'   *Psychological Methods* 11, 54-71.
#' @seealso [ilm_model()] with `family = "beta"`, and `ziformula` for the case
#'   where the boundary values are a separate process rather than rounding.
#' @examples
#' y <- c(0, 0.2, 0.5, 0.9, 1)
#' ilm_squeeze(y)
#' @export
ilm_squeeze <- function(y, n = length(y), quiet = FALSE) {
  if (!is.numeric(y)) stop("`y` must be numeric", call. = FALSE)
  ok <- !is.na(y)
  if (any(y[ok] < 0 | y[ok] > 1))
    stop("`y` has values outside [0, 1]; this transformation is for ",
         "proportions.", call. = FALSE)
  if (!is.numeric(n) || length(n) != 1L || n < 2)
    stop("`n` must be a single number of at least 2", call. = FALSE)
  out <- (y * (n - 1) + 0.5) / n
  if (!quiet) {
    nb <- sum(y[ok] <= 0 | y[ok] >= 1)
    message(sprintf(
      "ilm_squeeze: %d value%s at a boundary; every value moved towards 1/2 ",
      nb, if (nb == 1L) "" else "s"),
      sprintf("by up to %.4g. An observed 0 is now %.4g and an observed 1 is %.4g.",
              max(abs(out - y), na.rm = TRUE), 0.5 / n, (n - 0.5) / n))
  }
  out
}

#' Score and weight for the beta family
#'
#' The derivative of a row's log-likelihood with respect to its own linear
#' predictor, which is what the sandwich estimator needs. For the beta
#' log-likelihood
#'
#'   l = log G(phi) - log G(mu phi) - log G((1-mu) phi)
#'       + (mu phi - 1) log y + ((1-mu) phi - 1) log(1-y)
#'
#' the derivative with respect to mu is phi times the difference between the
#' observed logit of y and the digamma pair, and dmu/deta is mu (1 - mu).
#'
#' @keywords internal
#' @noRd
ilm_beta_u <- function(y, eta, phi) {
  mu <- stats::plogis(eta)
  ystar <- log(y) - log1p(-y)
  mustar <- digamma(mu * phi) - digamma((1 - mu) * phi)
  phi * (ystar - mustar) * mu * (1 - mu)
}

#' @keywords internal
#' @noRd
ilm_beta_w <- function(eta, phi) {
  mu <- stats::plogis(eta)
  phi * (trigamma(mu * phi) + trigamma((1 - mu) * phi)) * (mu * (1 - mu))^2
}
