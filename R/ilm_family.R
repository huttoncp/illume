#' Response distributions supported by illume
#'
#' Describes how a response distribution enters the model. The random-effect
#' machinery -- covariance structures, random slopes, smooths, AR(1) -- is
#' identical for every family and never sees the response; a family supplies only
#' the piece that turns a linear predictor into a log-likelihood.
#'
#' @details
#' Supported families:
#'
#' \describe{
#'   \item{`"gaussian"`}{continuous response, identity link, one dispersion
#'     parameter (the residual standard deviation).}
#'   \item{`"binomial"`}{two-category response, logit link. The response may be
#'     0/1, a two-level factor, `TRUE`/`FALSE`, or a character column -- all
#'     four are treated identically. For a factor or character response the
#'     **second** level is the one being modelled, as in [stats::glm()], so
#'     `factor(c("no", "yes"))` models the probability of `"yes"`. A proportion
#'     with `weights` giving the number of trials also works.}
#'   \item{`"poisson"`}{non-negative counts, log link, no dispersion parameter.}
#'   \item{`"nbinom"`}{counts with more variability than Poisson allows, log
#'     link, one dispersion parameter. Variance is `mu + mu^2/k`, so smaller `k`
#'     means more overdispersion; as `k` grows it approaches the Poisson.}
#'   \item{`"multinomial"`}{three or more unordered categories. The only family
#'     needing more than one linear predictor: with `J` categories it uses
#'     `J - 1` dimensions with sum-to-zero coding.}
#' }
#'
#' Zero-inflation, hurdle models, Tweedie and other specialised families are
#' deliberately out of scope. `glmmTMB` covers those well and there is nothing to
#' gain from a weaker reimplementation.
#'
#' @section On numerical stability:
#' The binomial log-likelihood uses `logspace_add()` rather than
#' `log(1 + exp(eta))`. They are mathematically identical, but the naive form
#' overflows to infinity once `eta` exceeds about 709, which an optimiser can
#' easily reach while exploring. The stable form is exact throughout.
#'
#' @param family Character: one of `"gaussian"`, `"binomial"`, `"poisson"`,
#'   `"nbinom"`, `"multinomial"`.
#' @return A list describing the family, with elements `name`, `link`,
#'   `n_disp` (number of dispersion parameters), `disp_names`, `C_of()` (number
#'   of linear predictor dimensions), `nll()`, `linkinv()` and `sim()`.
#' @references
#' McCullagh, P., & Nelder, J. A. (1989). *Generalized Linear Models*, 2nd ed.
#' Chapman & Hall.
#'
#' Hilbe, J. M. (2011). *Negative Binomial Regression*, 2nd ed. Cambridge
#' University Press.
#' @examples
#' ilm_family("poisson")$n_disp      # Poisson has no dispersion parameter
#' ilm_family("nbinom")$disp_names
#' @export
ilm_family <- function(family = c("gaussian", "binomial", "poisson",
                                  "nbinom", "multinomial")) {
  family <- match.arg(family)

  ## Each nll() receives eta (N x C), the response, the frequency weights and
  ## the dispersion parameters on the log scale, and returns the negative
  ## log-likelihood.  Constants that do not involve parameters may be dropped.
  switch(family,

    gaussian = list(
      name = "gaussian", link = "identity", n_disp = 1L,
      disp_names = "log_sigma", C_of = function(J) 1L,
      nll = function(eta, y, w, disp, ...) {
        -sum(w * dnorm(y, eta[, 1], exp(disp[1]), log = TRUE))
      },
      linkinv = function(e) e,
      sim = function(eta, w, disp) stats::rnorm(nrow(eta), eta[, 1], exp(disp[1]))),

    binomial = list(
      name = "binomial", link = "logit", n_disp = 0L,
      disp_names = character(0), C_of = function(J) 1L,
      nll = function(eta, y, w, disp, ...) {
        e <- eta[, 1]
        ## y may be 0/1, or a proportion with w = number of trials
        -sum(w * (y * e - logspace_add(0 * e, e)))
      },
      linkinv = function(e) 1 / (1 + exp(-e)),
      sim = function(eta, w, disp)
        stats::rbinom(nrow(eta), size = pmax(1, round(w)),
                      prob = 1 / (1 + exp(-eta[, 1]))) / pmax(1, round(w))),

    poisson = list(
      name = "poisson", link = "log", n_disp = 0L,
      disp_names = character(0), C_of = function(J) 1L,
      nll = function(eta, y, w, disp, ...) {
        -sum(w * dpois(y, exp(eta[, 1]), log = TRUE))
      },
      linkinv = function(e) exp(e),
      sim = function(eta, w, disp) stats::rpois(nrow(eta), exp(eta[, 1]))),

    nbinom = list(
      name = "nbinom", link = "log", n_disp = 1L,
      disp_names = "log_k", C_of = function(J) 1L,
      nll = function(eta, y, w, disp, ...) {
        mu <- exp(eta[, 1]); k <- exp(disp[1])
        -sum(w * dnbinom2(y, mu = mu, var = mu + mu * mu / k, log = TRUE))
      },
      linkinv = function(e) exp(e),
      sim = function(eta, w, disp)
        stats::rnbinom(nrow(eta), size = exp(disp[1]), mu = exp(eta[, 1]))),

    multinomial = list(
      name = "multinomial", link = "logit", n_disp = 0L,
      disp_names = character(0), C_of = function(J) J - 1L,
      ## y arrives as a count matrix (N x J) and Tct projects the C free
      ## dimensions onto all J categories under sum-to-zero coding
      nll = function(eta, y, w, disp, Tct, ...) {
        etaJ <- eta %*% Tct
        -(sum(y * etaJ) - sum(w * log(rowSums(exp(etaJ)))))
      },
      linkinv = function(e) e,
      sim = function(eta, w, disp, Tct) {
        P <- exp(eta %*% Tct); P <- P / rowSums(P)
        cp <- t(apply(P, 1, cumsum))
        as.integer(rowSums(stats::runif(nrow(P)) > cp)) + 1L
      })
  )
}

#' Check that a response is valid for its family
#'
#' Catches the common mistakes early and explains them, rather than letting the
#' optimiser fail obscurely later: counts that are negative or fractional,
#' binomial responses outside 0 to 1, and so on.
#'
#' @param y The response, after any conversion.
#' @param family A family object from [ilm_family()].
#' @return `TRUE`, invisibly; called for its error messages.
#' @keywords internal
#' @noRd
ilm_check_response <- function(y, family) {
  nm <- family$name
  if (nm == "gaussian" && !is.numeric(y))
    stop("gaussian family needs a numeric response", call. = FALSE)
  if (nm %in% c("poisson", "nbinom")) {
    if (!is.numeric(y) || any(y < 0))
      stop(nm, " family needs non-negative counts", call. = FALSE)
    if (any(abs(y - round(y)) > 1e-8))
      stop(nm, " family needs whole-number counts; the response has fractional ",
           "values. If these are rates, use an offset instead.", call. = FALSE)
  }
  if (nm == "binomial" && (any(y < 0) || any(y > 1)))
    stop("binomial family needs a 0/1 response, or a proportion between 0 and 1 ",
         "with weights giving the number of trials", call. = FALSE)
  invisible(TRUE)
}
