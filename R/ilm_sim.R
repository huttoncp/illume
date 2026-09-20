## ---------------------------------------------------------------------------
## A simulated mixed-type dataset for testing and examples.
##
## The point is not to look like real data but to EXERCISE THE DIAGNOSTICS:
## every column below is there because some check should have something to say
## about it. A fixture where nothing is wrong tests nothing.
##
##   id         panel identifier, repeated across periods
##   date       monthly, regular, repeated across units (long format)
##   grp        factor with a deliberately rare level and an unused level
##   site       character with case/whitespace variants and empty strings
##   flag       balanced logical
##   consented  near-constant logical (separation risk)
##   score      gaussian
##   income     lognormal: right-skewed, bounded at zero
##   visits     Poisson counts (dispersion ~ 1)
##   claims     negative binomial counts (dispersion > 1)
##   downtime   zero-inflated counts
##   cohort     constant
##   lab_value  gaussian with missing values
## ---------------------------------------------------------------------------

#' A simulated mixed-type dataset for testing and examples
#'
#' Entirely synthetic, so it carries no third-party copyright. Every column is
#' present because some diagnostic should have something to say about it: a
#' fixture where nothing is wrong tests nothing.
#'
#' \describe{
#'   \item{id}{panel identifier, repeated across periods}
#'   \item{date}{monthly, regular, repeated across units (long format)}
#'   \item{grp}{factor with a deliberately rare level and an unused level}
#'   \item{site}{character with case and whitespace variants, and empty strings}
#'   \item{flag}{balanced logical}
#'   \item{consented}{near-constant logical, a separation risk}
#'   \item{score}{gaussian}
#'   \item{income}{lognormal: right-skewed and bounded at zero}
#'   \item{visits}{Poisson counts, dispersion near 1}
#'   \item{claims}{negative binomial counts, dispersion above 1}
#'   \item{downtime}{zero-inflated counts}
#'   \item{cohort}{constant}
#'   \item{lab_value}{gaussian with missing values}
#' }
#'
#' @param n_id Number of units.
#' @param n_period Observations per unit.
#' @param seed Random seed, so the fixture is reproducible.
#' @return A data frame with `n_id * n_period` rows and 13 columns.
#' @examples
#' d <- ilm_sim()
#' str(d)
#' @export
ilm_sim <- function(n_id = 75L, n_period = 12L, seed = 2026L) {
  stopifnot(n_id >= 2L, n_period >= 1L)
  set.seed(seed)
  n <- n_id * n_period
  id <- rep(seq_len(n_id), each = n_period)

  ## regular monthly spacing, repeated across units: the shape panel data
  ## actually takes, and the case where "duplicate dates" is normal
  date <- rep(seq(as.Date("2024-01-01"), by = "month", length.out = n_period),
              times = n_id)

  ## one level is rare on purpose -- thin levels are the usual reason a factor
  ## model will not fit -- and one level is declared but never observed
  ## the rare level is assigned explicitly rather than sampled, so it stays
  ## below the default rare_n threshold at any n and the check always fires
  gv <- sample(c("alpha", "beta", "gamma"), n, TRUE, prob = c(0.45, 0.35, 0.20))
  gv[sample.int(n, 3L)] <- "delta"
  grp <- factor(gv, levels = c("alpha", "beta", "gamma", "delta", "epsilon"))

  site <- sample(c("North", "north", "North ", "South", "", "East"), n, TRUE,
                 prob = c(0.30, 0.06, 0.04, 0.30, 0.05, 0.25))

  flag      <- sample(c(TRUE, FALSE), n, TRUE)
  consented <- sample(c(TRUE, FALSE), n, TRUE, prob = c(0.995, 0.005))

  u <- stats::rnorm(n_id, 0, 0.8)[id]          # a per-unit random effect
  score  <- round(50 + 6 * u + stats::rnorm(n, 0, 5), 2)
  income <- round(stats::rlnorm(n, meanlog = 10.4 + 0.25 * u, sdlog = 0.55), 2)
  visits <- stats::rpois(n, exp(1.1 + 0.3 * u))
  claims <- stats::rnbinom(n, size = 1.2, mu = exp(1.4 + 0.3 * u))
  downtime <- ifelse(stats::runif(n) < 0.65, 0L, stats::rpois(n, 4))

  lab_value <- round(stats::rnorm(n, 7.4, 1.1), 2)
  lab_value[sample.int(n, floor(0.12 * n))] <- NA_real_

  data.frame(id = id, date = date, grp = grp, site = site,
             flag = flag, consented = consented,
             score = score, income = income, visits = visits,
             claims = claims, downtime = downtime,
             cohort = "2024", lab_value = lab_value,
             stringsAsFactors = FALSE)
}
