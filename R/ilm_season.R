## ---------------------------------------------------------------------------
## Cycles: the remedy that has to exist for the diagnosis to be worth making.
##
## ilm_check_ar() can now tell a cycle apart from autoregression and say so.
## That is only useful if the package can also FIT the cycle, otherwise the
## verdict sends the user somewhere else. A cycle of known period is ordinary
## fixed effects -- there is no new likelihood here, only the right columns.
##
## Two ways to write those columns:
##
##   ilm_fourier(t, 12)          a smooth cycle, 2 columns per harmonic
##   ilm_cyclic(t, 12, df = 5)   a smooth cycle with local flexibility
##   factor(month)               an arbitrary shape, period - 1 columns
##
## All three are smooth-to-jumpy in that order, and cheap-to-expensive in that
## order too. A 12-level factor spends 11 degrees of freedom on a shape that two
## or four columns usually capture, and the spare degrees of freedom come
## straight out of the precision of everything else in the model.
##
## Measured at MATCHED degrees of freedom, the two smooth bases are close.
## Across four shapes on a 12-phase cycle (sinusoid, asymmetric, narrow spike,
## two peaks) the better of the two won by 0 to 4 AIC on the smooth shapes, and
## on the two shapes that genuinely jumped the FACTOR beat both by 100 AIC or
## more. On a long cycle -- day of year, one narrow summer feature -- Fourier
## led below about six degrees of freedom and the spline led above it, by at
## most 23 AIC.
##
## So the spline is not the basis that wins; it is the basis whose flexibility
## is local. It earns its place on long cycles, where a factor is unaffordable
## and adding resolution to one part of the cycle should not disturb the rest.
##
## Both close on themselves, which an ordinary bs() or ns() on a phase variable
## does not -- that leaves a step between December and January and spends a
## degree of freedom estimating a jump that is not there.
## ---------------------------------------------------------------------------

#' Fourier terms for a cycle of known period
#'
#' Sine and cosine columns for use directly in a model formula, when the
#' response has a repeating cycle -- a month-of-year effect in monthly data, a
#' day-of-week effect in daily data, a tide, a shift pattern.
#'
#' @details
#' `K` sets how much shape the cycle is allowed to have. `K = 1` is a single
#' smooth peak and trough per cycle; higher `K` adds harmonics and can bend into
#' sharper or double-peaked shapes, at two more columns each. A cycle of length
#' `period` supports at most `floor(period / 2)` harmonics, beyond which the
#' columns are not distinguishable at the observed times.
#'
#' Unlike [stats::poly()] this depends on nothing but its arguments, so a model
#' containing it predicts correctly on new data, including times outside the
#' range that was fitted.
#'
#' @section Choosing a basis:
#' A factor for the phase -- `factor(month)` -- fits any shape at all, but
#' spends `period - 1` degrees of freedom doing it, and gives no reason for
#' neighbouring phases to resemble each other. Fourier terms spend `2K` and are
#' smooth by construction. For a 12-month cycle, `K = 2` (four columns) covers
#' most real seasonal shapes.
#'
#' Go to a factor when the cycle has a genuine step in it, such as a policy that
#' starts in April: measured on a narrow two-month spike, a factor beat both
#' smooth bases by more than 140 AIC at matched sample size, because no smooth
#' basis can make a step cheaply. See [ilm_cyclic()] for the local alternative.
#'
#' @param x Time, as a number, integer or `Date`. One unit is one step, so the
#'   period is expressed in the same units.
#' @param period Length of one full cycle: 12 for months of a year, 7 for days
#'   of a week, 365.25 for days of a year.
#' @param K Number of harmonics.
#' @return A numeric matrix with one row per element of `x` and `2 * K` columns
#'   (one fewer when `2 * K` equals `period`, where the last sine column is
#'   identically zero).
#' @seealso [ilm_check_ar()], which names the period when the residuals contain
#'   a cycle.
#' @examples
#' head(ilm_fourier(1:24, period = 12))
#'
#' set.seed(1)
#' d <- data.frame(id = factor(rep(1:20, each = 24)), t = rep(1:24, 20))
#' d$y <- 2 * sin(2 * pi * d$t / 12) + rnorm(480)
#' f <- ilm_model(y ~ ilm_fourier(t, 12) + (1 | id), data = d,
#'                family = "gaussian", verbose = FALSE)
#' coef(f)
#' @export
ilm_fourier <- function(x, period, K = 1L) {
  if (missing(period))
    stop("`period` is required: it is the length of one full cycle, in the ",
         "same units as `x` (12 for months of a year, 7 for days of a week).",
         call. = FALSE)
  if (is.factor(x))
    stop("`x` is a factor. Converting it here would use the level order rather ",
         "than the labels, and silently change the period. Pass as.integer(x) ",
         "if the levels are equally spaced steps.", call. = FALSE)
  xn <- suppressWarnings(as.numeric(x))
  if (any(is.na(xn) & !is.na(x)))
    stop("`x` must be numeric, integer or Date, not ", class(x)[1],
         call. = FALSE)
  period <- suppressWarnings(as.numeric(period)[1])
  if (!is.finite(period) || period <= 1)
    stop("`period` must be a number greater than 1, not ", period, call. = FALSE)
  K <- suppressWarnings(as.integer(K)[1])
  if (is.na(K) || K < 1L)
    stop("`K` must be a whole number of at least 1", call. = FALSE)
  kmax <- floor(period / 2)
  if (K > kmax)
    stop("`K` is ", K, " but a period of ", period, " supports at most ", kmax,
         " harmonic", if (kmax == 1L) "" else "s",
         ": beyond that the columns repeat the ones already there. Use K <= ",
         kmax, ".", call. = FALSE)

  ang <- 2 * pi * xn / period
  cols <- list(); nms <- character(0)
  for (k in seq_len(K)) {
    ## at 2k == period the sine is zero at every whole-numbered time, so the
    ## column is empty rather than merely collinear
    if (2 * k != period) {
      cols[[length(cols) + 1L]] <- sin(k * ang); nms <- c(nms, paste0("sin", k))
    }
    cols[[length(cols) + 1L]] <- cos(k * ang); nms <- c(nms, paste0("cos", k))
  }
  out <- matrix(unlist(cols), nrow = length(xn),
                dimnames = list(NULL, nms))
  attr(out, "period") <- period
  attr(out, "K") <- K
  out
}

#' Cyclic cubic spline basis for a cycle of known period
#'
#' A smooth periodic basis that joins up at the wrap-around: the value and the
#' first two derivatives match at the end of one cycle and the start of the
#' next, so December runs into January without a step.
#'
#' @details
#' The alternative to [ilm_fourier()], and the one whose flexibility is
#' *local*: changing the curve near one knot leaves the rest of the cycle
#' alone, where adding a Fourier harmonic changes the shape everywhere.
#'
#' Do not expect a large gain from that. Compared against [ilm_fourier()] at
#' matched degrees of freedom on a twelve-phase cycle, the two were within a few
#' AIC of each other for every smooth shape tried, and where the shape really
#' jumped a factor beat both by a hundred or more. On a long cycle -- day of
#' year with one narrow summer feature -- Fourier led below about six degrees of
#' freedom and this led above it, by at most 23 AIC.
#'
#' Reach for it when the cycle is long, when you want more resolution in one
#' part of it without disturbing the rest, or when the alternative is a factor
#' with more levels than the data can pay for. Reach for [ilm_fourier()] on a
#' short cycle or a plain rise and fall, and for a factor when the shape has a
#' genuine step in it.
#'
#' `df` is the number of columns returned, and so the degrees of freedom the
#' term spends -- the same convention as [splines::ns()]. One basis function is
#' dropped, because a cyclic basis sums to one at every point and would
#' otherwise be collinear with the intercept; dropping it leaves the fitted
#' curve and the predictions unchanged, since the remaining columns and the
#' intercept span exactly the same space.
#'
#' Do not use an ordinary [splines::bs()] or [splines::ns()] on a phase variable
#' for this. Neither knows the two ends of the cycle are the same point, so both
#' leave a discontinuity at the wrap-around and waste degrees of freedom
#' estimating a jump that is not there.
#'
#' @param x Time, as a number, integer or `Date`. Values are reduced modulo
#'   `period`, so the raw time works as well as the phase.
#' @param period Length of one full cycle, in the units of `x`.
#' @param df Number of columns, and so the degrees of freedom spent. Must be at
#'   least 3.
#' @return A numeric matrix with `df` columns.
#' @seealso [ilm_fourier()] for the cheaper smooth alternative,
#'   [ilm_check_ar()], which names the period when the residuals contain a cycle.
#' @examples
#' b <- ilm_cyclic(1:12, period = 12, df = 4)
#' dim(b)
#' # the basis is periodic: one cycle on is the same point
#' all.equal(ilm_cyclic(1:6, 12, 4), ilm_cyclic(13:18, 12, 4),
#'           check.attributes = FALSE)
#' @export
ilm_cyclic <- function(x, period, df = 5L) {
  if (missing(period))
    stop("`period` is required: it is the length of one full cycle, in the ",
         "same units as `x`.", call. = FALSE)
  if (is.factor(x))
    stop("`x` is a factor. Converting it here would use the level order rather ",
         "than the labels, and silently change the period. Pass as.integer(x) ",
         "if the levels are equally spaced steps.", call. = FALSE)
  xn <- suppressWarnings(as.numeric(x))
  if (any(is.na(xn) & !is.na(x)))
    stop("`x` must be numeric, integer or Date, not ", class(x)[1],
         call. = FALSE)
  period <- suppressWarnings(as.numeric(period)[1])
  if (!is.finite(period) || period <= 1)
    stop("`period` must be a number greater than 1, not ", period, call. = FALSE)
  df <- suppressWarnings(as.integer(df)[1])
  if (is.na(df) || df < 3L)
    stop("`df` must be at least 3: a cubic spline that closes on itself needs ",
         "four basis functions, one of which is dropped for the intercept.",
         call. = FALSE)
  if (df + 1L > period && abs(period - round(period)) < 1e-8)
    stop("`df` is ", df, " but a period of ", round(period), " has only ",
         round(period), " distinct phases, so at most ", round(period) - 1L,
         " columns are identifiable. Use df <= ", round(period) - 1L,
         ", or factor(x %% ", round(period), ") for a free shape.",
         call. = FALSE)

  ## nk periodic basis functions over nk equal intervals, built by extending the
  ## knots three beyond each end (cubic, order 4) and folding the three columns
  ## that wrap back onto the three they duplicate
  nk <- df + 1L
  h <- period / nk
  kn <- seq(0, period, length.out = nk + 1L)
  ext <- c(kn[1L] - h * 3:1, kn, kn[nk + 1L] + h * 1:3)
  xm <- xn %% period
  ok <- is.finite(xm)
  B <- matrix(NA_real_, length(xm), nk + 3L)
  if (any(ok))
    B[ok, ] <- splines::splineDesign(ext, xm[ok], ord = 4L, outer.ok = TRUE)
  out <- B[, seq_len(nk), drop = FALSE]
  out[, 1:3] <- out[, 1:3] + B[, nk + 1:3, drop = FALSE]

  ## the columns sum to 1 at every x, so one of them is the intercept written
  ## a long way round; dropping it changes nothing the model can see
  out <- out[, -nk, drop = FALSE]
  dimnames(out) <- list(NULL, paste0("c", seq_len(df)))
  attr(out, "period") <- period
  attr(out, "df") <- df
  out
}
