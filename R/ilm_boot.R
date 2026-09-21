## ---------------------------------------------------------------------------
## Bootstrap intervals, group differences, and missingness summaries.
##
## BCa intervals are computed directly -- bias correction from the bootstrap
## distribution, acceleration from the jackknife -- rather than delegated to
## `boot`. That is a speed decision rather than a dependency one: `boot`
## reaches its influence values through a generic path that re-enters the
## statistic with weight vectors, which costs 7x at n = 500 and 17x at
## n = 2000 against a direct jackknife.
## ---------------------------------------------------------------------------

ILM_CI_TYPES <- c("percentile", "bca", "normal", "basic")

## ---- one interval ----------------------------------------------------------

#' @keywords internal
#' @noRd
ilm_boot_stat <- function(y, stat_fun, R, conf, ci_type, stat_label = NULL,
                          progress = NULL) {
  y <- y[!is.na(y)]
  n <- length(y)
  obs <- stat_fun(y)
  if (n < 2L || !is.finite(obs))
    return(list(observed = obs, lower = NA_real_, upper = NA_real_, n = n))
  ## One replicate at a time. Drawing all the indices at once into an n x R
  ## matrix reads well and does not scale: at n = 20,000 and R = 2000 that is a
  ## 160 MB allocation, and at a million rows it is 8 GB, so the function stops
  ## working rather than merely slowing down. Per replicate the same n indices
  ## are reused, and it measured 12% faster besides.
  pb <- ilm_progress(R, progress)
  th  <- vapply(seq_len(R), function(r) {
    pb$tick(r); stat_fun(y[sample.int(n, n, replace = TRUE)]) }, 1)
  pb$done()
  th  <- th[is.finite(th)]
  if (!length(th))
    return(list(observed = obs, lower = NA_real_, upper = NA_real_, n = n))
  a2 <- (1 - conf) / 2
  lu <- switch(ci_type,
    percentile = unname(stats::quantile(th, c(a2, 1 - a2), names = FALSE)),
    basic      = c(2 * obs - stats::quantile(th, 1 - a2, names = FALSE),
                   2 * obs - stats::quantile(th, a2, names = FALSE)),
    normal     = obs + c(-1, 1) * stats::qnorm(1 - a2) * stats::sd(th),
    bca        = ilm_bca(y, th, obs, stat_fun, conf, stat_label))
  list(observed = obs, lower = lu[1], upper = lu[2], n = n)
}

## Bias-corrected and accelerated percentiles (Efron & Tibshirani 1993, ch. 14).
## z0 measures how far the bootstrap distribution sits from the observed value;
## the acceleration `a` comes from the jackknife and corrects for a statistic
## whose variance changes with its value.
#' @keywords internal
#' @noRd
ilm_bca <- function(y, th, obs, stat_fun, conf, stat_label = NULL) {
  n <- length(y)
  prop <- mean(th < obs)
  ## with no bootstrap replicate on one side the correction is undefined, so
  ## fall back to percentile rather than returning an infinite endpoint
  if (prop <= 0 || prop >= 1)
    return(unname(stats::quantile(th, c((1 - conf) / 2, 1 - (1 - conf) / 2),
                                  names = FALSE)))
  z0 <- stats::qnorm(prop)
  jk <- ilm_jackknife(y, stat_fun, stat_label)
  jm <- mean(jk); dv <- jm - jk
  den <- 6 * (sum(dv^2))^1.5
  a <- if (den == 0) 0 else sum(dv^3) / den
  z <- stats::qnorm(c((1 - conf) / 2, 1 - (1 - conf) / 2))
  adj <- stats::pnorm(z0 + (z0 + z) / (1 - a * (z0 + z)))
  unname(stats::quantile(th, adj, names = FALSE))
}

## Every value's leave-one-out statistic.
##
## Done literally this is O(n^2) -- n calls, each on a vector of n - 1 -- and it
## dominates a BCa interval: 4.31s of a 6.07s call at n = 20,000, where the
## resampling itself took 2.55s. The mean, variance and standard deviation all
## have exact leave-one-out forms in terms of the running sums, which is O(n)
## and measured at 0.00s on the same data. Anything else, including a
## user-supplied function, still takes the loop, because nothing can be assumed
## about it.
#' @keywords internal
#' @noRd
ilm_jackknife <- function(y, stat_fun, stat_label = NULL) {
  n <- length(y)
  if (!is.null(stat_label) && stat_label %in% c("mean", "var", "sd")) {
    S <- sum(y)
    if (stat_label == "mean") return((S - y) / (n - 1))
    if (n < 3L) return(vapply(seq_len(n), function(i) stat_fun(y[-i]), 1))
    ## var of the leave-one-out sample, from the sums rather than the sample:
    ## SS_{-i} = SS - y_i^2 and mean_{-i} = (S - y_i)/(n - 1)
    SS <- sum(y^2)
    v <- (SS - y^2 - (S - y)^2 / (n - 1)) / (n - 2)
    v[v < 0] <- 0                     # only ever a rounding artefact
    return(if (stat_label == "sd") sqrt(v) else v)
  }
  ## A custom statistic leaves no choice but the loop, and at large n that is
  ## minutes rather than seconds. Say so rather than appearing to hang.
  if (n > 5000L)
    message("BCa with a custom statistic needs ", format(n, big.mark = ","),
            " leave-one-out evaluations, which will take a while. ",
            "ci_type = \"percentile\" avoids the jackknife entirely.")
  vapply(seq_len(n), function(i) stat_fun(y[-i]), 1)
}

#' @keywords internal
#' @noRd
ilm_stat_fun <- function(stat) {
  if (is.function(stat)) return(stat)
  switch(stat, mean = base::mean, median = stats::median,
         sd = stats::sd, var = stats::var,
         stop("unknown `stat`: ", sQuote(stat),
              ". Options are 'mean', 'median', 'sd', 'var', or a function ",
              "taking a numeric vector.", call. = FALSE))
}

## ---- user-facing -----------------------------------------------------------

#' Bootstrap confidence interval for a statistic
#'
#' @param data A data frame, or a numeric vector when `y` is `NULL`.
#' @param y Name of the numeric column to summarise.
#' @param by Optional grouping columns.
#' @param stat `"mean"`, `"median"`, `"sd"`, `"var"`, or a function taking a
#'   numeric vector.
#' @param R Bootstrap replicates.
#' @param conf Confidence level.
#' @param ci_type `"percentile"`, `"bca"`, `"normal"` or `"basic"`. On skewed
#'   data percentile and BCa hold their nominal coverage better than the other
#'   two.
#' @param seed Random seed.
#' @param progress Show a progress bar. Defaults to [interactive()], so a
#'   bar appears when someone is watching and nothing is written in a
#'   script or a knitted document. See [ilm_progress_arg].
#' @return A one-row data frame per group, with `observed`, `lower`, `upper`
#'   and the settings used.
#' @references
#' Efron, B. and Tibshirani, R. J. (1993). An Introduction to the Bootstrap.
#' Chapman and Hall.
#' @seealso [ilm_boot_diff()] for a difference between two groups.
#' @examples
#' d <- ilm_sim()
#' ilm_boot_ci(d, "score", R = 200, seed = 1)
#' ilm_boot_ci(d, "score", by = "grp", R = 200, seed = 1)
#' @export
ilm_boot_ci <- function(data, y = NULL, by = NULL, stat = "mean",
                        R = 2000L, conf = 0.95, ci_type = "percentile",
                        seed = NULL, progress = NULL) {
  if (length(ci_type) != 1L || !ci_type %in% ILM_CI_TYPES)
    stop("unknown `ci_type`: ", paste(sQuote(ci_type), collapse = ", "),
         ". Options are ", paste(sQuote(ILM_CI_TYPES), collapse = ", "), ".",
         call. = FALSE)
  if (!is.numeric(conf) || length(conf) != 1L || conf <= 0 || conf >= 1)
    stop("`conf` must be a single number strictly between 0 and 1", call. = FALSE)
  if (!is.numeric(R) || length(R) != 1L || R < 2)
    stop("`R` must be a single number of at least 2", call. = FALSE)
  sf <- ilm_stat_fun(stat)
  lab <- if (is.function(stat)) "custom" else stat

  if (is.numeric(data) && is.null(y)) { v <- data; data <- NULL } else {
    if (is.null(y)) stop("`y` must name the numeric column to summarise",
                         call. = FALSE)
    miss <- setdiff(c(y, by), names(data))
    if (length(miss))
      stop("column(s) not found in the data: ", paste(miss, collapse = ", "),
           ". Available: ", paste(utils::head(names(data), 12), collapse = ", "),
           call. = FALSE)
    if (!is.numeric(data[[y]]))
      stop("`y` (", y, ") must be numeric; it is ", class(data[[y]])[1],
           call. = FALSE)
    v <- data[[y]]
  }
  if (!is.null(seed)) set.seed(seed)

  one <- function(vec) {
    s <- ilm_boot_stat(vec, sf, R, conf, ci_type,
                       if (is.function(stat)) NULL else stat, progress)
    data.frame(stat = lab, observed = s$observed, lower = s$lower,
               upper = s$upper, conf = conf, R = as.integer(R),
               ci_type = ci_type, n = s$n, stringsAsFactors = FALSE)
  }
  if (is.null(by) || is.null(data)) return(one(v))
  g <- interaction(data[by], drop = TRUE)
  parts <- lapply(split(v, g), one)
  cbind(setNames(data.frame(names(parts), stringsAsFactors = FALSE),
                 paste(by, collapse = ".")),
        do.call(rbind, parts), row.names = NULL)
}

## ---- differences between groups ---------------------------------------------
##
## The difference itself with its own interval: the quantity a reader actually
## wants, rather than a set of intervals to eyeball for overlap (judging a
## difference by whether separate intervals overlap is conservative and lossy).
##
## Past two groups the multiplicity IS the problem. Five levels give ten
## comparisons, and ten intervals each nominally 95% do not jointly cover at
## 95%. So the default reports SIMULTANEOUS intervals, calibrated by the
## bootstrap rather than read off a distributional table:
##
##   1. resample within each group once per replicate, so every comparison in a
##      replicate sees the same resampled groups and the comparisons keep the
##      correlation they actually have -- which is what makes a joint statement
##      possible, and is why the groups cannot be resampled one comparison at a
##      time;
##   2. standardise each comparison by its own bootstrap standard error, so a
##      wide comparison and a narrow one are on one scale;
##   3. take the largest absolute standardised value across comparisons within
##      that replicate, and use the `conf` quantile of that maximum as a single
##      critical value for all of them.
##
## That is the studentized-maximum, or single-step "max-t", construction.
## Against Tukey's range test it assumes neither normality nor a common
## variance, which is the reason to be bootstrapping in the first place.

ILM_ADJUST <- c("max_t", "bonferroni", "none")

## One interval from one bootstrap series. Factored out so that the adjusted and
## unadjusted paths below cannot drift apart: with several comparisons the
## switch would otherwise be written once per path and applied m times.
#' @keywords internal
#' @noRd
ilm_boot_lu <- function(th, obs, conf, ci_type) {
  a2 <- (1 - conf) / 2
  switch(ci_type,
    percentile = unname(stats::quantile(th, c(a2, 1 - a2), names = FALSE)),
    basic      = c(2 * obs - stats::quantile(th, 1 - a2, names = FALSE),
                   2 * obs - stats::quantile(th, a2, names = FALSE)),
    normal     = obs + c(-1, 1) * stats::qnorm(1 - a2) * stats::sd(th),
    ## a difference between two samples has no single jackknife series, so BCa
    ## uses the bias correction alone; stated rather than silently approximated
    bca = {
      prop <- mean(th < obs)
      if (prop <= 0 || prop >= 1)
        unname(stats::quantile(th, c(a2, 1 - a2), names = FALSE))
      else {
        z0 <- stats::qnorm(prop); z <- stats::qnorm(c(a2, 1 - a2))
        unname(stats::quantile(th, stats::pnorm(2 * z0 + z), names = FALSE))
      }
    })
}

## Which pairs to compare: every pair, or every level against one reference.
#' @keywords internal
#' @noRd
ilm_boot_pairs <- function(lv, ref = NULL) {
  if (is.null(ref)) {
    ij <- utils::combn(length(lv), 2L)
    return(list(i = ij[1, ], j = ij[2, ]))
  }
  if (length(ref) != 1L || is.na(ref))
    stop("`ref` must be a single level of the grouping variable", call. = FALSE)
  k <- match(as.character(ref), lv)
  if (is.na(k))
    stop("`ref` (", ref, ") is not a level of the grouping variable. ",
         "Levels are: ", paste(lv, collapse = ", "), ".", call. = FALSE)
  ## the reference is the "from" side, so a positive difference means the other
  ## level sits above it
  list(i = rep(k, length(lv) - 1L), j = setdiff(seq_along(lv), k))
}

#' Bootstrap intervals for differences between groups
#'
#' Reports each difference with its own interval, rather than a set of intervals
#' to compare by eye: judging a difference by whether separate intervals overlap
#' is conservative and lossy.
#'
#' With more than two groups every pair is compared, and the intervals are
#' **simultaneous by default**. Ten comparisons each at a nominal 95% do not
#' jointly cover at 95%, and reporting them as though they did is the usual way
#' a pairwise table misleads. `adjust` chooses how that is handled:
#'
#' * `"max_t"` (default) resamples all groups together, standardises each
#'   comparison by its own bootstrap standard error, and takes the `conf`
#'   quantile of the largest absolute standardised value across comparisons as
#'   one critical value for all of them. This is the single-step studentized
#'   maximum. It assumes neither normality nor a common variance, which is what
#'   separates it from Tukey's range test. Its intervals are symmetric about the
#'   observed difference, so `ci_type` does not apply to them.
#' * `"bonferroni"` builds each interval at level `1 - (1 - conf) / m` in the
#'   requested `ci_type`. More conservative than `"max_t"`, but it keeps the
#'   shape of a percentile or BCa interval, which matters on skewed data.
#' * `"none"` builds each interval at `conf`. Correct for one comparison,
#'   optimistic for several; use it when the comparisons were chosen in advance.
#'
#' With a single comparison there is nothing to adjust, and the `adjust` column
#' of the result reads `"none"` whatever was asked for.
#'
#' `p_value` and `p_adj` are built the same way the interval is, so the two
#' cannot contradict each other: under `"max_t"` both come off the studentized
#' maximum, and otherwise both come off the bootstrap distribution directly.
#' They are read off the same `R` replicates from opposite directions -- a
#' quantile and a tail proportion -- so a comparison sitting within one
#' replicate of the critical value can still land either side. Raise `R` if a
#' result is that close to the line. With `ci_type` of `"basic"`, `"normal"` or
#' `"bca"` the correspondence is close rather than exact, since those reshape
#' the same replicates while the p-value does not.
#'
#' @section What this is calibrated for:
#'
#' Against `stats::TukeyHSD()` on normal, equal-variance, balanced data -- the
#' case Tukey is exactly right for -- the endpoints agree to 0.006 and the
#' adjusted p-values to 0.011. Family-wise error over six comparisons at group
#' sizes of 45 to 80 came to 0.068 on normal data with a common variance and
#' 0.072 with variances differing fourfold, against a nominal 0.05.
#'
#' Heavy skew at small n is the case to know about. On lognormal data with a
#' spread parameter up to 1.2 and those same group sizes, a SINGLE unadjusted
#' comparison already erred 0.090 of the time, and no `ci_type` moved it
#' (percentile 0.090, BCa 0.089, basic 0.093, normal 0.084). That is the
#' bootstrapped mean of a heavily skewed small sample, not the multiplicity
#' adjustment, which is still doing its work: those six comparisons reject
#' 0.302 of the time unadjusted against 0.126 under `"max_t"`. It converges --
#' at group sizes of 300 to 500 the same design gives 0.049 per comparison and
#' 0.062 family-wise.
#'
#' `stat = "median"` is the remedy when the data look like that. On the same
#' lognormal design, under a null placing the medians rather than the means
#' together, it gave 0.043 per comparison and 0.050 family-wise.
#'
#' @param x A data frame, or a formula `y ~ g` with `data` supplied.
#' @param y Name of the numeric column. Ignored when `x` is a formula.
#' @param group Name of the grouping column, or several to cross with
#'   [interaction()]. Ignored when `x` is a formula.
#' @param data The data frame, when `x` is a formula.
#' @param stat `"mean"`, `"median"`, `"sd"`, `"var"`, or a function taking a
#'   numeric vector.
#' @param R Bootstrap replicates.
#' @param conf Confidence level. Under `"max_t"` and `"bonferroni"` this is the
#'   level of the whole family of comparisons, not of one interval.
#' @param ci_type `"percentile"`, `"bca"`, `"normal"` or `"basic"`. Applies to
#'   `"bonferroni"` and `"none"`; `"max_t"` has its own construction.
#' @param adjust `"max_t"`, `"bonferroni"` or `"none"`; see Details.
#' @param ref Optional level to compare every other level against, giving
#'   `J - 1` comparisons instead of `J * (J - 1) / 2`.
#' @param seed Random seed.
#' @param ... Passed between methods.
#' @return One row per comparison, with `from`, `to`, `observed`, `lower`,
#'   `upper`, `p_value`, `p_adj`, `excludes_zero` and the settings used.
#'   Differences are `to` minus `from`.
#' @references
#' Efron, B. and Tibshirani, R. J. (1993). An Introduction to the Bootstrap.
#' Chapman and Hall.
#'
#' Westfall, P. H. and Young, S. S. (1993). Resampling-Based Multiple Testing.
#' Wiley.
#' @seealso [ilm_boot_ci()] for a single group.
#' @examples
#' d <- ilm_sim()
#' ## all pairs, simultaneous by default
#' ilm_boot_diff(d, "score", "grp", R = 300, seed = 1)
#' ## the same thing through a formula
#' ilm_boot_diff(score ~ grp, data = d, R = 300, seed = 1)
#' ## every level against one reference, on medians
#' ilm_boot_diff(score ~ grp, data = d, stat = "median", ref = "alpha",
#'               R = 300, seed = 1)
#' @export
ilm_boot_diff <- function(x, ...) UseMethod("ilm_boot_diff")

#' @rdname ilm_boot_diff
#' @export
ilm_boot_diff.formula <- function(x, data = NULL, ...) {
  if (length(x) != 3L)
    stop("the formula needs a left and a right hand side, as in y ~ g",
         call. = FALSE)
  if (is.null(data))
    stop("`data` must be supplied with a formula, as in ",
         "ilm_boot_diff(y ~ g, data = d)", call. = FALSE)
  if (!is.data.frame(data))
    stop("`data` must be a data frame; it is ", class(data)[1], call. = FALSE)
  ## na.pass so the per-group missingness reporting below stays the caller's to
  ## see, rather than rows vanishing here
  mf <- stats::model.frame(x, data = data, na.action = stats::na.pass)
  if (ncol(mf) < 2L)
    stop("the right hand side must name at least one grouping variable",
         call. = FALSE)
  ilm_boot_diff.data.frame(mf, y = names(mf)[1L], group = names(mf)[-1L], ...)
}

#' @rdname ilm_boot_diff
#' @export
#' @param progress Show a progress bar. Defaults to [interactive()], so a
#'   bar appears when someone is watching and nothing is written in a
#'   script or a knitted document. See [ilm_progress_arg].
ilm_boot_diff.data.frame <- function(x, y = NULL, group = NULL, stat = "mean",
                                     R = 2000L, conf = 0.95,
                                     ci_type = "percentile",
                                     adjust = "max_t", ref = NULL,
                                     seed = NULL, progress = NULL, ...) {
  data <- x
  ## `...` exists for the generic and for the formula method to pass through;
  ## anything still sitting in it here is a typo, and silence would hide it
  dots <- list(...)
  if (length(dots))
    stop("unused argument(s): ",
         paste(if (is.null(names(dots))) "" else names(dots), collapse = ", "),
         call. = FALSE)
  if (is.null(y) || is.null(group))
    stop("`y` and `group` must name columns, or pass a formula as in ",
         "ilm_boot_diff(y ~ g, data = d)", call. = FALSE)
  if (length(ci_type) != 1L || !ci_type %in% ILM_CI_TYPES)
    stop("unknown `ci_type`: ", paste(sQuote(ci_type), collapse = ", "),
         ". Options are ", paste(sQuote(ILM_CI_TYPES), collapse = ", "), ".",
         call. = FALSE)
  if (length(adjust) != 1L || !adjust %in% ILM_ADJUST)
    stop("unknown `adjust`: ", paste(sQuote(adjust), collapse = ", "),
         ". Options are ", paste(sQuote(ILM_ADJUST), collapse = ", "), ".",
         call. = FALSE)
  if (!is.numeric(conf) || length(conf) != 1L || conf <= 0 || conf >= 1)
    stop("`conf` must be a single number strictly between 0 and 1", call. = FALSE)
  if (!is.numeric(R) || length(R) != 1L || R < 2)
    stop("`R` must be a single number of at least 2", call. = FALSE)
  if (length(y) != 1L)
    stop("`y` must name exactly one column", call. = FALSE)
  miss <- setdiff(c(y, group), names(data))
  if (length(miss))
    stop("column(s) not found in the data: ", paste(miss, collapse = ", "),
         ". Available: ", paste(utils::head(names(data), 12), collapse = ", "),
         call. = FALSE)
  if (!is.numeric(data[[y]]))
    stop("`y` (", y, ") must be numeric; it is ", class(data[[y]])[1],
         call. = FALSE)
  R <- as.integer(R)
  gname <- paste(group, collapse = ".")
  gf <- if (length(group) == 1L) factor(data[[group]])
        else interaction(data[group], drop = TRUE, sep = ".")
  ## an unused level would otherwise become an empty group and a column of NaN
  gf <- droplevels(gf)
  lv <- levels(gf)
  if (length(lv) < 2L)
    stop("`group` (", gname, ") needs at least 2 levels to compare; it has ",
         length(lv),
         if (length(lv) == 1L) paste0(" (", lv, ")") else "", ".", call. = FALSE)

  sf <- ilm_stat_fun(stat)
  vals <- data[[y]]
  keep_row <- !is.na(gf) & !is.na(vals)
  grps <- split(vals[keep_row], droplevels(gf[keep_row]))
  ## a level present only on rows with a missing y is gone from `grps`, so ask
  ## for the counts by level rather than by position
  ng <- vapply(lv, function(l) length(grps[[l]]), 1L)
  thin <- lv[ng < 2L]
  if (length(thin))
    stop("each group needs at least 2 non-missing observations; ",
         paste(sprintf("%s: %d", thin, ng[thin]), collapse = ", "),
         ". Drop those levels, or use ilm_boot_ci() to describe them.",
         call. = FALSE)

  pr <- ilm_boot_pairs(lv, ref)
  m  <- length(pr$i)
  ## with one comparison there is no family, so nothing is adjusted whatever
  ## was asked for -- and the reported `adjust` says what was actually done
  if (m == 1L) adjust <- "none"

  if (!is.null(seed)) set.seed(seed)
  ## Every group resampled once per replicate: the rows of TH are joint draws,
  ## which is what lets the maximum below be a statement about all comparisons
  ## at once.
  pb <- ilm_progress(length(lv) * R, progress); .tk <- 0L
  TH <- vapply(lv, function(l) {
    v <- grps[[l]]; n <- length(v)
    ## per replicate, for the reason given in ilm_boot_stat(): the all-at-once
    ## index matrix is an allocation proportional to n * R
    z <- vapply(seq_len(R), function(r) {
      .tk <<- .tk + 1L; pb$tick(.tk)
      sf(v[sample.int(n, n, replace = TRUE)]) }, 1)
    z
  }, numeric(R))
  pb$done()
  obs_g <- vapply(lv, function(l) sf(grps[[l]]), 1)

  D  <- TH[, pr$j, drop = FALSE] - TH[, pr$i, drop = FALSE]   # R x m
  dh <- unname(obs_g[pr$j] - obs_g[pr$i])
  ## a replicate counts only if every comparison in it is finite, since a
  ## maximum taken over a shifting set of comparisons is not a maximum
  usable <- apply(is.finite(D), 1L, all)
  D <- D[usable, , drop = FALSE]
  if (!nrow(D))
    stop("no bootstrap replicate produced a finite value for every comparison",
         call. = FALSE)

  se <- apply(D, 2L, stats::sd)
  ## a degenerate comparison has no scale to standardise by and cannot enter
  ## the maximum; its own interval still comes back, as a point
  scal <- is.finite(se) & se > 0
  Z <- sweep(sweep(D, 2L, dh, "-"), 2L, ifelse(scal, se, 1), "/")
  Mr <- if (any(scal)) apply(abs(Z[, scal, drop = FALSE]), 1L, max)
        else rep(0, nrow(D))
  zobs <- ifelse(scal, abs(dh) / se, 0)

  lo <- up <- numeric(m)
  if (adjust == "max_t") {
    q <- stats::quantile(Mr, conf, names = FALSE)
    lo <- dh - q * se; up <- dh + q * se
  } else {
    cf <- if (adjust == "bonferroni") 1 - (1 - conf) / m else conf
    for (k in seq_len(m)) {
      lu <- ilm_boot_lu(D[, k], dh[k], cf, ci_type)
      lo[k] <- lu[1]; up[k] <- lu[2]
    }
  }
  ## The p-value is built the same way the interval is, so the two cannot
  ## contradict each other -- a table saying p = 0.03 beside an interval
  ## covering zero is the kind of disagreement this package treats as a bug.
  ## Under max_t both come off the studentized maximum. Otherwise both come off
  ## the bootstrap distribution itself: a percentile interval at level 1 - a
  ## excludes zero exactly when less than a / 2 of the replicates lie on the far
  ## side of it, which is this p-value. For `basic`, `normal` and `bca` the
  ## duality is close rather than exact, since those reshape the same
  ## replicates.
  pv <- if (adjust == "max_t")
    vapply(seq_len(m), function(k) mean(abs(Z[, k]) >= zobs[k]), 1)
  else
    vapply(seq_len(m), function(k)
      min(1, 2 * min(mean(D[, k] <= 0), mean(D[, k] >= 0))), 1)
  padj <- switch(adjust,
    max_t      = vapply(zobs, function(z) mean(Mr >= z), 1),
    bonferroni = pmin(1, m * pv),
    none       = pv)

  ## The proportion of replicates on each side of zero. Not a p-value and not
  ## a posterior probability -- it is a description of where the bootstrap
  ## distribution sits, which is what ilm_plot_boot_diff() draws.
  p_sup <- vapply(seq_len(m), function(k) mean(D[, k] > 0), 1)

  out <- data.frame(stat = if (is.function(stat)) "custom" else stat,
             group = gname, from = lv[pr$i], to = lv[pr$j],
             observed = dh, lower = lo, upper = up,
             conf = conf, R = R, ci_type = if (adjust == "max_t") "max_t" else ci_type,
             adjust = adjust, n_comparisons = m,
             n_from = unname(ng[pr$i]), n_to = unname(ng[pr$j]),
             p_value = pv, p_adj = padj, p_superiority = round(p_sup, 4),
             excludes_zero = is.finite(lo) & (lo > 0 | up < 0),
             stringsAsFactors = FALSE, row.names = NULL)
  ## the replicate differences themselves, one column per comparison, so the
  ## distribution behind a row can be drawn rather than only summarised
  colnames(D) <- paste(lv[pr$i], lv[pr$j], sep = " -> ")
  attr(out, "draws") <- D
  out
}

## ---- missingness -----------------------------------------------------------

#' Missingness in one variable
#'
#' @param data A data frame, or a vector when `y` is `NULL`.
#' @param y Name of the column.
#' @param by Optional grouping columns.
#' @param digits Rounding for `p_na`.
#' @return A data frame with `obs`, `n`, `na` and `p_na`.
#' @seealso [ilm_describe_na_all()], [ilm_plot_missing()].
#' @examples
#' ilm_describe_na(ilm_sim(), "lab_value")
#' @export
ilm_describe_na <- function(data, y = NULL, by = NULL, digits = 4) {
  if (is.null(y) && !is.data.frame(data)) { v <- data; data <- NULL } else {
    if (is.null(y)) stop("`y` must name a column, or pass a vector", call. = FALSE)
    miss <- setdiff(c(y, by), names(data))
    if (length(miss))
      stop("column(s) not found in the data: ", paste(miss, collapse = ", "),
           ". Available: ", paste(utils::head(names(data), 12), collapse = ", "),
           call. = FALSE)
    v <- data[[y]]
  }
  one <- function(vec) {
    nn <- fnobs(vec)
    data.frame(obs = length(vec), n = nn, na = length(vec) - nn,
               p_na = round((length(vec) - nn) / length(vec), digits),
               stringsAsFactors = FALSE)
  }
  if (is.null(by) || is.null(data)) return(one(v))
  g <- interaction(data[by], drop = TRUE)
  parts <- lapply(split(v, g), one)
  cbind(setNames(data.frame(names(parts), stringsAsFactors = FALSE),
                 paste(by, collapse = ".")),
        do.call(rbind, parts), row.names = NULL)
}

#' Missingness in every variable
#'
#' Sorted with the most missing first, since those are the variables worth
#' looking at.
#'
#' @inheritParams ilm_describe_na
#' @param sort Sort by proportion missing, descending.
#' @return A data frame with `variable`, `obs`, `n`, `na` and `p_na`.
#' @examples
#' ilm_describe_na_all(ilm_sim())
#' @export
ilm_describe_na_all <- function(data, by = NULL, digits = 4, sort = TRUE) {
  if (!is.data.frame(data)) stop("`data` must be a data frame", call. = FALSE)
  miss <- setdiff(by, names(data))
  if (length(miss))
    stop("`by` variable(s) not found in the data: ", paste(miss, collapse = ", "),
         call. = FALSE)
  cols <- setdiff(names(data), by)
  res <- do.call(rbind, lapply(cols, function(cn) {
    r <- ilm_describe_na(data, cn, by = by, digits = digits)
    cbind(variable = cn, r, stringsAsFactors = FALSE)
  }))
  ## the variables with the most missingness are the ones worth looking at, so
  ## they go first unless the caller wants the original column order
  if (sort) res <- res[order(-res$p_na, res$variable), , drop = FALSE]
  rownames(res) <- NULL
  res
}
