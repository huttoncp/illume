## The interval of a typical group's survival curve in ilm_survival(), for every
## family. The linear predictor of a typical group is linear in the fixed
## effects and in the penalised coefficients of the smooths that are not
## factor-smooth interactions, so its interval is the delta method on the
## linear predictor, carried through the survivor function, with the standard
## error read from the fit's covariance. No draws are used, so the interval is
## the same on every run. Each interval here is checked against a reference
## built from the whole dense inverse of the joint precision, by hand.

ti_data <- function(seed = 41, n = 250, gsd = 0.4) {
  set.seed(seed)
  d <- data.frame(x = stats::runif(n), z = stats::rnorm(n),
                  h = factor(sample(letters[1:8], n, TRUE)))
  hl <- stats::rnorm(8, 0, gsd)
  eta <- 1 + 0.5 * d$z + 1.2 * sin(2 * pi * d$x) + hl[as.integer(d$h)]
  tt <- exp(eta + 0.5 * log(stats::rexp(n)))
  ct <- stats::rexp(n, 1 / (2 * stats::median(tt)))
  d$time <- pmin(tt, ct); d$event <- as.integer(tt <= ct)
  d
}
ti_fit <- function(d, family, form, reml = FALSE, ...)
  suppressMessages(suppressWarnings(
    ilm_model(form, data = d, family = family, censor = ilm_surv(d$time, d$event),
              reml = reml, verbose = FALSE, ...)))
ti_nd <- function(d, x = c(0.25, 0.5, 0.75), z = c(-1, 0, 1))
  data.frame(x = x, z = z, h = factor("a", levels = levels(d$h)))
ti_has_smooth <- function(form) grepl("s(", paste(deparse(form), collapse = ""), fixed = TRUE)
ti_surv <- function(nm, z)
  switch(nm, weibull = exp(-exp(z)), lognormal = stats::pnorm(-z), loglogistic = 1 / (1 + exp(z)))

## the smooth terms a typical curve carries: not a factor-smooth interaction's
ti_used <- function(f)
  which(vapply(f$re, function(r) identical(r$kind, "basis"), TRUE) & !illume:::ilm_fs_terms(f))

## The interval by hand. A row's linear predictor is a' theta, with theta the
## fit's whole parameter vector and `a` zero except at the row's fixed design
## and at the row's smooth design (the baseline spline's too, for a flexible
## model, at each time). Its variance is a' Q^-1 a with the DENSE inverse of the
## whole joint precision, where the package solves for the columns it needs;
## with no smooth in the curve it is vcov(), the fixed effects' covariance.
ti_ref <- function(f, nd, times, conf = 0.95) {
  crit <- stats::qnorm(1 - (1 - conf) / 2)
  ndx <- illume:::ilm_newX(f, nd, fs_fill = TRUE)
  used <- ti_used(f)
  rp <- f$rp; p <- ncol(f$X)
  keep <- if (is.null(rp)) seq_len(p) else setdiff(seq_len(p), rp$cols)
  rn <- names(f$obj$env$par); pb <- which(rn == "beta"); pv <- which(rn == "bvec")
  Vfull <- if (length(used)) solve(as.matrix(illume:::ilm_joint_prec(f))) else NULL
  Vb <- if (!length(used)) as.matrix(vcov(f))
  se_of <- function(a_beta, a_smooth) {
    if (length(used)) {
      a <- numeric(nrow(Vfull)); a[pb] <- a_beta
      for (k in used) a[pv[f$b_idx[[k]]]] <- a_smooth[[as.character(k)]]
      sqrt(as.numeric(t(a) %*% Vfull %*% a))
    } else sqrt(as.numeric(t(a_beta) %*% Vb %*% a_beta))
  }
  out <- list()
  for (i in seq_len(nrow(nd))) {
    sm <- lapply(used, function(k) as.matrix(illume:::ilm_basis_new(f, ndx, k))[i, ])
    names(sm) <- as.character(used)
    off <- sum(ndx$X[i, ] * f$beta[keep, 1L]) +
      sum(vapply(used, function(k) sum(sm[[as.character(k)]] *
                    as.vector(illume:::ilm_Bhat_term(f, k))), 0))
    if (is.null(rp)) {
      ab <- ndx$X[i, ]
      se <- se_of(ab, sm)
      sc <- unname(f$dispersion[[1]]); nm <- f$family$name
      out[[i]] <- data.frame(row = i, time = times, se = se, scale = sc,
        surv = ti_surv(nm, (log(times) - off) / sc),
        lower = ti_surv(nm, (log(times) - (off - crit * se)) / sc),
        upper = ti_surv(nm, (log(times) - (off + crit * se)) / sc))
    } else {
      Bs <- illume:::ilm_rcs(log(times), rp$knots)
      e <- as.vector(Bs %*% f$beta[rp$cols, 1L]) + off
      se <- vapply(seq_along(times), function(j) {
        ab <- numeric(p); ab[rp$cols] <- Bs[j, ]; ab[keep] <- ndx$X[i, ]
        se_of(ab, sm)
      }, 0)
      Sf <- f$family$surv
      out[[i]] <- data.frame(row = i, time = times, surv = Sf(e),
                             lower = Sf(e + crit * se), upper = Sf(e - crit * se))
    }
  }
  do.call(rbind, out)
}
ti_typical <- function(f, nd, times, ...)
  suppressWarnings(ilm_survival(f, newdata = nd, times = times, groups = "typical", ...))

test_that("the covariance block formed from the sparse factor is the dense inverse's block", {
  d <- ti_data()
  for (reml in c(TRUE, FALSE)) {
    f <- ti_fit(d, "rp", time ~ s(x, k = 8) + z + (1 | h), reml = reml)
    rn <- names(f$obj$env$par)
    bidx <- unlist(f$b_idx[ti_used(f)])
    pos <- c(which(rn == "beta"), which(rn == "bvec")[bidx])
    V <- illume:::ilm_joint_cov(f, pos)
    Vd <- solve(as.matrix(illume:::ilm_joint_prec(f)))[pos, pos]
    expect_equal(unname(V), unname(Vd), tolerance = 1e-9, label = paste("reml =", reml))
    ## and one coordinate alone, in the middle of the vector
    expect_equal(unname(illume:::ilm_joint_cov(f, pos[5])), unname(Vd[5, 5, drop = FALSE]),
                 tolerance = 1e-9)
  }
})

test_that("by maximum likelihood the fixed effects' block of the joint covariance is vcov()", {
  ## so a curve with no smooth in it, which reads vcov(), reads the same block
  d <- ti_data(42)
  f <- ti_fit(d, "weibull", time ~ z + (1 | h), joint = TRUE)
  pb <- which(names(f$obj$env$par) == "beta")
  expect_equal(unname(illume:::ilm_joint_cov(f, pb)), unname(as.matrix(vcov(f))),
               tolerance = 1e-8)
})

test_that("a direction held within the joint distribution is taken out of the covariance as it is out of a draw", {
  ## no group variation in the data, a group term fitted anyway: its variance is
  ## held at the boundary, and ilm_draws() conditions on the held direction
  d <- ti_data(45, gsd = 0)
  f <- ti_fit(d, "weibull", time ~ s(x, k = 8) + z + (1 | h), reml = FALSE)
  skip_if(is.null(f$hessian_dirs) || !length(f$hessian_held),
          "this fit did not end at the boundary on this platform")
  rn <- names(f$obj$env$par)
  bidx <- unlist(f$b_idx[ti_used(f)])
  V <- illume:::ilm_joint_cov(f, c(which(rn == "beta"), which(rn == "bvec")[bidx]))
  dr <- ilm_draws(f, nsim = 20000L, seed = 3L, blocks = c("beta", "bvec"), natural = FALSE)
  w <- rownames(dr$draws)
  Ve <- stats::cov(t(dr$draws[c(which(w == "beta"), which(w == "bvec")[bidx]), ]))
  s <- sqrt(diag(V))
  ## in units of the coordinates' own SDs; 20000 draws put the noise near 0.01
  expect_lt(max(abs(Ve - V) / outer(s, s)), 0.05)
})

test_that("the design the interval uses reproduces the typical linear predictor's smooth part", {
  ## the smooths' design at the rows, times a draw of their coefficients, is
  ## ilm_smooth_part(typical = TRUE): a factor-smooth's blocks are not in it
  d <- ti_data(43)
  for (fam in c("weibull", "rp")) {
    f <- ti_fit(d, fam, time ~ s(x, k = 5) + z + s(x, h, bs = "fs", k = 4), reml = TRUE)
    expect_true(any(illume:::ilm_fs_terms(f)))
    nd <- ti_nd(d)
    ndx <- illume:::ilm_newX(f, nd, fs_fill = TRUE)
    tc <- illume:::ilm_typical_cov(f, ndx)
    bv <- stats::rnorm(length(f$sdr$par.random[names(f$sdr$par.random) == "bvec"]))
    bidx <- unlist(f$b_idx[ti_used(f)])
    expect_equal(as.vector(tc$Zs %*% bv[bidx]),
                 illume:::ilm_smooth_part(f, ndx, bv, typical = TRUE)[, 1L],
                 tolerance = 1e-10, label = fam)
    expect_equal(ncol(tc$Zs), length(bidx))
    expect_equal(dim(tc$V), rep(ncol(f$X) + length(bidx), 2L))
  }
})

test_that("an accelerated failure time typical interval is the delta method on the linear predictor, for every family", {
  d <- ti_data(); nd <- ti_nd(d); tg <- c(0.5, 2, 6)
  for (fam in c("weibull", "lognormal", "loglogistic")) {
    for (form in list(time ~ z, time ~ z + (1 | h), time ~ s(x, k = 8) + z + (1 | h))) {
      f <- ti_fit(d, fam, form, reml = ti_has_smooth(form))
      s <- ti_typical(f, nd, tg)
      r <- ti_ref(f, nd, tg)
      lab <- paste(fam, deparse(form))
      expect_equal(s$surv, r$surv, tolerance = 1e-10, label = lab)
      expect_equal(s$lower, r$lower, tolerance = 1e-7, label = paste("lower,", lab))
      expect_equal(s$upper, r$upper, tolerance = 1e-7, label = paste("upper,", lab))
    }
  }
})

test_that("for the Weibull the interval is today's complementary log-log formula, given the same standard error", {
  ## log(-log S) is linear in the linear predictor there, so shifting it by
  ## crit * se / scale and evaluating S at eta -/+ crit * se are the same
  ## thing; the lognormal and log-logistic are not, which is why they move
  d <- ti_data(); nd <- ti_nd(d); tg <- c(0.5, 2, 6)
  for (form in list(time ~ z, time ~ s(x, k = 8) + z + (1 | h))) {
    f <- ti_fit(d, "weibull", form, reml = ti_has_smooth(form))
    s <- ti_typical(f, nd, tg)
    r <- ti_ref(f, nd, tg)
    cll <- log(-log(pmin(pmax(r$surv, 1e-12), 1 - 1e-12)))
    half <- stats::qnorm(0.975) * r$se / r$scale
    keep <- r$surv > 1e-12 & r$surv < 1 - 1e-12            # where the old clamp did nothing
    expect_true(any(keep))
    expect_equal(s$lower[keep], exp(-exp(cll + half))[keep], tolerance = 1e-9)
    expect_equal(s$upper[keep], exp(-exp(cll - half))[keep], tolerance = 1e-9)
  }
})

test_that("a flexible-parametric typical interval is the delta method on the linear predictor, with and without a smooth", {
  d <- ti_data(44); nd <- ti_nd(d); tg <- c(0.5, 2, 6)
  for (reml in c(FALSE, TRUE)) {
    for (form in list(time ~ z, time ~ z + (1 | h), time ~ s(x, k = 8) + z + (1 | h))) {
      ## REML has nothing to estimate in a model with no random part at all
      f <- ti_fit(d, "rp", form, reml = reml &&
                    (ti_has_smooth(form) || grepl("|", paste(deparse(form), collapse = ""), fixed = TRUE)))
      s <- ti_typical(f, nd, tg)
      r <- ti_ref(f, nd, tg)
      lab <- paste(deparse(form), "reml =", reml)
      expect_equal(s$surv, r$surv, tolerance = 1e-10, label = lab)
      expect_equal(s$lower, r$lower, tolerance = 1e-7, label = paste("lower,", lab))
      expect_equal(s$upper, r$upper, tolerance = 1e-7, label = paste("upper,", lab))
    }
  }
})

test_that("with only a factor-smooth interaction the fixed effects' covariance gives the interval, and the term is at zero", {
  d <- ti_data(46)
  for (fam in c("weibull", "rp")) {
    f <- ti_fit(d, fam, time ~ z + s(x, h, bs = "fs", k = 4), reml = TRUE)
    expect_equal(length(ti_used(f)), 0L)
    nd <- ti_nd(d); tg <- c(1, 4)
    s <- ti_typical(f, nd, tg)
    r <- ti_ref(f, nd, tg)
    expect_equal(s$lower, r$lower, tolerance = 1e-7, label = fam)
    expect_equal(s$upper, r$upper, tolerance = 1e-7, label = fam)
  }
})

test_that("the flexible interval with no smooth is exactly what it was", {
  ## numbers from the previous release's code (maximum likelihood and REML),
  ## on this data, to the digits printed: the path is the same delta method on
  ## vcov(), so nothing moves
  d <- ti_data(41)
  nd <- data.frame(x = c(0.3, 0.7), z = c(-0.5, 1), h = factor("a", levels = levels(d$h)))
  want <- list(
    ml_re = list(
      surv  = c(0.634646471761, 0.338946545109, 0.795435861692, 0.580087263869),
      lower = c(0.557979729831, 0.259131524736, 0.723238767675, 0.476209604294),
      upper = c(0.701628070292, 0.420296891565, 0.850733524873, 0.670496304484)),
    reml_re = list(
      surv  = c(0.632905944259, 0.335831622620, 0.793230736956, 0.575481351922),
      lower = c(0.552719354496, 0.251522680199, 0.718444689081, 0.468084832113),
      upper = c(0.702636521340, 0.422057596642, 0.850210711393, 0.668849193550)),
    ml_plain = list(
      surv  = c(0.631694539030, 0.337651824848, 0.796603146455, 0.584212591446),
      lower = c(0.563023726150, 0.267605678200, 0.728271934447, 0.487761016972),
      upper = c(0.692586912427, 0.408915840311, 0.849521565924, 0.668712326717)))
  fits <- list(ml_re = ti_fit(d, "rp", time ~ z + (1 | h)),
               reml_re = ti_fit(d, "rp", time ~ z + (1 | h), reml = TRUE),
               ml_plain = ti_fit(d, "rp", time ~ z))
  for (k in names(want)) {
    s <- ti_typical(fits[[k]], nd, c(1, 3))
    for (col in c("surv", "lower", "upper"))
      expect_equal(s[[col]], want[[k]][[col]], tolerance = 1e-6, label = paste(k, col))
  }
})

test_that("nsim and seed do not reach a typical interval, in any family", {
  d <- ti_data(47, n = 200); nd <- ti_nd(d); tg <- c(1, 4)
  forms <- list(time ~ z, time ~ z + (1 | h), time ~ s(x, k = 6) + z + (1 | h))
  for (fam in c("weibull", "lognormal", "loglogistic", "rp")) for (form in forms) {
    f <- ti_fit(d, fam, form, reml = ti_has_smooth(form))
    a <- ti_typical(f, nd, tg, nsim = 10L, seed = 1L)
    b <- ti_typical(f, nd, tg, nsim = 500L, seed = 99L)
    expect_identical(a, b, label = paste(fam, deparse(form)))
  }
  ## and the user's random stream is not touched by it
  f <- ti_fit(d, "weibull", time ~ s(x, k = 6) + z)
  set.seed(5); u <- stats::runif(1)
  set.seed(5); invisible(ti_typical(f, nd, tg)); expect_equal(stats::runif(1), u)
})

test_that("a typical interval agrees with many joint draws on the linear predictor's scale, at two seeds", {
  ## The draws' own SD and 2.5% and 97.5% points, on the scale on which the
  ## curve is linear in the linear predictor (log(-log S), Weibull and the
  ## flexible hazard scale), in units of the interval's half-width over 1.96.
  ## On the curve's own scale the far tail magnifies the draws' noise, and a
  ## tolerance there is either loose or flaky.
  d <- ti_data(48); nd <- ti_nd(d, x = c(0.2, 0.5, 0.8), z = c(-1, 0, 1)); tg <- c(1, 4)
  cll <- function(S) log(-log(S))
  for (fam in c("weibull", "rp")) {
    f <- ti_fit(d, fam, time ~ s(x, k = 8) + z + (1 | h), reml = TRUE)
    if (!is.null(f$rp)) expect_identical(f$family$rp_scale, "hazard")
    s <- ti_typical(f, nd, tg)
    ndx <- illume:::ilm_newX(f, nd, fs_fill = TRUE)
    rp <- f$rp; p <- ncol(f$X)
    keep <- if (is.null(rp)) seq_len(p) else setdiff(seq_len(p), rp$cols)
    Bs <- if (!is.null(rp)) illume:::ilm_rcs(log(tg), rp$knots)
    sc <- if (is.null(rp)) unname(f$dispersion[[1]])
    ## the linear predictor (log(-log S)) at each row and time, for one draw
    lp <- function(b, bv) {
      o <- as.vector(ndx$X %*% b[keep]) + illume:::ilm_smooth_part(f, ndx, bv, typical = TRUE)[, 1L]
      unlist(lapply(o, function(oi) if (is.null(rp)) (log(tg) - oi) / sc
                                    else as.vector(Bs %*% b[rp$cols]) + oi))
    }
    half <- (cll(s$lower) - cll(s$upper)) / 2          # log(-log) falls as S rises
    est <- cll(s$surv)
    for (seed in c(11L, 12L)) {
      dr <- ilm_draws(f, nsim = 4000L, seed = seed, blocks = c("beta", "bvec"), natural = FALSE)
      w <- rownames(dr$draws)
      L <- sapply(seq_len(ncol(dr$draws)), function(j) lp(dr$draws[w == "beta", j], dr$draws[w == "bvec", j]))
      sdv <- apply(L, 1, stats::sd)
      lab <- paste(fam, "seed", seed)
      ## the SD against the interval's (half-width / 1.96): relative, 4000 draws put its noise near 1%
      expect_lt(max(abs(sdv / (half / stats::qnorm(0.975)) - 1)), 0.05, label = lab)
      ## the 2.5% and 97.5% points, from the estimate, in units of the SD
      q1 <- apply(L, 1, stats::quantile, probs = 0.025, names = FALSE)
      q2 <- apply(L, 1, stats::quantile, probs = 0.975, names = FALSE)
      expect_lt(max(abs((q1 - est) / sdv + stats::qnorm(0.975))), 0.2, label = lab)
      expect_lt(max(abs((q2 - est) / sdv - stats::qnorm(0.975))), 0.2, label = lab)
    }
  }
})

test_that("the limits enclose the estimate at every time, at survival near 1e-18 too", {
  ## the interval used to be formed from the estimate clamped to [1e-12, 1 - 1e-12],
  ## and at S near 1e-18 its lower limit sat above the estimate
  d <- ti_data(49); nd <- ti_nd(d); tg <- c(0.01, 1, 8, 30, 60)
  for (fam in c("weibull", "lognormal", "loglogistic", "rp")) {
    f <- ti_fit(d, fam, time ~ s(x, k = 6) + z + (1 | h), reml = TRUE)
    s <- ti_typical(f, nd, tg)
    expect_true(all(s$lower <= s$surv & s$surv <= s$upper), label = fam)
    expect_true(all(is.finite(s$lower) & is.finite(s$upper)), label = fam)
  }
  f <- ti_fit(d, "weibull", time ~ z)
  s <- ti_typical(f, ti_nd(d, z = c(-3, 0, 3)), c(40, 200))
  expect_true(any(s$surv < 1e-12 & s$surv > 0))
  expect_true(all(s$lower <= s$surv & s$surv <= s$upper))
})

test_that("an unusable covariance gives no interval and a warning in every family, never one of no width", {
  d <- ti_data(40, n = 200); nd <- ti_nd(d, x = c(0.3, 0.7), z = c(0, 1)); tg <- c(1, 3)
  for (fam in c("weibull", "lognormal", "loglogistic", "rp")) {
    ## a missing entry in the fixed effects' covariance, with no smooth in the curve
    f0 <- ti_fit(d, fam, time ~ z + (1 | h))
    f0$sdr$cov.fixed[1, 1] <- NA_real_
    expect_warning(s0 <- ilm_survival(f0, newdata = nd, times = tg, groups = "typical"),
                   "covariance is unusable")
    expect_true(all(is.na(s0$lower) & is.na(s0$upper)), label = fam)
    expect_true(all(is.finite(s0$surv)), label = fam)
    ## a joint precision of zeros, which has no Cholesky factor, with a smooth in it
    f1 <- ti_fit(d, fam, time ~ s(x, k = 6) + z + (1 | h), reml = TRUE)
    f1$jointPrecision <- f1$jointPrecision * 0
    ## the Matrix package also warns, as the Cholesky factor fails
    expect_warning(expect_warning(
      s1 <- ilm_survival(f1, newdata = nd, times = tg, groups = "typical"),
      "no intervals"), "not positive definite")
    expect_true(all(is.na(s1$lower) & is.na(s1$upper)), label = fam)
    expect_true(all(is.finite(s1$surv)), label = fam)
    ## a covariance that is not one: a negative variance beyond rounding
    f2 <- ti_fit(d, fam, time ~ z)
    f2$sdr$cov.fixed <- -f2$sdr$cov.fixed
    expect_warning(s2 <- ilm_survival(f2, newdata = nd, times = tg, groups = "typical"),
                   "negative or not finite")
    expect_true(all(is.na(s2$lower) & is.na(s2$upper)), label = fam)
    expect_true(all(is.finite(s2$surv)), label = fam)
    ## a Hessian that is not positive definite, whose covariance TMB hands back
    ## finite all the same: vcov() warns, and the curve has no interval and one
    ## warning of its own, whatever the covariance would give
    f3 <- ti_fit(d, fam, time ~ z + (1 | h))
    expect_true(all(is.finite(f3$sdr$cov.fixed)))
    f3$sdr$pdHess <- FALSE
    f3$hessian_held <- NULL
    w3 <- NULL
    s3 <- withCallingHandlers(ilm_survival(f3, newdata = nd, times = tg, groups = "typical"),
                              warning = function(w) { w3 <<- c(w3, conditionMessage(w))
                                                      invokeRestart("muffleWarning") })
    expect_length(w3, 1L)
    expect_match(w3, "covariance is unusable")
    expect_true(all(is.na(s3$lower) & is.na(s3$upper)), label = fam)
    expect_true(all(is.finite(s3$surv)), label = fam)
  }
})
