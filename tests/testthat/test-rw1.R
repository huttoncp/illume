# A random walk over time, ilm_rw1(), and the cells accessor, ilm_cells().
#
# The load-bearing tests hold a gaussian walk to the exact likelihood of the
# local-level model. With a gaussian response the Laplace approximation is
# exact, so the maximum likelihood and REML estimates, the log-likelihood and a
# forecast from the fitted walk all have closed forms to agree with. Those are
# computed here from V = s2_eta * K + s2_eps * I, where K[s, t] is the time
# both observations share since their group's first -- min(t_s, t_t) - t_first
# -- and the reference estimates were maximised from it independently, to a
# relative tolerance of 1e-15, by another agent.

## a local-level series: zero at each group's first time, then Brownian steps
rw_data <- function(time, group, beta, s2n, s2e, seed) {
  set.seed(seed)
  u <- numeric(length(time))
  for (g in unique(group)) {
    i <- which(group == g); i <- i[order(time[i])]
    u[i] <- c(0, cumsum(stats::rnorm(length(i) - 1, 0, sqrt(s2n * diff(time[i])))))
  }
  data.frame(time = time, group = group,
             y = beta + u + stats::rnorm(length(time), 0, sqrt(s2e)))
}

rw_cases <- function() list(
  uni_a = rw_data(1:60, rep("A", 60), 2, 0.5, 1, 101),
  uni_b = rw_data(1:60, rep("A", 60), 2, 0.1, 1, 102),
  uni_c = rw_data(1:60, rep("A", 60), 2, 1, 0.25, 103),
  panel = rw_data(rep(1:30, 3), rep(c("A", "B", "C"), each = 30), 2, 0.5, 1, 104),
  gaps  = local({
    set.seed(105); tt <- cumsum(c(1, sample(1:3, 49, TRUE)))
    rw_data(tt, rep("A", 50), 2, 0.5, 1, 106)
  }))

## the exact covariance of the observations, and the exact forecast of a
## group's last cell: its mean, and the variance of beta + u there
rw_dense <- function(d, s2n, s2e, g, beta = NULL) {
  n <- nrow(d); K <- matrix(0, n, n)
  for (h in unique(d$group)) {
    i <- which(d$group == h)
    K[i, i] <- outer(d$time[i], d$time[i], pmin) - min(d$time[i])
  }
  Vi <- solve(s2n * K + diag(s2e, n))
  a <- sum(Vi)
  reml <- is.null(beta)
  ## REML integrates beta out, so its mode is the GLS estimate; under ML the
  ## forecast conditions on the fitted beta
  if (reml) beta <- sum(Vi %*% d$y) / a
  i <- which(d$group == g); last <- i[which.max(d$time[i])]
  k <- s2n * K[, last]
  v <- s2n * K[last, last] - drop(crossprod(k, Vi %*% k))
  if (reml) v <- v + (1 - sum(Vi %*% k))^2 / a
  c(mean = beta + drop(crossprod(k, Vi %*% (d$y - beta))), var = v)
}

## ---- the index arithmetic --------------------------------------------------

test_that("the cells are the distinct (group, time) pairs, anchored at each group's first", {
  d <- data.frame(id = c(2, 1, 1, 1, 2, 2, 1),
                  t  = c(4, 0, 3, 3, 1, 9, 10))
  s <- ilm_rw1(d$t, d$id)
  ## group 1 has times 0, 3, 10 and group 2 has 1, 4, 9: six cells, in group
  ## then time order, and the two observations at time 3 share one
  expect_equal(s$n_cell, 6L)
  expect_equal(s$idx, c(5L, 1L, 2L, 2L, 4L, 6L, 3L))
  expect_equal(s$first, c(1L, 4L))
  expect_equal(s$rest, c(2L, 3L, 5L, 6L))
  expect_equal(s$prev, c(1L, 2L, 4L, 5L))
  expect_equal(s$gap, c(3, 7, 3, 5))
  ## the anchors are held at zero, so they are not latent values
  expect_equal(s$n_latent, 4L)
  expect_equal(s$obs_per_latent, 7 / 4)
  expect_s3_class(s, c("ilm_rw1", "ilm_cor"))
  expect_error(ilm_rw1(c(1, 1, 2), c("a", "a", "b")), "only one distinct time")
})

## ---- the exact local-level likelihood ---------------------------------------

test_that("a gaussian walk is the local-level model, by ML and by REML", {
  ## case, criterion, s2_eta, s2_eps, beta, logLik: the exact likelihood,
  ## maximised independently
  ref <- read.table(header = TRUE, text = "
    case  crit s2_eta         s2_eps         beta          loglik
    uni_a ML   0.621719903563 0.48167575929  1.49037374973 -95.0090561949
    uni_a REML 0.608067079823 0.507997380182 1.48687664258 -94.6535158536
    uni_b ML   0.221090606605 0.890272514677 2.14893296407 -95.9641210607
    uni_b REML 0.230706499874 0.898023014166 2.13864609247 -95.56949265
    uni_c ML   0.540958278012 0.56527845609  1.85326564303 -95.675155569
    uni_c REML 0.52985449874  0.591344051493 1.84052306365 -95.2814298699
    panel ML   0.335631203344 0.973451506428 1.98833702322 -151.015291567
    panel REML 0.341627891935 0.981942025225 1.98848139311 -151.067110137
    gaps  ML   0.273443205951 1.26641004399  1.4803057495  -92.1724459792
    gaps  REML 0.280428136598 1.28947641958  1.48060745017 -91.5791266499")
  cs <- rw_cases()
  for (i in seq_len(nrow(ref))) {
    r <- ref[i, ]
    f <- ilm_model(y ~ 1, data = cs[[r$case]], family = "gaussian",
                   ar = ilm_rw1(~ time | group), reml = r$crit == "REML",
                   verbose = FALSE)
    lab <- paste(r$case, r$crit)
    expect_equal(f$opt$convergence, 0L, label = lab)
    expect_equal(f$Sigma$ar[1, 1], r$s2_eta, tolerance = 2e-5, label = lab)
    expect_equal(unname(f$dispersion)^2, r$s2_eps, tolerance = 2e-5, label = lab)
    expect_equal(unname(coef(f))[1], r$beta, tolerance = 2e-5, label = lab)
    ## the likelihood itself, constants and REML's included
    expect_equal(as.numeric(logLik(f)), r$loglik, tolerance = 1e-8, label = lab)
    expect_true(is.na(f$rho))
  }
})

test_that("a forecast from the fitted walk is the exact one", {
  ## what a forecast of a group's last cell needs is its mode and the
  ## variance of the level there, and both come out of the fit: the mode from
  ## par.random, the variance from the joint precision -- under ML the
  ## random block alone, beta held at its estimate; under REML beta is in the
  ## random block, and its uncertainty is part of the level's
  cs <- rw_cases()
  for (case in c("uni_a", "panel", "gaps")) for (reml in c(FALSE, TRUE)) {
    d <- cs[[case]]
    f <- ilm_model(y ~ 1, data = d, family = "gaussian",
                   ar = ilm_rw1(~ time | group), reml = reml, joint = TRUE,
                   verbose = FALSE)
    cl <- ilm_cells(f)
    g <- if (case == "panel") "B" else "A"
    k <- cl$index[cl$last & cl$group == g]
    Q <- f$jointPrecision; rn <- rownames(Q)
    ir <- which(rn %in% c("B_ar", if (reml) "beta"))
    e <- numeric(length(ir))
    e[which(rn[ir] == "B_ar")[k]] <- 1
    if (reml) e[rn[ir] == "beta"] <- 1
    v <- drop(crossprod(e, solve(as.matrix(Q[ir, ir]), e)))
    pr <- f$sdr$par.random
    m <- unname(coef(f))[1] + unname(pr[names(pr) == "B_ar"][k])
    want <- rw_dense(d, f$Sigma$ar[1, 1], unname(f$dispersion)^2, g,
                     beta = if (reml) NULL else unname(coef(f))[1])
    lab <- paste(case, if (reml) "REML" else "ML")
    expect_equal(m, unname(want["mean"]), tolerance = 1e-8, label = lab)
    expect_equal(v, unname(want["var"]), tolerance = 1e-8, label = lab)
  }
})

## ---- the structure given by name ------------------------------------------

test_that("a structure given by name is read from the model frame", {
  set.seed(3)
  d <- data.frame(id = rep(c("a", "b", "c"), each = 20), day = rep(1:20, 3))
  d$y <- 1 + stats::rnorm(60) + rep(cumsum(stats::rnorm(20, 0, 0.5)), 3)
  fn <- ilm_model(y ~ 1, data = d, family = "gaussian",
                  ar = ilm_rw1(~ day | id), verbose = FALSE)
  fv <- ilm_model(y ~ 1, data = d, family = "gaussian",
                  ar = ilm_rw1(d$day, d$id), verbose = FALSE)
  expect_equal(fn$opt$par, fv$opt$par, tolerance = 1e-8)
  expect_equal(fn$ar$vars, c(time = "day", group = "id"))
  expect_null(fv$ar$vars)
  ## a row dropped for a missing value leaves the structure with it, where
  ## vectors taken from the whole data frame fall out of step
  dn <- d; dn$y[c(5, 41)] <- NA
  fd <- ilm_model(y ~ 1, data = dn, family = "gaussian",
                  ar = ilm_rw1(~ day | id), verbose = FALSE)
  ok <- !is.na(dn$y)
  fc <- ilm_model(y ~ 1, data = dn[ok, ], family = "gaussian",
                  ar = ilm_rw1(dn$day[ok], dn$id[ok]), verbose = FALSE)
  expect_equal(fd$opt$par, fc$opt$par, tolerance = 1e-8)
  expect_error(ilm_model(y ~ 1, data = dn, family = "gaussian",
                         ar = ilm_rw1(dn$day, dn$id), verbose = FALSE),
               "model matrix has 58 rows")
  ## the same for the two older structures
  expect_equal(ilm_ar1(~ day | id)$vars, c(time = "day", group = "id"))
  expect_s3_class(ilm_car1(~ day | id), "ilm_cor_named")
  ## what cannot be read
  expect_error(ilm_rw1(~ day | id, d$id), "and a `group` as well")
  expect_error(ilm_rw1(~ day), "one-sided formula")
  expect_error(ilm_car1(~ day + id), "one-sided formula")
  expect_error(ilm_model(y ~ 1, data = d, family = "gaussian",
                         ar = ilm_rw1(~ dy | id), verbose = FALSE),
               "'dy', which is not a column")
  X <- matrix(1, 60, 1, dimnames = list(NULL, "(Intercept)"))
  expect_error(ilm_fit(X, d$y, family = "gaussian",
                       ar = ilm_rw1(~ day | id), verbose = FALSE),
               "ilm_fit\\(\\) does not see")
})

## ---- ilm_cells() ------------------------------------------------------------

test_that("ilm_cells() lays out every structure's cells", {
  set.seed(4)
  ## twelve units seen weekly for eight weeks; one misses week 3 and one
  ## stops after week 5
  d <- expand.grid(wk = 0:7, id = sprintf("s%02d", 1:12), stringsAsFactors = FALSE)
  d <- d[!(d$id == "s01" & d$wk == 3) & !(d$id == "s02" & d$wk > 5), ]
  d$day <- as.Date("2024-01-01") + 7 * d$wk
  d$y <- stats::rnorm(nrow(d)) + stats::rnorm(12)[as.integer(factor(d$id))]
  fit <- function(ar) suppressMessages(suppressWarnings(
    ilm_model(y ~ 1, data = d, family = "gaussian", ar = ar, verbose = FALSE)))

  ## AR(1): the whole grid for every unit, empty steps included
  ca <- ilm_cells(fit(ilm_ar1(~ day | id)))
  expect_equal(nrow(ca), 12L * 8L)
  expect_equal(names(ca), c("term", "type", "group", "time", "index", "n_obs",
                            "anchor", "last"))
  expect_s3_class(ca$time, "Date")
  expect_equal(levels(ca$group), sprintf("s%02d", 1:12))
  expect_equal(ca$index, seq_len(96))
  expect_equal(ca$n_obs[ca$group == "s01"], c(1, 1, 1, 0, 1, 1, 1, 1))
  ## a unit seen only early still runs to the end of the grid
  s02 <- ca[ca$group == "s02", ]
  expect_equal(s02$n_obs, c(rep(1, 6), 0, 0))
  expect_equal(which(s02$last), 8L)
  expect_equal(attr(ca, "step"), 7)
  expect_equal(attr(ca, "origin"), as.Date("2024-01-01"))
  expect_false(any(ca$anchor))

  ## CAR(1): one cell per observed time
  cc <- ilm_cells(fit(ilm_car1(~ day | id)))
  expect_equal(nrow(cc), nrow(d))
  expect_true(all(cc$n_obs == 1L))
  expect_equal(sum(cc$last), 12L)
  expect_null(attr(cc, "step"))

  ## a random walk: its anchors have no index
  fr <- fit(ilm_rw1(~ day | id))
  cr <- ilm_cells(fr)
  expect_equal(cr$anchor, !duplicated(cr$group))
  expect_true(all(is.na(cr$index[cr$anchor])))
  expect_equal(cr$index[!cr$anchor], seq_len(nrow(d) - 12L))
  expect_equal(attr(cr, "n_latent"), nrow(d) - 12L)
  expect_match(attr(cr, "parameterisation"), "variance per unit")
  expect_equal(attr(cr, "vars"), c(time = "day", group = "id"))

  ## every fitted row sits on the cell with its own group and time
  oc <- attr(cr, "obs_cell")
  expect_equal(cr$time[oc], fr$model$day)
  expect_equal(as.character(cr$group[oc]), as.character(fr$model$id))
  ## and each cell's index finds its value in the fit
  pr <- fr$sdr$par.random
  Ba <- illume:::ilm_Bar_hat(fr)
  expect_equal(Ba[!cr$anchor, 1], unname(pr[names(pr) == "B_ar"]))
  expect_true(all(Ba[cr$anchor, ] == 0))

  expect_null(ilm_cells(ilm_model(y ~ 1, data = d, family = "gaussian",
                                  verbose = FALSE)))
  expect_error(ilm_cells(d), "must be a fitted ilm_model")
})

## ---- the walk away from a gaussian response --------------------------------

test_that("a population average over a walk spreads with the time since it began", {
  skip_on_cran()
  set.seed(11)
  ng <- 8; nt <- 25
  d <- data.frame(id = factor(rep(sprintf("u%02d", 1:ng), each = nt)),
                  day = rep(seq(0, by = 2, length.out = nt), ng),
                  x = stats::rnorm(ng * nt))
  u <- unlist(lapply(1:ng, function(i)
    c(0, cumsum(stats::rnorm(nt - 1, 0, sqrt(0.02 * 2))))))
  d$y <- stats::rpois(nrow(d), exp(1 + 0.3 * d$x + u))
  f <- ilm_model(y ~ x, data = d, family = "poisson",
                 ar = ilm_rw1(~ day | id), verbose = FALSE)
  expect_equal(f$opt$convergence, 0L)
  s2 <- f$Sigma$ar[1, 1]
  ## through a log link the average is exact: exp(eta + s2 * elapsed / 2)
  el <- d$day - stats::ave(d$day, d$id, FUN = min)
  eta <- as.vector(f$X %*% f$beta)
  expect_equal(unname(predict(f, marginal = TRUE)[, 1]),
               exp(eta + 0.5 * s2 * el), tolerance = 1e-10)
  ## new rows: a fitted group from its own first time, and one the fit has
  ## not seen from its earliest row
  nd <- data.frame(x = 0, day = c(10, 48, 100, 104), id = c("u01", "u01", "zz", "zz"))
  expect_equal(unname(predict(f, newdata = nd, marginal = TRUE)[, 1]),
               exp(f$beta[1] + 0.5 * s2 * c(10, 48, 0, 4)), tolerance = 1e-10)
  ## the typical group is unaffected
  expect_equal(unname(predict(f, newdata = nd)[, 1]), rep(exp(f$beta[1]), 4),
               tolerance = 1e-12)
  ## given as vectors, the walk cannot place new rows
  fv <- ilm_model(y ~ x, data = d, family = "poisson",
                  ar = ilm_rw1(d$day, d$id), verbose = FALSE)
  expect_error(predict(fv, newdata = nd, marginal = TRUE), "Give it by name")
  expect_equal(predict(fv, marginal = TRUE), predict(f, marginal = TRUE),
               tolerance = 1e-6)
  expect_output(print(summary(f)), "random walk over time")
})

test_that("the variance is per unit of time, and the fit is the same in any unit", {
  skip_on_cran()
  ## the same walk timed in days and in seconds: the variance per second is
  ## the variance per day over 86400, and the likelihood does not move -- a
  ## step's variance, gap * Sigma, is the same number either way
  set.seed(21)
  ng <- 6; nt <- 20
  d <- data.frame(id = rep(sprintf("g%d", 1:ng), each = nt),
                  day = rep(cumsum(c(0, sample(1:3, nt - 1, TRUE))), ng))
  u <- unlist(lapply(1:ng, function(i) {
    dd <- d$day[d$id == sprintf("g%d", i)]
    c(0, cumsum(stats::rnorm(nt - 1, 0, sqrt(0.05 * diff(dd)))))
  }))
  d$y <- stats::rpois(nrow(d), exp(1.2 + u))
  d$when <- as.POSIXct("2024-03-01", tz = "UTC") + 86400 * d$day
  fd <- ilm_model(y ~ 1, data = d, family = "poisson",
                  ar = ilm_rw1(~ day | id), verbose = FALSE)
  fs <- ilm_model(y ~ 1, data = d, family = "poisson",
                  ar = ilm_rw1(~ when | id), verbose = FALSE)
  expect_equal(as.numeric(logLik(fs)), as.numeric(logLik(fd)), tolerance = 1e-7)
  expect_equal(fs$Sigma$ar[1, 1] * 86400, fd$Sigma$ar[1, 1], tolerance = 1e-5)
  expect_equal(coef(fs), coef(fd), tolerance = 1e-5)
  cs <- ilm_cells(fs)
  expect_s3_class(cs$time, "POSIXct")
  expect_equal(cs$time[attr(cs, "obs_cell")], d$when)
  ## and for a multinomial walk, whose covariance has a floor on its diagonal
  set.seed(22)
  w <- unlist(lapply(1:ng, function(i) {
    dd <- d$day[d$id == sprintf("g%d", i)]
    c(0, cumsum(stats::rnorm(nt - 1, 0, sqrt(0.08 * diff(dd)))))
  }))
  eta <- cbind(0, 0.4 + u, -0.3 + w)
  pr <- exp(eta) / rowSums(exp(eta))
  d$k <- factor(apply(pr, 1, function(q) sample(c("a", "b", "c"), 1, prob = q)))
  md <- ilm_model(k ~ 1, data = d, family = "multinomial",
                  ar = ilm_rw1(~ day | id), verbose = FALSE)
  ms <- ilm_model(k ~ 1, data = d, family = "multinomial",
                  ar = ilm_rw1(~ when | id), verbose = FALSE)
  expect_equal(as.numeric(logLik(ms)), as.numeric(logLik(md)), tolerance = 1e-6)
  expect_equal(ms$Sigma$ar * 86400, md$Sigma$ar, tolerance = 1e-3)
})

test_that("a simulated walk starts at zero and spreads with elapsed time", {
  ar <- ilm_rw1(rep(c(0, 1, 3, 7), 2), rep(1:2, each = 4))
  set.seed(5)
  dr <- replicate(6000, illume:::ilm_ar_draw(ar, NA_real_, matrix(0.5), 1L)[, 1])
  expect_true(all(dr[c(1, 5), ] == 0))
  v <- apply(dr, 1, stats::var)
  expect_equal(v[c(2, 3, 4)], 0.5 * c(1, 3, 7), tolerance = 0.1)
  ## the increments are independent: the step from 3 to 7 is uncorrelated
  ## with the level at 3
  expect_lt(abs(stats::cor(dr[3, ], dr[4, ] - dr[3, ])), 0.05)
})

test_that("a multinomial walk holds its anchors at zero in every category", {
  skip_on_cran()
  set.seed(12)
  ng <- 30; nt <- 12
  d <- data.frame(id = rep(1:ng, each = nt), t = rep(1:nt, ng))
  w1 <- unlist(lapply(1:ng, function(i) c(0, cumsum(stats::rnorm(nt - 1, 0, 0.3)))))
  w2 <- unlist(lapply(1:ng, function(i) c(0, cumsum(stats::rnorm(nt - 1, 0, 0.3)))))
  eta <- cbind(0, 0.3 + w1, -0.2 + w2)
  pr <- exp(eta) / rowSums(exp(eta))
  d$y <- factor(apply(pr, 1, function(p) sample(c("a", "b", "c"), 1, prob = p)))
  f <- ilm_model(y ~ 1, data = d, family = "multinomial",
                 ar = ilm_rw1(~ t | id), verbose = FALSE)
  expect_equal(f$opt$convergence, 0L)
  cl <- ilm_cells(f); nl <- attr(cl, "n_latent")
  v <- f$sdr$par.random[names(f$sdr$par.random) == "B_ar"]
  expect_length(v, nl * 2L)
  Ba <- illume:::ilm_Bar_hat(f)
  expect_true(all(Ba[cl$anchor, ] == 0))
  ## the second linear predictor's values follow the first's
  free <- !cl$anchor
  expect_equal(Ba[free, 2], unname(v[cl$index[free] + nl]))
  P <- predict(f, marginal = TRUE, ndraw = 50)
  expect_equal(unname(rowSums(P)), rep(1, nrow(d)), tolerance = 1e-12)
})

test_that("the walk survives a rebuild, a resample and a print", {
  set.seed(7)
  d <- data.frame(id = rep(c("a", "b"), each = 15), t = rep(1:15, 2))
  d$y <- stats::rnorm(30, 0, 0.5) +
    c(cumsum(c(0, stats::rnorm(14, 0, 0.6))), cumsum(c(0, stats::rnorm(14, 0, 0.6))))
  f <- ilm_model(y ~ 1, data = d, family = "gaussian",
                 ar = ilm_rw1(~ t | id), verbose = FALSE)
  ## a new parameter vector moves the variance and leaves no correlation
  p <- f$opt$par; p[names(p) == "lchol_ar"] <- log(0.3)
  rb <- illume:::ilm_rebuild(f, p)
  expect_equal(rb$Sigma$ar[1, 1], 0.09)
  expect_true(is.na(rb$rho))
  ## the rows a simulated study draws keep their gaps
  s <- illume:::ilm_ar_rows(f$ar, c(16:30, 1:15), rep(1:2, each = 15))
  expect_s3_class(s, "ilm_rw1")
  expect_equal(s$gap, f$ar$gap)
  expect_equal(illume:::ilm_ar_group(f$ar), rep(1:2, each = 15))
  out <- capture.output(print(f$ar))
  expect_true(any(grepl("Random walk, continuous time", out, fixed = TRUE)))
  expect_true(any(grepl("after the 2 held at zero", out, fixed = TRUE)))
  expect_output(print(ilm_rw1(~ t | id)), "built from the data")
})
