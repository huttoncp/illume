## A random-effect covariance at the edge of its range -- a variance of zero,
## a correlation of +/-1 -- leaves the likelihood flat in that direction. The
## fixed effects are still identified, and ilm_model() reports them, with the
## verdict saying which parts of the fit can be trusted.
##
## Measured on the messy-data regime with a true random-effect SD of 0.05
## (studies/scripts/messy_compare.R, its own seeds, 400 fits): usable fits
## rose from 63.2% to 97.5%, the newly usable ones covered at 0.944, and their
## standard errors were 1.011 times those of the same model without the random
## term -- which is what a variance of zero says they should be.
##
## At the boundary only the direction the data cannot resolve is held (see
## ilm_hess_recover()), so a boundary fit's standard errors are those of the
## reduced model -- a covariance of lower rank, or the term dropped -- and no
## longer depend on what TMB made of its own Hessian. The measurement behind
## that choice is studies/scripts/boundary_se.R.

bnd_data <- function(seed) {
  set.seed(seed); n <- 300
  d <- data.frame(x = rnorm(n),
                  g = factor(sample(c("north", "central", "south"), n, TRUE)),
                  site = factor(sample(30, n, TRUE)))
  ## no site effect at all
  eta <- cbind(0.4 * d$x, -0.3 * d$x + 0.5 * (d$g == "south"))
  P <- exp(cbind(eta, 0)); P <- P / rowSums(P)
  d$y <- factor(apply(P, 1, function(p) sample(c("lo", "mid", "hi"), 1, prob = p)),
                levels = c("lo", "mid", "hi"))
  d
}

test_that("a covariance at its boundary gives the fixed effects of the reduced model", {
  d <- bnd_data(4)
  f <- suppressMessages(ilm_model(y ~ x + g + (1 | site), data = d,
                                  family = "multinomial", verbose = FALSE))
  f0 <- ilm_model(y ~ x + g, data = d, family = "multinomial", verbose = FALSE)
  ## The optimum is a correlation of 1 (logLik -311.3): a covariance of rank
  ## one. Only the direction that reaches the boundary is held, whatever TMB
  ## made of its own Hessian -- which R 4.4 and R 4.6.1 decided differently,
  ## and which used to decide the route -- so the standard errors are those
  ## of the rank-one model itself, on every platform.
  expect_equal(f$hessian_how, "boundary")
  expect_equal(f$hessian_held, "site")
  expect_equal(f$hessian_flat, 1L)
  expect_true(f$ok)
  expect_true("site" %in% f$boundary_terms)
  fr <- ilm_model(y ~ x + g + (1 | site), data = d, family = "multinomial",
                  verbose = FALSE, re_struct = list(site = list(type = "rr", rank = 1L)))
  expect_silent(V <- vcov(f))
  expect_true(all(is.finite(V)))
  expect_equal(sqrt(diag(V)), sqrt(diag(vcov(fr))), tolerance = 1e-3)
  ## and, with no site effect in the data, close to the model without the term
  r <- sqrt(diag(V)) / sqrt(diag(vcov(f0)))
  expect_true(all(r > 0.97 & r < 1.05))
  expect_lt(max(abs(coef(f) - coef(f0))), 0.05)
  expect_silent(suppressMessages(ilm_anova(f)))
  expect_false(isTRUE(f$sdr$pdHess))          # the record stays honest
  expect_equal(f$checks$status[f$checks$check == "hessian"], "BOUNDARY")
  ## the parameters along the held direction carry no uncertainty, and
  ## asking for them says so
  expect_warning(Vf <- vcov(f, full = TRUE), "held at its estimate")
  expect_true(all(is.finite(Vf)))
  ## and the verdict says, where the numbers are read, what to trust
  out <- paste(utils::capture.output(summary(f)), collapse = "\n")
  expect_match(out, "What can be trusted", fixed = TRUE)
  expect_match(out, "held at its estimate", fixed = TRUE)
  expect_match(paste(utils::capture.output(print(f)), collapse = "\n"),
               "[BOUNDARY: site at a boundary; fixed effects usable]",
               fixed = TRUE)
})

test_that("a boundary is held the same way whatever TMB made of its Hessian", {
  ## the property that makes a boundary fit's standard errors the same on
  ## every platform: TMB's verdict on its own Hessian, which moves with the
  ## arithmetic, no longer decides anything
  d <- bnd_data(4)
  f <- suppressMessages(ilm_model(y ~ x + g + (1 | site), data = d,
                                  family = "multinomial", verbose = FALSE))
  pn <- names(f$opt$par)
  cb <- ilm_cov_blocks(f$re, f$Sigma, f$Sigma_d, f$ty, f$rk, f$dk, f$toff,
                       f$ar, f$opt$par, pn)
  expect_true("site" %in% cb$flagged)
  hs <- lapply(c(TRUE, FALSE), function(v) {
    sdr <- f$sdr; sdr$pdHess <- v
    ilm_hess_recover(f$obj, f$opt, sdr, cb, joint = FALSE)
  })
  for (h in hs) {
    expect_equal(h$how, "boundary")
    expect_equal(h$held, "site")
  }
  ib <- pn == "beta"
  expect_equal(hs[[1]]$sdr$cov.fixed[ib, ib], hs[[2]]$sdr$cov.fixed[ib, ib],
               tolerance = 1e-10)
})

test_that("a boundary term curved in every direction gets the accurate Hessian, whatever TMB said", {
  ## A term can be called a boundary by its estimate and still have no
  ## direction flat enough to hold. Its standard errors then come from the
  ## accurate Hessian, never TMB's cruder one -- so they, too, are the same
  ## on every platform. The 0.0.8.9000 run of studies/scripts/boundary_se.R
  ## found one such fit, where TMB's Hessian put them 2 to 10% off.
  f <- fit_basic(1)
  pn <- names(f$opt$par)
  cb <- ilm_cov_blocks(f$re, f$Sigma, f$Sigma_d, f$ty, f$rk, f$dk, f$toff,
                       f$ar, f$opt$par, pn)
  cb$flagged <- "subj"                  # as if at a boundary, curved everywhere
  hs <- lapply(c(TRUE, FALSE), function(v) {
    sdr <- f$sdr; sdr$pdHess <- v
    ilm_hess_recover(f$obj, f$opt, sdr, cb, joint = FALSE)
  })
  for (h in hs) expect_equal(h$how, "recomputed")
  ib <- pn == "beta"
  expect_equal(hs[[1]]$sdr$cov.fixed[ib, ib], hs[[2]]$sdr$cov.fixed[ib, ib],
               tolerance = 1e-10)
})

test_that("a covariance at a correlation of -1 is a boundary, not a failure", {
  ## the rank-deficient covariance used to be graded FAIL, with "the standard
  ## errors above are not usable" beneath standard errors within 4% of the
  ## model without the term
  d <- bnd_data(1)
  f <- ilm_model(y ~ x + g + (1 | site), data = d, family = "multinomial",
                 verbose = FALSE)
  ## wherever on the ridge the optimiser stops -- RTMB 1.x at a log-Cholesky
  ## diagonal of -9.9, RTMB 2.0 at the floor -- the direction along it is held
  expect_equal(f$hessian_how, "boundary")
  expect_true(f$ok)
  expect_equal(f$checks$status[f$checks$check == "sigma_rank[site]"], "BOUNDARY")
  expect_true("site" %in% f$boundary_terms)
  out <- paste(utils::capture.output(summary(f)), collapse = "\n")
  expect_match(out, "What can be trusted", fixed = TRUE)
  expect_match(out, "held at its estimate", fixed = TRUE)
  ## A correlation of -1 is a rank-one covariance, so no fit can beat the
  ## rank-one model by more than rounding. RTMB 2.0 once "beat" it by 0.3,
  ## deep in the region where the Laplace arithmetic is noise, and printed
  ## standard errors of 0.
  fr <- ilm_model(y ~ x + g + (1 | site), data = d, family = "multinomial",
                  verbose = FALSE, re_struct = list(site = list(type = "rr", rank = 1L)))
  expect_lt(as.numeric(logLik(f)) - as.numeric(logLik(fr)), 1e-3)
  ## and it IS the rank-one model, standard errors included
  expect_equal(sqrt(diag(vcov(f))), sqrt(diag(vcov(fr))), tolerance = 1e-3)
  f0 <- ilm_model(y ~ x + g, data = d, family = "multinomial", verbose = FALSE)
  r <- sqrt(diag(vcov(f))) / sqrt(diag(vcov(f0)))
  expect_true(all(r > 0.95 & r < 1.1))
})

test_that("a random slope at a variance of zero keeps the uncertainty of what is still estimated", {
  set.seed(21); n <- 480
  d <- data.frame(g = factor(rep(1:60, each = 8)), x = rnorm(n))
  d$y <- rbinom(n, 1, plogis(-0.3 + 0.5 * d$x + rnorm(60)[d$g]))  # no slope variance
  f <- suppressMessages(ilm_model(y ~ x + (1 + x | g), data = d, verbose = FALSE))
  expect_equal(f$hessian_how, "boundary")
  expect_equal(f$hessian_held, "g")
  expect_true(f$ok)
  f0 <- suppressMessages(ilm_model(y ~ x + (1 | g), data = d, verbose = FALSE))
  ## Only the direction that reaches zero is held. The intercept-slope term
  ## beside it sits at zero but is curved, so its uncertainty is kept, and
  ## the standard errors come out a little above the random-intercept
  ## model's rather than below: holding the whole term would have put them
  ## below, which is how it undercovered in studies/scripts/boundary_se.R.
  r <- sqrt(diag(vcov(f))) / sqrt(diag(vcov(f0)))
  expect_true(all(r >= 0.999 & r < 1.05))
})

test_that("a correlation taken past the floor is brought back to it", {
  ## put the optimiser where RTMB 2.0 took it: a log-Cholesky diagonal of -18,
  ## where the objective is noise
  d <- bnd_data(1)
  f <- ilm_model(y ~ x + g + (1 | site), data = d, family = "multinomial",
                 verbose = FALSE)
  pn <- names(f$opt$par); it <- which(pn == "theta")
  pos <- ilm_floor_pos(f$re, f$ty, f$toff, f$C, pn)
  expect_equal(pos, it[c(1, 3)])                # both diagonals of a 2 x 2 factor
  ctl <- list(iter.max = 3000, eval.max = 3000)
  bad <- f$opt; bad$par[it] <- c(log(0.048), -0.191, -18)
  o <- ilm_floor_refit(f$obj, bad, pos, ilm_logsd_floor, ctl)
  expect_true(all(o$par[pos] >= ilm_logsd_floor))
  expect_lt(abs(o$objective - f$opt$objective), 1e-3)
  ## a fit that never went past the floor is left exactly as it was
  expect_identical(ilm_floor_refit(f$obj, f$opt, pos, ilm_logsd_floor, ctl), f$opt)
  ## nothing to floor with one predictor, or in a reduced-rank or diagonal term
  expect_length(ilm_floor_pos(f$re, "rr", f$toff, f$C, pn), 0L)
  expect_length(ilm_floor_pos(f$re, "diag", f$toff, f$C, pn), 0L)
  expect_length(ilm_floor_pos(f$re, f$ty, f$toff, 1L, pn), 0L)
  ## end to end: unconstrained, this one goes to -13 with RTMB 1.x; the site
  ## effect is null, so the fixed effects are those of the model without it
  d8 <- bnd_data(8)
  f8 <- ilm_model(y ~ x + g + (1 | site), data = d8, family = "multinomial",
                  verbose = FALSE)
  p8 <- ilm_floor_pos(f8$re, f8$ty, f8$toff, f8$C, names(f8$opt$par))
  expect_gte(min(f8$opt$par[p8]), ilm_logsd_floor)
  expect_true(f8$ok)
  f80 <- ilm_model(y ~ x + g, data = d8, family = "multinomial", verbose = FALSE)
  r <- sqrt(diag(vcov(f8))) / sqrt(diag(vcov(f80)))
  expect_true(all(r > 0.97 & r < 1.05))
  expect_lt(max(abs(coef(f8) - coef(f80))), 0.05)
})

test_that("fixed effects that are not identified are not rescued", {
  ## aliased columns: the singular direction is in the fixed effects
  ## themselves, and no holding of a covariance can make them estimable
  set.seed(2); n <- 200
  d <- data.frame(x1 = rnorm(n), site = factor(sample(10, n, TRUE)))
  d$x2 <- 2 * d$x1
  d$y <- rpois(n, exp(0.3 + 0.2 * d$x1))
  f <- suppressWarnings(ilm_model(y ~ x1 + x2 + (1 | site), data = d,
                                  family = "poisson", verbose = FALSE))
  expect_false(f$hessian_how %in% c("recomputed", "boundary"))
  expect_false(f$ok)
  expect_warning(vcov(f), "not positive definite")
})

test_that("the recomputed Hessian agrees with the exact one where both exist", {
  ## a fixed-effects fit has an exact Hessian through automatic
  ## differentiation, so the finite-difference one can be checked against it
  set.seed(5); n <- 200
  d <- data.frame(x = rnorm(n)); d$y <- rpois(n, exp(0.2 + 0.4 * d$x))
  f <- ilm_model(y ~ x, data = d, family = "poisson", verbose = FALSE)
  H <- ilm_hessian(function(p) as.numeric(f$obj$gr(p)), f$opt$par)
  expect_equal(H, unname(f$obj$he(f$opt$par)), tolerance = 1e-6)
})

test_that("the interpretation says which parts of a boundary fit stand", {
  d <- bnd_data(4)
  f <- ilm_model(y ~ x + g + (1 | site), data = d, family = "multinomial",
                 verbose = FALSE)
  it <- ilm_interpret(f, ame = FALSE)
  dg <- paste(it$sections$diagnostics, collapse = " ")
  ## not "all checks passed": the covariance of site is not to be read
  expect_false(grepl("fitting checks passed", dg, fixed = TRUE))
  expect_match(dg, "BOUNDARY --", fixed = TRUE)
  expect_match(dg, "`site`", fixed = TRUE)
  expect_match(dg, "still usable", fixed = TRUE)
  ## and a caveat about few groups, with the calibrated alternative
  expect_match(paste(it$sections$caveats, collapse = " "), "ilm_pb_lrt()",
               fixed = TRUE)
})
