## The df of differences between slopes, and a fallback that says so.
##
## A difference between two slopes tested on Satterthwaite's df is tested on
## its own Satterthwaite df, as emmeans and lmerTest test it -- it was once a
## z test beside the t tests of the slopes it compared. And a row whose df
## cannot be formed is a z test that is named as one, never one that passes
## for the t test the rest of the table has.

## random slopes, so a difference between arms' slopes has a df set by the
## number of subjects -- far from the z test's infinity
slope_data <- function(seed = 5) {
  set.seed(seed)
  ns <- 18
  d <- expand.grid(id = factor(seq_len(ns)), week = 0:5)
  d$arm <- factor(rep(c("ctl", "drugA", "drugB"),
                      each = ns / 3))[as.integer(d$id)]
  sl <- c(ctl = 0.05, drugA = 0.45, drugB = 0.20)
  u0 <- rnorm(ns, 0, 1.0); u1 <- rnorm(ns, 0, 0.25)
  d$y <- 3 + (sl[as.character(d$arm)] + u1[as.integer(d$id)]) * d$week +
    u0[as.integer(d$id)] + rnorm(nrow(d), 0, 0.6)
  d
}

## emmeans writes "drugA - ctl" where illume may write "ctl - drugA": match
## the pairs by their levels, whatever the order
pair_key <- function(s) vapply(strsplit(s, " - ", fixed = TRUE),
                               function(p) paste(sort(trimws(p)), collapse = "|"), "")

test_that("differences between Satterthwaite slopes match emmeans and lmerTest", {
  skip_if_not_installed("lmerTest")
  skip_if_not_installed("emmeans")
  d <- slope_data()
  f <- ilm_model(y ~ arm * week + (1 + week | id), data = d,
                 family = "gaussian", reml = TRUE, verbose = FALSE)
  m <- lmerTest::lmer(y ~ arm * week + (1 + week | id), data = d, REML = TRUE)
  ct <- ilm_contrast(ilm_trends(f, "arm", var = "week"), adjust = "none")
  ref <- as.data.frame(emmeans::contrast(
    emmeans::emtrends(m, "arm", var = "week", lmer.df = "satterthwaite"),
    "pairwise", adjust = "none"))
  i <- match(pair_key(ref$contrast), pair_key(ct$contrast))
  expect_false(anyNA(i))
  expect_equal(abs(ct$estimate[i]), abs(ref$estimate), tolerance = 1e-5)
  expect_equal(ct$se[i], ref$SE, tolerance = 1e-5)
  expect_equal(ct$df[i], ref$df, tolerance = 1e-3)
  expect_equal(ct$p_value[i], ref$p.value, tolerance = 1e-4)
  ## a t test on about the number of subjects, not a z test
  expect_true(all(ct$df < 30))
  expect_identical(attr(ct, "df_method"), "satterthwaite")

  ## and a difference that is one coefficient has that coefficient's df in
  ## lmerTest's own summary
  lt <- coef(summary(m))
  k <- which(pair_key(ct$contrast) == pair_key("drugA - ctl"))
  expect_equal(ct$df[k], unname(lt["armdrugA:week", "df"]), tolerance = 1e-3)
})

test_that("differences between Kenward-Roger slopes match emmeans", {
  skip_if_not_installed("lmerTest")
  skip_if_not_installed("emmeans")
  skip_if_not_installed("pbkrtest")
  d <- slope_data()
  f <- ilm_model(y ~ arm * week + (1 + week | id), data = d,
                 family = "gaussian", reml = TRUE, verbose = FALSE)
  m <- lmerTest::lmer(y ~ arm * week + (1 + week | id), data = d, REML = TRUE)
  ct <- ilm_contrast(ilm_trends(f, "arm", var = "week", df = "kenward-roger"),
                     adjust = "none")
  ref <- as.data.frame(emmeans::contrast(
    emmeans::emtrends(m, "arm", var = "week", lmer.df = "kenward-roger"),
    "pairwise", adjust = "none"))
  i <- match(pair_key(ref$contrast), pair_key(ct$contrast))
  ## the adjusted covariance, and Kenward-Roger's df for each difference
  expect_equal(ct$se[i], ref$SE, tolerance = 1e-5)
  expect_equal(ct$df[i], ref$df, tolerance = 1e-3)
  expect_identical(attr(ct, "df_method"), "kenward-roger")
})

test_that("the joint max-t reference uses the smallest df and says t", {
  d <- slope_data()
  f <- ilm_model(y ~ arm * week + (1 + week | id), data = d,
                 family = "gaussian", reml = TRUE, verbose = FALSE)
  tr <- ilm_trends(f, "arm", var = "week")
  cn <- ilm_contrast(tr, adjust = "none")
  cm <- ilm_contrast(tr)
  expect_equal(cm$df, cn$df)
  ## joint intervals are wider than one-at-a-time ones on the same df
  expect_true(all(cm$upper - cm$lower > cn$upper - cn$lower))
  out <- gsub("[[:space:]]+", " ", paste(capture.output(print(cm)), collapse = " "))
  expect_match(out, "t tests on satterthwaite df", fixed = TRUE)
})

test_that("contrasts of marginal means take the reference of their means", {
  ## exact t with nothing integrated out, Satterthwaite's t otherwise
  d <- slope_data()
  fl <- ilm_model(y ~ arm + week, data = d, family = "gaussian",
                  verbose = FALSE)
  cl <- ilm_contrast(ilm_emmeans(fl, "arm"), adjust = "none")
  expect_equal(cl$df, rep(fl$resid_df, 3))
  fm <- ilm_model(y ~ arm + week + (1 | id), data = d, family = "gaussian",
                  verbose = FALSE)
  cm <- ilm_contrast(ilm_emmeans(fm, "arm"), adjust = "none")
  expect_true(all(is.finite(cm$df) & cm$df > 0))
  expect_identical(attr(cm, "df_method"), "satterthwaite")
  ## and z on request
  cz <- ilm_contrast(ilm_emmeans(fm, "arm", df = "asymptotic"), adjust = "none")
  expect_equal(cz$df, rep(Inf, 3))
  expect_true(any(grepl("z tests", capture.output(print(cz)), fixed = TRUE)))
})

test_that("a difference whose df cannot be formed is a z test, and says so", {
  ## the slope along week differs between arms and not between sites, so a
  ## pair of cells in the same arm differs by nothing: its variance is zero
  ## and no Satterthwaite df exists for it, while the others have theirs
  d <- slope_data()
  d$site <- factor(rep(c("north", "south"), length.out = nrow(d)))
  f <- ilm_model(y ~ arm * week + site + (1 + week | id), data = d,
                 family = "gaussian", reml = TRUE, verbose = FALSE)
  tr <- ilm_trends(f, c("arm", "site"), var = "week")
  expect_warning(ct <- ilm_contrast(tr, adjust = "none"),
                 "could not be formed for contrasts .*the contrast has no variance")
  fb <- attr(ct, "df_fallback")
  expect_false(is.null(fb))
  same_arm <- vapply(strsplit(ct$contrast, " - ", fixed = TRUE), function(p)
    identical(sub(" .*", "", p[1]), sub(" .*", "", p[2])), TRUE)
  expect_identical(fb$rows, which(same_arm))
  expect_true(all(is.infinite(ct$df[same_arm])))
  expect_true(all(is.finite(ct$df[!same_arm])))
  out <- gsub("[[:space:]]+", " ", paste(capture.output(print(ct)), collapse = " "))
  expect_match(out, "t tests on satterthwaite df; z tests for contrasts")
  expect_match(out, "the contrast has no variance")
})

test_that("slopes whose df cannot be formed are z tests, named in the header", {
  d <- slope_data()
  f <- ilm_model(y ~ arm * week + (1 + week | id), data = d,
                 family = "gaussian", reml = TRUE, verbose = FALSE)
  ## a fit whose variance components have no usable sampling covariance --
  ## as when sdreport could not invert their block of the Hessian -- while
  ## the fixed effects' own covariance, and so every slope's SE, is intact
  th <- ilm_par_blocks(f)$theta
  f$sdr$cov.fixed[th, th] <- NaN
  expect_warning(tr <- ilm_trends(f, "arm", var = "week"),
                 "could not be formed for any slope: the variance components")
  expect_identical(attr(tr, "df_method"), "asymptotic")
  expect_true(all(is.infinite(tr$df)))
  expect_true(all(is.finite(tr$se) & tr$se > 0))
  expect_equal(tr$p.value, 2 * stats::pnorm(-abs(tr$statistic)))
  out <- gsub("[[:space:]]+", " ", paste(capture.output(print(tr)), collapse = " "))
  expect_match(out, "z tests: Satterthwaite's df could not be formed")
  ## their differences are tried again, and fall back out loud too
  expect_warning(ilm_contrast(tr, adjust = "none"),
                 "could not be formed for any contrast")
})

test_that("one row falling back leaves the other rows their t", {
  d <- slope_data()
  f <- ilm_model(y ~ arm * week + (1 + week | id), data = d,
                 family = "gaussian", reml = TRUE, verbose = FALSE)
  L <- rbind(replace(numeric(length(f$beta)), 2, 1), 0)
  expect_warning(dd <- ilm_trend_df(f, L, "auto"),
                 "could not be formed for slope 2: the contrast has no variance")
  expect_identical(attr(dd, "method"), "satterthwaite")
  expect_true(is.finite(dd[1]) && is.infinite(dd[2]))
  expect_identical(attr(dd, "fallback")$rows, 2L)
})

test_that("a method named for a fit it does not apply to stops", {
  set.seed(3)
  d <- data.frame(id = factor(rep(1:12, each = 5)), x = rnorm(60))
  d$y <- rpois(60, exp(0.2 + 0.3 * d$x + rnorm(12, 0, 0.4)[d$id]))
  f <- ilm_model(y ~ x + (1 | id), data = d, family = "poisson",
                 verbose = FALSE)
  ## once a quiet z test; "auto" is still z for this family, by design
  expect_error(ilm_trends(f, specs = NULL, var = "x", df = "satterthwaite"),
               "LINEAR mixed models")
  tr <- expect_silent(ilm_trends(f, specs = NULL, var = "x"))
  expect_identical(attr(tr, "df_method"), "asymptotic")
})
