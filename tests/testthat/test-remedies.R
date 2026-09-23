## ilm_remedies() and ilm_apply_remedy(): every check that can be not OK has a
## rule, and the remedies a refit can make are made through the fit's own call.

## x and y are centred within the levels of both g and h, which are crossed
## and balanced, so neither grouping differs by anything and a random
## intercept over either is estimated at a variance of zero -- exactly, on
## any platform, which a simulated null effect would not promise
zero_data <- function(n_g = 20, per = 10, seed = 1) {
  set.seed(seed)
  d <- data.frame(g = factor(rep(seq_len(n_g), each = per)),
                  h = factor(rep(seq_len(per), n_g)),
                  x = rnorm(n_g * per), e = rnorm(n_g * per))
  dc <- function(v) v - ave(v, d$g) - ave(v, d$h) + mean(v)
  d$x <- dc(d$x)
  d$y <- 1 + 0.5 * d$x + dc(d$e)
  d
}

test_that("every check the fit can report has a remedy rule, and no rule is stale", {
  ## the check names, read from the code that makes them: every
  ## ilm_add_check() call in the package, by its name before any `[term]`
  ns <- asNamespace("illume"); found <- character(0)
  walk <- function(e) {
    if (!is.call(e)) return(invisible())
    if (identical(e[[1L]], as.name("ilm_add_check"))) {
      a <- as.list(e)[-1L]
      nm <- names(a); if (is.null(nm)) nm <- rep("", length(a))
      ck <- if ("check" %in% nm) a[["check"]] else a[nm == ""][[2L]]
      s <- if (is.character(ck)) ck else Filter(is.character, as.list(ck))[[1L]]
      found <<- c(found, sub("\\[$", "", s))
    }
    for (x in as.list(e)[-1L]) if (!missing(x)) walk(x)
  }
  for (f in ls(ns, all.names = TRUE)) {
    fn <- get(f, envir = ns)
    if (is.function(fn) && !is.primitive(fn)) walk(body(fn))
  }
  found <- unique(found)
  expect_gte(length(found), 15L)
  expect_setequal(found, illume:::ilm_rem_known)
})

test_that("a random intercept at a variance of zero is dropped, and nothing else moves", {
  d <- zero_data()
  f <- suppressMessages(ilm_model(y ~ x + (1 | g), data = d, verbose = FALSE))
  expect_true("variance_boundary" %in% f$checks$check[f$checks$status != "OK"])
  rem <- ilm_remedies(f)
  expect_s3_class(rem, "ilm_remedies")
  hit <- which(rem$change == "formula = y ~ x")
  expect_length(hit, 1L)
  expect_identical(rem$tier[hit], "structural")
  expect_true(any(rem$change == "boundary = \"avoid\""))
  expect_message(f2 <- ilm_apply_remedy(f, rem, rem$id[hit]), "no longer checked")
  expect_equal(as.numeric(fixef(f2)), unname(coef(lm(y ~ x, d))), tolerance = 1e-5)
  expect_equal(as.numeric(fixef(f2)), as.numeric(fixef(f)), tolerance = 1e-5)
  ## the call is the one a person would have written, and it refits
  expect_identical(getCall(f2)$data, quote(d))
  expect_identical(getCall(f2)$verbose, FALSE)
  expect_identical(nrow(f2$remedy_log), 1L)
  expect_identical(f2$remedy_log$tier, "structural")
  u <- update(f2, . ~ . + I(x^2))
  expect_s3_class(u, "ilm_model")
  expect_output(print(rem), "structural -- ")
  expect_output(print(ilm_remedies(f2)), "No remedies")
})

test_that("remedies chain, and each fit logs how it was reached", {
  d <- zero_data()
  f <- suppressMessages(ilm_model(y ~ x + (1 | g) + (1 | h), data = d, verbose = FALSE))
  r1 <- ilm_remedies(f)
  i1 <- which(r1$change == "formula = y ~ x + (1 | h)")
  expect_length(i1, 1L)
  f2 <- suppressMessages(ilm_apply_remedy(f, r1, r1$id[i1]))
  r2 <- ilm_remedies(f2)
  i2 <- which(r2$change == "formula = y ~ x")
  expect_length(i2, 1L)
  f3 <- suppressMessages(ilm_apply_remedy(f2, r2, r2$id[i2]))
  expect_identical(nrow(f3$remedy_log), 2L)
  expect_identical(f3$remedy_log$change,
                   c("formula = y ~ x + (1 | h)", "formula = y ~ x"))
  expect_identical(getCall(f3)$data, quote(d))
})

test_that("a nested bar is written out when one of its terms is dropped", {
  set.seed(9)
  d <- data.frame(a = factor(rep(1:6, each = 40)), b = factor(rep(rep(1:4, each = 10), 6)),
                  x = rnorm(240), e = rnorm(240))
  d$x <- d$x - ave(d$x, d$a, d$b)
  d$y <- 1 + 0.5 * d$x + d$e - ave(d$e, d$a, d$b)
  f <- suppressMessages(ilm_model(y ~ x + (1 | a/b), data = d, verbose = FALSE))
  rem <- ilm_remedies(f)
  ## dropping one of the two terms (1 | a/b) stands for keeps the other
  keep_a <- grepl("^formula = y ~ x \\+ \\(1 \\| a\\)$", rem$change)
  expect_true(any(keep_a))
  f2 <- suppressMessages(ilm_apply_remedy(f, rem, rem$id[which(keep_a)[1]]))
  expect_identical(names(f2$re), "a")
})

test_that("a smooth shrunk to a straight line is replaced by the line, at the same likelihood", {
  set.seed(3); n <- 300
  d <- data.frame(x = runif(n), site = factor(sample(20, n, TRUE)))
  d$y <- rnorm(n, 1 + 0.5 * d$x + rnorm(20, 0, 0.5)[d$site])
  f <- ilm_model(y ~ s(x, k = 6) + (1 | site), data = d, family = "gaussian",
                 verbose = FALSE)
  row <- f$checks$check == "smooth_shrinkage[s(x)]"
  expect_true(f$checks$status[row] %in% c("FAIL", "BOUNDARY"))
  ## the check no longer says to drop the term, which would drop the line too
  expect_match(f$checks$suggestion[row], "unpenalised part")
  rem <- ilm_remedies(f)
  hit <- grep("^formula = ", rem$change)
  expect_length(hit, 1L)
  expect_identical(rem$tier[hit], "structural")
  f2 <- suppressMessages(ilm_apply_remedy(f, rem, rem$id[hit]))
  expect_false(any(grepl("s(x", deparse1(formula(f2)), fixed = TRUE)))
  expect_equal(as.numeric(logLik(f2)), as.numeric(logLik(f)), tolerance = 1e-6)
  expect_equal(attr(logLik(f2), "df"), attr(logLik(f), "df") - 1)
})

test_that("a category covariance at its edge is given the rank the data support", {
  set.seed(1); n <- 300
  d <- data.frame(x = rnorm(n),
                  g = factor(sample(c("north", "central", "south"), n, TRUE)),
                  site = factor(sample(30, n, TRUE)))
  eta <- cbind(0.4 * d$x, -0.3 * d$x + 0.5 * (d$g == "south"))   # no site effect
  P <- exp(cbind(eta, 0)); P <- P / rowSums(P)
  d$y <- factor(apply(P, 1, function(p) sample(c("lo", "mid", "hi"), 1, prob = p)),
                levels = c("lo", "mid", "hi"))
  f <- suppressMessages(ilm_model(y ~ x + g + (1 | site), data = d,
                                  family = "multinomial", verbose = FALSE))
  rem <- ilm_remedies(f)
  hit <- which(rem$change == "re_struct = list(site = list(type = \"rr\", rank = 1L))")
  expect_length(hit, 1L)
  expect_identical(rem$tier[hit], "structural")
  f2 <- suppressMessages(ilm_apply_remedy(f, rem, rem$id[hit]))
  expect_identical(f2$re_struct$site$type, "rr")
  expect_identical(f2$re_struct$site$rank, 1L)
  expect_identical(f2$checks$status[f2$checks$check == "sigma_rank[site]"], "OK")
})

test_that("with one category dimension a level check names what can come down, not a rank", {
  set.seed(4); n <- 240
  d <- data.frame(x = rnorm(n), g = factor(rep(1:8, each = 30)))
  d$y <- rbinom(n, 1, plogis(0.3 + 0.5 * d$x + rnorm(8, 0, 0.7)[d$g]))
  f <- suppressMessages(ilm_model(y ~ x + (1 + x | g), data = d, verbose = FALSE))
  row <- f$checks$check == "re_levels[g]"
  expect_true(f$checks$status[row] %in% c("WARN", "FAIL"))
  ## rr(1) IS the structure fitted when there is one category dimension
  expect_false(grepl("rr(", f$checks$suggestion[row], fixed = TRUE))
  expect_match(f$checks$suggestion[row], "d_cor = FALSE")
  rem <- ilm_remedies(f)
  expect_false(any(grepl("type = \"rr\"", rem$change, fixed = TRUE)))
  expect_true(any(rem$change == "re_struct = list(g = list(type = \"us\", d_cor = FALSE))"))
  expect_true(any(rem$change == "formula = y ~ x + (1 | g)"))
  expect_true(all(rem$tier %in% c("numerical", "structural", "estimand")))
  ## pooling levels is a remedy only a person can make, and says so
  man <- rem$id[!nzchar(rem$change)]
  expect_gte(length(man), 1L)
  expect_error(ilm_apply_remedy(f, rem, man[1]), "made by hand")
  ## a subset keeps each remedy's own change
  sub <- rem[rem$id != 1L, ]
  k <- sub$id[nzchar(sub$change)][1]
  f2 <- suppressMessages(ilm_apply_remedy(f, sub, k))
  expect_identical(f2$remedy_log$change, rem$change[rem$id == k])
})

test_that("re_struct may name one term, and a name that is not a term stops", {
  set.seed(5); n <- 300
  d <- data.frame(x = runif(n), site = factor(sample(20, n, TRUE)),
                  grp = factor(sample(6, n, TRUE)))
  d$y <- rnorm(n, 1 + sin(3 * d$x) + rnorm(20, 0, 0.5)[d$site] + rnorm(6, 0, 0.5)[d$grp])
  f <- ilm_model(y ~ s(x) + (1 | site) + (1 | grp), data = d, family = "gaussian",
                 verbose = FALSE, re_struct = list(site = list(type = "diag")))
  expect_identical(f$re_struct$site$type, "diag")
  expect_identical(f$re_struct$grp$type, "us")
  expect_error(ilm_model(y ~ x + (1 | site), data = d, family = "gaussian",
                         verbose = FALSE, re_struct = list(sit = list(type = "us"))),
               "'sit' is not a term of this model. The terms are: 'site'")
})

test_that("the standalone checks' remedies are written out and applied", {
  set.seed(6); n <- 300
  d <- data.frame(x = rnorm(n), g = factor(rep(1:30, each = 10)))
  d$y <- rnbinom(n, mu = exp(1 + 0.4 * d$x + rnorm(30, 0, 0.3)[d$g]), size = 1.5)
  f <- ilm_model(y ~ x + (1 | g), data = d, family = "poisson", verbose = FALSE)
  dsp <- suppressMessages(ilm_check_dispersion(f, B = 200))
  expect_identical(dsp$status, "FAIL")
  rem <- ilm_remedies(f, dispersion = dsp)
  hit <- which(rem$change == "family = \"nbinom\"")
  expect_length(hit, 1L)
  expect_output(print(rem), "ilm_apply_remedy(fit, <this list>, id)", fixed = TRUE)
  f2 <- suppressMessages(ilm_apply_remedy(f, rem, rem$id[hit]))
  expect_identical(f2$family$name, "nbinom")

  set.seed(7); n <- 400
  d <- data.frame(x = rnorm(n), g = factor(rep(1:40, each = 10)))
  mu <- exp(1 + 0.4 * d$x + rnorm(40, 0, 0.3)[d$g])
  d$y <- ifelse(runif(n) < 0.3, 0L, rpois(n, mu))
  f <- ilm_model(y ~ x + (1 | g), data = d, family = "poisson", verbose = FALSE)
  zz <- suppressMessages(ilm_check_zeros(f, B = 200))
  expect_identical(zz$status, "FAIL")
  rem <- ilm_remedies(f, zeros = zz)
  expect_true(all(c("ziformula = ~1", "ziformula = ~1, zi_type = \"hurdle\"") %in% rem$change))
  f2 <- suppressMessages(ilm_apply_remedy(f, rem, rem$id[rem$change == "ziformula = ~1"]))
  expect_false(is.null(f2$Zzi))

  ## a spread that differs by group: the dispersion formula is written in the
  ## column the check was given
  set.seed(8); n <- 300
  d <- data.frame(x = rnorm(n), grp = factor(sample(c("a", "b", "c"), n, TRUE)))
  d$y <- 1 + 0.5 * d$x + rnorm(n, sd = c(a = 0.5, b = 1, c = 2)[as.character(d$grp)])
  f <- ilm_model(y ~ x + grp, data = d, family = "gaussian", verbose = FALSE)
  vres <- list(trend = 0.01, p_trend = 0.6, ratio = 9, p_ratio = 0.001,
               n_refits = 100L, B = 100L, status = "FAIL", suggestion = "", by = "grp")
  rem <- ilm_remedies(f, variance = vres)
  expect_identical(rem$change, "dispformula = ~grp")
  f2 <- suppressMessages(ilm_apply_remedy(f, rem, 1L))
  expect_identical(deparse1(getCall(f2)$dispformula), "~grp")
  expect_error(ilm_remedies(f, variance = list(status = "FAIL")), "ilm_check_variance")
})

test_that("remedies belong to their fit, and the data are found or asked for", {
  d <- zero_data()
  f <- suppressMessages(ilm_model(y ~ x + (1 | g), data = d, verbose = FALSE))
  other <- suppressMessages(ilm_model(y ~ x + (1 | h), data = d, verbose = FALSE))
  expect_error(ilm_apply_remedy(f, ilm_remedies(other), 1L), "different fit")
  expect_error(ilm_apply_remedy(f, ilm_remedies(f), 99L), "must be one of the remedy ids")
  ## a model fitted inside a function from a formula written outside it:
  ## the formula's environment does not hold the function's data
  fo <- y ~ x + (1 | g)
  mk <- function() { dd <- d; ilm_model(fo, data = dd, verbose = FALSE) }
  fk <- suppressMessages(mk())
  rk <- ilm_remedies(fk)
  expect_error(ilm_apply_remedy(fk, rk, 1L), "cannot be found .* Pass them as `data`")
  fk2 <- suppressMessages(ilm_apply_remedy(fk, rk, 1L, data = d))
  expect_s3_class(fk2, "ilm_model")
  ## the rules for a check that is not OK read the fit, whatever status it
  ## came with: more restarts for a gradient that is not at zero
  r <- illume:::ilm_rem_rules(f, "gradient", "FAIL")
  expect_identical(r[[1L]]$tier, "numerical")
  expect_identical(r[[1L]]$args, list(restarts = 10L))
  expect_length(illume:::ilm_rem_rules(f, "gradient", "OK"), 0L)
  expect_length(illume:::ilm_rem_rules(f, "parameter_aliasing", "INCONCLUSIVE"), 0L)
  ## the aliased pair is named whole, though a name holds a parenthesis
  al <- illume:::ilm_rem_rules(f, "parameter_aliasing", "WARN",
    detail = "largest |parameter correlation| = 0.970 (beta.(Intercept) <-> theta)")
  expect_match(al[[1L]]$remedy, "(beta.(Intercept) <-> theta)", fixed = TRUE)
})
