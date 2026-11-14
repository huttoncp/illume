## `groups` on a factor-smooth interaction, s(x, id, bs = "fs"). Its blocks
## are stored as basis terms, which predict() used to evaluate whatever
## `groups` said: an unseen subject got a curve of zero without a word under
## "fitted", and "typical" and "population" gave each row its own subject's
## curve. A subject's curve is a draw from a population of subjects, so the
## term now follows `groups` as a grouping term does: the row's own subject
## under "fitted", zero under "typical", and the average over the subjects'
## distribution under "population".

fsg_data <- function(seed, n_id = 30, n_per = 10, pois = FALSE) {
  set.seed(seed)
  id <- factor(rep(seq_len(n_id), each = n_per))
  x <- stats::runif(length(id))
  a <- stats::rnorm(n_id, 0, 0.6)[id]; b <- stats::rnorm(n_id, 0, 0.8)[id]
  eta <- sin(2 * pi * x) + a + b * (x - 0.5)
  y <- if (pois) stats::rpois(length(x), exp(0.5 + 0.5 * eta))
       else eta + stats::rnorm(length(x), 0, 0.4)
  data.frame(y = y, x = x, id = id)
}
fsg_fit <- function(form, d, family) suppressMessages(suppressWarnings(
  ilm_model(form, data = d, family = family, reml = TRUE, verbose = FALSE)))
fsg_p <- function(...) suppressWarnings(as.numeric(predict(...)))

test_that("the fs blocks are found among the random terms", {
  d <- fsg_data(1)
  f <- fsg_fit(y ~ s(x, k = 5) + s(x, id, bs = "fs", k = 5), d, "gaussian")
  fs <- illume:::ilm_fs_terms(f)
  expect_length(fs, length(f$re))
  expect_true(all(grepl("^s\\(x,id\\)", names(f$re)[fs])))
  expect_false(any(fs[names(f$re) == "s(x)"]))
})

test_that("groups = 'fitted': each row's own subject, and an unseen or missing one is an error", {
  d <- fsg_data(2, pois = TRUE)
  f <- fsg_fit(y ~ s(x, id, bs = "fs", k = 5), d, "poisson")
  ## the fit's own rows: its fitted values, with or without newdata
  fv <- as.numeric(ilm_fitted(f))
  expect_equal(fsg_p(f), fv, tolerance = 1e-8)
  expect_equal(fsg_p(f, newdata = d), fv, tolerance = 1e-8)
  ## an unseen subject, as a factor level or as text
  un <- data.frame(x = 0.5, id = factor("31", levels = c(levels(d$id), "31")))
  expect_error(predict(f, newdata = un), "not in the fit.*s\\(x,id\\).*typical")
  expect_error(predict(f, newdata = data.frame(x = 0.5, id = "31")),
               "1 subject of `id`.*\\(31\\) was not in the fit")
  ## no subject column
  expect_error(predict(f, newdata = data.frame(x = 0.5)),
               "does not have 'id', the subject factor.*typical")
})

test_that("each random term follows its own rule under 'fitted'", {
  d <- fsg_data(3, pois = TRUE)
  set.seed(3); d$h <- factor(sample(letters[1:6], nrow(d), TRUE))
  f <- fsg_fit(y ~ s(x, id, bs = "fs", k = 5) + (1 | h), d, "poisson")
  lv <- function(v, new) factor(v, levels = c(levels(d[[new]]), "z"))
  ok <- data.frame(x = 0.5, id = lv("1", "id"), h = lv("a", "h"))
  expect_true(is.finite(fsg_p(f, newdata = ok)))
  expect_error(predict(f, newdata = transform(ok, id = lv("z", "id"))),
               "subject of `id`")
  expect_error(predict(f, newdata = transform(ok, h = lv("z", "h"))),
               "group of `h`")
})

test_that("groups = 'typical' puts the fs term at zero, in the estimate and in the draws", {
  d <- fsg_data(4, pois = TRUE)
  f <- fsg_fit(y ~ s(x, k = 5) + s(x, id, bs = "fs", k = 5), d, "poisson")
  xs <- c(0.1, 0.5, 0.9)
  ## the same for every subject, seen, unseen or not given
  at <- function(id) data.frame(x = xs, id = factor(id, levels = c(levels(d$id), "new")))
  p1 <- fsg_p(f, newdata = at("1"), groups = "typical")
  expect_equal(fsg_p(f, newdata = at("7"), groups = "typical"), p1)
  expect_equal(fsg_p(f, newdata = at("new"), groups = "typical"), p1)
  expect_equal(fsg_p(f, newdata = data.frame(x = xs), groups = "typical"), p1)
  ## and equal to the fixed part with the common smooth alone
  nd <- illume:::ilm_newX(f, at("1"))
  eta <- nd$X %*% f$beta
  for (k in which(!illume:::ilm_fs_terms(f) &
                  vapply(f$re, function(e) identical(e$kind, "basis"), TRUE)))
    eta <- eta + illume:::ilm_basis_new(f, nd, k) %*% illume:::ilm_Bhat_term(f, k)
  expect_equal(p1, exp(as.numeric(eta)), tolerance = 1e-10)
  ## the draws behind an interval leave it out as well: subjects 1 and 7
  ## have the same standard error, as they would not if their curves were drawn
  s1 <- suppressWarnings(predict(f, newdata = at("1"), groups = "typical", se.fit = TRUE, nsim = 50L))
  s7 <- suppressWarnings(predict(f, newdata = at("7"), groups = "typical", se.fit = TRUE, nsim = 50L))
  expect_true(s1$joint)
  expect_equal(as.numeric(s1$se.fit), as.numeric(s7$se.fit))
})

test_that("ilm_eta() keeps the fs term by default, as ilm_power() needs", {
  d <- fsg_data(5)
  f <- fsg_fit(y ~ s(x, id, bs = "fs", k = 5), d, "gaussian")
  nd <- illume:::ilm_newX(f, d)
  expect_equal(as.numeric(illume:::ilm_eta(f, nd)), as.numeric(ilm_fitted(f)),
               tolerance = 1e-8)
  expect_equal(as.numeric(illume:::ilm_eta(f, nd, fs = FALSE)),
               rep(as.numeric(f$beta[1]), nrow(d)), tolerance = 1e-10)
})

test_that("groups = 'population' integrates the fs term out: the log link in closed form", {
  d <- fsg_data(6, pois = TRUE)
  f <- fsg_fit(y ~ s(x, id, bs = "fs", k = 5), d, "poisson")
  xs <- c(0.05, 0.3, 0.6, 0.95)
  nd1 <- data.frame(x = xs, id = factor("1", levels = levels(d$id)))
  pp <- fsg_p(f, newdata = nd1, groups = "population")
  ## E[exp(eta + Z)] = exp(eta + v / 2) for Z ~ N(0, v), v the sum over the
  ## blocks of the row's sum of squares times the block's variance
  sob <- f$smooths[["s(x,id)"]]
  D <- illume:::ilm_smooth_design(sob, nd1)
  fsk <- which(illume:::ilm_fs_terms(f))
  v <- Reduce(`+`, lapply(seq_along(fsk), function(j)
    rowSums(D$Xr[[j]]^2) * f$Sigma[[fsk[j]]][1, 1]))
  expect_true(all(v > 0))
  expect_equal(pp, exp(as.numeric(f$beta[1]) + v / 2), tolerance = 1e-8)
  ## the same whatever subject, seen, unseen or not given
  un <- data.frame(x = xs, id = factor("new", levels = c(levels(d$id), "new")))
  expect_equal(fsg_p(f, newdata = un, groups = "population"), pp, tolerance = 1e-12)
  expect_equal(fsg_p(f, newdata = data.frame(x = xs), groups = "population"), pp,
               tolerance = 1e-12)
  expect_equal(fsg_p(f, newdata = transform(nd1, id = factor("9", levels = levels(d$id))),
                     groups = "population"), pp, tolerance = 1e-12)
  ## without newdata it is the fit's own rows
  expect_equal(fsg_p(f, groups = "population"),
               fsg_p(f, newdata = d, groups = "population"), tolerance = 1e-12)
})

test_that("groups = 'population' under an identity link is the typical subject", {
  d <- fsg_data(7)
  f <- fsg_fit(y ~ s(x, k = 5) + s(x, id, bs = "fs", k = 5), d, "gaussian")
  nd <- data.frame(x = c(0.2, 0.7))
  expect_equal(fsg_p(f, newdata = nd, groups = "population"),
               fsg_p(f, newdata = nd, groups = "typical"))
})

test_that("a multinomial fit with an fs term is refused under 'population', naming 'typical'", {
  d <- fsg_data(8, n_id = 20, n_per = 15)
  set.seed(8)
  d$y3 <- cut(d$y + stats::rnorm(nrow(d), 0, 0.5), c(-Inf, -0.5, 0.5, Inf),
              labels = c("lo", "mid", "hi"))
  f <- fsg_fit(y3 ~ s(x, id, bs = "fs", k = 4), d, "multinomial")
  nd <- data.frame(x = 0.5, id = factor("1", levels = levels(d$id)))
  expect_error(predict(f, newdata = nd, groups = "population"),
               "multinomial fit with a factor-smooth.*groups = \"typical\"")
  pt <- suppressWarnings(predict(f, newdata = nd, groups = "typical"))
  expect_equal(sum(pt), 1, tolerance = 1e-12)
})

test_that("ilm_scenario() averages over the fs term, as over every other random term", {
  d <- fsg_data(9, pois = TRUE)
  d$g <- factor(ifelse(as.integer(d$id) %% 2 == 0, "a", "b"))
  f <- fsg_fit(y ~ g + s(x, id, bs = "fs", k = 5), d, "poisson")
  sc <- suppressWarnings(ilm_scenario(f, g = c("a", "b"), sims = 20L, progress = FALSE))
  expect_true(attr(sc, "mixed"))
  pop <- vapply(c("a", "b"), function(v) {
    nd <- d; nd$g <- factor(v, levels = levels(d$g))
    mean(fsg_p(f, newdata = nd, groups = "population"))
  }, 0)
  expect_equal(as.data.frame(sc)$estimate, unname(pop), tolerance = 1e-10)
})

test_that("ilm_ame() does not depend on which subject each row has", {
  d <- fsg_data(10, pois = TRUE)
  d$g <- factor(ifelse(as.integer(d$id) %% 2 == 0, "a", "b"))
  f <- fsg_fit(y ~ g + s(x, id, bs = "fs", k = 5), d, "poisson")
  a <- suppressWarnings(as.data.frame(ilm_ame(f, groups = "population")))$estimate
  ## the same average from predict(), with the subjects shuffled
  set.seed(10); ds <- d; ds$id <- ds$id[sample(nrow(ds))]
  m <- function(dd, v) { dd$g <- factor(v, levels = levels(d$g))
    mean(fsg_p(f, newdata = dd, groups = "population")) }
  expect_equal(a, m(ds, "b") - m(ds, "a"), tolerance = 1e-8)
})

test_that("ilm_survival() gives a typical subject's curve, and refuses the other two", {
  d <- fsg_data(11)
  set.seed(11)
  tt <- exp(1 + 0.5 * d$y + 0.6 * log(stats::rexp(nrow(d))))
  ct <- stats::rexp(nrow(d), 1 / (2 * stats::median(tt)))
  d$time <- pmin(tt, ct); d$event <- as.integer(tt <= ct)
  f <- suppressMessages(suppressWarnings(
    ilm_model(time ~ s(x, id, bs = "fs", k = 5), data = d, family = "weibull",
              censor = ilm_surv(d$time, d$event), reml = TRUE, verbose = FALSE)))
  nd <- data.frame(x = 0.5, id = factor("1", levels = levels(d$id)))
  expect_error(ilm_survival(f, newdata = nd, times = 2), "groups = \"fitted\".*typical")
  expect_error(ilm_survival(f, newdata = nd, times = 2, groups = "population"),
               "groups = \"population\".*typical")
  s1 <- suppressWarnings(ilm_survival(f, newdata = nd, times = c(1, 2), groups = "typical"))
  s3 <- suppressWarnings(ilm_survival(f, newdata = transform(nd, id = factor("3", levels = levels(d$id))),
                                      times = c(1, 2), groups = "typical"))
  expect_equal(s1$surv, s3$surv)
  ## the Weibull survivor function at the typical subject's linear predictor
  sc <- unname(f$dispersion[[1]])
  expect_equal(s1$surv, exp(-exp((log(c(1, 2)) - as.numeric(f$beta[1])) / sc)),
               tolerance = 1e-8)
})

## ---- what "typical" and "population" rest on, and the in-sample functions --

test_that("a row's sum of squares in each fs block is the same whichever fitted subject it is given", {
  ## groups = "typical" and "population" build the design with the first
  ## fitted subject put in, and "population" uses each row's sum of squares in
  ## each block; that is right only if the sums do not depend on the subject
  cases <- list(
    list(seed = 21, n_id = 30, form = y ~ s(x, id, bs = "fs", k = 5)),
    list(seed = 22, n_id = 6,  form = y ~ s(x, id, bs = "fs", k = 3, m = 1)),
    list(seed = 23, n_id = 12, form = y ~ s(x, k = 5) + s(x, id, bs = "fs", k = 8)),
    list(seed = 24, n_id = 10, form = y ~ s(x, id, bs = "fs", k = 6, xt = list(bs = "cr"))))
  xs <- c(0, 0.13, 0.5, 0.77, 1)
  for (z in cases) {
    d <- fsg_data(z$seed, n_id = z$n_id)
    f <- fsg_fit(z$form, d, "gaussian")
    sob <- f$smooths[["s(x,id)"]]
    ss <- lapply(levels(d$id), function(l) {
      D <- illume:::ilm_smooth_design(sob, data.frame(x = xs, id = factor(l, levels = levels(d$id))))
      sapply(D$Xr, function(M) rowSums(M^2))
    })
    expect_true(all(ss[[1]] > 0))
    for (j in seq_along(ss)[-1]) expect_equal(ss[[j]], ss[[1]], tolerance = 1e-12)
    ## and the filled-in design gives the same sums as any subject's own
    D <- illume:::ilm_smooth_design(sob, illume:::ilm_fs_fill(sob, data.frame(x = xs)))
    expect_equal(sapply(D$Xr, function(M) rowSums(M^2)), ss[[length(ss)]], tolerance = 1e-12)
  }
})

fsg_multi <- function(seed) {
  d <- fsg_data(seed, n_id = 20, n_per = 15)
  set.seed(seed)
  d$y3 <- cut(d$y + stats::rnorm(nrow(d), 0, 0.5), c(-Inf, -0.5, 0.5, Inf),
              labels = c("lo", "mid", "hi"))
  d
}

test_that("ilm_fitted(groups = 'typical') follows predict(): the fs term at zero", {
  d <- fsg_data(25, pois = TRUE)
  f <- fsg_fit(y ~ s(x, k = 5) + s(x, id, bs = "fs", k = 5), d, "poisson")
  expect_equal(as.numeric(ilm_fitted(f, groups = "typical")),
               fsg_p(f, groups = "typical"), tolerance = 1e-10)
  ## "fitted" keeps each row's own subject, as predict() does
  expect_equal(as.numeric(ilm_fitted(f, groups = "fitted")), fsg_p(f), tolerance = 1e-8)
  dm <- fsg_multi(26)
  fm <- fsg_fit(y3 ~ s(x, id, bs = "fs", k = 4), dm, "multinomial")
  expect_equal(unname(ilm_fitted(fm, groups = "typical")),
               unname(suppressWarnings(predict(fm, groups = "typical"))), tolerance = 1e-10)
})

test_that("ilm_rqr(groups = 'typical') takes each row's distribution with the fs term at zero", {
  d <- fsg_data(27, pois = TRUE)
  f <- fsg_fit(y ~ s(x, id, bs = "fs", k = 5), d, "poisson")
  mu <- fsg_p(f, groups = "typical")
  set.seed(1)
  lo <- stats::ppois(d$y - 1, mu); hi <- stats::ppois(d$y, mu)
  u <- lo + stats::runif(nrow(d)) * (hi - lo)
  r <- ilm_rqr(f, groups = "typical", seed = 1L)
  expect_equal(as.numeric(r), pmin(pmax(u, 1e-10), 1 - 1e-10), tolerance = 1e-10)
  expect_false(isTRUE(all.equal(as.numeric(r), as.numeric(ilm_rqr(f, groups = "fitted", seed = 1L)))))
})

test_that("ilm_scores(groups = 'typical') scores predict()'s typical probabilities", {
  score <- function(P, y) c(log_score = mean(-log(P[cbind(seq_along(y), y)])),
                            accuracy = mean(max.col(P, ties.method = "first") == y))
  dm <- fsg_multi(28)
  fm <- fsg_fit(y3 ~ s(x, id, bs = "fs", k = 4), dm, "multinomial")
  P <- suppressWarnings(predict(fm, groups = "typical"))
  expect_equal(ilm_scores(fm, groups = "typical")[c("log_score", "accuracy")],
               score(P, as.integer(dm$y3)), tolerance = 1e-10)
  ## an ordinal fit, whose probabilities come from its latent predictor
  dm$yo <- factor(dm$y3, ordered = TRUE)
  ## (an ordinal fit needs a predictor besides the smooth)
  fo <- fsg_fit(yo ~ x + s(x, id, bs = "fs", k = 4), dm, "ordinal")
  Po <- suppressWarnings(predict(fo, groups = "typical"))
  expect_equal(ilm_scores(fo, groups = "typical")[c("log_score", "accuracy")],
               score(Po, as.integer(dm$yo)), tolerance = 1e-10)
})

test_that("with no fs term the three are unchanged, and ilm_eta_hat() keeps fs by default", {
  ## the flag reaches only a factor-smooth's blocks: on fits without one the
  ## typical linear predictor is the same with it as without it
  d <- fsg_data(29, pois = TRUE)
  set.seed(29); d$h <- factor(sample(letters[1:6], nrow(d), TRUE))
  for (form in list(y ~ s(x, k = 5) + (1 | h), y ~ x + (1 | id))) {
    f <- fsg_fit(form, d, "poisson")
    expect_identical(illume:::ilm_eta_hat(f, FALSE, fs = FALSE), illume:::ilm_eta_hat(f, FALSE))
    expect_equal(as.numeric(ilm_fitted(f, groups = "typical")),
                 fsg_p(f, groups = "typical"), tolerance = 1e-10)
  }
  ## ilm_robust() and ilm_rp_offset() call ilm_eta_hat(conditional = FALSE)
  ## without the flag, and on an fs fit that still includes the term
  f <- fsg_fit(y ~ s(x, id, bs = "fs", k = 5), d, "poisson")
  expect_equal(illume:::ilm_eta_hat(f, FALSE), illume:::ilm_eta_hat(f, TRUE))
})
