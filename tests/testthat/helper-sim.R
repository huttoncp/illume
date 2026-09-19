# Small simulated datasets shared by the tests.  Deliberately tiny: the tests
# check behaviour and contracts, not statistical performance, and the whole
# suite has to run quickly.

sim_mlmm <- function(seed = 1, n_subj = 25, per = 12, J = 3, with_time = FALSE,
                     with_smooth = FALSE, rho = 0) {
  set.seed(seed)
  C <- J - 1L
  N <- n_subj * per
  Tc <- stats::contr.sum(J)
  dd <- data.frame(
    subj = factor(rep(seq_len(n_subj), each = per)),
    time = rep(seq_len(per), n_subj),
    x1   = stats::rnorm(N),
    grp  = factor(sample(c("a", "b", "c"), N, TRUE)),
    xs   = stats::runif(N)
  )
  X <- stats::model.matrix(~ x1 + grp, dd)
  bt <- matrix(stats::runif(ncol(X) * C, -0.7, 0.7), ncol(X), C)
  b <- matrix(stats::rnorm(n_subj * C), n_subj, C) %*%
       diag(sqrt(c(0.6, 0.35, rep(0.3, max(0, C - 2)))[seq_len(C)]), C)
  eta <- X %*% bt + b[as.integer(dd$subj), , drop = FALSE]
  if (with_time) eta <- eta + 0.4 * dd$time / per
  if (with_smooth) eta <- eta + cbind(sin(2 * pi * dd$xs),
                                      matrix(0, N, C - 1L))[, seq_len(C), drop = FALSE]
  P <- exp(eta %*% t(Tc)); P <- P / rowSums(P)
  labs <- paste0("c", seq_len(J))
  dd$y <- factor(labs[apply(P, 1, function(pr) sample.int(J, 1L, prob = pr))],
                 levels = labs)
  attr(dd, "beta_true") <- bt
  dd
}

# a fitted model most tests can reuse
fit_basic <- function(seed = 1, ...) {
  dd <- sim_mlmm(seed = seed, ...)
  illume::lum_model(y ~ x1 + grp + (1 | subj), data = dd,
                    family = "multinomial", verbose = FALSE)
}
