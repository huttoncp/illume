## Subject curves: a common curve plus each subject's intercept, slope and a
## small wiggle. Used by identity.R.
sim <- function(seed, n_id = 30, n_per = 10, pois = FALSE, unbal = FALSE) {
  set.seed(seed)
  sizes <- if (unbal) sample(3:18, n_id, replace = TRUE) else rep(n_per, n_id)
  id <- factor(rep(seq_len(n_id), times = sizes))
  x <- runif(length(id))
  a <- rnorm(n_id, 0, 0.6)[id]; b <- rnorm(n_id, 0, 0.8)[id]
  w <- rnorm(n_id, 0, 0.3)[id]
  eta <- sin(2 * pi * x) + a + b * (x - 0.5) + w * cos(3 * pi * x)
  y <- if (pois) rpois(length(x), exp(0.5 + 0.5 * eta)) else eta + rnorm(length(x), 0, 0.4)
  data.frame(y = y, x = x, id = id)
}
