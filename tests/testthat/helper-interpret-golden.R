## The interpreter's golden cases: fits that reach every branch of
## ilm_interpret() for a model and for contrasts -- each family the effects
## are said for, numbers and factors of two and more levels, association and
## causal wording, the coefficients' fallback with its ratios and average
## marginal effects, separation, an interaction, grouping, an offset, a family
## read off the response, and a p-value on each side of 0.001. The English
## each gives is kept in tests/testthat/golden/interpret/, written by
## dev/make_interpret_golden.R, and test-interpret-golden.R holds the
## interpreter to it word for word.

golden_interpret_cases <- function() {
  q <- function(e) suppressMessages(suppressWarnings(e))
  set.seed(20260929)
  n <- 240
  d <- data.frame(x = stats::rnorm(n), z = stats::rnorm(n),
                  g2 = factor(sample(c("control", "treated"), n, TRUE)),
                  h3 = factor(sample(c("north", "south", "west"), n, TRUE)),
                  site = factor(rep(seq_len(24), each = 10)),
                  income = round(stats::runif(n, 20000, 90000)),
                  t = stats::runif(n, 0.5, 3))
  u <- stats::rnorm(24, 0, 0.5)[d$site]
  d$y <- 1 + 0.6 * d$x + 0.05 * d$z + 0.8 * (d$g2 == "treated") +
    c(0, 0.3, -0.2)[d$h3] + u + stats::rnorm(n)
  d$yi <- 20 + 0.0002 * d$income + stats::rnorm(n, 0, 4)
  d$yb <- stats::rbinom(n, 1, stats::plogis(-0.3 + 0.9 * d$x + 0.5 * (d$g2 == "treated") + u))
  d$k <- stats::rpois(n, exp(0.4 + 0.3 * d$x + 0.2 * (d$g2 == "treated") + u))
  d$kt <- stats::rpois(n, d$t * exp(0.2 + 0.3 * d$x))
  d$kn <- stats::rnbinom(n, mu = exp(1 + 0.4 * d$x), size = 2)
  d$m3 <- factor(c("a", "b", "c")[1 + (d$x + stats::rnorm(n) > 0) + (d$z + stats::rnorm(n) > 0.8)])
  d$o <- factor(cut(0.8 * d$x + stats::rlogis(n), c(-Inf, -0.5, 0.7, Inf),
                    labels = c("low", "mid", "high")), ordered = TRUE)
  d$pr <- stats::plogis(0.2 + 0.5 * d$x + stats::rnorm(n, 0, 0.3))
  ## a level whose outcome never varies, for the separation sentence
  ds <- d; ds$ys <- ds$yb; ds$ys[ds$h3 == "west"] <- 0L
  m <- function(f, data = d, ...) q(ilm_model(f, data = data, verbose = FALSE, ...))
  fits <- list(
    gaussian_fixed = m(y ~ x + g2 + h3, family = "gaussian"),
    ## maximum likelihood, as these sentences were written for (a gaussian
    ## mixed model is fitted by REML by default, Craig's item 249)
    gaussian_mixed = m(y ~ x + z + h3 + (1 | site), family = "gaussian", reml = FALSE),
    gaussian_interaction = m(y ~ x * g2 + z, family = "gaussian"),
    gaussian_poly = m(y ~ poly(x, 2) + g2, family = "gaussian"),
    gaussian_whole = m(yi ~ income, family = "gaussian"),
    gaussian_inferred = m(y ~ x + g2),
    binomial = m(yb ~ x + g2 + z, family = "binomial"),
    binomial_mixed = m(yb ~ x + g2 + (1 | site), family = "binomial"),
    binomial_interaction = m(yb ~ x * g2, family = "binomial"),
    binomial_separation = m(ys ~ x + h3, data = ds, family = "binomial"),
    poisson_mixed = m(k ~ x + g2 + (1 | site), family = "poisson"),
    poisson_offset = m(kt ~ x + offset(log(t)), family = "poisson"),
    nbinom = m(kn ~ x + g2, family = "nbinom"),
    multinomial = m(m3 ~ x + g2, family = "multinomial"),
    ordinal = m(o ~ x + g2, family = "ordinal"),
    beta = m(pr ~ x, family = "beta"))
  out <- lapply(fits, function(f) list(object = f, args = list()))
  out$gaussian_causal <- list(object = fits$gaussian_fixed, args = list(causal = TRUE))
  out$binomial_no_ame <- list(object = fits$binomial, args = list(ame = FALSE))
  out$contrast_gaussian <- list(object = q(ilm_contrast(q(ilm_emmeans(fits$gaussian_fixed, "h3")))),
                                args = list())
  out$contrast_poisson <- list(object = q(ilm_contrast(q(ilm_emmeans(fits$poisson_mixed, "g2")))),
                               args = list())
  out
}

## what a case says, section by section, as printed text
golden_interpret_text <- function(case) {
  it <- suppressMessages(suppressWarnings(do.call(ilm_interpret, c(list(case$object), case$args))))
  s <- unclass(it)$sections
  unlist(lapply(names(s), function(k) c(paste0("## ", k), s[[k]], "")), use.names = FALSE)
}
