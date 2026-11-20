## ilm_survival() against the survivor function evaluated by hand, on models
## with a smooth beside a random intercept: does the curve carry the smooth's
## penalised part, under each `groups`? And does a single time work?
##
## For each fit, at three covariate rows, the survival at t = 3 from
## ilm_survival() is compared with the same quantity built from the fit's own
## parts: the fixed effects, the smooth's penalised part at its estimated
## coefficients, and for "fitted" the row's group's estimated effect. For
## "population" the survivor function is averaged over the group effect's
## normal distribution by integrate(). Three families of fit:
##   1. Weibull, time ~ s(x) + (1 | h), by REML;
##   2. the same without the smooth, time ~ x + (1 | h): the control;
##   3. flexible parametric (family = "rp"), time ~ s(x) + (1 | h).
## Then ilm_survival() with one time against two.
##
## Run from the repository root, with the build to be measured installed:
##   Rscript dev/measurements/survival-smooth/measure.R <label for the build>
suppressPackageStartupMessages({ library(illume); library(mgcv) })
lab <- commandArgs(trailingOnly = TRUE)
cat("build:", if (length(lab)) lab[1] else "(unlabelled)",
    "| illume", as.character(packageVersion("illume")),
    "| mgcv", as.character(packageVersion("mgcv")),
    "| R", paste0(R.version$major, ".", R.version$minor),
    "| date", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"), "\n")
q <- function(e) suppressMessages(suppressWarnings(e))

set.seed(31)
n <- 400
d <- data.frame(x = runif(n), h = factor(sample(letters[1:8], n, TRUE)))
hl <- rnorm(8, 0, 0.4)
eta <- 1 + 1.2 * sin(2 * pi * d$x) + hl[as.integer(d$h)]
tt <- exp(eta + 0.5 * log(rexp(n)))
ct <- rexp(n, 1 / (2 * median(tt)))
d$time <- pmin(tt, ct); d$event <- as.integer(tt <= ct)
nd <- data.frame(x = c(0.25, 0.5, 0.75), h = factor("a", levels = levels(d$h)))
tm <- 3

fit <- function(form, family) q(ilm_model(form, data = d, family = family,
                                         censor = ilm_surv(d$time, d$event),
                                         reml = TRUE, verbose = FALSE))

## The linear predictor at `nd` from the fit's parts, without the baseline:
## the fixed effects, every smooth's penalised part, and for "fitted" the
## estimated effect of the row's level of h.
eta_parts <- function(f, groups) {
  ndx <- illume:::ilm_newX(f, nd)
  keep <- if (is.null(f$rp)) seq_len(ncol(f$X)) else setdiff(seq_len(ncol(f$X)), f$rp$cols)
  e <- as.vector(ndx$X %*% f$beta[keep, 1L])
  pen <- numeric(nrow(nd))
  for (k in which(vapply(f$re, function(r) identical(r$kind, "basis"), TRUE)))
    pen <- pen + as.vector(illume:::ilm_basis_new(f, ndx, k) %*% illume:::ilm_Bhat_term(f, k))
  kh <- which(names(f$re) == "h")
  bh <- as.vector(illume:::ilm_Bhat_term(f, kh))[match(as.character(nd$h), f$re[[kh]]$levels)]
  list(fixed = e, pen = pen, h = if (groups == "fitted") bh else 0,
       vh = f$Sigma[["h"]][1, 1])
}
## the survivor function at t for a linear predictor e (baseline included)
S_of <- function(f) {
  if (is.null(f$rp)) {
    sc <- unname(f$dispersion[[1]])
    function(t, e) exp(-exp((log(t) - e) / sc))
  } else {
    base <- function(t) as.vector(illume:::ilm_rcs(log(t), f$rp$knots) %*% f$beta[f$rp$cols, 1L])
    function(t, e) f$family$surv(base(t) + e)
  }
}
by_hand <- function(f, groups, drop_pen = FALSE) {
  p <- eta_parts(f, groups); S <- S_of(f)
  e <- p$fixed + (if (drop_pen) 0 else p$pen) + p$h
  if (groups != "population") return(S(tm, e))
  vapply(e, function(ei) integrate(function(z) S(tm, ei + sqrt(p$vh) * z) * dnorm(z),
                                   -Inf, Inf, rel.tol = 1e-10)$value, 0)
}
S3 <- function(r) r$surv[r$time == tm]
from_ilm <- function(f, groups, times = c(2, tm))
  tryCatch(S3(q(ilm_survival(f, newdata = nd, times = times, groups = groups))),
           error = function(e) structure(rep(NA_real_, nrow(nd)), msg = conditionMessage(e)))

compare <- function(title, f) {
  cat("\n==", title, "==\n")
  if (any(vapply(f$re, function(r) identical(r$kind, "basis"), TRUE))) {
    p <- eta_parts(f, "fitted")
    cat("smooth's penalised part at x = 0.25, 0.5, 0.75 (link scale):", signif(p$pen, 5), "\n")
  }
  cat(sprintf("%-11s %-34s %-34s %s\n", "groups", "ilm_survival() S(3)", "by hand S(3)", "max abs diff"))
  for (gr in c("typical", "fitted", "population")) {
    a <- from_ilm(f, gr); b <- by_hand(f, gr)
    if (!is.null(attr(a, "msg"))) cat("  ilm_survival() ERROR:", attr(a, "msg"), "\n")
    cat(sprintf("%-11s %-34s %-34s %s\n", gr, paste(signif(a, 6), collapse = " "),
                paste(signif(b, 6), collapse = " "), signif(max(abs(a - b)), 3)))
  }
  invisible(NULL)
}

f1 <- fit(time ~ s(x, k = 8) + (1 | h), "weibull")
compare("1. Weibull, time ~ s(x, k = 8) + (1 | h)", f1)
cat("by hand without the penalised part, \"fitted\":",
    signif(by_hand(f1, "fitted", drop_pen = TRUE), 6), "\n")
cat("by hand without the penalised part, \"population\":",
    signif(by_hand(f1, "population", drop_pen = TRUE), 6), "\n")
f2 <- fit(time ~ x + (1 | h), "weibull")
compare("2. control: Weibull, time ~ x + (1 | h)", f2)
f3 <- fit(time ~ s(x, k = 8) + (1 | h), "rp")
compare("3. flexible parametric, time ~ s(x, k = 8) + (1 | h)", f3)
cat("by hand without the penalised part, \"typical\":",
    signif(by_hand(f3, "typical", drop_pen = TRUE), 6), "\n")
cat("by hand without the penalised part, \"fitted\":",
    signif(by_hand(f3, "fitted", drop_pen = TRUE), 6), "\n")

cat("\n== one time against two (fit 1, then fit 3) ==\n")
for (f in list(f1, f3)) for (gr in c("typical", "fitted", "population")) for (t_ in list(tm, c(2, tm))) {
  a <- from_ilm(f, gr, t_)
  cat(sprintf("  %-12s groups = %-12s times = %-9s %s\n",
              if (is.null(f$rp)) "Weibull" else "rp", dQuote(gr, FALSE), deparse(t_),
              if (!is.null(attr(a, "msg"))) paste("ERROR:", attr(a, "msg"))
              else paste(signif(a, 6), collapse = " ")))
}
