## The measured identity between a factor-smooth interaction's random blocks,
## as illume stores them from mgcv::smooth2random(), and mgcv::PredictMat()
## on the fitting data:
##   cbind(rand) == (PredictMat(sm, data) %*% diag(trans.D))[, order(rind)]
## Found by measurement on the objects mgcv returns; mgcv's source was not
## read, nor its help pages. Evidence for the entry in inst/PROVENANCE.md.
## Run from the repository root, with illume installed:
##   Rscript dev/provenance/fs-predict/identity.R
## The output kept beside it was made with illume at 62e3e6a (0.0.8.9003),
## illumex at 02e118b, mgcv 1.9.1 and R 4.4.3.
suppressPackageStartupMessages({ library(illume); library(mgcv) })
source(file.path("dev", "provenance", "fs-predict", "sim.R"))
cat("illume", as.character(packageVersion("illume")), "| mgcv",
    as.character(packageVersion("mgcv")), "|", R.version.string, "\n")
cat("date", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"), "\n\n")

## the mapping, and the one tried first that does not hold
ident <- function(sob, data) {
  X0 <- PredictMat(sob$sm, data); R <- do.call(cbind, sob$rand)
  D <- sob$re$trans.D; rind <- sob$re$rind
  c(order_rind = max(abs(R - (X0 %*% diag(D, length(D)))[, order(rind), drop = FALSE])),
    rind_first_guess = max(abs(R - X0[, rind, drop = FALSE] %*% diag(D[rind], length(D)))))
}
## illume's own stored smooth, from a fit
from_fit <- function(form, d) {
  f <- suppressMessages(suppressWarnings(
    ilm_model(form, data = d, family = "gaussian", reml = TRUE, verbose = FALSE)))
  i <- grep("bs = \"fs\"|id\\)", names(f$smooths))
  sob <- f$smooths[[if (length(i)) i[1] else length(f$smooths)]]
  list(sob = sob, data = d)
}
## the same construction illume makes (ilm_smooth(): smoothCon with
## absorb.cons and scale.penalty, then smooth2random type 2), for set-ups
## illume may not fit
direct <- function(spec, d) {
  sm <- smoothCon(spec, data = d, absorb.cons = TRUE, scale.penalty = TRUE)[[1]]
  re <- smooth2random(sm, "", type = 2)
  list(sob = list(sm = sm, re = re, rand = re$rand), data = d)
}

d30 <- sim(1); d5 <- sim(8, n_id = 5, n_per = 12); d60 <- sim(9, n_id = 60, n_per = 6)
du <- sim(5, unbal = TRUE); dz <- sim(10); dz$z <- runif(nrow(dz), 0.5, 1.5)
cases <- list(
  "fit: 30 levels, k=5, m=2"               = from_fit(y ~ s(x, id, bs = "fs", k = 5), d30),
  "fit: 30 levels, k=3, m=2"               = from_fit(y ~ s(x, id, bs = "fs", k = 3), d30),
  "fit: 30 levels, k=8, m=1, with s(x)"    = from_fit(y ~ s(x) + s(x, id, bs = "fs", k = 8, m = 1), d30),
  "fit: 5 levels, k=4"                     = from_fit(y ~ s(x, id, bs = "fs", k = 4), d5),
  "fit: 60 levels, k=5"                    = from_fit(y ~ s(x, id, bs = "fs", k = 5), d60),
  "fit: unbalanced, k=6, cr basis"         = from_fit(y ~ s(x, id, bs = "fs", k = 6, xt = list(bs = "cr")), du),
  "direct: 30 levels, k=5, m=3"            = direct(s(x, id, bs = "fs", k = 6, m = 3), d30),
  "direct: numeric by z, k=5"              = direct(s(x, id, bs = "fs", k = 5, by = z), dz),
  "direct: 2-d smooth, tp, k=10"           = { dd <- d30; dd$w <- runif(nrow(dd)); direct(s(x, w, id, bs = "fs", k = 10), dd) })
cat(sprintf("%-38s %22s %22s\n", "case", "order(rind) mapping", "rind (first guess)"))
for (n in names(cases)) {
  r <- tryCatch(ident(cases[[n]]$sob, cases[[n]]$data), error = function(e) conditionMessage(e))
  if (is.character(r)) cat(sprintf("%-38s ERROR: %s\n", n, r)) else
    cat(sprintf("%-38s %22.3g %22.3g\n", n, r[1], r[2]))
}
