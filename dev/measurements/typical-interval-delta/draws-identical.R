## ilm_draws() before and after its set-up was moved into a function of its own
## (so that the typical curve's covariance is read from the same distribution
## the draws come from): the draws, their modes and their layout, for fits of
## several kinds and each `given`, should be identical() between two builds.
##   Rscript dev/measurements/typical-interval-delta/draws-identical.R save <file.rds> <label for the build>
##   Rscript dev/measurements/typical-interval-delta/draws-identical.R compare <a.rds> <b.rds>
args <- commandArgs(trailingOnly = TRUE)
if (args[1] == "compare") {
  a <- readRDS(args[2]); b <- readRDS(args[3])
  stopifnot(identical(names(a), names(b)))
  cat(sprintf("%-46s %s\n", "case", "identical"))
  for (k in names(a)) cat(sprintf("%-46s %s\n", k, identical(a[[k]], b[[k]])))
  cat("all identical:", all(vapply(names(a), function(k) identical(a[[k]], b[[k]]), TRUE)), "\n")
  quit(save = "no")
}
suppressPackageStartupMessages(library(illume))
cat("build:", args[3], "| illume", as.character(packageVersion("illume")), "| R",
    paste0(R.version$major, ".", R.version$minor), "| date",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"), "\n")
q <- function(e) suppressMessages(suppressWarnings(e))
mk <- function(seed, gsd, n = 250) {
  set.seed(seed)
  d <- data.frame(x = runif(n), z = rnorm(n), h = factor(sample(letters[1:8], n, TRUE)))
  hl <- rnorm(8, 0, gsd)
  eta <- 1 + 0.5 * d$z + 1.2 * sin(2 * pi * d$x) + hl[as.integer(d$h)]
  tt <- exp(eta + 0.5 * log(rexp(n)))
  ct <- rexp(n, 1 / (2 * median(tt)))
  d$time <- pmin(tt, ct); d$event <- as.integer(tt <= ct)
  d
}
fit <- function(d, form, family, reml = FALSE)
  q(ilm_model(form, data = d, family = family, censor = ilm_surv(d$time, d$event),
              reml = reml, verbose = FALSE))
d1 <- mk(41, 0.4); d0 <- mk(45, 0)           # the second ends at the boundary: a held direction
fits <- list(
  "weibull, smooth + (1 | h), REML"   = fit(d1, time ~ s(x, k = 8) + z + (1 | h), "weibull", TRUE),
  "weibull, (1 | h), ML, joint"       = q(ilm_model(time ~ z + (1 | h), data = d1, family = "weibull",
                                                    censor = ilm_surv(d1$time, d1$event), joint = TRUE, verbose = FALSE)),
  "rp, smooth + (1 | h), ML"          = fit(d1, time ~ s(x, k = 8) + z + (1 | h), "rp"),
  "weibull, smooth + (1 | h), boundary" = fit(d0, time ~ s(x, k = 8) + z + (1 | h), "weibull"),
  "weibull, factor-smooth, REML"      = fit(d1, time ~ z + s(x, k = 5) + s(x, h, bs = "fs", k = 4), "weibull", TRUE))
out <- list()
for (k in names(fits)) for (g in c("none", "theta", "parameters")) {
  r <- ilm_draws(fits[[k]], nsim = 50L, seed = 1L, given = g)
  out[[paste(k, "| given =", g)]] <- list(draws = r$draws, mode = r$mode, map = r$map,
                                          natural = r$natural, held = r$held)
}
r <- ilm_draws(fits[[1]], nsim = 50L, seed = 1L, blocks = c("beta", "bvec"), natural = FALSE)
out[["weibull, smooth + (1 | h), REML | blocks"]] <- list(draws = r$draws, mode = r$mode, map = r$map)
saveRDS(out, args[2]); cat("saved", length(out), "cases to", args[2], "\n")
