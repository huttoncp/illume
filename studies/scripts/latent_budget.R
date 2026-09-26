## B6: how much information per latent value does the Laplace approximation
## need, for families with one linear predictor?
##
## The latent_budget and obs_per_ar_latent checks count ROWS per latent value,
## with thresholds calibrated on multinomial data, where one categorical
## observation says little about its latent. A Poisson count with a large mean
## says a great deal; a binary outcome at p = 0.1 very little. So the proposal
## is to judge the expected Fisher information per latent value instead:
##   pre-fit   rows per latent x the Fisher weight at the response's marginal
##             mean -- binomial p(1 - p), Poisson the mean, negative binomial
##             mu / (1 + mu / k) with k by the method of moments
##   post-fit  each latent value's own information, the Fisher weights at the
##             fitted conditional means summed over its rows (median across
##             latent values), recorded to validate the pre-fit version
## and to find where on it the approximation starts to fail: the latent SD
## attenuated, a slope's Wald interval under-covering, fits unusable.
##
## Usage: Rscript latent_budget.R <nrep> <ncore> <outdir>
## ---------------------------------------------------------------------------
args  <- commandArgs(trailingOnly = TRUE)
NREP  <- if (length(args) >= 1) as.integer(args[1]) else 20L
NCORE <- if (length(args) >= 2) as.integer(args[2]) else 4L
sp    <- if (length(args) >= 3) args[3] else "."
dir.create(sp, showWarnings = FALSE, recursive = TRUE)

responses <- data.frame(
  resp   = c("bin50", "bin10", "pois05", "pois3", "nb3"),
  family = c("binomial", "binomial", "poisson", "poisson", "nbinom"),
  b0     = c(0, stats::qlogis(0.1), log(0.5), log(3), log(3)),
  k      = c(NA, NA, NA, NA, 2), stringsAsFactors = FALSE)
structures <- data.frame(
  structure = c("ar1", "ar1", "ar1", "rw1", "rw1", "ri", "ri", "ri"),
  rows      = c(1, 2, 4, 1, 2, 2, 4, 8), stringsAsFactors = FALSE)
cells <- merge(merge(responses, structures), data.frame(sd = c(0.5, 1)))
cells$cell <- seq_len(nrow(cells))
jobs <- merge(cells, data.frame(rep = seq_len(NREP)))

one <- function(job) {
  suppressMessages(library(illume))
  seed <- 104729L * job$cell + job$rep
  set.seed(seed)
  G <- 20L; Tn <- 15L; L <- 60L
  if (job$structure %in% c("ar1", "rw1")) {
    cl <- expand.grid(t = seq_len(Tn), g = seq_len(G))
    lat <- unlist(lapply(seq_len(G), function(i) {
      if (job$structure == "ar1")
        as.numeric(stats::arima.sim(list(ar = 0.7), Tn,
                                    sd = job$sd * sqrt(1 - 0.49)))
      else c(0, cumsum(stats::rnorm(Tn - 1L, 0, job$sd / 2)))
    }))
    idx <- rep(seq_len(nrow(cl)), each = job$rows)
    d <- data.frame(t = cl$t[idx], g = factor(cl$g[idx]))
    true_sd <- if (job$structure == "ar1") job$sd else job$sd / 2
  } else {
    idx <- rep(seq_len(L), each = job$rows)
    lat <- stats::rnorm(L, 0, job$sd)
    d <- data.frame(g = factor(idx))
    true_sd <- job$sd
  }
  n <- nrow(d)
  d$x <- stats::rnorm(n)
  eta <- job$b0 + 0.3 * d$x + lat[idx]
  d$y <- switch(job$family,
    binomial = stats::rbinom(n, 1, stats::plogis(eta)),
    poisson  = stats::rpois(n, exp(eta)),
    nbinom   = stats::rnbinom(n, size = job$k, mu = exp(eta)))
  ar <- switch(job$structure, ar1 = ilm_ar1(~ t | g), rw1 = ilm_rw1(~ t | g),
               ri = NULL)
  fml <- if (job$structure == "ri") y ~ x + (1 | g) else y ~ x
  nlat <- length(lat)

  ## pre-fit information per latent value
  yb <- mean(d$y)
  w_pre <- switch(job$family,
    binomial = yb * (1 - yb),
    poisson  = yb,
    nbinom   = {
      v <- stats::var(d$y)
      k_mm <- if (v > yb) yb^2 / (v - yb) else Inf
      yb / (1 + yb / k_mm)
    })
  info_pre <- n / nlat * w_pre

  out <- data.frame(job[c("cell", "rep", "resp", "family", "structure",
                          "rows", "sd")],
                    seed = seed, n = n, nlat = nlat, true_sd = true_sd,
                    info_pre = info_pre, ok_fit = FALSE, usable = FALSE,
                    info_post = NA_real_, sd_hat = NA_real_, rho_hat = NA_real_,
                    slope = NA_real_, slope_se = NA_real_, cover = NA,
                    lb_status = NA_character_, ar_status = NA_character_,
                    seconds = NA_real_, stringsAsFactors = FALSE)
  t0 <- proc.time()[[3]]
  f <- tryCatch(suppressMessages(suppressWarnings(
    ilm_model(fml, data = d, family = job$family, ar = ar, verbose = FALSE))),
    error = function(e) NULL)
  out$seconds <- proc.time()[[3]] - t0
  if (is.null(f)) return(out)
  out$ok_fit <- TRUE
  out$usable <- isTRUE(illume:::ilm_fixed_usable(f)) &&
    !any(f$checks$check %in% c("gradient", "optimizer") &
           f$checks$status == "FAIL")
  ck <- f$checks
  out$lb_status <- ck$status[ck$check == "latent_budget"][1]
  if (any(ck$check == "obs_per_ar_latent"))
    out$ar_status <- ck$status[ck$check == "obs_per_ar_latent"][1]
  v <- ilm_varcorr(f)
  out$sd_hat <- if (job$structure == "ar1") sqrt(v$latent$Sigma[1, 1])
                else if (job$structure == "rw1") sqrt(v$latent$var_per_time[1, 1])
                else sqrt(v$re$g[1, 1])
  if (job$structure == "ar1") out$rho_hat <- v$latent$rho
  b <- coef(f)[["x"]]; se <- sqrt(diag(vcov(f)))[["x"]]
  out$slope <- b; out$slope_se <- se
  out$cover <- is.finite(se) && abs(b - 0.3) <= stats::qnorm(0.975) * se
  ## post-fit: the Fisher weight at each row's fitted conditional mean,
  ## summed over each latent value's rows
  ## (the conditional linear predictor: fixed effects, each group's effect
  ## and each cell's latent value, as the fit has them)
  mu <- tryCatch(as.numeric(f$family$linkinv(illume:::ilm_eta_hat(f, TRUE)[, 1])),
                 error = function(e) NULL)
  if (!is.null(mu) && length(mu) == n) {
    w <- switch(job$family,
      binomial = mu * (1 - mu),
      poisson  = mu,
      nbinom   = { k <- f$dispersion[[1]]; mu / (1 + mu / k) })
    out$info_post <- stats::median(tapply(w, idx, sum))
  }
  out
}

t0 <- Sys.time()
cl <- parallel::makeCluster(NCORE)
res <- do.call(rbind, parallel::parLapplyLB(cl, split(jobs, seq_len(nrow(jobs))), one))
parallel::stopCluster(cl)
utils::write.csv(res, file.path(sp, "latent_budget.csv"), row.names = FALSE)
cat("fits:", nrow(res), " errors:", sum(!res$ok_fit), " minutes:",
    round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "\n")

## by cell: the approximation's behaviour against the information per latent
agg <- do.call(rbind, lapply(split(res, res$cell), function(x) {
  u <- x[x$usable, ]
  data.frame(resp = x$resp[1], structure = x$structure[1], rows = x$rows[1],
             sd = x$sd[1], info_pre = stats::median(x$info_pre),
             info_post = stats::median(x$info_post, na.rm = TRUE),
             usable = mean(x$usable),
             sd_ratio = stats::median(u$sd_hat / u$true_sd),
             cover = mean(u$cover, na.rm = TRUE),
             rho = if (x$structure[1] == "ar1") stats::median(u$rho_hat) else NA,
             lb = paste(names(table(x$lb_status)), table(x$lb_status), collapse = " "),
             stringsAsFactors = FALSE)
}))
agg <- agg[order(agg$info_pre), ]
utils::write.csv(agg, file.path(sp, "latent_budget_cells.csv"), row.names = FALSE)
print(agg, digits = 3, row.names = FALSE)
