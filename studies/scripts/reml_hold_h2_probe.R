## A reading made AFTER the hold-rule study ran (2026-09-29), not part of its
## pre-registration: why rule C's SE differs from the intercept-only refit's
## in arm H2. Cell 3's data for its first 40 replicates, fitted by REML and by
## ML, full and reduced; for each held fit, the ratio and the boundary it
## reached. Usage: Rscript reml_hold_h2_probe.R <run dir with reml_hold.R> <out csv>
## what was held.
suppressPackageStartupMessages(library(illume))
src <- readLines(file.path(commandArgs(TRUE)[1], "reml_hold.R"))
## the study's own data generator, and nothing else from the script
i <- grep("^gen <- function", src); j <- grep("^## ---- one replicate", src)
cells <- rbind(data.frame(arm = "H1", G = c(8, 20), m = 5),
               data.frame(arm = "H2", expand.grid(G = c(8, 20), m = c(5, 10))))
cells$cell <- seq_len(nrow(cells)); TRUTH <- 0.5
eval(parse(text = src[i:(j - 1L)]))
q <- function(e) tryCatch(suppressMessages(suppressWarnings(e)), error = function(e) NULL)
se_x <- function(f) if (is.null(f)) NA_real_ else ilm_coef_table(f)["x", "Std. Error"]
rows <- list()
for (rep in 1:40) {
  ce <- cells[cells$cell == 3, ]
  d <- gen(ce, 104729L * 3L + rep)
  out <- list(rep = rep)
  for (est in c("reml", "ml")) {
    r <- est == "reml"
    f <- q(ilm_model(y ~ x + (1 + x | g), data = d, family = "gaussian", reml = r, verbose = FALSE))
    fx <- q(ilm_model(y ~ x + (1 | g), data = d, family = "gaussian", reml = r, verbose = FALSE))
    out[[paste0(est, "_held")]] <- !is.null(f) && length(f$hessian_held) > 0L
    out[[paste0(est, "_how")]] <- if (is.null(f)) NA else f$hessian_how
    out[[paste0(est, "_ratio")]] <- se_x(f) / se_x(fx)
    S <- if (!is.null(f)) as.matrix(f$Sigma_d[["g"]]) else NULL
    out[[paste0(est, "_sd_slope")]] <- if (is.null(S)) NA else sqrt(max(S[2, 2], 0))
    out[[paste0(est, "_cor")]] <- if (is.null(S) || S[2, 2] <= 0) NA else S[1, 2] / sqrt(S[1, 1] * S[2, 2])
    out[[paste0(est, "_sd_int")]] <- if (is.null(S)) NA else sqrt(max(S[1, 1], 0))
    out[[paste0(est, "_sd_int_red")]] <- if (is.null(fx)) NA else sqrt(as.matrix(if (is.null(fx$Sigma_d[["g"]])) fx$Sigma[["g"]] else fx$Sigma_d[["g"]])[1, 1])
  }
  rows[[rep]] <- as.data.frame(out)
}
r <- do.call(rbind, rows)
print(r, digits = 3, row.names = FALSE)
h <- r[r$reml_held & r$ml_held, ]
cat("\nheld by both:", nrow(h), "\n")
cat("median ratio  REML", round(median(h$reml_ratio), 4), " ML", round(median(h$ml_ratio), 4), "\n")
cat("share within [0.97, 1.03]  REML", mean(abs(h$reml_ratio - 1) <= 0.03),
    " ML", mean(abs(h$ml_ratio - 1) <= 0.03), "\n")
utils::write.csv(r, commandArgs(TRUE)[2], row.names = FALSE)
## which boundary each held fit reached: the slope variance at zero (its SD
## below 1e-3 of the intercept's), or a correlation of +/-1 with a slope
for (est in c("reml", "ml")) {
  held <- r[[paste0(est, "_held")]] %in% TRUE
  zero <- r[[paste0(est, "_sd_slope")]] < 1e-3 * pmax(r[[paste0(est, "_sd_int")]], 1e-8)
  rat <- r[[paste0(est, "_ratio")]]
  for (k in c("slope at zero", "correlation at one")) {
    s <- held & (if (k == "slope at zero") zero else !zero)
    cat(sprintf("%-4s %-19s n %2d  median ratio %.4f  within 3%% %.2f\n", est, k, sum(s),
                if (any(s)) median(rat[s]) else NA, if (any(s)) mean(abs(rat[s] - 1) <= 0.03) else NA))
  }
}
