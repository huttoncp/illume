## The four V2 "false holds": is the dropped model's fit (the label) at its
## optimum? Refit it with more restarts, and from a different start.
suppressPackageStartupMessages(library(illume))
src <- readLines(commandArgs(TRUE)[1]); n0 <- grep("^## ---- the three studies", src); n1 <- grep("cells <- switch(STUDY", src, fixed = TRUE) - 1
eval(parse(text = src[n0:n1]))
q <- function(e) tryCatch(suppressMessages(suppressWarnings(e)), error = function(e) NULL)
for (cr in list(c(12, 32), c(12, 56), c(12, 85), c(3, 89))) {
  job <- merge(re_cells[re_cells$cell == cr[1], ], data.frame(rep = cr[2])); s <- gen_re(job)
  f  <- q(ilm_model(s$fml, data = s$d, family = s$family, ar = s$ar, verbose = FALSE))
  f0 <- q(ilm_model(s$drop_fml, data = s$d, family = s$family, ar = s$ar, verbose = FALSE))
  f0r <- q(ilm_model(s$drop_fml, data = s$d, family = s$family, ar = s$ar, verbose = FALSE, restarts = 10))
  cat(sprintf("cell %d rep %d (%s, %s): ll full %.4f | dropped %.4f (conv %d, %s) | dropped, 10 restarts %.4f (conv %d)\n",
      cr[1], cr[2], job$family, job$design, logLik(f), logLik(f0), f0$opt$convergence, f0$hessian_how,
      logLik(f0r), f0r$opt$convergence))
  if (!is.null(s$ar)) cat("   full rho", round(f$rho, 3), " dropped rho", round(f0$rho, 3), "\n")
}
