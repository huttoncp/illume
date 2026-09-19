## ---------------------------------------------------------------------------
## Enforce the retention policy: keep raw per-replicate data for the three most
## recent runs of each study, and for older runs keep only the generated
## findings.
##
## "Three most recent" is per STUDY, not overall.  Counting globally would mean
## three re-runs of the coverage study silently discarded every benchmark and
## comparison result, which is not what a retention policy is for.
##
## A run is only ever deleted AFTER its findings have been written and verified
## to exist on disk.  If summarising fails, the run is kept and reported.
##
## IMPORTANT LIMIT: this prunes the WORKING TREE, not git history.  Blobs that
## have been committed stay in .git forever, so the repository keeps growing
## even as the checkout stays small.  This policy keeps the working tree and
## fresh clones manageable and the directory comprehensible; it is not a way to
## cap repository size.  For that you would need to avoid committing raw .rds
## at all, use git-lfs, or rewrite history -- all bigger decisions.
##
## Usage: Rscript prune_runs.R <studies_dir> [keep] [--dry-run]
## ---------------------------------------------------------------------------

a <- commandArgs(trailingOnly = TRUE)
studies <- if (length(a) >= 1 && !startsWith(a[1], "--")) a[1] else "."
keep    <- suppressWarnings(as.integer(a[2]))
if (is.na(keep)) keep <- 3L
dry     <- any(a == "--dry-run")

src <- file.path(studies, "scripts", "summarise_run.R")
if (!file.exists(src)) stop("cannot find summarise_run.R next to this script")
source(src, local = TRUE)

runs_root <- file.path(studies, "runs")
if (!dir.exists(runs_root)) stop("no runs/ directory under ", studies)

dir_bytes <- function(d) {
  f <- list.files(d, recursive = TRUE, full.names = TRUE)
  if (!length(f)) return(0)
  sum(file.info(f)$size, na.rm = TRUE)
}

cat(sprintf("Retention: %d most recent run(s) per study%s\n\n", keep,
            if (dry) "   [DRY RUN - nothing will be deleted]" else ""))

freed <- 0
for (sdir in list.dirs(runs_root, recursive = FALSE)) {
  study <- basename(sdir)
  rs <- sort(basename(list.dirs(sdir, recursive = FALSE)), decreasing = TRUE)
  if (!length(rs)) next
  keepers <- utils::head(rs, keep)
  drops   <- setdiff(rs, keepers)

  ## refresh findings for every run, so the permanent record is current even
  ## for runs that are being retained
  for (r in rs) write_findings(studies, study, r)

  cat(sprintf("%-10s %d run(s): keeping %s\n", study, length(rs),
              paste(keepers, collapse = ", ")))
  if (!length(drops)) next

  fp <- file.path(studies, "findings", paste0(study, ".md"))
  rec <- if (file.exists(fp)) readLines(fp, warn = FALSE) else character(0)
  for (r in drops) {
    ## never discard raw data whose findings are not actually on disk
    if (!any(startsWith(rec, paste0("## ", r)))) {
      cat(sprintf("           KEEPING %s - no findings recorded, refusing to delete\n", r))
      next
    }
    sz <- dir_bytes(file.path(sdir, r))
    cat(sprintf("           pruning %s (%.1f KB, findings retained)\n", r, sz / 1024))
    if (!dry) unlink(file.path(sdir, r), recursive = TRUE, force = TRUE)
    freed <- freed + sz
  }
}
cat(sprintf("\n%s %.1f KB\n", if (dry) "Would free" else "Freed", freed / 1024))
if (!dry && freed > 0)
  cat("Note: committed data remains in git history; this shrinks the working tree only.\n")
