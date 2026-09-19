## ---------------------------------------------------------------------------
## Enforce the retention policy: keep the raw per-replicate data from the three
## most recent PACKAGE VERSIONS, and for older versions keep only the generated
## findings.
##
## The unit of retention is a version's whole SET of studies, not each study
## separately. Results are evidence about a particular build of the package, so
## coverage from version 5 and benchmarks from version 3 are not a coherent
## picture; sets are kept or dropped together.
##
## A set is only ever deleted AFTER every study in it has findings written and
## verified present on disk. If any is missing, the set is kept and reported.
##
## IMPORTANT LIMIT: this prunes the WORKING TREE, not git history. Blobs that
## have been committed stay in .git forever, so the repository keeps growing
## even as the checkout stays small. This policy keeps the working tree and
## fresh clones manageable and the directory comprehensible; it is not a way to
## cap repository size. For that you would need to avoid committing raw .rds at
## all, use git-lfs, or rewrite history -- all bigger decisions.
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

has_findings <- function(study, version) {
  fp <- file.path(studies, "findings", paste0(study, ".md"))
  if (!file.exists(fp)) return(FALSE)
  any(startsWith(readLines(fp, warn = FALSE), paste0("## ", version)))
}

versions <- ver_sort(basename(list.dirs(runs_root, recursive = FALSE)))
if (!length(versions)) { cat("no run sets found\n"); quit(save = "no") }

cat(sprintf("Retention: %d most recent version(s)%s\n\n", keep,
            if (dry) "   [DRY RUN - nothing will be deleted]" else ""))

## refresh findings for every set, so the permanent record is current even for
## sets being retained
for (v in versions)
  for (s in list.dirs(file.path(runs_root, v), recursive = FALSE))
    write_findings(studies, basename(s), v)

keepers <- utils::head(versions, keep)
drops   <- setdiff(versions, keepers)

for (v in versions) {
  studs <- basename(list.dirs(file.path(runs_root, v), recursive = FALSE))
  cat(sprintf("  %-14s %2d stud%s  %7.1f KB  %s\n", v, length(studs),
              if (length(studs) == 1) "y" else "ies",
              dir_bytes(file.path(runs_root, v)) / 1024,
              if (v %in% keepers) "KEEP" else "prune"))
}
cat("\n")

freed <- 0
for (v in drops) {
  studs <- basename(list.dirs(file.path(runs_root, v), recursive = FALSE))
  missing <- studs[!vapply(studs, has_findings, TRUE, version = v)]
  if (length(missing)) {
    cat(sprintf("%-14s KEEPING - no findings recorded for: %s\n", v,
                paste(missing, collapse = ", ")))
    next
  }
  sz <- dir_bytes(file.path(runs_root, v))
  cat(sprintf("%-14s pruning (%.1f KB; findings for %d stud%s retained)\n", v,
              sz / 1024, length(studs), if (length(studs) == 1) "y" else "ies"))
  if (!dry) unlink(file.path(runs_root, v), recursive = TRUE, force = TRUE)
  freed <- freed + sz
}
cat(sprintf("\n%s %.1f KB\n", if (dry) "Would free" else "Freed", freed / 1024))
if (!dry && freed > 0)
  cat("Note: committed data remains in git history; this shrinks the working tree only.\n")
