## Generate a grouped pkgdown reference index from the actual Rd topics, so
## that no topic can be left out of the index by hand.
## Run from the package root.

rd <- list.files("man", pattern = "[.]Rd$")
topic <- sub("[.]Rd$", "", rd)
## pkgdown indexes by the Rd \name, which is what the filename encodes except
## where roxygen escaped a character; read the name to be certain
nm <- vapply(file.path("man", rd), function(f) {
  x <- readLines(f, warn = FALSE)
  i <- grep("^\\\\name\\{", x)[1]
  if (is.na(i)) NA_character_ else sub("^\\\\name\\{(.*)\\}.*$", "\\1", x[i])
}, character(1))
topic <- ifelse(is.na(nm), topic, nm)
topic <- setdiff(topic, "illume-package")

## ---- groups, in the order they should appear -------------------------------
## Each entry: title, a one-line description, a predicate over topic names,
## and optionally the topics to list first, in that order -- the front door
## of a group before its parts. Whatever else matches follows alphabetically.
##
## Edit THIS file, not _pkgdown.yml: the yml is regenerated from it, and a
## hand edit there is lost the next time an export is added.
grp <- list(
  list("Fitting models",
       "One engine, every family. The formula is lme4's and the object that comes back is the same whatever was fitted.",
       function(x) x %in% c("ilm_model", "ilm_fit", "ilm_family",
         "ilm_model_formula", "ilm_surv", "ilm_censor", "ilm_ar1", "ilm_car1",
         "ilm_fourier", "ilm_cyclic", "ilm_squeeze", "ilm_thresholds")),

  list("Other designs",
       "Identification strategies and sampling designs that are not a single regression.",
       function(x) grepl("^ilm_(iv|did|rdd|design|svy|mediate)", x)),

  list("Diagnostics",
       "Each check reports whether an assumption is consistent with the data, and names a remedy that exists in this package when it is not.",
       function(x) grepl("^ilm_(check_|appraise|rqr|binned|calibration|consistency|variogram|re_mahalanobis|scores)", x)),

  list("ANOVA",
       "Factorial and repeated-measures designs, specified by naming columns rather than by writing a formula with an error term.",
       function(x) grepl("^ilm_aov", x)),

  list("Inference and interpretation",
       "What the model says, on a scale someone can read.",
       function(x) grepl("^ilm_(anova|effects|emmeans|contrast|trends|ame|robust|vcov_cluster|denom_df|pb_lrt|rp_lrt|coef_table|se_fixef|zi_|scenario|interpret|moderation)", x) ||
                   x == "ilm_plot_moderation"),

  list("Design and power",
       "Before the data exist.",
       function(x) grepl("^ilm_power|^plot[.]ilm_power$|^ilm_scaffold$", x),
       c("ilm_scaffold", "ilm_power_design", "ilm_power", "ilm_power_n",
         "plot.ilm_power")),

  ## Describing data, outliers and anomalies, and structure (profile,
  ## cluster, reduce) are illumex's, and indexed on its site.

  list("Missing data",
       "Filling it in honestly, and pooling across the imputations. Describing it is illumex's: ilm_check_missing().",
       function(x) grepl("^ilm_(impute|mi_pool)", x)),

  list("Causal models",
       "A DAG, what it implies, and what it licenses you to say.",
       function(x) grepl("^ilm_(dag|adjust_sets|dsep)", x)),

  list("Plots",
       "Plots of a fitted model, built on tinyplot. Plots of the data themselves are illumex's.",
       function(x) grepl("^ilm_plot", x)),

  list("Simulation",
       "Draws from a fitted model. Example data with a known structure, ilm_sim(), is illumex's.",
       function(x) grepl("^ilm_(simulate|fitted|survival)$", x)),

  list("Working with other packages",
       "Registration shims and the methods that let the wider ecosystem dispatch on an illume fit.",
       function(x) grepl("^ilm_register|[.]ilm_model$", x))
)

assigned <- character(0)
lines <- c("reference:")
for (g in grp) {
  hit <- Filter(function(x) isTRUE(tryCatch(g[[3]](x), error = function(e) FALSE)),
                setdiff(topic, assigned))
  hit <- sort(hit)
  if (length(g) >= 4L) hit <- c(intersect(g[[4]], hit), setdiff(hit, g[[4]]))
  if (!length(hit)) next
  assigned <- c(assigned, hit)
  lines <- c(lines,
             paste0("- title: \"", g[[1]], "\""),   # quoted: titles may hold a colon
             paste0("  desc: >"),
             paste0("    ", g[[2]]),
             "  contents:",
             paste0("  - ", hit))
}

left <- sort(setdiff(topic, assigned))
if (length(left)) {
  lines <- c(lines, "- title: \"Shared parameters and print methods\"",
             "  desc: >",
             "    Documentation shared across functions, and methods you call by printing rather than by name.",
             "  contents:",
             paste0("  - ", left))
}

head <- c(
  "url: https://huttoncp.github.io/illume/",
  "template:",
  "  bootstrap: 5",
  "  bslib:",
  "    primary: \"#2a6f97\"",
  "",
  "home:",
  "  title: A Unified Engine for Exploration and Frequentist Inference",
  "",
  "navbar:",
  "  structure:",
  "    left:  [intro, reference, articles, news, illumex]",
  "    right: [search, github]",
  "  components:",
  "    illumex:",
  "      text: Exploration (illumex)",
  "      href: https://huttoncp.github.io/illumex/",
  "",
  "articles:",
  "- title: Start here",
  "  navbar: ~",
  "  contents:",
  "  - illume",
  "  - workflow",
  "  - benchmarking",
  "- title: Modelling",
  "  navbar: Modelling",
  "  contents:",
  "  - regression-models",
  "  - causal-models",
  "  - effect-size-and-power",
  "  - anova",
  "  - moderation",
  "  - missing-data",
  "")

writeLines(c(head, lines), "_pkgdown.yml")
cat("topics:", length(topic), " assigned:", length(assigned),
    " unassigned:", length(left), "\n")
if (length(left)) cat("UNASSIGNED:", paste(left, collapse = ", "), "\n")
